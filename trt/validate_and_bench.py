"""
Validate and benchmark a ProstT5 3Di TensorRT engine against the FP32 PyTorch
reference on the same GPU.

For each test sequence we compute the 3Di string two ways:
  * reference: fused FP32 encoder + CNN head in PyTorch,
  * engine:    the FP16 TensorRT engine,
then report per-residue agreement and per-sequence exact-match rate, plus a
latency benchmark (PyTorch FP16 vs TensorRT) so we can quote a speedup.

Usage (inside a job on the target GPU):
    uv run python trt/validate_and_bench.py \
        --engine $SCRATCH/prostt5-trt/engines/prostt5_3di.RTX2080Ti.fp16.engine \
        --cnn-ckpt $SCRATCH/prostt5-trt/cnn_head.pt \
        --cache-dir $SCRATCH/prostt5-trt/hf-cache \
        --fasta $SCRATCH/prostt5-trt/test_seqs.fasta \
        --runs 5
"""

import argparse
import time

import torch

from prostt5_3di import (
    load_model, load_tokenizer, tokenize_one, logits_to_3di,
)
from trt_runtime import ProstT5Engine

# Diverse real proteins used when no FASTA is supplied (human hemoglobin alpha,
# lysozyme C, ubiquitin, a short zinc finger).
DEFAULT_SEQS = {
    "HBA_HUMAN": (
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD"
        "ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    ),
    "LYSC_HUMAN": (
        "KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACH"
        "LSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV"
    ),
    "UBIQ_HUMAN": (
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    ),
    "ZNF_SHORT": "MERPYACPVESCDRRFSRSDELTRHIRIHTGQKP",
}


def read_fasta(path):
    seqs = {}
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def reference_3di(model, tokenizer, seq):
    ids, mask, real_len = tokenize_one(tokenizer, seq, device="cuda")
    with torch.no_grad():
        logits = model(ids, mask)  # (1, 20, L)
    return logits_to_3di(logits.float().cpu(), real_len)


def engine_3di(engine, tokenizer, seq):
    ids, mask, real_len = tokenize_one(tokenizer, seq, device="cuda")
    logits = engine.infer_logits(ids, mask)  # (20, L)
    return logits_to_3di(logits.float().cpu(), real_len)


def agreement(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 0.0, 0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / n, n


def bench(fn, seqs, runs):
    # warmup
    for s in seqs.values():
        fn(s)
    torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(runs):
        for s in seqs.values():
            fn(s)
    torch.cuda.synchronize()
    total = time.perf_counter() - t0
    n = runs * len(seqs)
    return total / n  # mean seconds per sequence


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--cnn-ckpt", required=True)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--fasta", default=None, help="Optional FASTA of test sequences")
    ap.add_argument("--max-seqs", type=int, default=20)
    ap.add_argument("--runs", type=int, default=5)
    args = ap.parse_args()

    assert torch.cuda.is_available(), "CUDA GPU required"
    gpu_name = torch.cuda.get_device_name(0)
    print(f"GPU: {gpu_name}")

    seqs = read_fasta(args.fasta) if args.fasta else dict(DEFAULT_SEQS)
    # Drop sequences longer than the engine's max profile length (2048) and cap count.
    seqs = {k: v for k, v in seqs.items() if len(v) <= 2048}
    seqs = dict(list(seqs.items())[: args.max_seqs])
    print(f"Validating on {len(seqs)} sequences (lengths "
          f"{min(len(v) for v in seqs.values())}-{max(len(v) for v in seqs.values())})")

    print("Loading FP32 reference model...")
    model = load_model(args.cnn_ckpt, cache_dir=args.cache_dir).to("cuda").eval()
    tokenizer = load_tokenizer(args.cache_dir)

    print("Loading TensorRT engine...")
    engine = ProstT5Engine(args.engine)

    # ---- correctness ----
    total_same, total_res, exact = 0, 0, 0
    print("\nPer-sequence 3Di agreement (TRT-FP16 vs PyTorch-FP32):")
    for name, seq in seqs.items():
        ref = reference_3di(model, tokenizer, seq)
        trt_out = engine_3di(engine, tokenizer, seq)
        frac, n = agreement(ref, trt_out)
        total_same += int(frac * n)
        total_res += n
        if ref == trt_out:
            exact += 1
        print(f"  {name:14s} len={len(seq):4d}  agreement={frac*100:6.2f}%"
              f"  {'EXACT' if ref == trt_out else ''}")
    overall = 100.0 * total_same / total_res if total_res else 0.0
    print(f"\nOverall per-residue agreement: {overall:.3f}%  "
          f"({total_same}/{total_res})")
    print(f"Exact-match sequences: {exact}/{len(seqs)}")

    # ---- benchmark ----
    print("\nBenchmarking...")
    model_fp16 = load_model(args.cnn_ckpt, cache_dir=args.cache_dir).half().to("cuda").eval()

    def pt_fp16(seq):
        ids, mask, real_len = tokenize_one(tokenizer, seq, device="cuda")
        with torch.no_grad():
            logits = model_fp16(ids, mask)
        return logits_to_3di(logits.float().cpu(), real_len)

    pt_t = bench(pt_fp16, seqs, args.runs)
    trt_t = bench(lambda s: engine_3di(engine, tokenizer, s), seqs, args.runs)
    print(f"PyTorch FP16 : {pt_t*1000:8.2f} ms/seq")
    print(f"TensorRT FP16: {trt_t*1000:8.2f} ms/seq")
    print(f"Speedup      : {pt_t/trt_t:6.2f}x")

    print("\nRESULT_SUMMARY "
          f"gpu='{gpu_name}' agreement={overall:.3f} exact={exact}/{len(seqs)} "
          f"pt_ms={pt_t*1000:.2f} trt_ms={trt_t*1000:.2f} speedup={pt_t/trt_t:.2f}")


if __name__ == "__main__":
    main()
