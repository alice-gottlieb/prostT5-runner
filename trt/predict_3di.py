"""
Standalone TensorRT 3Di predictor.

Reads an amino-acid FASTA, runs the ProstT5 3Di TensorRT engine for the current
GPU, and writes a 3Di FASTA (one 3Di string per input sequence). The output can
be fed to foldseek for structural search.

Usage (inside a job on the target GPU):
    uv run python trt/predict_3di.py \
        --engine $SCRATCH/prostt5-trt/engines/prostt5_3di.H200.fp16.engine \
        --cache-dir $SCRATCH/prostt5-trt/hf-cache \
        --in proteins.fasta --out proteins_3di.fasta
"""

import argparse
import time

import torch

from prostt5_3di import load_tokenizer, tokenize_one, logits_to_3di
from trt_runtime import ProstT5Engine


def read_fasta(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--in", dest="infile", required=True, help="Input AA FASTA")
    ap.add_argument("--out", required=True, help="Output 3Di FASTA")
    ap.add_argument("--lowercase", action="store_true",
                    help="Emit lowercase 3Di (ProstT5 convention) instead of upper")
    ap.add_argument("--max-len", type=int, default=2048,
                    help="Skip sequences longer than the engine's max profile length")
    args = ap.parse_args()

    assert torch.cuda.is_available(), "CUDA GPU required"
    print(f"GPU: {torch.cuda.get_device_name(0)}")

    tokenizer = load_tokenizer(args.cache_dir)
    engine = ProstT5Engine(args.engine)

    n, skipped, t0 = 0, 0, time.perf_counter()
    with open(args.out, "w") as out:
        for name, seq in read_fasta(args.infile):
            if len(seq) > args.max_len:
                skipped += 1
                continue
            ids, mask, real_len = tokenize_one(tokenizer, seq, device="cuda")
            logits = engine.infer_logits(ids, mask)
            three_di = logits_to_3di(logits.float().cpu(), real_len)
            if args.lowercase:
                three_di = three_di.lower()
            out.write(f">{name}\n{three_di}\n")
            n += 1
    dt = time.perf_counter() - t0
    print(f"Wrote {n} 3Di sequences to {args.out} "
          f"({skipped} skipped > {args.max_len} aa) in {dt:.1f}s "
          f"({1000*dt/max(n,1):.1f} ms/seq)")


if __name__ == "__main__":
    main()
