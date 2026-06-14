"""
Diagnostic: compare 3Di logits from three paths on a single sequence to locate
where the TensorRT pipeline diverges from the PyTorch reference.

  1. PyTorch FP32 fused model      (ground truth)
  2. onnxruntime on the ONNX file  (isolates the exported graph)
  3. the TensorRT engine           (isolates the TRT build)

If (2) already disagrees with (1), the ONNX export is wrong. If only (3)
disagrees, the TRT build/runtime is wrong (e.g. FP16 overflow, I/O dtype).

Usage (GPU job):
    uv run python trt/diag_compare.py \
        --engine $SCRATCH/prostt5-trt/engines/prostt5_3di.A100.fp16.engine \
        --onnx   $SCRATCH/prostt5-trt/prostt5_3di.onnx \
        --cnn-ckpt $SCRATCH/prostt5-trt/cnn_head.pt \
        --cache-dir $SCRATCH/prostt5-trt/hf-cache
"""

import argparse

import numpy as np
import torch

from prostt5_3di import load_model, load_tokenizer, tokenize_one, logits_to_3di
from trt_runtime import ProstT5Engine

SEQ = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVAD"
       "ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")


def stats(name, t):
    t = t.float()
    print(f"  {name}: shape={tuple(t.shape)} min={t.min():.3f} max={t.max():.3f} "
          f"mean={t.mean():.3f} nan={int(torch.isnan(t).sum())} inf={int(torch.isinf(t).sum())}")


def agree(a, b):
    n = min(len(a), len(b))
    return sum(a[i] == b[i] for i in range(n)) / n * 100 if n else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--onnx", required=True)
    ap.add_argument("--cnn-ckpt", required=True)
    ap.add_argument("--cache-dir", default=None)
    args = ap.parse_args()

    tokenizer = load_tokenizer(args.cache_dir)
    ids, mask, real_len = tokenize_one(tokenizer, SEQ, device="cpu")
    print(f"seq len={len(SEQ)} real_len={real_len} token_len={ids.shape[1]}")
    print("input_ids[:8] =", ids[0, :8].tolist())

    # 1. PyTorch FP32
    model = load_model(args.cnn_ckpt, cache_dir=args.cache_dir).to("cuda").eval()
    with torch.no_grad():
        ref = model(ids.cuda(), mask.cuda()).cpu()  # (1,20,L)
    print("\n[1] PyTorch FP32")
    stats("logits", ref)
    ref_3di = logits_to_3di(ref, real_len)
    print("  argmax[:40]:", ref.argmax(1)[0, :40].tolist())
    print("  3Di:", ref_3di[:60])

    # 2. onnxruntime (CPU) on the exported ONNX
    print("\n[2] onnxruntime (ONNX graph)")
    try:
        import onnxruntime as ort

        sess = ort.InferenceSession(args.onnx, providers=["CPUExecutionProvider"])
        ort_out = sess.run(
            ["logits"],
            {"input_ids": ids.numpy().astype(np.int64),
             "attention_mask": mask.numpy().astype(np.int64)},
        )[0]
        ort_t = torch.from_numpy(ort_out)
        stats("logits", ort_t)
        ort_3di = logits_to_3di(ort_t, real_len)
        print("  argmax[:40]:", ort_t.argmax(1)[0, :40].tolist())
        print("  3Di:", ort_3di[:60])
        print(f"  agreement vs PyTorch: {agree(ref_3di, ort_3di):.2f}%")
    except Exception as e:
        print("  onnxruntime check failed:", e)

    # 3. TensorRT engine
    print("\n[3] TensorRT engine")
    engine = ProstT5Engine(args.engine)
    print("  engine input dtypes:", {n: str(engine.dtypes[n]) for n in engine.input_names})
    print("  engine output dtype:", str(engine.dtypes[engine.output_names[0]]))
    trt_logits = engine.infer_logits(ids.cuda(), mask.cuda()).cpu()  # (20,L)
    stats("logits", trt_logits)
    trt_3di = logits_to_3di(trt_logits, real_len)
    print("  argmax[:40]:", trt_logits.argmax(0)[:40].tolist())
    print("  3Di:", trt_3di[:60])
    print(f"  agreement vs PyTorch: {agree(ref_3di, trt_3di):.2f}%")


if __name__ == "__main__":
    main()
