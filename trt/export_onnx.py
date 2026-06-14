"""
Export the fused ProstT5 encoder + CNN 3Di head to ONNX.

The exported graph is architecture-independent, so this runs once. TensorRT
engines are then built per-GPU from this single ONNX file.

Inputs  : input_ids (B, L) int64, attention_mask (B, L) int64
Output  : logits (B, 20, L) float32  (3Di-class scores per residue)

Batch and length are dynamic axes; downstream we build batch-1 engines with a
dynamic length profile.

Usage:
    uv run python trt/export_onnx.py \
        --cnn-ckpt $SCRATCH/prostt5-trt/cnn_head.pt \
        --cache-dir $SCRATCH/prostt5-trt/hf-cache \
        --out $SCRATCH/prostt5-trt/prostt5_3di.onnx
"""

import argparse

import torch

from prostt5_3di import load_model, load_tokenizer, tokenize_one


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnn-ckpt", required=True, help="Path to CNN head checkpoint")
    ap.add_argument("--cache-dir", default=None, help="HuggingFace cache dir")
    ap.add_argument("--out", required=True, help="Output .onnx path")
    ap.add_argument("--opset", type=int, default=17)
    args = ap.parse_args()

    print("Loading FP32 reference model (CPU)...")
    model = load_model(args.cnn_ckpt, cache_dir=args.cache_dir)
    tokenizer = load_tokenizer(args.cache_dir)

    # A real short sequence gives the exporter representative shapes.
    sample = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF"
    input_ids, attention_mask, _ = tokenize_one(tokenizer, sample, device="cpu")
    print(f"Tracing with sample shape {tuple(input_ids.shape)}")

    dynamic_axes = {
        "input_ids": {0: "batch", 1: "length"},
        "attention_mask": {0: "batch", 1: "length"},
        "logits": {0: "batch", 2: "length"},
    }

    with torch.no_grad():
        torch.onnx.export(
            model,
            (input_ids, attention_mask),
            args.out,
            input_names=["input_ids", "attention_mask"],
            output_names=["logits"],
            dynamic_axes=dynamic_axes,
            opset_version=args.opset,
            do_constant_folding=True,
        )
    print(f"Wrote ONNX to {args.out}")

    # Sanity-check the exported graph loads and has the expected I/O.
    import onnx

    m = onnx.load(args.out)
    onnx.checker.check_model(m)
    print("ONNX check passed. Inputs/outputs:")
    for t in list(m.graph.input) + list(m.graph.output):
        print("  ", t.name)


if __name__ == "__main__":
    main()
