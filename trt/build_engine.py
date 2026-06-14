"""
Build a TensorRT engine from the ProstT5 3Di ONNX file.

TensorRT engines are specific to the GPU architecture and TensorRT version they
are built on, so this script is run once per GPU type (RTX2080Ti / A100 / H200),
each time on a node of that type.

We build a batch-1 engine with a dynamic sequence-length optimization profile
and FP16 enabled (the big win for the 1.2B-param encoder).

Usage (inside a job on the target GPU):
    uv run python trt/build_engine.py \
        --onnx $SCRATCH/prostt5-trt/prostt5_3di.onnx \
        --engine $SCRATCH/prostt5-trt/engines/prostt5_3di.RTX2080Ti.fp16.engine \
        --min-len 16 --opt-len 384 --max-len 2048 \
        --fp16
"""

import argparse
import os

import tensorrt as trt


def build(onnx_path, engine_path, min_len, opt_len, max_len, fp16, workspace_gb,
          int8=False, calib_fasta=None, cache_dir=None, calib_cache=None):
    logger = trt.Logger(trt.Logger.INFO)
    builder = trt.Builder(logger)
    # TensorRT 10: networks are strongly-typed/explicit-batch by default (flag 0).
    network = builder.create_network(0)
    parser = trt.OnnxParser(network, logger)

    print(f"Parsing ONNX: {onnx_path}")
    # parse_from_file (not parse(bytes)) so TensorRT resolves the external-data
    # weight files that sit next to the graph (the model is >2GB).
    if not parser.parse_from_file(onnx_path):
        for i in range(parser.num_errors):
            print("  ONNX parse error:", parser.get_error(i))
        raise RuntimeError("Failed to parse ONNX")

    config = builder.create_builder_config()
    config.set_memory_pool_limit(
        trt.MemoryPoolType.WORKSPACE, int(workspace_gb * (1 << 30))
    )
    if fp16:
        if not builder.platform_has_fast_fp16:
            print("WARNING: platform reports no fast FP16; building FP16 anyway")
        config.set_flag(trt.BuilderFlag.FP16)
        print("FP16 enabled")

    # Batch-1, dynamic length optimization profile.
    profile = builder.create_optimization_profile()
    for name in ("input_ids", "attention_mask"):
        profile.set_shape(name, (1, min_len), (1, opt_len), (1, max_len))
    config.add_optimization_profile(profile)

    if int8:
        from int8_calibrator import ProstT5Int8Calibrator

        config.set_flag(trt.BuilderFlag.INT8)
        # FP16 stays on too, so TensorRT can keep precision-sensitive layers in
        # FP16 and only drop robust layers to INT8.
        # Calibration runs at a single fixed shape (the calib profile's opt dims).
        calib_profile = builder.create_optimization_profile()
        for name in ("input_ids", "attention_mask"):
            calib_profile.set_shape(name, (1, opt_len), (1, opt_len), (1, opt_len))
        config.set_calibration_profile(calib_profile)
        config.int8_calibrator = ProstT5Int8Calibrator(
            calib_fasta, cache_dir, opt_len, calib_cache
        )
        print("INT8 enabled (with FP16 fallback) using calibrator")

    print(
        f"Building engine (min/opt/max length = {min_len}/{opt_len}/{max_len})... "
        "this can take several minutes."
    )
    serialized = builder.build_serialized_network(network, config)
    if serialized is None:
        raise RuntimeError("Engine build failed (build_serialized_network returned None)")

    os.makedirs(os.path.dirname(engine_path), exist_ok=True)
    with open(engine_path, "wb") as f:
        f.write(serialized)
    size_mb = os.path.getsize(engine_path) / (1 << 20)
    print(f"Wrote engine: {engine_path} ({size_mb:.1f} MB)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--onnx", required=True)
    ap.add_argument("--engine", required=True)
    ap.add_argument("--min-len", type=int, default=16)
    ap.add_argument("--opt-len", type=int, default=384)
    ap.add_argument("--max-len", type=int, default=2048)
    ap.add_argument("--fp16", action="store_true")
    ap.add_argument("--int8", action="store_true",
                    help="Build an INT8 engine (requires --calib-fasta)")
    ap.add_argument("--calib-fasta", default=None, help="FASTA for INT8 calibration")
    ap.add_argument("--cache-dir", default=None, help="HF cache dir (for tokenizer)")
    ap.add_argument("--calib-cache", default=None,
                    help="Path to read/write the INT8 calibration table")
    ap.add_argument("--workspace-gb", type=float, default=8.0)
    args = ap.parse_args()

    if args.int8 and not args.calib_fasta:
        ap.error("--int8 requires --calib-fasta")

    print("TensorRT version:", trt.__version__)
    build(
        args.onnx, args.engine, args.min_len, args.opt_len, args.max_len,
        args.fp16, args.workspace_gb,
        int8=args.int8, calib_fasta=args.calib_fasta, cache_dir=args.cache_dir,
        calib_cache=args.calib_cache,
    )


if __name__ == "__main__":
    main()
