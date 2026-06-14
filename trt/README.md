# ProstT5 → 3Di, optimized with TensorRT

TensorRT-optimized version of the ProstT5 amino-acid → 3Di predictor (the same
encoder-only + CNN-head model foldseek uses for `createdb --prostt5-model`),
built per GPU architecture on Hoffman (RTX2080Ti, A100, H200).

The optimized engine produces the foldseek **3Di** structural alphabet directly
and can feed `foldseek search`, at a multiple of PyTorch throughput.

## What it is

- **Model**: `Rostlab/ProstT5` T5 **encoder** → per-residue 1024-d embeddings →
  the official 2-layer **CNN head** → 20-class 3Di logits → argmax. Fused into a
  single graph ([prostt5_3di.py](prostt5_3di.py)).
- **Export**: one ONNX, weights as external data (>2 GB)
  ([export_onnx.py](export_onnx.py)). `--fp16` exports an FP16 graph in which
  `T5LayerNorm` keeps its variance in FP32.
- **Engines**: one per GPU arch (engines are arch + TRT-version specific). FP16
  uses a **strongly-typed** network; INT8 uses INT8 matmuls with a BF16 fallback
  ([build_engine.py](build_engine.py)).
- **Runtime**: TensorRT 10 with PyTorch CUDA tensors as the device buffers — no
  pycuda/cuda-python ([trt_runtime.py](trt_runtime.py)).
- **CLI**: [predict_3di.py](predict_3di.py) — AA FASTA → 3Di FASTA.

Batch-1, dynamic sequence length (the foldseek-style per-protein pattern).

## The key gotcha: T5 + FP16

A naive FP16 build collapses the output to a constant ("DDDD…", ~16% agreement).
T5's RMSNorm computes `mean(x²)`; in FP16 the large activations overflow
(`x² > 65504 → inf → rsqrt(inf)=0`), zeroing the embeddings. PyTorch avoids this
because `T5LayerNorm` upcasts the variance to FP32.

Fix: export an **FP16 ONNX** (variance stays FP32) and build a **strongly-typed**
engine that obeys those dtypes — FP16 everywhere except the FP32 variance. This
works on all archs including Turing (no BF16 needed). For INT8, pair it with
**BF16** (FP32 exponent range) as the high-precision fallback (Ampere+).

A plain FP32 engine matches PyTorch FP32 at 100.000%, confirming the ONNX graph
is correct and the issue was purely FP16 numerics.

## Results (per-residue 3Di agreement vs PyTorch FP32; batch-1)

| GPU        | Precision | Agreement | Speedup vs PyTorch FP16 |
| ---------- | --------- | --------- | ----------------------- |
| RTX2080Ti  | FP16      | 99.86%    | 1.91×                   |
| A100-40GB  | FP16      | 99.80%    | 2.04×                   |
| A100-40GB  | INT8+BF16 | 99.13%    | 1.90×                   |
| H200       | FP16      | 99.80%    | 1.66×                   |
| H200       | INT8+BF16 | 99.03%    | 2.39×                   |

Remaining ~0.1–1% are near-decision-boundary residues flipping under reduced
precision. INT8 is built only on A100/H200 (the BF16 fallback needs Ampere+;
Turing's RTX2080Ti gets FP16). Speedups are batch-1 vs PyTorch FP16 on the same
GPU; absolute latency is much lower on H200 (≈7 ms/seq) than RTX2080Ti
(≈39 ms/seq).

## Run it on Hoffman

All artifacts live in `$SCRATCH/prostt5-trt`; everything runs via `qsub`.

```bash
# 0. one-time env (login node): trt/setup_and_prep.sh creates the venv
# 1. bootstrap (downloads + FASTAs + FP32 ONNX) — CPU job
qsub -pe shared 4 -l h_data=8G,h_rt=2:00:00 trt/job_bootstrap.sh

# 2. per-arch FP16 build + validate (first run exports the FP16 ONNX: DO_EXPORT=1)
qsub -pe shared 4 -l gpu,A100,cuda=1,h_data=8G,h_rt=2:00:00 \
     -v GPU_LABEL=A100,MAX_LEN=2048,WORKSPACE_GB=8,DO_EXPORT=1 trt/job_build_validate.sh
qsub -pe shared 4 -l gpu,H200,cuda=1,h_data=8G,h_rt=2:00:00 \
     -v GPU_LABEL=H200,MAX_LEN=2048,WORKSPACE_GB=8 trt/job_build_validate.sh
qsub -pe shared 4 -l gpu,RTX2080Ti,cuda=1,h_data=6G,h_rt=5:00:00 \
     -v GPU_LABEL=RTX2080Ti,MAX_LEN=1024,WORKSPACE_GB=4 trt/job_build_validate.sh

# 3. INT8 (A100/H200 only)
qsub -pe shared 4 -l gpu,A100,cuda=1,h_data=8G,h_rt=4:00:00 \
     -v GPU_LABEL=A100,MAX_LEN=2048,WORKSPACE_GB=8 trt/job_build_int8.sh

# 4. predict (deliverable)
uv run python trt/predict_3di.py \
    --engine $SCRATCH/prostt5-trt/engines/prostt5_3di.A100.fp16.engine \
    --cache-dir $SCRATCH/prostt5-trt/hf-cache \
    --in proteins.fasta --out proteins_3di.fasta
```

## Environment notes (Hoffman / CentOS 7)

- GCC 4.8.5 can't compile modern packages → install wheels only
  (`uv pip install --only-binary :all:`).
- Use `tensorrt-cu12<11` from `https://pypi.nvidia.com` (the `tensorrt`
  meta-package pulls CUDA-13 / TRT-11). Installed: TensorRT 10.7.
- `module load cuda` isn't needed in batch jobs (the torch + tensorrt-cu12
  wheels bundle CUDA); `module` isn't even initialised there.
- Don't `import tensorrt` on the login node (memory ulimit can't mmap its .so);
  it works on the GPU compute nodes.
- The cluster Python's urllib has no CA bundle → downloads use a certifi context.
- int32 (not int64) ids/mask: TensorRT's int64→float cast is unreliable.
