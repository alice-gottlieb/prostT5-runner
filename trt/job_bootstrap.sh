#!/bin/bash
# One-time bootstrap for the ProstT5 TensorRT work, run as a qsub CPU job
# (compute nodes have internet). Reuses the venv already created in
# $SCRATCH/prostt5-trt/venv and downloads the model weights + builds the
# validation and INT8-calibration FASTAs. Nothing is written to $HOME or to the
# distil-prostt5 directory.
#
# Submit:
#   qsub -pe shared 4 -l h_data=8G,h_rt=2:00:00 trt/job_bootstrap.sh
#$ -cwd
#$ -N prostt5_trt_boot
#$ -o /u/scratch/a/aliceg/logs/misc/$JOB_ID.out
#$ -j y
#$ -M $USER@ucla.edu
#$ -m ea
set -euo pipefail

WORK="$SCRATCH/prostt5-trt"
REPO="$HOME/prostT5-runner"
export UV_CACHE_DIR="$WORK/uv-cache"
export HF_HOME="$WORK/hf-cache"

mkdir -p "$WORK/engines" "$SCRATCH/logs/misc"
cp "$REPO/trt/job_bootstrap.sh" "$WORK/job_bootstrap.sh"   # audit copy

source "$WORK/venv/bin/activate"
cd "$REPO/trt"

echo "=== Node ==="; hostname
python -c "import torch, transformers; print('torch', torch.__version__, '| transformers', transformers.__version__)"

echo "=== Download CNN head checkpoint ==="
python - <<PY
from prostt5_3di import download_cnn_checkpoint
download_cnn_checkpoint("$WORK/cnn_head.pt")
print("CNN head ready")
PY

echo "=== Pre-cache ProstT5 weights ==="
python - <<PY
from transformers import T5EncoderModel, T5Tokenizer
T5Tokenizer.from_pretrained("Rostlab/ProstT5", do_lower_case=False, legacy=True)
T5EncoderModel.from_pretrained("Rostlab/ProstT5")
print("ProstT5 cached under $HF_HOME")
PY

echo "=== Build validation + calibration FASTAs from a genome .faa ==="
python - <<PY
import glob, random

def read_faa(path):
    recs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name and buf:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name and buf:
        recs.append((name, "".join(buf)))
    return recs

faas = sorted(glob.glob("$SCRATCH/ncbi_genomes/*.faa"))
assert faas, "no genome .faa files found in \$SCRATCH/ncbi_genomes"
# Pull from a couple of genomes for diversity.
recs = []
for p in faas[:3]:
    recs.extend(read_faa(p))
recs = [(n, s) for n, s in recs if 40 <= len(s) <= 1500]
random.seed(0)
random.shuffle(recs)

def write(path, items):
    with open(path, "w") as fh:
        for n, s in items:
            fh.write(f">{n}\n{s}\n")
    print(f"wrote {len(items)} seqs -> {path}")

write("$WORK/test_seqs.fasta", recs[:30])
write("$WORK/calib_seqs.fasta", recs[30:30+256])
PY

# Export the (architecture-independent) ONNX here on CPU so the per-arch GPU
# build jobs are fully independent and can run in parallel.
echo "=== Exporting ONNX (CPU) ==="
python export_onnx.py \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME" \
    --out "$WORK/prostt5_3di.onnx"

echo "=== DONE bootstrap ==="
ls -la "$WORK"
