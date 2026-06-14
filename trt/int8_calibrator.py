"""
INT8 entropy calibrator for the ProstT5 3Di engine.

TensorRT calibrates a dynamic-shape network at a single fixed shape (the opt
dims of the calibration profile), so every calibration sample is tokenized and
padded/truncated to `opt_len`. We feed one real protein per batch (batch size 1)
and stage the int64 input buffers as persistent torch CUDA tensors, handing their
device pointers to TensorRT.
"""

import os

import numpy as np
import tensorrt as trt
import torch

from prostt5_3di import load_tokenizer, preprocess_sequences


def _read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name and buf:
                    seqs.append("".join(buf))
                name, buf = line[1:], []
            elif line:
                buf.append(line)
    if name and buf:
        seqs.append("".join(buf))
    return seqs


class ProstT5Int8Calibrator(trt.IInt8EntropyCalibrator2):
    def __init__(self, calib_fasta, cache_dir, opt_len, cache_path, max_samples=256):
        super().__init__()
        self.cache_path = cache_path
        self.opt_len = opt_len

        tokenizer = load_tokenizer(cache_dir)
        pad_id = tokenizer.pad_token_id or 0
        seqs = _read_fasta(calib_fasta)[:max_samples]
        print(f"INT8 calibrator: {len(seqs)} sequences, fixed length {opt_len}")

        # Pre-tokenize each sequence to exactly opt_len (pad with pad_id / mask 0,
        # or truncate). Kept on CPU; copied to GPU one at a time in get_batch.
        self.samples = []
        for s in seqs:
            enc = tokenizer(preprocess_sequences([s]), add_special_tokens=True,
                            return_tensors="np")
            ids = enc["input_ids"][0]
            mask = enc["attention_mask"][0]
            ids = self._fit(ids, opt_len, pad_id)
            mask = self._fit(mask, opt_len, 0)
            # int32 to match the engine's int32 inputs.
            self.samples.append((ids.astype(np.int32), mask.astype(np.int32)))

        # Persistent GPU buffers (stable device pointers for TensorRT).
        self.d_ids = torch.zeros(1, opt_len, dtype=torch.int32, device="cuda")
        self.d_mask = torch.zeros(1, opt_len, dtype=torch.int32, device="cuda")
        self.idx = 0

    @staticmethod
    def _fit(arr, length, pad_value):
        if len(arr) >= length:
            return arr[:length]
        out = np.full(length, pad_value, dtype=arr.dtype)
        out[: len(arr)] = arr
        return out

    def get_batch_size(self):
        return 1

    def get_batch(self, names):
        if self.idx >= len(self.samples):
            return None
        ids, mask = self.samples[self.idx]
        self.idx += 1
        self.d_ids.copy_(torch.from_numpy(ids).view(1, -1).cuda())
        self.d_mask.copy_(torch.from_numpy(mask).view(1, -1).cuda())
        ptr = {"input_ids": int(self.d_ids.data_ptr()),
               "attention_mask": int(self.d_mask.data_ptr())}
        return [ptr[n] for n in names]

    def read_calibration_cache(self):
        if os.path.exists(self.cache_path):
            with open(self.cache_path, "rb") as f:
                return f.read()
        return None

    def write_calibration_cache(self, cache):
        with open(self.cache_path, "wb") as f:
            f.write(cache)
