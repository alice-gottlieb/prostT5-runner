"""
Shared ProstT5 -> 3Di model definition and helpers.

This reproduces the official encoder-only 3Di predictor from
https://github.com/mheinzinger/ProstT5 (scripts/predict_3Di_encoderOnly.py):
the ProstT5 T5 encoder produces per-residue embeddings, and a small 2-layer
CNN head maps each residue's embedding to one of 20 3Di states.

The same module is used three ways:
  * as the FP32 PyTorch reference (ground truth for validation),
  * as the thing we export to ONNX, and
  * (indirectly) as the network TensorRT compiles.

Keeping a single definition guarantees the reference and the engine compute the
same function.
"""

import os
import re

import torch
import torch.nn as nn
from transformers import T5EncoderModel, T5Tokenizer

# ProstT5 model id. The fp16 repo holds the same weights stored in half; we load
# the fp32 master so the PyTorch reference is the true ground truth and let
# TensorRT do the FP16 conversion itself.
PROSTT5_MODEL_ID = "Rostlab/ProstT5"

# Official CNN head checkpoint (maps 1024-d embeddings -> 20 3Di classes).
CNN_CKPT_URL = "https://github.com/mheinzinger/ProstT5/raw/main/cnn_chkpnt/model.pt"

# Class index -> 3Di letter, exactly as in the reference script.
SS_MAPPING = {
    0: "A", 1: "C", 2: "D", 3: "E", 4: "F", 5: "G", 6: "H", 7: "I", 8: "K",
    9: "L", 10: "M", 11: "N", 12: "P", 13: "Q", 14: "R", 15: "S", 16: "T",
    17: "V", 18: "W", 19: "Y",
}
SS_LETTERS = "".join(SS_MAPPING[i] for i in range(len(SS_MAPPING)))


class CNNHead(nn.Module):
    """Two-layer CNN that turns per-residue embeddings into 3Di-class logits.

    Input  x: (batch, length, 1024)
    Output  : (batch, 20, length)   -- logits over the 20 3Di classes
    """

    def __init__(self):
        super().__init__()
        # Named `classifier` to match the published checkpoint's state_dict keys.
        self.classifier = nn.Sequential(
            nn.Conv2d(1024, 32, kernel_size=(7, 1), padding=(3, 0)),
            nn.ReLU(),
            nn.Dropout(0.0),
            nn.Conv2d(32, 20, kernel_size=(7, 1), padding=(3, 0)),
        )

    def forward(self, x):
        # (B, L, 1024) -> (B, 1024, L, 1) so width-1 2D convs act along length.
        x = x.permute(0, 2, 1).unsqueeze(dim=-1)
        logits = self.classifier(x)        # (B, 20, L, 1)
        return logits.squeeze(dim=-1)      # (B, 20, L)


class ProstT5For3Di(nn.Module):
    """Encoder + CNN head, fused so the whole AA->3Di-logits path is one graph."""

    def __init__(self, encoder: T5EncoderModel, head: CNNHead):
        super().__init__()
        self.encoder = encoder
        self.head = head

    def forward(self, input_ids, attention_mask):
        emb = self.encoder(
            input_ids=input_ids, attention_mask=attention_mask
        ).last_hidden_state                       # (B, L, 1024)
        # Zero out padding/special-token positions before the CNN, matching the
        # reference (it multiplies embeddings by the attention mask).
        emb = emb * attention_mask.unsqueeze(-1).to(emb.dtype)
        return self.head(emb)                      # (B, 20, L)


def download_cnn_checkpoint(dest_path: str) -> str:
    """Download the CNN head checkpoint to dest_path if not already present.

    Must run where there is internet (the Hoffman login node, not compute nodes).
    """
    if os.path.exists(dest_path):
        return dest_path
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    print(f"Downloading CNN head checkpoint -> {dest_path}")
    # The cluster Python's urllib has no system CA bundle, so point it at
    # certifi's (certifi ships with httpx in this env). urllib follows GitHub's
    # raw -> githubusercontent redirect automatically.
    import ssl
    import urllib.request

    import certifi

    ctx = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(CNN_CKPT_URL, context=ctx, timeout=120) as r:
        data = r.read()
    with open(dest_path, "wb") as f:
        f.write(data)
    return dest_path


def load_cnn_head(ckpt_path: str) -> CNNHead:
    head = CNNHead()
    state = torch.load(ckpt_path, map_location="cpu")
    # The published file wraps weights under a "state_dict" key.
    state_dict = state.get("state_dict", state)
    head.load_state_dict(state_dict)
    head.eval()
    return head


def load_tokenizer(cache_dir: str | None = None) -> T5Tokenizer:
    return T5Tokenizer.from_pretrained(
        PROSTT5_MODEL_ID, do_lower_case=False, cache_dir=cache_dir, legacy=True
    )


def load_model(cnn_ckpt_path: str, cache_dir: str | None = None) -> ProstT5For3Di:
    """Load the fused FP32 encoder+head module in eval mode (on CPU)."""
    encoder = T5EncoderModel.from_pretrained(
        PROSTT5_MODEL_ID, cache_dir=cache_dir, torch_dtype=torch.float32
    )
    encoder.eval()
    head = load_cnn_head(cnn_ckpt_path)
    model = ProstT5For3Di(encoder, head)
    model.eval()
    return model


def preprocess_sequences(sequences: list[str]) -> list[str]:
    """Apply the ProstT5 AA->3Di input convention.

    Upper-case, replace rare/ambiguous residues U/Z/O/B with X, space out the
    residues, and prepend the <AA2fold> direction token.
    """
    out = []
    for seq in sequences:
        seq = re.sub(r"[UZOB]", "X", seq.upper())
        out.append("<AA2fold> " + " ".join(list(seq)))
    return out


def logits_to_3di(logits: torch.Tensor, real_len: int) -> str:
    """Argmax a (20, L) or (1, 20, L) logits tensor and map the first
    `real_len` positions to a 3Di string."""
    if logits.dim() == 3:
        logits = logits[0]
    preds = logits.argmax(dim=0).tolist()       # length L
    return "".join(SS_MAPPING[p] for p in preds[:real_len])


def tokenize_one(tokenizer: T5Tokenizer, sequence: str, device: str = "cpu"):
    """Tokenize a single raw AA sequence for batch-1 inference.

    Returns (input_ids, attention_mask, real_len) where real_len is the number
    of amino acids (so trailing special tokens can be dropped from the output).
    """
    proc = preprocess_sequences([sequence])
    enc = tokenizer(
        proc, add_special_tokens=True, padding="longest", return_tensors="pt"
    )
    input_ids = enc.input_ids.to(device)
    attention_mask = enc.attention_mask.to(device)
    real_len = len(re.sub(r"[UZOB]", "X", sequence.upper()))
    return input_ids, attention_mask, real_len
