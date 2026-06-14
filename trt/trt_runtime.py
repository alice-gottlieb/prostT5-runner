"""
Minimal TensorRT runtime for the ProstT5 3Di engine (batch-1, dynamic length).

Device memory and the CUDA stream come from PyTorch CUDA tensors, so the only
extra dependency beyond torch is the `tensorrt` Python package (no pycuda /
cuda-python needed).
"""

import numpy as np
import tensorrt as trt
import torch

# Map TensorRT dtypes to the torch dtypes we allocate I/O buffers with.
_TRT_TO_TORCH = {
    trt.DataType.FLOAT: torch.float32,
    trt.DataType.HALF: torch.float16,
    trt.DataType.INT32: torch.int32,
    trt.DataType.INT64: torch.int64,
    trt.DataType.BOOL: torch.bool,
}


class ProstT5Engine:
    def __init__(self, engine_path: str, device: str = "cuda"):
        self.device = device
        self.logger = trt.Logger(trt.Logger.WARNING)
        runtime = trt.Runtime(self.logger)
        with open(engine_path, "rb") as f:
            self.engine = runtime.deserialize_cuda_engine(f.read())
        if self.engine is None:
            raise RuntimeError(f"Failed to deserialize engine: {engine_path}")
        self.context = self.engine.create_execution_context()

        # Classify tensors as input/output and remember their dtypes.
        self.input_names, self.output_names = [], []
        self.dtypes = {}
        for i in range(self.engine.num_io_tensors):
            name = self.engine.get_tensor_name(i)
            self.dtypes[name] = _TRT_TO_TORCH[self.engine.get_tensor_dtype(name)]
            if self.engine.get_tensor_mode(name) == trt.TensorIOMode.INPUT:
                self.input_names.append(name)
            else:
                self.output_names.append(name)

    def infer_logits(self, input_ids: torch.Tensor, attention_mask: torch.Tensor):
        """Run the engine on one (1, L) example, returning (20, L) logits on GPU."""
        ids = input_ids.to(self.device, self.dtypes["input_ids"]).contiguous()
        mask = attention_mask.to(self.device, self.dtypes["attention_mask"]).contiguous()

        self.context.set_input_shape("input_ids", tuple(ids.shape))
        self.context.set_input_shape("attention_mask", tuple(mask.shape))

        out_shape = tuple(self.context.get_tensor_shape("logits"))
        out = torch.empty(out_shape, dtype=self.dtypes["logits"], device=self.device)

        self.context.set_tensor_address("input_ids", ids.data_ptr())
        self.context.set_tensor_address("attention_mask", mask.data_ptr())
        self.context.set_tensor_address("logits", out.data_ptr())

        stream = torch.cuda.current_stream(self.device)
        ok = self.context.execute_async_v3(stream.cuda_stream)
        if not ok:
            raise RuntimeError("execute_async_v3 returned False")
        stream.synchronize()
        return out[0]  # (20, L)
