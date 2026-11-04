def get_accelerator() -> str:
    """
    Return the best available PyTorch Lightning accelerator string.

    Probes PyTorch at runtime in priority order: CUDA GPU → Apple MPS → CPU.

    Returns
    -------
    str
        "gpu" if a CUDA device is available, "mps" if Apple Silicon GPU is
        available, or "cpu" as the fallback.
    """
    try:
        import torch
        if torch.cuda.is_available():
            return "gpu"
        if torch.backends.mps.is_available():
            return "mps"
    except ImportError:
        pass
    return "cpu"
