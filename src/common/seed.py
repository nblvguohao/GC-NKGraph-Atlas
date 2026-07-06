"""
GC-NKGraph-Atlas seed utilities.

Ensure reproducible random seeds across all frameworks.
"""

import random
import os

import numpy as np


def set_seed(seed: int = 0, deterministic: bool = True):
    """Set Python, NumPy, and PyTorch random seeds.

    Args:
        seed: Integer seed value.
        deterministic: If True, configure PyTorch for deterministic operations
                      (may reduce performance).
    """
    random.seed(seed)
    np.random.seed(seed)

    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)
        if deterministic:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass

    os.environ["PYTHONHASHSEED"] = str(seed)
