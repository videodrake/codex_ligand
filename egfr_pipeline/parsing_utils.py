# -*- coding: utf-8 -*-
"""Shared safe-parsing utilities used across all pipeline phases.

Consolidates the 20+ independent _safe_float/_safe_int implementations
into a single module with two calling conventions:

    # Pattern A: Returns None on failure (phase1/phase2/phase3 style)
    safe_float("3.14")        -> 3.14
    safe_float("bad")         -> None

    # Pattern B: Returns a default on failure (verdict/report/phase4 style)
    safe_float("bad", 0.0)    -> 0.0
    safe_int("bad", 0)        -> 0
"""
from typing import Optional, overload, Union


def safe_float(val: object, default: Optional[float] = None) -> Optional[float]:
    """Convert *val* to float, returning *default* on failure.

    Handles None, empty strings, and non-numeric values gracefully.
    """
    if val is None or (isinstance(val, str) and val.strip() == ""):
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val: object, default: Optional[int] = None) -> Optional[int]:
    """Convert *val* to int, returning *default* on failure.

    Tries int() first; falls back to int(float()) for strings like "3.0".
    """
    if val is None or (isinstance(val, str) and val.strip() == ""):
        return default
    try:
        return int(val)
    except (ValueError, TypeError):
        try:
            return int(float(val))
        except (ValueError, TypeError):
            return default
