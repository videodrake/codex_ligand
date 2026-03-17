# -*- coding: utf-8 -*-
"""Centralized logging configuration for the PPI docking pipeline.

All format strings, level constants, and setup functions live here.
Other modules import ``get_logger`` / ``get_worker_logger`` instead of
creating loggers with ad-hoc names.

Logger hierarchy::

    pipeline            Main process (INFO→stdout, DEBUG→file)
    pipeline.worker     Worker processes (DEBUG→workers.log)
"""

import logging
import os
import sys

# ---- Format constants ---- #
LOG_FMT_MAIN = "%(asctime)s [%(levelname)s] %(message)s"
LOG_FMT_WORKER = (
    "%(asctime)s [%(name)s] (%(levelname)s) [PID:%(process)d] %(message)s"
)
DATE_FMT = "%Y-%m-%d %H:%M:%S"

# ---- Sentinel for idempotent worker setup ---- #
_worker_initialized_path: str = ""


# ---- Main process ---- #

def setup_main_logging(log_dir: str) -> logging.Logger:
    """Configure and return the ``pipeline`` logger.

    * StreamHandler on *stdout* at INFO level
    * FileHandler ``pipeline.log`` at DEBUG level

    Called once from ``PipelineManager.__init__``.
    """
    logger = logging.getLogger("pipeline")
    logger.setLevel(logging.DEBUG)

    # Remove stale handlers (re-entry safety)
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        try:
            handler.close()
        except Exception:
            pass

    main_fmt = logging.Formatter(LOG_FMT_MAIN, datefmt=DATE_FMT)

    # stdout — INFO and above
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(main_fmt)
    logger.addHandler(sh)

    # file — DEBUG and above
    fh = logging.FileHandler(
        os.path.join(log_dir, "pipeline.log"), encoding="utf-8"
    )
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(main_fmt)
    logger.addHandler(fh)

    logger.propagate = False
    return logger


# ---- Worker processes ---- #

def setup_worker_logging() -> None:
    """Configure the ``pipeline.worker`` logger for the current process.

    Reads ``PYROSETTA_WORKER_LOG_DIR`` to locate ``workers.log``.
    Idempotent — only attaches handlers once per (process, path) pair.
    """
    global _worker_initialized_path

    log_dir = os.environ.get("PYROSETTA_WORKER_LOG_DIR", "").strip()
    if not log_dir:
        log_dir = os.getcwd()
    os.makedirs(log_dir, exist_ok=True)

    desired_path = os.path.abspath(os.path.join(log_dir, "workers.log"))

    if _worker_initialized_path == desired_path:
        return

    logger = logging.getLogger("pipeline.worker")
    logger.setLevel(logging.DEBUG)

    # Remove handlers pointing elsewhere (directory change across seeds)
    for handler in list(logger.handlers):
        base = getattr(handler, "baseFilename", "")
        if base and os.path.abspath(base) != desired_path:
            logger.removeHandler(handler)
            try:
                handler.close()
            except Exception:
                pass

    has_handler = any(
        os.path.abspath(getattr(h, "baseFilename", "")) == desired_path
        for h in logger.handlers
        if getattr(h, "baseFilename", "")
    )
    if not has_handler:
        fmt = logging.Formatter(LOG_FMT_WORKER, datefmt=DATE_FMT)
        fh = logging.FileHandler(desired_path, encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    logger.propagate = False
    _worker_initialized_path = desired_path


# ---- Convenience accessors ---- #

def get_logger() -> logging.Logger:
    """Return the main ``pipeline`` logger."""
    return logging.getLogger("pipeline")


def get_worker_logger() -> logging.Logger:
    """Return the ``pipeline.worker`` logger."""
    return logging.getLogger("pipeline.worker")
