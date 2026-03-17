# -*- coding: utf-8 -*-
import logging
import os
import sys
from typing import TYPE_CHECKING, Dict, Optional, Set

import pyrosetta
import pyrosetta.rosetta.core.import_pose
import pyrosetta.rosetta.std

if TYPE_CHECKING:
    from pyrosetta.rosetta.core.pose import Pose

TOP_LEVEL_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_worker_logging_initialized = False
_worker_logging_path = ""
_master_seed = None


def _worker_log_dir() -> str:
    log_dir = os.environ.get("PYROSETTA_WORKER_LOG_DIR", "").strip()
    if not log_dir:
        log_dir = TOP_LEVEL_DIR
    os.makedirs(log_dir, exist_ok=True)
    return log_dir


def _worker_log_path() -> str:
    return os.path.join(_worker_log_dir(), f"worker_{os.getpid()}.log")


def set_master_seed(seed: Optional[int]) -> None:
    """Set the master random seed for reproducibility across workers."""
    global _master_seed
    _master_seed = seed


def get_master_seed() -> Optional[int]:
    """Return the current master seed (None if not set)."""
    return _master_seed

def setup_worker_logging() -> None:
    """
    Configure logging handlers for worker processes.
    Idempotent: only attaches handlers once per process.
    Must be called at the start of every worker task function.
    """
    global _worker_logging_initialized
    global _worker_logging_path

    logger = logging.getLogger('python_internal')
    logger.setLevel(logging.DEBUG)

    desired_path = os.path.abspath(_worker_log_path())
    if _worker_logging_initialized and _worker_logging_path == desired_path:
        return
    for handler in list(logger.handlers):
        base_filename = getattr(handler, "baseFilename", "")
        if base_filename and os.path.abspath(base_filename) != desired_path:
            logger.removeHandler(handler)
            try:
                handler.close()
            except Exception:
                pass

    has_file_handler = any(
        os.path.abspath(getattr(handler, "baseFilename", "")) == desired_path
        for handler in logger.handlers
        if getattr(handler, "baseFilename", "")
    )
    if not has_file_handler:
        fmt = logging.Formatter(
            '%(asctime)s [%(name)s] (%(levelname)s) [PID:%(process)d] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        fh = logging.FileHandler(desired_path, encoding='utf-8')
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    logger.propagate = False

    _worker_logging_initialized = True
    _worker_logging_path = desired_path


internal_logger = logging.getLogger('python_internal')


def parse_residue_ranges(range_str: Optional[str]) -> Set[int]:
    """
    Parse residue range string into a set of integers.
    Example: '700-807,840-854,863-875' -> {700,701,...,807,840,...,854,863,...,875}
    Returns empty set if input is empty or None.
    """
    residues = set()
    if not range_str or not isinstance(range_str, str) or range_str.strip() == '':
        return residues
    for part in range_str.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            tokens = part.split('-', 1)
            try:
                start, end = int(tokens[0].strip()), int(tokens[1].strip())
                for r in range(start, end + 1):
                    residues.add(r)
            except ValueError:
                internal_logger.warning(f"Invalid residue range: '{part}', skipping.")
        else:
            try:
                residues.add(int(part))
            except ValueError:
                internal_logger.warning(f"Invalid residue number: '{part}', skipping.")
    return residues


def parse_residue_weights(weight_str: Optional[str]) -> Dict[int, float]:
    """Parse residue weight string into {residue_number: weight} dict.

    Format: "961:2.0,962-964:1.5,968" (no weight = 1.0)
    Returns empty dict if input is empty or None.
    """
    weights: Dict[int, float] = {}
    if not weight_str or not isinstance(weight_str, str) or weight_str.strip() == '':
        return weights
    for part in weight_str.split(','):
        part = part.strip()
        if not part:
            continue
        # Split off optional weight suffix
        if ':' in part:
            res_part, w_part = part.rsplit(':', 1)
            try:
                w = float(w_part.strip())
            except ValueError:
                internal_logger.warning(f"Invalid weight in '{part}', using 1.0")
                w = 1.0
        else:
            res_part = part
            w = 1.0
        # Parse residue range
        res_part = res_part.strip()
        if '-' in res_part:
            tokens = res_part.split('-', 1)
            try:
                start, end = int(tokens[0].strip()), int(tokens[1].strip())
                for r in range(start, end + 1):
                    weights[r] = w
            except ValueError:
                internal_logger.warning(f"Invalid residue range in '{part}', skipping.")
        else:
            try:
                weights[int(res_part)] = w
            except ValueError:
                internal_logger.warning(f"Invalid residue number in '{part}', skipping.")
    return weights


def init_rosetta() -> None:
    """
    Initialize PyRosetta safely.
    Ensures idempotency and proper flag settings for HPC environments.
    """
    setup_worker_logging()

    if _master_seed is not None:
        seed = _master_seed + os.getpid()
    else:
        seed = os.getpid()
    cmd_flags = (
        "-mute all "
        f"-jran {seed} "
        "-multithreading:total_threads 1 "
    )

    # Temporarily redirect stdout/stderr to suppress PyRosetta banner
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')

    try:
        pyrosetta.init(cmd_flags)
    except RuntimeError as e:
        if "already initialized" in str(e) or "Glut" in str(e):
            pass
        else:
            internal_logger.critical(f"PyRosetta Init RuntimeError: {e}")
            raise
    except Exception as e:
        internal_logger.critical(f"Unexpected PyRosetta Init Error: {e}")
        raise
    finally:
        # Restore stdout/stderr (only close if redirected)
        if sys.stdout is not original_stdout:
            sys.stdout.close()
        if sys.stderr is not original_stderr:
            sys.stderr.close()
        sys.stdout = original_stdout
        sys.stderr = original_stderr


def string_to_pose(pdb_string: str) -> Optional["Pose"]:
    """
    Converts PDB string to Pose object with validation.
    Returns None if input is invalid.
    """
    if not pdb_string or not isinstance(pdb_string, str):
        internal_logger.error("string_to_pose received invalid input.")
        return None

    init_rosetta()

    try:
        pose = pyrosetta.Pose()
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdb_string)

        if pose.total_residue() == 0:
            internal_logger.warning("Converted Pose has 0 residues.")

        return pose
    except Exception as e:
        internal_logger.error(f"string_to_pose failed: {e}")
        return None


def count_excluded_contacts(pose: "Pose", excluded_residues_A: Set[int], contact_dist: float = 10.0) -> int:
    """
    Count how many excluded Chain A residues are contacted by Chain B within contact_dist.
    Uses CA-CA distance. Residues not present in the PDB are silently skipped.
    Returns integer count.
    """
    if not excluded_residues_A or pose.num_chains() < 2:
        return 0

    conf = pose.conformation()
    chain_a_end = conf.chain_end(1)
    chain_b_start = conf.chain_begin(2)
    chain_b_end = conf.chain_end(2)

    b_cas = []
    for i in range(chain_b_start, chain_b_end + 1):
        if pose.residue(i).has("CA"):
            b_cas.append(pose.residue(i).xyz("CA"))

    if not b_cas:
        return 0

    dist_sq = contact_dist * contact_dist
    count = 0

    for i in range(1, chain_a_end + 1):
        pdb_num = pose.pdb_info().number(i)
        if pdb_num not in excluded_residues_A:
            continue
        if not pose.residue(i).has("CA"):
            continue
        a_ca = pose.residue(i).xyz("CA")
        for b_ca in b_cas:
            dx = a_ca.x - b_ca.x
            dy = a_ca.y - b_ca.y
            dz = a_ca.z - b_ca.z
            if dx * dx + dy * dy + dz * dz <= dist_sq:
                count += 1
                break

    return count


def pose_to_string(pose: Optional["Pose"]) -> str:
    """
    Converts Pose object to PDB string safely.
    """
    if pose is None:
        return ""

    try:
        s = pyrosetta.rosetta.std.ostringstream()
        pose.dump_pdb(s)
        return s.str()
    except Exception as e:
        internal_logger.error(f"pose_to_string failed: {e}")
        return ""
