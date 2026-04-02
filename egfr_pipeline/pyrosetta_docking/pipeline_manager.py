# -*- coding: utf-8 -*-
import base64
import configparser
import csv
import gc
import html as html_mod
import json
import logging
import math
import multiprocessing
import os
import re
import shutil
import sys
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose

from . import scoring, pyrosetta_init, movers
from .pyrosetta_init import TOP_LEVEL_DIR
from ..runtime import resolve_runtime_resources
from .run_metadata import (
    build_decoy_score_rows,
    build_input_validation_report,
    build_input_validation_markdown,
    build_output_root_name,
    build_run_metadata,
)

# ---- Module Constants ---- #
REPORT_WIDTH = 72
BOLTZMANN_KT = 1.0
SELECTION_POOL_MULTIPLIER = 3
COM_THRESHOLD_MIN = 30.0
RMSD_ERROR_VALUE = 999.9
UNSAT_FALLBACK_VALUE = 99

QUALITY_THRESHOLDS = {
    "energy_gap": {"strong": 10, "moderate": 5, "weak": 2},
    "dSASA": {"excellent": 1000, "good": 500, "marginal": 300},
    "sc_value": {"excellent": 0.70, "good": 0.50, "marginal": 0.35},
    "packstat": {"good": 0.55, "acceptable": 0.40},
    "dG_density": {"excellent": -1.5, "good": -1.0, "marginal": -0.5},
    "unsatHb": {"excellent": 5, "acceptable": 10},
    "nres_int": {"large": 20, "adequate": 15, "small": 10},
    "hbonds_int": {"good": 3, "minimum": 1},
    "binding_energy": {"good": -10, "marginal": -3},
    "convergence": {"strong": 0.40, "moderate": 0.20},
    "sampling": {"good": 10000, "marginal": 3000},
    "pipeline": {"good": 90, "marginal": 50},
    "interface_size": {"good": 500, "marginal": 300},
    "shape_fit": {"good": 0.50, "marginal": 0.35},
    "dG_density_check": {"good": -1.0, "marginal": -0.5},
    "nres_check": {"good": 15, "marginal": 8},
}


class PipelineManager:
    def __init__(
        self,
        config_file: str,
        override_input_pdb: Optional[str] = None,
        *,
        allocated_cpus: Optional[int] = None,
        runtime_dir: Optional[str] = None,
    ) -> None:
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file not found: {config_file}")

        self.config_file = os.path.abspath(config_file)
        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        req_cpus = self.config.getint('System', 'n_cpus', fallback=1)
        self.runtime_resources = resolve_runtime_resources(
            requested_cpus=req_cpus,
            allocated_cpus=allocated_cpus,
        )
        self.requested_cpus = req_cpus
        self.n_cpus = self.runtime_resources.effective_cpus
        self.chunksize = max(1, self.n_cpus * 2)  # batch tasks to reduce IPC overhead

        if override_input_pdb:
            self.input_pdb = os.path.abspath(override_input_pdb)
        else:
            self.input_pdb = os.path.abspath(self.config['Path']['input_pdb_name'])

        self.filename = os.path.basename(self.input_pdb)
        self.root_dir = os.path.abspath(
            build_output_root_name(
                self.config,
                self.input_pdb,
                os.path.splitext(self.filename)[0],
            )
        )
        self.runtime_dir = os.path.abspath(runtime_dir or os.path.join(self.root_dir, "runtime"))
        self.log_dir = os.path.join(self.root_dir, "logs")
        os.makedirs(self.log_dir, exist_ok=True)
        os.environ["PYROSETTA_WORKER_LOG_DIR"] = self.log_dir

        self._setup_logging()
        self._validate_config()

        self.dir_cache = os.path.join(TOP_LEVEL_DIR, "relaxed_cache")
        os.makedirs(self.dir_cache, exist_ok=True)

        self.dir_filter = os.path.join(self.root_dir, "filter_passed")
        self.dir_cluster = os.path.join(self.root_dir, "cluster_results")
        self.dir_final = os.path.join(self.root_dir, "final_result")

        os.makedirs(self.root_dir, exist_ok=True)
        for d in [self.dir_filter, self.dir_cluster, self.dir_final]:
            if os.path.exists(d):
                shutil.rmtree(d)
            os.makedirs(d)

        self.stats = {}
        self.stage_times = {}
        self.run_metadata = build_run_metadata(
            config=self.config,
            config_file=self.config_file,
            input_pdb=self.input_pdb,
            root_dir=self.root_dir,
            filename=self.filename,
            requested_cpus=self.requested_cpus,
            used_cpus=self.n_cpus,
            master_seed=self.master_seed,
            dir_filter=self.dir_filter,
            dir_cluster=self.dir_cluster,
            dir_final=self.dir_final,
        )
        self.input_validation_report = build_input_validation_report(
            config=self.config,
            input_pdb=self.input_pdb,
            run_metadata=self.run_metadata,
        )

    # ------------------------------------------------------------------ #
    #  Logging
    # ------------------------------------------------------------------ #
    def _setup_logging(self) -> None:
        from .logging_config import setup_main_logging
        self.logger = setup_main_logging(self.log_dir, root_dir=self.root_dir)

    # ------------------------------------------------------------------ #
    #  Config
    # ------------------------------------------------------------------ #
    def _validate_config(self) -> None:
        try:
            self.total_global = self.config.getint('Docking', 'total_global_models')

            self.filter_percentile = self.config.getfloat('Filter', 'filter_percentile')
            self.min_dSASA = self.config.getfloat('Filter', 'min_dSASA')
            self.contact_distance = self.config.getfloat('Filter', 'contact_distance', fallback=6.0)
            self.min_sc_value = self.config.getfloat('Filter', 'min_sc_value', fallback=0.3)
            self.min_survivors = self.config.getint('Filter', 'min_survivors')

            self.refine_per_struct = self.config.getint('Refinement', 'refine_per_struct', fallback=0)
            self.perturb_t = self.config.getfloat('Refinement', 'perturb_trans', fallback=0.1)
            self.perturb_r = self.config.getfloat('Refinement', 'perturb_rot', fallback=1.0)

            # threshold — auto or float
            _thr_raw = self.config.get('Cluster', 'threshold').strip().lower()
            if _thr_raw == 'auto':
                self.cluster_threshold_auto = True
                self.cluster_threshold = 4.0  # provisional
            else:
                self.cluster_threshold_auto = False
                self.cluster_threshold = float(_thr_raw)

            # cluster_top_n — auto or int
            _topn_raw = self.config.get('Cluster', 'cluster_top_n').strip().lower()
            if _topn_raw == 'auto':
                self.cluster_top_n_auto = True
                self.cluster_top_n = 20  # provisional
            else:
                self.cluster_top_n_auto = False
                self.cluster_top_n = int(_topn_raw)

            # members_per_cluster (no auto)
            self.members_per_cluster = self.config.getint('Cluster', 'members_per_cluster')

            # member_diversity_dist — auto or float
            _div_raw = self.config.get('Cluster', 'member_diversity_dist').strip().lower()
            if _div_raw == 'auto':
                self.member_diversity_dist_auto = True
                self.member_diversity_dist = 2.0  # provisional
            else:
                self.member_diversity_dist_auto = False
                self.member_diversity_dist = float(_div_raw)

            self.site_group_jaccard = self.config.getfloat(
                'Cluster', 'site_group_jaccard', fallback=0.3)

            self.save_top_n = self.config.getint('Output', 'save_top_n')
            self.save_filter_max = self.config.getint('Output', 'save_filter_max', fallback=1000)
            _archive_raw = self.config.getint('Output', 'archive_top_n', fallback=0)
            self.archive_top_n = _archive_raw if _archive_raw > 0 else self.save_top_n * 3

            # --- Constraints (optional) ---
            excl_str = self.config.get('Constraints', 'excluded_residues_A', fallback='')
            key_str = self.config.get('Constraints', 'key_residues_B', fallback='')
            self.excluded_residues_A = pyrosetta_init.parse_residue_ranges(excl_str)
            self.key_residues_B = pyrosetta_init.parse_residue_ranges(key_str)
            self.exclusion_contact_dist = self.config.getfloat(
                'Constraints', 'exclusion_contact_dist', fallback=10.0)
            self.enable_early_rejection = self.config.getboolean(
                'Constraints', 'enable_early_rejection', fallback=False)
            self.max_excluded_contacts = self.config.getint(
                'Constraints', 'max_excluded_contacts', fallback=0)
            self.key_residue_bonus_weight = self.config.getfloat(
                'Constraints', 'key_residue_bonus_weight', fallback=0.0)
            weight_str = self.config.get('Constraints', 'key_residue_weights', fallback='')
            self.key_residue_weights = pyrosetta_init.parse_residue_weights(weight_str)
            # If weights not specified but key_residues_B exists, fill uniform 1.0
            if not self.key_residue_weights and self.key_residues_B:
                self.key_residue_weights = {r: 1.0 for r in self.key_residues_B}

            if self.total_global <= 0: raise ValueError("total_global_models must be > 0")
            if self.n_cpus <= 0: raise ValueError("n_cpus must be > 0")
            if not (0 < self.filter_percentile < 100): raise ValueError("filter_percentile must be between 0 and 100")

            # --- Random seed (optional) ---
            _seed_raw = self.config.get('System', 'random_seed', fallback='').strip()
            if _seed_raw.lower() == 'auto':
                self.master_seed = int(time.time())
            elif _seed_raw and _seed_raw.isdigit():
                self.master_seed = int(_seed_raw)
            else:
                self.master_seed = None

            # --- Mini Refinement (optional) ---
            self.mini_refine_enabled = self.config.getboolean(
                'MiniRefinement', 'enabled', fallback=False)
            self.mini_refine_mode = self.config.get(
                'MiniRefinement', 'mode', fallback='pack_only')
            self.mini_refine_n_rounds = self.config.getint(
                'MiniRefinement', 'n_rounds', fallback=3)

            # --- v2.0 Filter Stages (optional, backward-compatible) ---
            self.filter_v2_enabled = self.config.has_section('FilterStage1')

            if self.filter_v2_enabled:
                # Stage 1: Coarse energy
                self.stage1_total_score_pct = self.config.getfloat(
                    'FilterStage1', 'total_score_percentile', fallback=10.0)
                self.stage1_max_dG = self.config.getfloat(
                    'FilterStage1', 'max_dG_separated', fallback=0.0)
                self.stage1_enabled = self.config.getboolean(
                    'FilterStage1', 'enabled', fallback=True)

                # Stage 2: Interface quality
                self.stage2_min_dSASA = self.config.getfloat(
                    'FilterStage2', 'min_dSASA', fallback=800.0)
                self.stage2_max_dG_density = self.config.getfloat(
                    'FilterStage2', 'max_dG_density', fallback=-1.5)
                self.stage2_min_sc = self.config.getfloat(
                    'FilterStage2', 'min_sc_value', fallback=0.65)
                self.stage2_min_packstat = self.config.getfloat(
                    'FilterStage2', 'min_packstat', fallback=0.65)
                self.stage2_max_unsat = self.config.getint(
                    'FilterStage2', 'max_delta_unsatHbonds', fallback=5)
                self.stage2_min_nres = self.config.getint(
                    'FilterStage2', 'min_nres_int', fallback=15)
                self.stage2_min_hbonds = self.config.getint(
                    'FilterStage2', 'min_hbonds_int', fallback=1)
                self.stage2_expensive = self.config.getboolean(
                    'FilterStage2', 'enable_expensive_metrics', fallback=True)
                self.stage2_min_survivors = self.config.getint(
                    'FilterStage2', 'min_survivors', fallback=50)

            # --- Experimental Data (optional) ---
            self.exp_data_configured = False
            if self.config.has_section('ExperimentalData'):
                _crit = self.config.get('ExperimentalData', 'critical_residues_B', fallback='').strip()
                _nonb = self.config.get('ExperimentalData', 'non_binding_residues_B', fallback='').strip()
                _kba = self.config.get('ExperimentalData', 'known_binding_region_A', fallback='').strip()
                self.exp_critical_B = pyrosetta_init.parse_residue_ranges(_crit)
                self.exp_non_binding_B = pyrosetta_init.parse_residue_ranges(_nonb)
                self.exp_binding_region_A = pyrosetta_init.parse_residue_ranges(_kba)
                self.exp_confidence = self.config.get('ExperimentalData', 'confidence', fallback='medium').strip()
                if self.exp_critical_B or self.exp_non_binding_B:
                    self.exp_data_configured = True

        except (KeyError, ValueError) as e:
            self.logger.critical(f"[Config Error] {e}")
            sys.exit(1)

    # ------------------------------------------------------------------ #
    #  Preflight: chain label / numbering / constraint validation
    # ------------------------------------------------------------------ #
    def _preflight_check(self, relaxed_pdb: str) -> bool:
        """Run preflight checks on the input structure before docking.

        Checks:
          1. Chain labels match expected A/B
          2. Chain count is exactly 2
          3. Excluded residue coverage (hard fail if < 50%)
          4. Log PDB numbering range per chain

        Returns True if all critical checks pass, False otherwise.
        """
        pose = pyrosetta_init.string_to_pose(relaxed_pdb)
        if pose is None:
            self.logger.critical("    [Preflight FAIL] Could not load relaxed pose")
            return False

        n_chains = pose.num_chains()
        self.logger.info(f"    [Preflight] Number of chains: {n_chains}")

        if n_chains < 2:
            self.logger.critical(
                f"    [Preflight FAIL] Expected 2 chains (A+B), found {n_chains}. "
                f"Input PDB must be a 2-chain dimer.")
            return False

        conf = pose.conformation()
        pdb_info = pose.pdb_info()

        # Check chain labels
        chain_labels = {}
        for chain_idx in range(1, n_chains + 1):
            first_res = conf.chain_begin(chain_idx)
            label = pdb_info.chain(first_res)
            chain_labels[chain_idx] = label

        self.logger.info(
            f"    [Preflight] Chain labels: {chain_labels}")

        expected_a = chain_labels.get(1, "?")
        expected_b = chain_labels.get(2, "?")

        if expected_a != "A" or expected_b != "B":
            self.logger.warning(
                f"    [Preflight WARN] Expected chain labels A/B, got "
                f"{expected_a}/{expected_b}. Docking uses hardcoded A_B partners. "
                f"This will cause errors if chain labels don't match.")

        # Log numbering ranges per chain
        for chain_idx in range(1, min(n_chains + 1, 4)):  # max 3 chains
            c_start = conf.chain_begin(chain_idx)
            c_end = conf.chain_end(chain_idx)
            pdb_start = pdb_info.number(c_start)
            pdb_end = pdb_info.number(c_end)
            chain_label = pdb_info.chain(c_start)
            n_res = c_end - c_start + 1
            self.logger.info(
                f"    [Preflight] Chain {chain_label} (internal {chain_idx}): "
                f"residues {pdb_start}-{pdb_end} ({n_res} res)")

        # Excluded residue coverage (hard check)
        if self.excluded_residues_A:
            chain_a_end = conf.chain_end(1)
            chain_a_pdb_nums = set()
            for i in range(1, chain_a_end + 1):
                chain_a_pdb_nums.add(pdb_info.number(i))

            excl_found = self.excluded_residues_A & chain_a_pdb_nums
            excl_missing = self.excluded_residues_A - chain_a_pdb_nums
            coverage_pct = len(excl_found) / len(self.excluded_residues_A) * 100

            self.logger.info(
                f"    [Preflight] Excluded residues: {len(excl_found)}/"
                f"{len(self.excluded_residues_A)} found ({coverage_pct:.0f}%)")

            if excl_missing and len(excl_missing) <= 20:
                self.logger.info(
                    f"    [Preflight] Missing excluded residues: "
                    f"{sorted(excl_missing)}")

            if coverage_pct < 50:
                self.logger.critical(
                    f"    [Preflight FAIL] Excluded residue coverage {coverage_pct:.0f}% < 50%. "
                    f"Numbering mismatch between config and input PDB. "
                    f"Check excluded_residues_A values against actual PDB numbering.")
                return False
            elif coverage_pct < 80:
                self.logger.warning(
                    f"    [Preflight WARN] Excluded residue coverage {coverage_pct:.0f}% < 80%. "
                    f"Some constraint residues may not be effective.")

        # Key residues B coverage
        if self.key_residues_B:
            chain_b_start = conf.chain_begin(2)
            chain_b_end = conf.chain_end(2)
            chain_b_pdb_nums = set()
            for i in range(chain_b_start, chain_b_end + 1):
                chain_b_pdb_nums.add(pdb_info.number(i))

            key_found = self.key_residues_B & chain_b_pdb_nums
            key_pct = len(key_found) / len(self.key_residues_B) * 100
            self.logger.info(
                f"    [Preflight] Key residues B: {len(key_found)}/"
                f"{len(self.key_residues_B)} found ({key_pct:.0f}%)")

        self.logger.info("    [Preflight] All checks passed.")
        return True

    # ------------------------------------------------------------------ #
    #  Filter threshold telemetry
    # ------------------------------------------------------------------ #
    def _log_filter_thresholds(self) -> None:
        """Log all active filter thresholds for reproducibility and debugging."""
        self.logger.info("    [Filter Thresholds] Active configuration:")

        if self.filter_v2_enabled:
            self.logger.info("    [Filter] Version: v2.0 (2-pass)")
            self.logger.info(
                f"    [Stage1] total_score_percentile={self.stage1_total_score_pct}, "
                f"max_dG={self.stage1_max_dG}, enabled={self.stage1_enabled}")
            self.logger.info(
                f"    [Stage2] min_dSASA={self.stage2_min_dSASA}, "
                f"max_dG_density={self.stage2_max_dG_density}, "
                f"min_sc={self.stage2_min_sc}")
            self.logger.info(
                f"    [Stage2] min_packstat={self.stage2_min_packstat}, "
                f"max_unsatHb={self.stage2_max_unsat}, "
                f"min_nres_int={self.stage2_min_nres}, "
                f"min_hbonds={self.stage2_min_hbonds}")
            self.logger.info(
                f"    [Stage2] expensive_metrics={self.stage2_expensive}, "
                f"min_survivors={self.stage2_min_survivors}")
            if self.mini_refine_enabled:
                self.logger.info(
                    f"    [MiniRefine] enabled, mode={self.mini_refine_mode}, "
                    f"n_rounds={self.mini_refine_n_rounds}")
        else:
            self.logger.info("    [Filter] Version: v1.0 (legacy)")
            self.logger.info(
                f"    [v1] filter_percentile={self.filter_percentile}, "
                f"min_dSASA={self.min_dSASA}, "
                f"min_sc={self.min_sc_value}, "
                f"min_survivors={self.min_survivors}")

        if self.exp_data_configured:
            self.logger.info(
                f"    [ExperimentalData] critical_B={len(self.exp_critical_B)} residues, "
                f"non_binding_B={len(self.exp_non_binding_B)} residues, "
                f"confidence={self.exp_confidence}")
        else:
            self.logger.info(
                "    [ExperimentalData] Not configured (sensitivity/specificity analysis skipped)")

    # ------------------------------------------------------------------ #
    #  L_RMSD helper (for clustering in main process)
    # ------------------------------------------------------------------ #
    def _compute_l_rmsd(self, mobile_pose: Any, ref_pose: Any) -> float:
        """Compute Ligand RMSD. Clones mobile internally — inputs are not modified."""
        mob = mobile_pose.clone()
        calpha_superimpose_pose(mob, ref_pose)

        chain_b_start = mob.conformation().chain_begin(2)
        chain_b_end = mob.conformation().chain_end(2)

        rmsd_sum = 0.0
        cnt = 0
        for i in range(chain_b_start, chain_b_end + 1):
            if mob.residue(i).has("CA"):
                d2 = mob.residue(i).xyz("CA").distance_squared(ref_pose.residue(i).xyz("CA"))
                rmsd_sum += d2
                cnt += 1

        if cnt == 0:
            return RMSD_ERROR_VALUE
        return math.sqrt(rmsd_sum / cnt)

    # ------------------------------------------------------------------ #
    #  Chain size resolver (for adaptive clustering)
    # ------------------------------------------------------------------ #
    def _resolve_chain_sizes(self, relaxed_pdb: str) -> None:
        """Compute chain sizes from relaxed pose and resolve 'auto' clustering parameters."""
        pose = pyrosetta_init.string_to_pose(relaxed_pdb)
        if pose is None or pose.num_chains() < 2:
            self.logger.warning(
                "    [Chain Size] Could not determine chain sizes. Using defaults.")
            self.chain_a_size = 0
            self.chain_b_size = 0
            return

        conf = pose.conformation()
        self.chain_a_size = conf.chain_end(1)
        self.chain_b_size = conf.chain_end(2) - conf.chain_begin(2) + 1

        self.logger.info(
            f"    [Chain Size] Chain A (receptor): {self.chain_a_size} residues, "
            f"Chain B (ligand): {self.chain_b_size} residues")

        # --- Resolve auto cluster_top_n ---
        if getattr(self, 'cluster_top_n_auto', False):
            if self.chain_a_size <= 200:
                self.cluster_top_n = 15
            elif self.chain_a_size <= 400:
                self.cluster_top_n = 25
            else:
                self.cluster_top_n = 35
            self.logger.info(
                f"    [Auto] cluster_top_n = {self.cluster_top_n} "
                f"(receptor {self.chain_a_size} residues)")

        # --- Resolve auto threshold ---
        if getattr(self, 'cluster_threshold_auto', False):
            BASE_THR = 4.0
            REF_RES = 50
            FLOOR, CAP = 3.0, 15.0
            raw = BASE_THR * math.sqrt(self.chain_b_size / REF_RES)
            self.cluster_threshold = max(FLOOR, min(CAP, round(raw, 1)))
            self.logger.info(
                f"    [Auto] cluster_threshold = {self.cluster_threshold:.1f} A "
                f"(ligand {self.chain_b_size} res, raw={raw:.1f})")

        # --- Resolve auto member_diversity_dist ---
        if getattr(self, 'member_diversity_dist_auto', False):
            self.member_diversity_dist = round(self.cluster_threshold * 0.5, 1)
            self.logger.info(
                f"    [Auto] member_diversity_dist = {self.member_diversity_dist:.1f} A "
                f"(50% of threshold)")

        # --- Excluded residue coverage ---
        if self.excluded_residues_A:
            chain_a_pdb_nums = set()
            for i in range(1, conf.chain_end(1) + 1):
                chain_a_pdb_nums.add(pose.pdb_info().number(i))
            excl_present = self.excluded_residues_A & chain_a_pdb_nums
            excl_coverage = len(excl_present) / self.chain_a_size * 100 if self.chain_a_size > 0 else 0
            self.stats['excl_zone_present'] = len(excl_present)
            self.stats['excl_zone_total'] = len(self.excluded_residues_A)
            self.stats['excl_zone_coverage'] = excl_coverage
            self.logger.info(
                f"    [Constraints] Excluded zone: {len(excl_present)}/{len(self.excluded_residues_A)} "
                f"residues present in Chain A ({excl_coverage:.1f}% of chain)")

        # Store in stats for Validation Report
        self.stats['chain_a_size'] = self.chain_a_size
        self.stats['chain_b_size'] = self.chain_b_size
        self.stats['cluster_top_n_auto'] = getattr(self, 'cluster_top_n_auto', False)
        self.stats['cluster_threshold_auto'] = getattr(self, 'cluster_threshold_auto', False)
        self.stats['member_diversity_dist_auto'] = getattr(self, 'member_diversity_dist_auto', False)
        self.stats['cluster_top_n_resolved'] = self.cluster_top_n
        self.stats['cluster_threshold_resolved'] = self.cluster_threshold
        self.stats['member_diversity_dist_resolved'] = self.member_diversity_dist

    # ------------------------------------------------------------------ #
    #  Config snapshot
    # ------------------------------------------------------------------ #
    def _save_config_snapshot(self) -> None:
        """Save a copy of the config file with runtime-resolved values appended."""
        snap_path = os.path.join(self.root_dir, "config_snapshot.ini")
        shutil.copy2(self.config_file, snap_path)
        with open(snap_path, 'a', encoding='utf-8') as f:
            f.write("\n\n# === Runtime Resolved Values ===\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# input_pdb (resolved) = {self.input_pdb}\n")
            f.write(f"# chain_a_size = {getattr(self, 'chain_a_size', '?')}\n")
            f.write(f"# chain_b_size = {getattr(self, 'chain_b_size', '?')}\n")
            f.write(f"# cluster_threshold (resolved) = {self.cluster_threshold}\n")
            f.write(f"# cluster_top_n (resolved) = {self.cluster_top_n}\n")
            f.write(f"# member_diversity_dist (resolved) = {self.member_diversity_dist}\n")
            f.write(f"# n_cpus (resolved) = {self.n_cpus}\n")
            f.write(f"# master_seed = {self.master_seed}\n")
            f.write(f"# filter_version = {'v2.0' if self.filter_v2_enabled else 'v1.0'}\n")
            # Append active filter thresholds for full reproducibility
            if self.filter_v2_enabled:
                f.write(f"# [Stage1] total_score_pct={self.stage1_total_score_pct}, "
                        f"max_dG={self.stage1_max_dG}\n")
                f.write(f"# [Stage2] min_dSASA={self.stage2_min_dSASA}, "
                        f"max_dG_density={self.stage2_max_dG_density}, "
                        f"min_sc={self.stage2_min_sc}\n")
                f.write(f"# [Stage2] min_packstat={self.stage2_min_packstat}, "
                        f"max_unsatHb={self.stage2_max_unsat}, "
                        f"min_nres={self.stage2_min_nres}, "
                        f"min_hbonds={self.stage2_min_hbonds}\n")
            else:
                f.write(f"# [v1] filter_pct={self.filter_percentile}, "
                        f"min_dSASA={self.min_dSASA}, "
                        f"min_sc={self.min_sc_value}\n")
        self.logger.info(f"    > [Config] Snapshot saved: {snap_path}")

    def _save_run_metadata(self, status: str) -> None:
        """Persist normalized run metadata for reproducibility and handoff."""
        metadata = dict(self.run_metadata)
        metadata["run_status"] = status
        metadata["input_validation_status"] = self.input_validation_report.get("status", "")
        metadata["stats"] = dict(self.stats)
        metadata["stage_times"] = dict(self.stage_times)
        metadata["resolved_cluster_threshold"] = getattr(self, "cluster_threshold", "")
        metadata["resolved_cluster_top_n"] = getattr(self, "cluster_top_n", "")

        out_path = os.path.join(self.root_dir, "pyrosetta_run_metadata.json")
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, ensure_ascii=False)
        self.stats["pyrosetta_run_metadata_path"] = out_path
        self.logger.info(f"    > [Output] Run metadata: {out_path}")

    def _save_input_validation_report(self) -> None:
        """Persist static Phase 1 input validation details before docking begins."""
        json_path = os.path.join(self.root_dir, "phase1_input_validation_report.json")
        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(self.input_validation_report, handle, indent=2, ensure_ascii=False)
        markdown_path = os.path.join(self.root_dir, "phase1_input_validation_summary.md")
        with open(markdown_path, "w", encoding="utf-8") as handle:
            handle.write(build_input_validation_markdown(self.input_validation_report))

        self.stats["phase1_input_validation_report_path"] = json_path
        self.stats["phase1_input_validation_summary_path"] = markdown_path

        status = self.input_validation_report.get("status", "unknown")
        self.logger.info(
            f"    > [Output] Input validation report ({status}): {json_path}"
        )
        self.logger.info(
            f"    > [Output] Input validation summary ({status}): {markdown_path}"
        )
        for warning in self.input_validation_report.get("warnings", []):
            self.logger.warning(f"    [Input Validation WARN] {warning}")

    def _export_decoy_scores(self, fully_scored: List[Dict]) -> None:
        """Write standardized full-scoring export for downstream review."""
        out_path = os.path.join(self.root_dir, "pyrosetta_decoy_scores.csv")
        rows = build_decoy_score_rows(fully_scored, self.run_metadata)
        pd.DataFrame(rows).to_csv(out_path, index=False)
        self.stats["pyrosetta_decoy_scores_path"] = out_path
        self.logger.info(
            f"    > [Output] Standardized decoy scores: {out_path} ({len(rows)} rows)"
        )

    # ------------------------------------------------------------------ #
    #  Comprehensive scored-model CSV exports (all models, pre-filter)
    # ------------------------------------------------------------------ #

    def _save_scored_all_models_csv(
        self, combined_data: List[Dict], filter_status_map: Optional[Dict[str, str]] = None
    ) -> str:
        """Write scored_all_models.csv with Pass 1 metrics for EVERY scored model.

        Called twice:
        1) Before filtering — filter_status_map is None, column written as empty.
        2) After filtering — filter_status_map maps id -> status string.

        Returns the CSV path.
        """
        csv_path = os.path.join(self.root_dir, "scored_all_models.csv")
        columns = [
            'id', 'dG_separated', 'dSASA', 'dSASA_polar', 'dSASA_hydrophobic',
            'sc_value', 'total_score',
            'center_x', 'center_y', 'center_z',
            'excluded_contacts', 'key_contact_ratio', 'filter_status',
        ]

        rows: List[Dict[str, Any]] = []
        for d in combined_data:
            row: Dict[str, Any] = {
                'id': d.get('id', ''),
                'dG_separated': d.get('dG_separated', ''),
                'dSASA': d.get('dSASA', ''),
                'dSASA_polar': d.get('dSASA_polar', 0.0),
                'dSASA_hydrophobic': d.get('dSASA_hydrophobic', 0.0),
                'sc_value': d.get('sc_value', ''),
                'total_score': d.get('total_score', ''),
                'center_x': d.get('center_x', ''),
                'center_y': d.get('center_y', ''),
                'center_z': d.get('center_z', ''),
                'excluded_contacts': d.get('excluded_contacts', ''),
                'key_contact_ratio': d.get('key_contact_ratio', ''),
                'filter_status': '',
            }
            if filter_status_map is not None:
                row['filter_status'] = filter_status_map.get(
                    str(d.get('id', '')), '')
            rows.append(row)

        df = pd.DataFrame(rows, columns=columns)
        df.to_csv(csv_path, index=False)
        self.logger.info(
            f"    > [Output] All scored models: {csv_path} ({len(rows)} rows)")
        return csv_path

    def _save_scored_stage2_models_csv(self, stage2_all: List[Dict]) -> str:
        """Write scored_stage2_models.csv for models that received expensive metrics.

        Includes both Stage 2 pass and Stage 2 fail models (all that had
        expensive metrics computed).
        """
        csv_path = os.path.join(self.root_dir, "scored_stage2_models.csv")
        columns = [
            'id', 'dG_separated', 'dSASA', 'sc_value', 'total_score',
            'packstat', 'delta_unsatHbonds', 'nres_int', 'hbonds_int',
            'filter_status',
        ]

        rows: List[Dict[str, Any]] = []
        for d in stage2_all:
            rows.append({
                'id': d.get('id', ''),
                'dG_separated': d.get('dG_separated', ''),
                'dSASA': d.get('dSASA', ''),
                'sc_value': d.get('sc_value', ''),
                'total_score': d.get('total_score', ''),
                'packstat': d.get('packstat', ''),
                'delta_unsatHbonds': d.get('delta_unsatHbonds', ''),
                'nres_int': d.get('nres_int', ''),
                'hbonds_int': d.get('hbonds_int', ''),
                'filter_status': d.get('_filter_status', ''),
            })

        df = pd.DataFrame(rows, columns=columns)
        df.to_csv(csv_path, index=False)
        self.logger.info(
            f"    > [Output] Stage 2 scored models: {csv_path} ({len(rows)} rows)")
        return csv_path

    def _update_scored_all_models_filter_status(
        self, filter_status_map: Dict[str, str]
    ) -> None:
        """Update the filter_status column in scored_all_models.csv.

        Reads the previously saved CSV, maps id -> filter_status, and re-writes.
        """
        csv_path = os.path.join(self.root_dir, "scored_all_models.csv")
        if not os.path.exists(csv_path):
            self.logger.warning(
                "scored_all_models.csv not found; skipping filter_status update")
            return

        df = pd.read_csv(csv_path)
        df['id'] = df['id'].astype(str)
        df['filter_status'] = df['id'].map(filter_status_map).fillna('')
        df.to_csv(csv_path, index=False)
        n_pass = (df['filter_status'] == 'pass').sum()
        n_reject = (df['filter_status'] != 'pass').sum()
        self.logger.info(
            f"    > [Output] Updated filter_status in scored_all_models.csv "
            f"(pass={n_pass}, reject={n_reject})")

    # ------------------------------------------------------------------ #
    #  Experimental data correlation
    # ------------------------------------------------------------------ #
    def _compute_exp_correlation(self, ranking_df) -> Dict:
        """Compute correlation between docking results and experimental data."""
        result = {}
        if not getattr(self, 'exp_data_configured', False):
            return result

        # Collect all observed binding residues from Chain B
        observed_B = set()
        per_cluster = defaultdict(set)
        for _, row in ranking_df.iterrows():
            br_b = row.get('Binding_Residues_B', '')
            nums = self._parse_binding_residue_nums(br_b)
            if not nums:
                continue
            observed_B.update(nums)
            parent = row.get('Parent', '')
            cid = parent.split('_')[0] if '_' in str(parent) else parent
            per_cluster[cid].update(nums)

        result['observed_B'] = observed_B
        result['per_cluster'] = dict(per_cluster)

        critical = getattr(self, 'exp_critical_B', set())
        non_binding = getattr(self, 'exp_non_binding_B', set())

        if critical:
            found = critical & observed_B
            result['sensitivity'] = len(found) / len(critical) if critical else 0
            result['critical_found'] = found
            result['critical_missing'] = critical - observed_B
            # Enrichment
            chain_b_size = getattr(self, 'chain_b_size', 0)
            if chain_b_size > 0 and observed_B:
                expected_rate = len(critical) / chain_b_size
                observed_rate = len(found) / len(observed_B) if observed_B else 0
                result['enrichment'] = observed_rate / expected_rate if expected_rate > 0 else 0
            else:
                result['enrichment'] = 0

        if non_binding:
            correctly_excluded = non_binding - observed_B
            result['specificity'] = len(correctly_excluded) / len(non_binding) if non_binding else 0
            result['false_positives'] = non_binding & observed_B

        if critical and non_binding:
            sens = result.get('sensitivity', 0)
            spec = result.get('specificity', 0)
            result['f1'] = (2 * sens * spec / (sens + spec)) if (sens + spec) > 0 else 0

        # Per-cluster sensitivity
        if critical:
            cluster_sens = {}
            for cid, res_set in per_cluster.items():
                found_c = critical & res_set
                cluster_sens[cid] = (len(found_c), len(critical), sorted(found_c))
            result['cluster_sensitivity'] = cluster_sens

        return result

    # ------------------------------------------------------------------ #
    #  Site grouping (Jaccard similarity)
    # ------------------------------------------------------------------ #
    @staticmethod
    def _parse_binding_residue_nums(br_str):
        """Extract residue numbers from 'A:LEU819,A:ALA822' -> {819, 822}."""
        nums = set()
        if not isinstance(br_str, str) or br_str in ('', 'None', 'No_Chain_2', 'Analysis_Failed'):
            return nums
        for token in br_str.split(','):
            m = re.search(r'(\d+)', token.strip())
            if m:
                nums.add(int(m.group(1)))
        return nums

    def _compute_site_groups(self, final_representatives):
        """Group clusters by Jaccard similarity of binding residues (union-find)."""
        threshold = self.site_group_jaccard
        if threshold <= 0:
            return None

        # Collect leader data (Member_ID == 1) per cluster
        leaders = {}
        for item in final_representatives:
            if item.get('Member_ID') == 1:
                cid = f"C{item['Cluster_ID']:02d}"
                br_a = self._parse_binding_residue_nums(
                    item.get('Binding_Residues_A', ''))
                br_b = self._parse_binding_residue_nums(
                    item.get('Binding_Residues_B', ''))
                leaders[cid] = {'res_a': br_a, 'res_b': br_b}

        cids = sorted(leaders.keys())
        if len(cids) < 2:
            return None

        # Union-Find
        parent = {c: c for c in cids}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        # Compute pairwise Jaccard and merge
        for i in range(len(cids)):
            for j in range(i + 1, len(cids)):
                ci, cj = cids[i], cids[j]
                a_i, a_j = leaders[ci]['res_a'], leaders[cj]['res_a']
                b_i, b_j = leaders[ci]['res_b'], leaders[cj]['res_b']

                # Jaccard per chain
                union_a = a_i | a_j
                j_a = len(a_i & a_j) / len(union_a) if union_a else 0
                union_b = b_i | b_j
                j_b = len(b_i & b_j) / len(union_b) if union_b else 0

                j_avg = j_a * 0.7 + j_b * 0.3  # receptor-weighted

                if j_avg >= threshold:
                    union(ci, cj)

        # Build groups
        groups = defaultdict(list)
        for c in cids:
            groups[find(c)].append(c)
        group_list = sorted(groups.values(), key=lambda g: g[0])

        cluster_to_group = {}
        for gid, members in enumerate(group_list, 1):
            for c in members:
                cluster_to_group[c] = gid

        return {
            'groups': group_list,
            'cluster_to_group': cluster_to_group,
            'n_groups': len(group_list),
        }

    # ------------------------------------------------------------------ #
    #  Funnel data computation
    # ------------------------------------------------------------------ #
    def _compute_funnel_data(self, survivor_pool: List[Dict]) -> None:
        """Compute dual L_RMSD for energy funnel plots and P_near metric.

        Saves all_scored_summary.csv and stores P_near in self.stats.
        Must be called after Stage 1 (v2) or dG filter (v1).
        """
        import tempfile

        t_start = time.time()
        count = len(survivor_pool)
        if count == 0:
            self.logger.warning("    > [Funnel Data] No survivors for funnel analysis.")
            return
        self.logger.info(
            f"    > [Funnel Data] Computing dual L_RMSD for {count} survivors...")

        # --- Find best dG structure and write to temp file ---
        best_item = min(survivor_pool, key=lambda d: d.get('dG_separated', 0))
        best_tmp = None
        try:
            fd, best_tmp = tempfile.mkstemp(suffix='.pdb', prefix='best_dG_')
            with os.fdopen(fd, 'w') as f:
                f.write(best_item['pdb_data'])
        except Exception as e:
            self.logger.warning(f"[Funnel Data] Failed to write best PDB: {e}")
            best_tmp = None

        relaxed_ref = getattr(self, 'cache_path', None)

        # --- Multiprocessing dual L_RMSD ---
        lrmsd_tasks = (
            (d['pdb_data'], relaxed_ref, best_tmp)
            for d in survivor_pool
        )

        lrmsd_best_map = {}
        with multiprocessing.Pool(self.n_cpus) as pool:
            lrmsd_iter = pool.imap(
                scoring.run_lrmsd_dual_task,
                lrmsd_tasks, chunksize=self.chunksize)

            for i, res in enumerate(lrmsd_iter):
                d = survivor_pool[i]
                d['L_RMSD'] = res.get('L_RMSD', RMSD_ERROR_VALUE)
                d['L_RMSD_best'] = res.get('L_RMSD_best', RMSD_ERROR_VALUE)
                lrmsd_best_map[d.get('id')] = d['L_RMSD_best']

        # Clean up temp file
        if best_tmp and os.path.exists(best_tmp):
            try:
                os.remove(best_tmp)
            except Exception:
                pass

        self.lrmsd_best_map = lrmsd_best_map

        # --- Save all_scored_summary.csv ---
        csv_path = os.path.join(self.root_dir, "all_scored_summary.csv")
        try:
            rows = []
            for d in survivor_pool:
                rows.append({
                    'id': d.get('id', ''),
                    'dG_separated': d.get('dG_separated', 0),
                    'L_RMSD': d.get('L_RMSD', RMSD_ERROR_VALUE),
                    'L_RMSD_best': d.get('L_RMSD_best', RMSD_ERROR_VALUE),
                    'dSASA': d.get('dSASA', 0),
                    'dSASA_polar': d.get('dSASA_polar', 0.0),
                    'dSASA_hydrophobic': d.get('dSASA_hydrophobic', 0.0),
                    'sc_value': d.get('sc_value', 0.0),
                    'total_score': d.get('total_score', 0.0),
                })
            pd.DataFrame(rows).to_csv(csv_path, index=False)
            self.logger.info(
                f"    > [Funnel Data] Saved {csv_path} ({len(rows)} rows)")
        except Exception as e:
            self.logger.warning(f"[Funnel Data] CSV save error: {e}")

        # --- Compute P_near ---
        try:
            LAMBDA = 2.0
            kT = BOLTZMANN_KT
            energies = np.array([d.get('dG_separated', 0) for d in survivor_pool])
            rmsds = np.array([d.get('L_RMSD_best', RMSD_ERROR_VALUE)
                              for d in survivor_pool])

            # Filter out error values
            valid = rmsds < RMSD_ERROR_VALUE - 1
            if valid.sum() > 0:
                e_valid = energies[valid]
                r_valid = rmsds[valid]
                e_min = e_valid.min()

                boltz = np.exp(-(e_valid - e_min) / kT)
                gaussian = np.exp(-(r_valid ** 2) / (LAMBDA ** 2))

                z = boltz.sum()
                if z > 0:
                    p_near = float((gaussian * boltz).sum() / z)
                else:
                    p_near = 0.0
            else:
                p_near = 0.0

            self.stats['p_near'] = p_near
            self.logger.info(
                f"    > [Funnel Data] P_near = {p_near:.4f} "
                f"(lambda={LAMBDA}A, kT={kT})")
        except Exception as e:
            self.logger.warning(f"[Funnel Data] P_near error: {e}")
            self.stats['p_near'] = None

        elapsed = time.time() - t_start
        self.logger.info(
            f"    > [Funnel Data] Done ({elapsed:.1f} sec)")

    # ------------------------------------------------------------------ #
    #  Visualization helpers
    # ------------------------------------------------------------------ #
    def generate_pymol_script(self, final_csv_path: str) -> None:
        try:
            df = pd.read_csv(final_csv_path)
            output_pml = os.path.join(self.dir_final, "view_results.pml")
            pml_dir = os.path.dirname(output_pml)
            with open(output_pml, "w") as f:
                f.write("# PyMOL Script\nreinitialize\nbg_color white\n")
                for _, row in df.iterrows():
                    rel_path = os.path.relpath(os.path.abspath(row['File_PDB']), pml_dir)
                    f.write(f"load {rel_path}\n")
                f.write("hide all\nshow cartoon\nspectrum b, blue_white_red, minimum=-5, maximum=5\n")
                f.write("python\n")
                f.write("first = cmd.get_object_list()[0]\n")
                f.write("for obj in cmd.get_object_list()[1:]:\n")
                f.write("    cmd.align(f'{obj} and chain A', f'{first} and chain A')\n")
                f.write("python end\n")
            self.logger.info("    > [Output] PyMOL script generated.")
        except Exception as e:
            self.logger.error(f"!!! PyMOL Script Error: {e}")

    def plot_energy_funnel(self, final_csv_path: str) -> None:
        try:
            df = pd.read_csv(final_csv_path)
            if 'dG_separated' not in df.columns or 'L_RMSD' not in df.columns:
                return

            # Try loading background data (all Stage 1 / dG-filtered survivors)
            bg_csv = os.path.join(self.root_dir, "all_scored_summary.csv")
            bg_df = None
            if os.path.exists(bg_csv):
                try:
                    bg_df = pd.read_csv(bg_csv)
                except Exception:
                    pass

            # --- Plot 1: Energy Funnel vs Relaxed ---
            output_png = os.path.join(self.root_dir, "energy_funnel.png")
            fig, ax = plt.subplots(figsize=(8, 6))

            if bg_df is not None and 'L_RMSD' in bg_df.columns:
                ax.scatter(bg_df['L_RMSD'], bg_df['dG_separated'],
                           c='gray', s=8, alpha=0.3, label='Stage 1 survivors')

            ax.scatter(df['L_RMSD'], df['dG_separated'],
                       c='red', s=50, alpha=0.8, edgecolors='darkred',
                       linewidths=0.5, label='Final ranking', zorder=5)

            ax.set_xlabel("L_RMSD vs Relaxed (A)")
            ax.set_ylabel("dG_separated (REU)")
            ax.set_title("Energy Funnel (vs Relaxed)")
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(output_png, dpi=150)
            plt.close(fig)
            self.logger.info("    > [Output] Energy funnel plot saved.")

            # --- Plot 2: Energy Funnel vs Best ---
            has_best_bg = (bg_df is not None and 'L_RMSD_best' in bg_df.columns)
            has_best_fg = 'L_RMSD_best' in df.columns
            if has_best_bg or has_best_fg:
                output_best = os.path.join(self.root_dir, "energy_funnel_best.png")
                fig2, ax2 = plt.subplots(figsize=(8, 6))

                if has_best_bg:
                    ax2.scatter(bg_df['L_RMSD_best'], bg_df['dG_separated'],
                                c='gray', s=8, alpha=0.3,
                                label='Stage 1 survivors')

                if has_best_fg:
                    ax2.scatter(df['L_RMSD_best'], df['dG_separated'],
                                c='red', s=50, alpha=0.8,
                                edgecolors='darkred', linewidths=0.5,
                                label='Final ranking', zorder=5)
                else:
                    ax2.scatter(df['L_RMSD'], df['dG_separated'],
                                c='red', s=50, alpha=0.8,
                                edgecolors='darkred', linewidths=0.5,
                                label='Final ranking', zorder=5)

                # Best model marker (star at x~0)
                best_dG = df['dG_separated'].min()
                ax2.scatter([0], [best_dG], marker='*', c='gold',
                            s=200, edgecolors='black', linewidths=0.8,
                            zorder=10, label='Best model')

                # P_near annotation
                p_near = self.stats.get('p_near')
                if p_near is not None:
                    ax2.text(0.97, 0.03,
                             f"P_near = {p_near:.3f}",
                             transform=ax2.transAxes, fontsize=9,
                             ha='right', va='bottom',
                             bbox=dict(boxstyle='round,pad=0.3',
                                       facecolor='lightyellow',
                                       edgecolor='gray', alpha=0.8))

                ax2.set_xlabel("L_RMSD vs Best (A)")
                ax2.set_ylabel("dG_separated (REU)")
                ax2.set_title("Energy Funnel (vs Best Model)")
                ax2.legend(fontsize=8)
                fig2.tight_layout()
                fig2.savefig(output_best, dpi=150)
                plt.close(fig2)
                self.logger.info(
                    "    > [Output] Energy funnel (vs best) plot saved.")

        except Exception as e:
            self.logger.error(f"!!! Plotting Error: {e}")

    def generate_cluster_overview(self) -> None:
        CLUSTER_COLORS = [
            'red', 'blue', 'green', 'orange', 'cyan',
            'magenta', 'yellow', 'purple', 'salmon', 'lime',
            'slate', 'hotpink', 'olive', 'teal', 'violet',
            'brown', 'pink', 'aquamarine', 'gold', 'wheat',
        ]
        try:
            files = sorted([os.path.join(self.dir_cluster, f)
                            for f in os.listdir(self.dir_cluster) if f.endswith(".pdb")])
            if files:
                output_pml = os.path.join(self.root_dir, "1_OVERVIEW_Clusters.pml")
                pml_dir = os.path.dirname(output_pml)
                # Extract unique cluster IDs (C01, C02, ...) from filenames
                cluster_ids = sorted(set(
                    os.path.basename(f).split('_M')[0] for f in files))
                with open(output_pml, "w") as f:
                    f.write("reinitialize\n")
                    for pdb in files:
                        rel_path = os.path.relpath(os.path.abspath(pdb), pml_dir)
                        f.write(f"load {rel_path}\n")
                    # Align by chain A
                    f.write("python\n")
                    f.write("first = cmd.get_object_list()[0]\n")
                    f.write("for obj in cmd.get_object_list()[1:]:\n")
                    f.write("    cmd.align(f'{obj} and chain A', f'{first} and chain A')\n")
                    f.write("python end\n")
                    # Display: chain A gray, chain B colored by cluster
                    f.write("hide all\nshow cartoon\n")
                    f.write("color gray80, chain A\n")
                    for i, cid in enumerate(cluster_ids):
                        color = CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
                        f.write(f"group {cid}, {cid}_*\n")
                        f.write(f"color {color}, {cid}_* and chain B\n")
        except Exception as e:
            self.logger.error(f"!!! Cluster Overview Error: {e}")

    def generate_per_cluster_views(self, final_csv_path: str) -> None:
        try:
            df = pd.read_csv(final_csv_path)
            df['Cluster_Group'] = df['Parent'].apply(lambda x: x.split('_')[0])
            for cid in df['Cluster_Group'].unique():
                cdf = df[df['Cluster_Group'] == cid]
                output_pml = os.path.join(self.dir_final, f"2_DETAIL_{cid}.pml")
                pml_dir = os.path.dirname(output_pml)
                with open(output_pml, "w") as f:
                    f.write("reinitialize\n")
                    for _, row in cdf.iterrows():
                        rel_path = os.path.relpath(os.path.abspath(row['File_PDB']), pml_dir)
                        f.write(f"load {rel_path}\n")
                    f.write("hide all\nshow cartoon\nspectrum b, blue_white_red, min=-5, max=5\n")
                    f.write("python\n")
                    f.write("first = cmd.get_object_list()[0]\n")
                    f.write("for obj in cmd.get_object_list()[1:]:\n")
                    f.write("    cmd.align(f'{obj} and chain A', f'{first} and chain A')\n")
                    f.write("python end\n")

                    # Interface residue sticks from Binding_Residues_A/B
                    resi_a = set()
                    resi_b = set()
                    _invalid = {'None', 'No_Chain_2', 'Analysis_Failed', ''}
                    for _, row in cdf.iterrows():
                        for col, rset in [('Binding_Residues_A', resi_a),
                                          ('Binding_Residues_B', resi_b)]:
                            br = str(row.get(col, ''))
                            if br in _invalid:
                                continue
                            for token in br.split(','):
                                token = token.strip()
                                # format: "A:LEU819" -> extract 819
                                if ':' in token:
                                    num_part = ''.join(c for c in token.split(':')[1] if c.isdigit() or c == '-')
                                else:
                                    num_part = ''.join(c for c in token if c.isdigit() or c == '-')
                                if num_part:
                                    rset.add(num_part)

                    if resi_a:
                        resi_str = '+'.join(sorted(resi_a, key=lambda x: int(x) if x.lstrip('-').isdigit() else 0))
                        f.write(f"select iface_A, chain A and resi {resi_str}\n")
                        f.write("show sticks, iface_A\n")
                        f.write("color cyan, iface_A\n")
                    if resi_b:
                        resi_str = '+'.join(sorted(resi_b, key=lambda x: int(x) if x.lstrip('-').isdigit() else 0))
                        f.write(f"select iface_B, chain B and resi {resi_str}\n")
                        f.write("show sticks, iface_B\n")
                        f.write("color yellow, iface_B\n")
                    if resi_a or resi_b:
                        f.write("deselect\n")
        except Exception as e:
            self.logger.error(f"!!! Per-Cluster View Error: {e}")

    # ================================================================== #
    #  Validation Report
    # ================================================================== #
    def generate_validation_report(self, ranking_df: Any) -> None:
        """Generate comprehensive docking validation report with quality assessment."""
        s = self.stats
        lines = []
        a = lines.append

        W = REPORT_WIDTH

        # ---- Header ----
        filter_ver = s.get('filter_version', 'v1.0')
        a("=" * W)
        a("        PPI 결합 사이트 탐색 리포트")
        a(f"        Global Blind Docking Pipeline (필터 {filter_ver})")
        a("=" * W)
        a("")
        a(f"  입력 구조          : {self.filename}")
        a(f"  샘플링 규모        : {self.total_global:,} 모델")
        a(f"  필터 버전          : {filter_ver}")
        auto_thr = " (auto)" if s.get('cluster_threshold_auto') else ""
        auto_top = " (auto)" if s.get('cluster_top_n_auto') else ""
        auto_div = " (auto)" if s.get('member_diversity_dist_auto') else ""
        a(f"  설정               : filter_pct={self.filter_percentile}%, "
          f"min_dSASA={self.min_dSASA}, cluster_thr={self.cluster_threshold}A{auto_thr}")
        a(f"  체인 크기          : Chain A={s.get('chain_a_size', '?')} res, "
          f"Chain B={s.get('chain_b_size', '?')} res")
        a(f"  클러스터 설정      : max={self.cluster_top_n}{auto_top}, "
          f"thr={self.cluster_threshold}A{auto_thr}, "
          f"diversity={self.member_diversity_dist}A{auto_div}")
        if self.master_seed is not None:
            a(f"  랜덤 시드          : {self.master_seed}")
        a("")

        # ---- 1. 파이프라인 실행 요약 ----
        a("=" * W)
        a("  1. 파이프라인 실행 요약")
        a("=" * W)
        a("")

        total_accepted = s['docking_success']
        total_rejected = s.get('docking_rejected', 0)
        total_errors = s.get('docking_errors', 0)
        dock_rate = total_accepted / self.total_global * 100
        a(f"  글로벌 도킹         : {total_accepted:,}/{self.total_global:,} "
          f"성공 ({dock_rate:.1f}%)")
        if total_rejected > 0 or total_errors > 0:
            a(f"  결과 분류:")
            a(f"    성공              : {total_accepted:,}개")
            if total_rejected > 0:
                rej_pct = total_rejected / self.total_global * 100
                a(f"    조기 거부         : {total_rejected:,}개 "
                  f"({rej_pct:.1f}%, 제외구역 접촉)")
            if total_errors > 0:
                err_pct = total_errors / self.total_global * 100
                a(f"    구조적 실패       : {total_errors:,}개 "
                  f"({err_pct:.1f}%, 도킹 에러)")

        # Excluded zone coverage
        excl_cov = s.get('excl_zone_coverage')
        if excl_cov is not None:
            excl_present = s.get('excl_zone_present', 0)
            excl_total = s.get('excl_zone_total', 0)
            a(f"  제외구역 커버리지   : {excl_present}/{excl_total}개 잔기 "
              f"(Chain A의 {excl_cov:.1f}%)")
            if excl_cov > 50:
                a("    ** 경고: 제외구역이 Chain A 표면의 과반 — 탐색 범위 제한")
            elif excl_cov > 30:
                a("    ** 참고: 제외구역이 넓음 — 탐색 가능 영역이 제한적일 수 있음")

        a(f"  스코어링            : {s['total_scored']:,}개 완료")
        excl_filt = s.get('excl_filter_rejected', 0)
        if excl_filt > 0:
            a(f"  제외구역 필터       : {excl_filt:,}개 제거됨 "
              f"(금지 구역)")

        if filter_ver == 'v2.0':
            # v2.0 단계별 보고
            stage1_surv = s.get('stage1_survivors', 0)
            stage2_surv = s.get('stage2_survivors', 0)
            total_scored = s['total_scored']
            stage1_pct = stage1_surv / total_scored * 100 if total_scored > 0 else 0
            a(f"  Stage 1 생존        : {stage1_surv:,}개 "
              f"(스코어링 대비 {stage1_pct:.1f}%)")
            if s.get('mini_refine_enabled', False):
                mr_ok = s.get('mini_refine_success', 0)
                mr_ct = s.get('mini_refine_count', 0)
                mr_mode = s.get('mini_refine_mode', 'pack_only')
                mr_nr = s.get('mini_refine_n_rounds', 3)
                a(f"  Mini Refinement     : {mr_ok:,}/{mr_ct:,}개 리파인됨 "
                  f"({mr_mode}, {mr_nr} rounds)")
            if stage1_surv > 0:
                a(f"  Stage 2 생존        : {stage2_surv:,}개 "
                  f"(Stage 1 대비 {stage2_surv/stage1_surv*100:.1f}%)")
            else:
                a(f"  Stage 2 생존        : {stage2_surv:,}개")
        else:
            # v1.0 보고
            if s['dG_cutoff'] is not None:
                a(f"  dG 필터 ({self.filter_percentile}%)      : "
                  f"{s['dG_passed']:,}개 통과 (기준값 = {s['dG_cutoff']:.2f} REU)")
            else:
                a(f"  dG 필터             : 건너뜀 (모든 점수 동일)")

            a(f"  dSASA 필터          : {s['dSASA_passed']:,}개 통과 "
              f"(>= {self.min_dSASA} A^2)")

        if s['fallback_used']:
            fb_level = s.get('fallback_level', 0)
            a(f"  ** 폴백 적용됨      : Level {fb_level} "
              f"(min_survivors={self.min_survivors})")

        a(f"  최종 후보           : {s['final_candidates']:,}개")
        a(f"  전체 스코어링       : {s.get('total_full_scored', 0):,}개")
        a(f"  클러스터 수         : {s['num_clusters']}")
        a(f"  대표 구조           : {s['num_representatives']}")
        a(f"  출력 모델           : {len(ranking_df)}개")
        a("")

        # ---- 1.1. 단계별 소요 시간 ----
        if self.stage_times:
            a("=" * W)
            a("  1.1. 단계별 소요 시간")
            a("=" * W)
            a("")
            stage_order = [
                'Relaxation', 'Global Docking', 'Scoring & Filtering',
                'Full Scoring', 'Clustering', 'Selection & Save'
            ]
            for stage_name in stage_order:
                elapsed = self.stage_times.get(stage_name)
                if elapsed is not None:
                    a(f"  {stage_name:<24s}: {timedelta(seconds=int(elapsed))}")
            total_rt = s.get('total_runtime', '')
            if total_rt:
                a(f"  {'총 실행 시간':<24s}: {total_rt}")
            a("")

        # ---- 1.5. Mini Refinement 통계 (v2.0, 활성 시) ----
        if s.get('mini_refine_enabled', False):
            a("=" * W)
            a("  1.5. Mini Refinement 통계")
            a("=" * W)
            a("")
            mr_mode = s.get('mini_refine_mode', 'pack_only')
            mr_nr = s.get('mini_refine_n_rounds', 3)
            mr_ct = s.get('mini_refine_count', 0)
            mr_ok = s.get('mini_refine_success', 0)
            mr_time = s.get('mini_refine_time', 0)
            a(f"  모드               : {mr_mode} ({mr_nr} rounds)")
            a(f"  대상/성공           : {mr_ok:,}/{mr_ct:,}")
            a(f"  소요 시간           : {mr_time:.1f} sec")
            a("")

            pre_dG = s.get('mini_refine_pre_dG_mean')
            post_dG = s.get('mini_refine_post_dG_mean')
            pre_dSASA = s.get('mini_refine_pre_dSASA_mean')
            post_dSASA = s.get('mini_refine_post_dSASA_mean')
            pre_sc = s.get('mini_refine_pre_sc_mean')
            post_sc = s.get('mini_refine_post_sc_mean')
            d_dG = s.get('mini_refine_delta_dG', 0)
            d_dSASA = s.get('mini_refine_delta_dSASA', 0)
            d_sc = s.get('mini_refine_delta_sc', 0)
            sc_deferred = s.get('mini_refine_sc_deferred', False)

            if pre_dG is not None:
                a("  리파인 전후 메트릭 비교 (Stage 1 생존자 평균):")
                a(f"  {'메트릭':>12s}  {'리파인 전':>10s}  {'리파인 후':>10s}  {'변화량':>10s}")
                a(f"  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*10}")
                a(f"  {'dG (REU)':>12s}  {pre_dG:>10.2f}  {post_dG:>10.2f}  {d_dG:>+10.2f}")
                a(f"  {'dSASA (A^2)':>12s}  {pre_dSASA:>10.1f}  {post_dSASA:>10.1f}  {d_dSASA:>+10.1f}")
                if sc_deferred:
                    a(f"  {'sc_value':>12s}  {'(deferred)':>10s}  {post_sc:>10.3f}  {'N/A':>10s}")
                else:
                    a(f"  {'sc_value':>12s}  {pre_sc:>10.3f}  {post_sc:>10.3f}  {d_sc:>+10.3f}")
                a("")

                if sc_deferred:
                    a("  --> sc는 Pass 1에서 스킵 → 리파인 후 최초 측정 (비용 절감)")
                elif d_sc > 0.01:
                    a("  --> Mini Refinement이 형상 상보성을 의미있게 개선했습니다.")
                elif d_sc > 0:
                    a("  --> Mini Refinement이 형상 상보성을 소폭 개선했습니다.")
                else:
                    a("  --> Mini Refinement 후 형상 상보성 변화가 미미합니다.")
                    a("      n_rounds 증가 또는 pack_min 모드를 고려하십시오.")

                # Post-refinement sc distribution (threshold calibration guide)
                sc_med = s.get('mini_refine_post_sc_median')
                sc_q25 = s.get('mini_refine_post_sc_q25')
                sc_q75 = s.get('mini_refine_post_sc_q75')
                if sc_med is not None and sc_q25 is not None and sc_q75 is not None:
                    a("")
                    a(f"  리파인 후 sc 분포: 중앙값={sc_med:.3f}, "
                      f"Q25={sc_q25:.3f}, Q75={sc_q75:.3f}")
                    stage2_sc = getattr(self, 'stage2_min_sc', None)
                    if stage2_sc is not None and sc_med > 0:
                        pct_above = sum(1 for x in [sc_q25, sc_med, sc_q75]
                                        if x is not None and x >= stage2_sc)
                        guide = "대부분" if pct_above >= 2 else "일부"
                        a(f"  현재 Stage 2 임계값: sc >= {stage2_sc}")
                        a(f"  → 생존자의 {guide}가 임계값 이상 (Q25 기준)")
            a("")

        # ---- 2. 에너지 분포 분석 ----
        a("=" * W)
        a("  2. 에너지 분포 분석")
        a("=" * W)
        a("")
        a("  전체 스코어링 모델 (모집단):")
        a(f"    dG 평균            : {s['mean_dG']:.2f} REU")
        a(f"    dG 중앙값          : {s.get('median_dG', s['mean_dG']):.2f} REU")
        a(f"    dG 표준편차        : {s['std_dG']:.2f}")
        a(f"    dG 최소 (최적)     : {s['min_dG']:.2f} REU")
        a(f"    dG 최대 (최악)     : {s['max_dG']:.2f} REU")
        a("")

        top_dG = ranking_df['dG_separated']
        a(f"  최종 랭킹 모델 (상위 {len(ranking_df)}개):")
        a(f"    dG 최적 (Rank 1)   : {top_dG.iloc[0]:.2f} REU")
        a(f"    dG 최악            : {top_dG.iloc[-1]:.2f} REU")
        a(f"    dG 평균            : {top_dG.mean():.2f} REU")
        a(f"    dG 중앙값          : {top_dG.median():.2f} REU")
        a("")

        energy_gap = s['mean_dG'] - top_dG.iloc[0]
        top_spread = top_dG.iloc[-1] - top_dG.iloc[0]
        a("  에너지 퍼널 분석:")
        a(f"    갭 (모집단 평균 - 최적) : {energy_gap:.2f} REU")
        a(f"    상위 모델 분산          : {top_spread:.2f} REU")

        eg = QUALITY_THRESHOLDS["energy_gap"]
        if energy_gap > eg["strong"]:
            a("    --> 강한 에너지 퍼널 감지됨")
        elif energy_gap > eg["moderate"]:
            a("    --> 보통 수준의 에너지 퍼널 감지됨")
        elif energy_gap > eg["weak"]:
            a("    --> 약한 에너지 퍼널")
        else:
            a("    --> 경고: 명확한 에너지 퍼널 없음 (평탄한 에너지 지형)")
        a("")

        p_near = s.get('p_near')
        if p_near is not None:
            a(f"    P_near (lambda=2.0A) : {p_near:.3f}")
            if p_near > 0.5:
                a("    --> 강한 에너지 퍼널 (best 근처에 저에너지 집중)")
            elif p_near > 0.1:
                a("    --> 보통 수준의 에너지 퍼널")
            else:
                a("    --> 퍼널 없음 (에너지가 공간적으로 분산)")
            a("")

        # ---- 3. 인터페이스 품질 ----
        a("=" * W)
        a(f"  3. 인터페이스 품질 (상위 {len(ranking_df)}개)")
        a("=" * W)
        a("")

        dsasa = ranking_df['dSASA']
        sc = ranking_df['sc_value']
        ps = ranking_df['packstat']

        a("  dSASA (매장 표면적):")
        a(f"    평균   : {dsasa.mean():.1f} A^2")
        a(f"    중앙값 : {dsasa.median():.1f} A^2")
        a(f"    범위   : {dsasa.min():.1f} ~ {dsasa.max():.1f} A^2")
        dt = QUALITY_THRESHOLDS["dSASA"]
        if dsasa.mean() > dt["excellent"]:
            a("    --> 우수 (넓고 잘 형성된 인터페이스)")
        elif dsasa.mean() > dt["good"]:
            a("    --> 양호 (의미있는 인터페이스)")
        elif dsasa.mean() > dt["marginal"]:
            a("    --> 경계 (작은 인터페이스)")
        else:
            a("    --> 경고: 매우 작은 인터페이스 - 비특이적 결합 가능성")
        a("")

        a("  형상 상보성 (sc_value):")
        a(f"    평균   : {sc.mean():.3f}")
        a(f"    중앙값 : {sc.median():.3f}")
        a(f"    범위   : {sc.min():.3f} ~ {sc.max():.3f}")
        st = QUALITY_THRESHOLDS["sc_value"]
        if sc.mean() > st["excellent"]:
            a("    --> 우수한 기하학적 맞물림")
        elif sc.mean() > st["good"]:
            a("    --> 양호한 기하학적 맞물림")
        elif sc.mean() > st["marginal"]:
            a("    --> 경계 수준의 기하학적 맞물림")
        else:
            a("    --> 경고: 불량한 기하학적 맞물림")
        a("")

        a("  원자 패킹 (packstat):")
        a(f"    평균   : {ps.mean():.3f}")
        a(f"    중앙값 : {ps.median():.3f}")
        a(f"    범위   : {ps.min():.3f} ~ {ps.max():.3f}")
        pt = QUALITY_THRESHOLDS["packstat"]
        if ps.mean() > pt["good"]:
            a("    --> 양호한 원자 패킹")
        elif ps.mean() > pt["acceptable"]:
            a("    --> 허용 가능한 원자 패킹")
        else:
            a("    --> 경고: 불량한 원자 패킹")
        a("")

        # v2.0 신규 메트릭 (dG_density, delta_unsatHbonds, nres_int)
        if 'dG_density' in ranking_df.columns:
            dg_dens = ranking_df['dG_density'].dropna()
            if len(dg_dens) > 0:
                a("  dG 밀도 (dG/dSASA*100):")
                a(f"    평균   : {dg_dens.mean():.3f}")
                a(f"    중앙값 : {dg_dens.median():.3f}")
                a(f"    범위   : {dg_dens.min():.3f} ~ {dg_dens.max():.3f}")
                dgt = QUALITY_THRESHOLDS["dG_density"]
                if dg_dens.mean() < dgt["excellent"]:
                    a("    --> 우수한 에너지 밀도")
                elif dg_dens.mean() < dgt["good"]:
                    a("    --> 양호한 에너지 밀도")
                elif dg_dens.mean() < dgt["marginal"]:
                    a("    --> 경계 수준의 에너지 밀도")
                else:
                    a("    --> 경고: 불량한 에너지 밀도")
                a("")

        if 'delta_unsatHbonds' in ranking_df.columns:
            unsat = pd.to_numeric(ranking_df['delta_unsatHbonds'], errors='coerce').dropna()
            if len(unsat) > 0 and unsat.mean() < UNSAT_FALLBACK_VALUE - 1:
                a("  미충족 수소결합 (delta_unsatHbonds):")
                a(f"    평균   : {unsat.mean():.1f}")
                a(f"    중앙값 : {unsat.median():.1f}")
                a(f"    범위   : {unsat.min():.0f} ~ {unsat.max():.0f}")
                ut = QUALITY_THRESHOLDS["unsatHb"]
                if unsat.mean() < ut["excellent"]:
                    a("    --> 우수한 수소결합 충족도")
                elif unsat.mean() < ut["acceptable"]:
                    a("    --> 양호한 수소결합 충족도")
                else:
                    a("    --> 경고: 미충족 수소결합이 많음")
                a("")

        if 'nres_int' in ranking_df.columns:
            nres = pd.to_numeric(ranking_df['nres_int'], errors='coerce').dropna()
            if len(nres) > 0 and nres.mean() > 0:
                a("  인터페이스 잔기 수 (nres_int):")
                a(f"    평균   : {nres.mean():.1f}")
                a(f"    중앙값 : {nres.median():.1f}")
                a(f"    범위   : {nres.min():.0f} ~ {nres.max():.0f}")
                nt = QUALITY_THRESHOLDS["nres_int"]
                if nres.mean() > nt["large"]:
                    a("    --> 넓은 인터페이스")
                elif nres.mean() > nt["adequate"]:
                    a("    --> 적절한 인터페이스 크기")
                else:
                    a("    --> 경고: 인터페이스 잔기 수 부족")
                a("")

        if 'hbonds_int' in ranking_df.columns:
            hb = pd.to_numeric(ranking_df['hbonds_int'], errors='coerce').dropna()
            if len(hb) > 0:
                a("  인터페이스 수소결합 (hbonds_int):")
                a(f"    평균   : {hb.mean():.1f}")
                a(f"    중앙값 : {hb.median():.1f}")
                a(f"    범위   : {hb.min():.0f} ~ {hb.max():.0f}")
                ht = QUALITY_THRESHOLDS["hbonds_int"]
                if hb.mean() >= ht["good"]:
                    a("    --> 양호한 수소결합 네트워크")
                elif hb.mean() >= ht["minimum"]:
                    a("    --> 최소한의 수소결합 존재")
                else:
                    a("    --> 경고: 인터페이스 수소결합 없음")
                a("")

        # ---- 4. 결합 사이트 탐색 분석 ----
        a("=" * W)
        a("  4. 결합 사이트 탐색 분석")
        a("=" * W)
        a("")
        a("  * Global Blind Docking은 단백질 표면 전체를 탐색하여")
        a("    다양한 후보 결합 사이트를 발견하는 것이 목적입니다.")
        a("    각 클러스터는 하나의 후보 결합 포켓을 나타냅니다.")
        a("")

        potential_extra = s.get('potential_extra_clusters', 0)
        if potential_extra > 0:
            total_possible = s.get('total_possible_clusters', s['num_clusters'])
            a(f"  발견 결합 사이트     : {s['num_clusters']}개 / {self.cluster_top_n} (최대), "
              f"실제 가능: {total_possible}개")
            if total_possible > self.cluster_top_n * 1.5:
                a(f"  --> 주의: 최대 클러스터 수({self.cluster_top_n})가 부족할 수 있습니다.")
            if s.get('cluster_top_n_escalated'):
                a(f"  클러스터 자동 확장   : {s['cluster_top_n_original']} -> {s['cluster_top_n_escalated_to']}")
        else:
            a(f"  발견 결합 사이트     : {s['num_clusters']}개 / {self.cluster_top_n} (최대)")
        cluster_sizes = s.get('cluster_sizes', [])
        if cluster_sizes:
            a(f"  최대 사이트 멤버     : {max(cluster_sizes)}개 포즈")
            a(f"  평균 사이트 멤버     : {np.mean(cluster_sizes):.1f}개 포즈")
        a("")

        # 최종 모델의 결합 사이트 분포
        ranking_copy = ranking_df.copy()
        ranking_copy['Cluster_Group'] = ranking_copy['Parent'].apply(
            lambda x: x.split('_')[0] if '_' in str(x) else x)
        cluster_counts = ranking_copy['Cluster_Group'].value_counts()

        a(f"  상위 {len(ranking_df)}개 모델 내 고유 결합 사이트 : {len(cluster_counts)}개")
        a("")
        a("  결합 사이트 분포 (최종 모델):")
        for cid, cnt in cluster_counts.head(5).items():
            pct = cnt / len(ranking_df) * 100
            bar = "#" * int(pct / 2)
            a(f"    {cid:>4s} : {cnt:2d}개 모델 ({pct:4.0f}%)  {bar}")
        a("")

        # 결합 사이트별 후보 수 (클러스터링 단계 기준)
        cluster_pops = s.get('cluster_populations', [])
        if cluster_pops:
            a("  결합 사이트별 후보 수 (클러스터링 단계 기준):")
            for i, pop in enumerate(cluster_pops):
                cid_str = f"C{i+1:02d}"
                bar = "#" * min(pop // 5, 30)
                a(f"    {cid_str} : {pop:4d}개 후보  {bar}")
            a("")

        # Boltzmann 가중 결합 사이트 확률
        boltz = s.get('boltzmann_probs', {})
        if boltz:
            a("  Boltzmann 가중 결합 사이트 확률 (kT=1.0 REU):")
            sorted_boltz = sorted(boltz.items(), key=lambda x: -x[1])
            for cid, prob in sorted_boltz[:10]:
                bar = "#" * int(prob * 50)
                a(f"    {cid:>4s} : {prob:6.1%}  {bar}")
            a("")

        # 중복 제거 통계
        n_dedup = s.get('n_dedup_removed', 0)
        if n_dedup > 0:
            a(f"  L_RMSD 중복 제거     : {n_dedup}개 유사 구조 제거됨")
            a("")

        # 사이트 탐색 결과 판정
        if cluster_pops:
            total_assigned = sum(cluster_pops)
            dominant_pop = max(cluster_pops)
            dominant_pop_pct = dominant_pop / total_assigned if total_assigned > 0 else 0.0
            dominant_pop_idx = cluster_pops.index(dominant_pop)
            dominant_pop_cid = f"C{dominant_pop_idx + 1:02d}"
            if dominant_pop_pct > 0.5:
                a(f"  --> 지배적 결합 사이트 발견: {dominant_pop_cid} "
                  f"({dominant_pop_pct:.0%}) — 최우선 후보")
            elif dominant_pop_pct > 0.3:
                a(f"  --> 유력 결합 사이트: {dominant_pop_cid} "
                  f"({dominant_pop_pct:.0%}) — 우세하나 대안 사이트도 존재")
            else:
                n_sites = len(cluster_pops)
                a(f"  --> {n_sites}개 후보 결합 사이트 발견 "
                  f"(최대 {dominant_pop_cid} {dominant_pop_pct:.0%})")
                a("      PyMOL에서 각 사이트의 생물학적 타당성을 개별 평가 권장")
        elif len(cluster_counts) > 0:
            dominant_pct = cluster_counts.iloc[0] / len(ranking_df)
            dominant_cid = cluster_counts.index[0]
            if dominant_pct > 0.5:
                a(f"  --> 지배적 결합 사이트: {dominant_cid} ({dominant_pct:.0%})")
            elif dominant_pct > 0.3:
                a(f"  --> 유력 결합 사이트: {dominant_cid} ({dominant_pct:.0%})")
            else:
                a(f"  --> {len(cluster_counts)}개 후보 사이트 발견 "
                  f"— 개별 평가 필요")
        else:
            a("  --> 결합 사이트 판정 불가 (클러스터 정보 없음)")
        a("")

        # L_RMSD 분포
        lrmsd = ranking_df['L_RMSD']
        a("  L_RMSD 분포 (최종 모델):")
        a(f"    평균   : {lrmsd.mean():.2f} A")
        a(f"    범위   : {lrmsd.min():.2f} ~ {lrmsd.max():.2f} A")
        a("    (Global Blind Docking에서는 표면 전체 탐색으로 인해")
        a("     높은 L_RMSD가 정상이며, 다양한 사이트 탐색을 반영합니다)")
        a("")

        # I_RMSD 분포 (인터페이스 RMSD)
        if 'I_RMSD' in ranking_df.columns:
            irmsd = ranking_df['I_RMSD']
            valid_irmsd = irmsd[(irmsd >= 0) & (irmsd.notna())]
            if len(valid_irmsd) > 0:
                a("  I_RMSD 분포 (Chain B 인터페이스 잔기만, 최종 모델):")
                a(f"    평균   : {valid_irmsd.mean():.2f} A")
                a(f"    범위   : {valid_irmsd.min():.2f} ~ {valid_irmsd.max():.2f} A")
                a(f"    유효   : {len(valid_irmsd)}/{len(irmsd)} 모델")
                a("    (I_RMSD는 결합 인터페이스 잔기만의 RMSD로,")
                a("     L_RMSD보다 결합 포즈 품질을 더 정확히 반영합니다)")
                a("")

        # 결합 핫스팟 분석 (Chain A / Chain B 분리)
        _invalid_tags = {'None', 'No_Chain_2', 'Analysis_Failed', ''}

        def _hotspot_section(chain_label, col_name, fallback_col=None):
            """Generate hotspot sub-section for one chain."""
            residues_all = []
            n_models = 0
            for _, row in ranking_df.iterrows():
                br = row.get(col_name, '')
                if (not isinstance(br, str) or br in _invalid_tags) and fallback_col:
                    br = row.get(fallback_col, '')
                if isinstance(br, str) and br not in _invalid_tags:
                    residues = [r.strip() for r in br.split(',') if r.strip()]
                    residues_all.extend(residues)
                    n_models += 1

            if not residues_all:
                a(f"  {chain_label} 결합 잔기: 접촉 감지되지 않음")
                a("")
                return
            res_counts = Counter(residues_all).most_common(10)
            a(f"  {chain_label} 핫스팟 잔기 ({n_models}개 접촉 모델 기준):")
            for res, cnt in res_counts:
                pct = cnt / n_models * 100 if n_models > 0 else 0
                a(f"    {res:>14s} : {cnt:2d}/{n_models} ({pct:.0f}%)")
            a("")

            n_consistent = sum(1 for _, cnt in res_counts
                               if n_models > 0
                               and cnt / n_models > 0.5)
            if n_consistent >= 5:
                a(f"  --> 강한 핫스팟: {n_consistent}개 잔기가 50%+ 모델에서 반복 출현")
            elif n_consistent >= 2:
                a(f"  --> 핫스팟 감지: {n_consistent}개 잔기가 50%+ 모델에서 반복 출현")
            else:
                a("  --> 다양한 사이트 탐색됨 — 사이트별 접촉 잔기 분석 필요")
            a("")

        _hotspot_section("Chain A (Receptor)", "Binding_Residues_A", "Binding_Residues")
        _hotspot_section("Chain B (Ligand)", "Binding_Residues_B")

        # 핵심 잔기 접촉 요약 (제약 조건 활성 시)
        if self.key_residues_B and 'key_contact_ratio' in ranking_df.columns:
            a("  핵심 잔기 접촉 (Chain B, 소프트 메트릭):")
            kcr = ranking_df['key_contact_ratio']
            a(f"    평균 접촉률        : {kcr.mean():.2%}")
            a(f"    범위               : {kcr.min():.2%} ~ {kcr.max():.2%}")
            n_with_key = (kcr > 0).sum()
            a(f"    핵심잔기 접촉 모델 : {n_with_key}/{len(ranking_df)}")
            if kcr.mean() > 0.5:
                a("    --> 핵심 잔기가 인터페이스에 자주 참여함")
            elif kcr.mean() > 0.2:
                a("    --> 핵심 잔기가 보통 수준으로 관여함")
            else:
                a("    --> 핵심 잔기가 인터페이스에 거의 없음 (다른 부위에 결합 가능)")
            a("")

        # 제약 조건 요약
        if self.excluded_residues_A or self.key_residues_B:
            a("  제약 조건 설정:")
            if self.excluded_residues_A:
                a(f"    제외 Chain A 잔기  : {len(self.excluded_residues_A)}개")
                a(f"    조기 거부 활성화   : {self.enable_early_rejection}")
                a(f"    제외 접촉 거리     : {self.exclusion_contact_dist} A")
            if self.key_residues_B:
                a(f"    핵심 Chain B 잔기  : {len(self.key_residues_B)}개")
                a(f"    보너스 가중치      : {self.key_residue_bonus_weight}")
            a("")

        # ---- 4.5. 사이트 그룹핑 (Jaccard) ----
        site_groups = s.get('site_groups')
        if site_groups:
            a("=" * W)
            a("  4.5. 사이트 그룹핑 (Jaccard 유사도)")
            a("=" * W)
            a("")
            a(f"  Jaccard 임계값      : {self.site_group_jaccard}")
            a(f"  그룹 수             : {site_groups.get('n_groups', 0)}개 "
              f"(클러스터 {len(site_groups.get('cluster_to_group', {}))}개)")
            a("")

            group_boltz = s.get('site_group_boltzmann', {})
            cluster_pops_list = s.get('cluster_populations', [])

            a("  그룹별 상세:")
            for gid, members in enumerate(site_groups.get('groups', []), 1):
                member_str = ', '.join(members)
                total_pop = 0
                for m in members:
                    try:
                        cidx = int(m[1:]) - 1
                        if 0 <= cidx < len(cluster_pops_list):
                            total_pop += cluster_pops_list[cidx]
                    except (ValueError, IndexError):
                        pass
                boltz_p = group_boltz.get(gid, 0)
                a(f"    Group {gid}: [{member_str}]  Pop={total_pop}  "
                  f"Boltzmann={boltz_p:.1%}")
            a("")

        # ---- 5. 종합 판정 ----
        a("=" * W)
        a("  5. 종합 판정")
        a("=" * W)
        a("")

        # 핵심 체크: C2(결합 에너지), C3(에너지 퍼널), C4(인터페이스 크기)
        CORE_CHECKS = {'C2_결합에너지', 'C3_에너지퍼널', 'C4_인터페이스크기'}
        checks = []

        # C1: 파이프라인 실행 (조기 거부 + 에러는 정상/예외 동작이므로 유효 성공률 사용)
        pip_t = QUALITY_THRESHOLDS["pipeline"]
        excluded_from_rate = total_rejected + total_errors
        if excluded_from_rate > 0:
            effective_attempted = self.total_global - excluded_from_rate
            effective_rate = (total_accepted / effective_attempted * 100
                              if effective_attempted > 0 else 0)
            c1_rate = effective_rate
            detail_parts = []
            if total_rejected > 0:
                detail_parts.append(f'조기거부 {total_rejected:,}개')
            if total_errors > 0:
                detail_parts.append(f'에러 {total_errors:,}개')
            c1_detail = (f'{dock_rate:.1f}% 전체 '
                         f'(유효 {effective_rate:.1f}%, '
                         f'{", ".join(detail_parts)} 제외)')
        else:
            c1_rate = dock_rate
            c1_detail = f'{dock_rate:.1f}% 성공'
        if c1_rate >= pip_t["good"]:
            checks.append(('C1_파이프라인', 'PASS', c1_detail))
        elif c1_rate >= pip_t["marginal"]:
            checks.append(('C1_파이프라인', 'WARNING', c1_detail))
        else:
            checks.append(('C1_파이프라인', 'FAIL', c1_detail))

        # C2: 결합 에너지
        best_dG = top_dG.iloc[0]
        be_t = QUALITY_THRESHOLDS["binding_energy"]
        if s['std_dG'] < 1e-6:
            checks.append(('C2_결합에너지', 'FAIL',
                           'dG 모두 동일 (분화 없음)'))
        elif best_dG < be_t["good"]:
            checks.append(('C2_결합에너지', 'PASS', f'{best_dG:.2f} REU'))
        elif best_dG < be_t["marginal"]:
            checks.append(('C2_결합에너지', 'WARNING', f'{best_dG:.2f} REU'))
        else:
            checks.append(('C2_결합에너지', 'FAIL',
                           f'{best_dG:.2f} REU (비호의적)'))

        # C3: 에너지 퍼널 (energy_gap + P_near 복합 판정)
        p_near = s.get('p_near')
        if p_near is not None:
            c3_detail = f'{energy_gap:.1f} REU 갭, P_near={p_near:.3f}'
            if energy_gap > eg["strong"] and p_near > 0.1:
                checks.append(('C3_에너지퍼널', 'PASS', c3_detail))
            elif energy_gap > eg["moderate"] or p_near > 0.3:
                checks.append(('C3_에너지퍼널', 'WARNING', c3_detail))
            else:
                checks.append(('C3_에너지퍼널', 'FAIL',
                               c3_detail + ' (평탄한 지형)'))
        else:
            # P_near 없으면 기존 energy_gap 단독 판정 유지 (하위호환)
            if energy_gap > eg["strong"]:
                checks.append(('C3_에너지퍼널', 'PASS',
                               f'{energy_gap:.1f} REU 갭'))
            elif energy_gap > eg["moderate"]:
                checks.append(('C3_에너지퍼널', 'WARNING',
                               f'{energy_gap:.1f} REU 갭'))
            else:
                checks.append(('C3_에너지퍼널', 'FAIL',
                               f'{energy_gap:.1f} REU 갭 (평탄한 지형)'))

        # C4: 인터페이스 크기
        mean_dsasa = dsasa.mean()
        is_t = QUALITY_THRESHOLDS["interface_size"]
        if mean_dsasa > is_t["good"]:
            checks.append(('C4_인터페이스크기', 'PASS',
                           f'{mean_dsasa:.0f} A^2'))
        elif mean_dsasa > is_t["marginal"]:
            checks.append(('C4_인터페이스크기', 'WARNING',
                           f'{mean_dsasa:.0f} A^2'))
        else:
            checks.append(('C4_인터페이스크기', 'FAIL',
                           f'{mean_dsasa:.0f} A^2 (너무 작음)'))

        # C5: 형상 상보성
        mean_sc = sc.mean()
        sf_t = QUALITY_THRESHOLDS["shape_fit"]
        if mean_sc > sf_t["good"]:
            checks.append(('C5_형상상보성', 'PASS', f'{mean_sc:.3f}'))
        elif mean_sc > sf_t["marginal"]:
            checks.append(('C5_형상상보성', 'WARNING', f'{mean_sc:.3f}'))
        else:
            checks.append(('C5_형상상보성', 'FAIL', f'{mean_sc:.3f}'))

        # C6: 사이트 탐색 다양성 (PPI 사이트 탐색에서는 다양한 사이트 발견이 정상)
        cv_t = QUALITY_THRESHOLDS["convergence"]
        if cluster_pops:
            total_assigned = sum(cluster_pops)
            c6_dominant = max(cluster_pops)
            c6_pct = c6_dominant / total_assigned if total_assigned > 0 else 0.0
            c6_idx = cluster_pops.index(c6_dominant)
            c6_cid = f"C{c6_idx + 1:02d}"
            c6_n_sites = len(cluster_pops)
        elif len(cluster_counts) > 0:
            c6_pct = cluster_counts.iloc[0] / len(ranking_df)
            c6_cid = cluster_counts.index[0]
            c6_n_sites = len(cluster_counts)
        else:
            c6_pct = 0.0
            c6_cid = "N/A"
            c6_n_sites = 0
        if c6_pct > cv_t["strong"]:
            checks.append(('C6_사이트탐색', 'PASS',
                           f'지배적 사이트 {c6_cid} ({c6_pct:.0%})'))
        elif c6_pct > cv_t["moderate"]:
            checks.append(('C6_사이트탐색', 'PASS',
                           f'유력 사이트 {c6_cid} ({c6_pct:.0%})'))
        else:
            checks.append(('C6_사이트탐색', 'PASS',
                           f'{c6_n_sites}개 후보 사이트 발견'))

        # C7: 샘플링 충분성
        samp_t = QUALITY_THRESHOLDS["sampling"]
        if self.total_global >= samp_t["good"]:
            checks.append(('C7_샘플링규모', 'PASS',
                           f'{self.total_global:,}개 모델'))
        elif self.total_global >= samp_t["marginal"]:
            checks.append(('C7_샘플링규모', 'WARNING',
                           f'{self.total_global:,}개 모델'))
        else:
            checks.append(('C7_샘플링규모', 'FAIL',
                           f'{self.total_global:,}개 모델 (부족)'))

        # C8: dG 밀도 (v2.0 메트릭, 가용 시에만)
        dgc_t = QUALITY_THRESHOLDS["dG_density_check"]
        if 'dG_density' in ranking_df.columns:
            dg_dens_vals = ranking_df['dG_density'].dropna()
            if len(dg_dens_vals) > 0:
                mean_dens = dg_dens_vals.mean()
                if mean_dens < dgc_t["good"]:
                    checks.append(('C8_dG밀도', 'PASS',
                                   f'{mean_dens:.3f}'))
                elif mean_dens < dgc_t["marginal"]:
                    checks.append(('C8_dG밀도', 'WARNING',
                                   f'{mean_dens:.3f}'))
                else:
                    checks.append(('C8_dG밀도', 'FAIL',
                                   f'{mean_dens:.3f} (에너지 밀도 낮음)'))

        # C9: 인터페이스 잔기 수 (v2.0 메트릭, 가용 시에만)
        nrc_t = QUALITY_THRESHOLDS["nres_check"]
        if 'nres_int' in ranking_df.columns:
            nres_vals = pd.to_numeric(ranking_df['nres_int'], errors='coerce').dropna()
            if len(nres_vals) > 0 and nres_vals.mean() > 0:
                mean_nres = nres_vals.mean()
                if mean_nres > nrc_t["good"]:
                    checks.append(('C9_인터페이스잔기', 'PASS',
                                   f'{mean_nres:.0f}개 잔기'))
                elif mean_nres > nrc_t["marginal"]:
                    checks.append(('C9_인터페이스잔기', 'WARNING',
                                   f'{mean_nres:.0f}개 잔기'))
                else:
                    checks.append(('C9_인터페이스잔기', 'FAIL',
                                   f'{mean_nres:.0f}개 잔기 (부족)'))

        # C10: 인터페이스 수소결합 (v2.0 메트릭, 가용 시에만)
        if 'hbonds_int' in ranking_df.columns:
            hb_vals = pd.to_numeric(ranking_df['hbonds_int'], errors='coerce').dropna()
            if len(hb_vals) > 0:
                mean_hb = hb_vals.mean()
                ht = QUALITY_THRESHOLDS["hbonds_int"]
                if mean_hb >= ht["good"]:
                    checks.append(('C10_수소결합', 'PASS',
                                   f'{mean_hb:.1f}개'))
                elif mean_hb >= ht["minimum"]:
                    checks.append(('C10_수소결합', 'WARNING',
                                   f'{mean_hb:.1f}개'))
                else:
                    checks.append(('C10_수소결합', 'FAIL',
                                   f'{mean_hb:.1f}개 (없음)'))

        for name, status, detail in checks:
            symbol = {'PASS': '[통과]   ', 'WARNING': '[주의]   ',
                      'FAIL': '[실패]   '}[status]
            a(f"  {symbol} {name:<22s} : {detail}")

        a("")

        n_pass = sum(1 for _, st, _ in checks if st == 'PASS')
        n_warn = sum(1 for _, st, _ in checks if st == 'WARNING')
        n_fail = sum(1 for _, st, _ in checks if st == 'FAIL')
        core_fail = sum(1 for nm, st, _ in checks
                        if st == 'FAIL' and nm in CORE_CHECKS)

        a(f"  점수: {n_pass} 통과 / {n_warn} 주의 / {n_fail} 실패"
          f"  (핵심 실패: {core_fail})")
        a("")

        # 판정 (핵심/보조 가중, PPI 사이트 탐색 관점)
        if core_fail == 0 and n_fail == 0 and n_warn <= 1:
            verdict = "높은 신뢰도"
            desc_lines = [
                "강한 결합 에너지와 잘 형성된 인터페이스를 가진",
                "후보 결합 사이트가 발견되었습니다.",
                "PyMOL로 상위 사이트의 생물학적 타당성을 확인하십시오.",
            ]
        elif core_fail == 0 and n_fail <= 1:
            verdict = "보통 신뢰도"
            desc_lines = [
                "후보 결합 사이트가 발견되었으나 일부 주의사항이 있습니다.",
                "주의/실패 항목을 검토하고, PyMOL에서 각 사이트를",
                "시각적으로 확인한 후 유력 후보를 선별하십시오.",
            ]
        elif core_fail <= 1 and n_fail <= 3:
            verdict = "낮은 신뢰도"
            desc_lines = [
                "후보 사이트가 발견되었으나 인터페이스 품질에 우려가 있습니다.",
                "실패 항목을 참조하고, 실험적 검증을 병행하십시오.",
                "입력 구조 품질 확인 및 제약 조건 재검토를 권장합니다.",
            ]
        else:
            verdict = "신뢰 불가"
            desc_lines = [
                "다수의 품질 기준이 실패했습니다.",
                "발견된 사이트가 실제 결합을 반영하지 않을 수 있습니다.",
                "입력 PDB, 체인 정의, 제약 조건을 확인하십시오.",
            ]

        a(f"  +{'=' * 68}+")
        a(f"  |{' ' * 68}|")
        a(f"  |  종합 판정: {verdict:<53s}|")
        a(f"  |{' ' * 68}|")
        for dl in desc_lines:
            a(f"  |  {dl:<65s}|")
        a(f"  |{' ' * 68}|")
        a(f"  +{'=' * 68}+")
        a("")

        # ---- 6. 실험 데이터 상관 분석 (활성 시) ----
        if getattr(self, 'exp_data_configured', False):
            exp_corr = self._compute_exp_correlation(ranking_df)
            if exp_corr:
                a("=" * W)
                a("  6. 실험 데이터 상관 분석")
                a("=" * W)
                a("")
                a(f"  신뢰도             : {getattr(self, 'exp_confidence', 'medium')}")

                critical = getattr(self, 'exp_critical_B', set())
                non_binding = getattr(self, 'exp_non_binding_B', set())
                if critical:
                    a(f"  핵심 잔기 (Chain B) : {len(critical)}개 ({sorted(critical)[:10]}{'...' if len(critical) > 10 else ''})")
                if non_binding:
                    a(f"  비결합 잔기 (Chain B): {len(non_binding)}개")
                a("")

                sens = exp_corr.get('sensitivity')
                spec = exp_corr.get('specificity')
                f1 = exp_corr.get('f1')
                enrich = exp_corr.get('enrichment')

                if sens is not None:
                    a(f"  민감도 (Sensitivity) : {sens:.2%}  (핵심 잔기 중 감지된 비율)")
                if spec is not None:
                    a(f"  특이도 (Specificity) : {spec:.2%}  (비결합 잔기 중 올바르게 제외)")
                if f1 is not None:
                    a(f"  F1 Score            : {f1:.3f}")
                if enrich is not None and enrich > 0:
                    a(f"  농축 배율 (Enrichment): {enrich:.2f}x")
                a("")

                # Per-cluster sensitivity
                cluster_sens = exp_corr.get('cluster_sensitivity', {})
                if cluster_sens:
                    a("  사이트별 핵심잔기 감지:")
                    sorted_cs = sorted(cluster_sens.items(),
                                       key=lambda x: -x[1][0])
                    for cid, (found_n, total_n, found_list) in sorted_cs[:15]:
                        pct = found_n / total_n * 100 if total_n > 0 else 0
                        residues_str = ','.join(str(r) for r in found_list[:5])
                        if len(found_list) > 5:
                            residues_str += '...'
                        a(f"    {cid:>4s} : {found_n:2d}/{total_n} "
                          f"({pct:4.0f}%)  [{residues_str}]")
                    a("")

                # Found / missing
                found = exp_corr.get('critical_found', set())
                missing = exp_corr.get('critical_missing', set())
                if found:
                    a(f"  감지된 핵심 잔기    : {sorted(found)}")
                if missing:
                    a(f"  미감지 핵심 잔기    : {sorted(missing)}")
                fp = exp_corr.get('false_positives', set())
                if fp:
                    a(f"  위양성 (비결합에 접촉): {sorted(fp)}")
                a("")

        # ---- 추천사항 ----
        a("  추천사항:")

        if self.total_global < 5000:
            a(f"  - 현재 샘플링 ({self.total_global:,}개 모델)이 비교적 적습니다.")
            a("    프로덕션 품질을 위해 10,000~20,000개 모델을 권장합니다.")

        if s['fallback_used']:
            fb_lv = s.get('fallback_level', 0)
            a(f"  - 필터 폴백 Level {fb_lv}이 작동했습니다. Stage 2 엄격 기준을")
            a("    통과한 포즈가 부족합니다. 필터 임계값 완화 또는 샘플링 증가를")
            a("    고려하십시오. (PPI 인터페이스는 패킹이 느슨할 수 있음)")

        potential_extra = s.get('potential_extra_clusters', 0)
        if potential_extra > 0:
            total_possible = s.get('total_possible_clusters', 0)
            a(f"  - 클러스터 최대값({self.cluster_top_n})에 도달 후 {potential_extra}개의 추가")
            a(f"    사이트가 감지되었습니다 (실제 가능: {total_possible}개).")
            dropped_path = s.get('dropped_candidates_path', '')
            n_dropped = s.get('n_dropped_candidates', 0)
            if dropped_path:
                a(f"    --> 드롭된 후보 {n_dropped}개의 메트릭: {os.path.basename(dropped_path)}")
            if s.get('cluster_top_n_escalated'):
                a(f"    --> 자동 확장 적용: {s['cluster_top_n_original']} -> {s['cluster_top_n_escalated_to']}")
            elif s.get('cluster_top_n_auto'):
                a("    auto 모드 결정값이지만, 수동으로 더 큰 값을 설정할 수 있습니다.")
            else:
                a("    cluster_top_n을 늘리거나 auto로 설정하는 것을 고려하십시오.")

        FAIL_RECOMMENDATIONS = {
            'C1_파이프라인': [
                "  - 파이프라인 성공률이 낮습니다. 입력 PDB 구조와",
                "    체인 정의를 확인하십시오. 제외구역 설정이 과도할 수 있습니다.",
            ],
            'C3_에너지퍼널': [
                "  - 명확한 에너지 퍼널이 없습니다. 샘플링을 늘리거나",
                "    이 단백질 쌍에 선호 결합 사이트가 존재하는지 확인하십시오.",
            ],
            'C4_인터페이스크기': [
                "  - 인터페이스가 매우 작습니다. 양쪽 체인이 완전하고",
                "    입력 PDB에서 적절히 분리되어 있는지 확인하십시오.",
            ],
            'C5_형상상보성': [
                "  - 형상 상보성이 불량합니다. 기하학적 불일치를 나타냅니다.",
                "    결합 표면의 구조적 변화가 필요할 수 있습니다.",
            ],
            'C7_샘플링규모': [
                "  - total_global_models를 최소 10,000개로 늘리십시오.",
                "    신뢰할 수 있는 사이트 탐색을 위해 20,000개 이상을 권장합니다.",
            ],
            'C8_dG밀도': [
                "  - dG 밀도가 낮습니다. 인터페이스 에너지가 너무 분산되어 있습니다.",
                "    크지만 최적화되지 않은 인터페이스일 수 있습니다.",
            ],
            'C9_인터페이스잔기': [
                "  - 인터페이스 잔기 수가 매우 적습니다.",
                "    결합 인터페이스가 너무 작아 유효한 결합 사이트가 아닐 수 있습니다.",
            ],
            'C10_수소결합': [
                "  - 인터페이스 수소결합이 감지되지 않았습니다.",
                "    생물학적 인터페이스는 거의 항상 교차 수소결합을 포함합니다.",
            ],
        }

        for name, status, detail in checks:
            if status != 'FAIL':
                continue
            # C2: 결합 에너지 — detail 기반 분기
            if name == 'C2_결합에너지':
                if '동일' in detail:
                    a("  - 모든 dG 점수가 동일합니다. 체인 간 접촉이 없을 수 있습니다.")
                    a("    입력 PDB 체인 정의와 jump 설정을 확인하십시오.")
                else:
                    a("  - 결합 에너지가 비호의적입니다. 안정적인 복합체를 형성하지")
                    a("    않거나 샘플링이 부족했을 수 있습니다.")
                continue
            if name in FAIL_RECOMMENDATIONS:
                for line in FAIL_RECOMMENDATIONS[name]:
                    a(line)

        # PPI 사이트 탐색 관점의 후속 조치 가이드
        a("")
        a("  후속 분석 가이드:")
        a("  - 1_OVERVIEW_Clusters.pml : 모든 결합 사이트를 한눈에 비교")
        a("  - 2_DETAIL_C##.pml       : 각 사이트의 멤버 포즈 상세 확인")
        a("  - view_results.pml       : 최종 랭킹 모델 (B-factor 컬러링)")
        a("  - cluster_summary.csv    : 사이트별 에너지/메트릭 수치 비교")
        a("  - 상위 사이트를 실험 데이터 (돌연변이, 교차결합 등)와 대조하여")
        a("    생물학적으로 가장 타당한 결합 사이트를 선별하십시오.")

        a("")

        # ---- 참고: 품질 기준표 ----
        a("=" * W)
        a("  참고: 품질 기준 임계값 (본 파이프라인 체크 기준)")
        a("-" * W)
        a("  메트릭               우수         양호        경계         불량")
        a("-" * W)
        a("  dG_separated (REU)   < -10        -10 ~ -3    -3 ~ 0       > 0")
        a("  dSASA (A^2)          > 1000       500~1000    300~500      < 300")
        a("  dG_density           < -1.5       -1.0~-1.5   -0.5~-1.0    > -0.5")
        a("  sc_value             > 0.70       0.50~0.70   0.35~0.50    < 0.35")
        a("  packstat             > 0.55       0.40~0.55   -            < 0.40")
        a("  delta_unsatHbonds    < 5          5~10        -            > 10")
        a("  hbonds_int           >= 3         1~3         -            0")
        a("  nres_int             > 20         15~20       8~15         < 8")
        a("  에너지 갭 (REU)      > 10         5~10        2~5          < 2")
        a("  수렴도 (지배 비율)   > 40%        20~40%      -            < 20%")
        a("-" * W)
        a("  * Rosetta 도킹 문헌 참고 (Gray 2003, Chaudhury 2011)")
        a("")

        # ---- 하단 정보 ----
        a("=" * W)
        a(f"  생성 시각  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        a(f"  실행 시간  : {s.get('total_runtime', 'N/A')}")
        a("=" * W)

        report_text = "\n".join(lines)

        # Save report
        report_path = os.path.join(self.dir_final, "docking_validation_report.txt")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(report_text)

        report_path_root = os.path.join(self.root_dir, "docking_validation_report.txt")
        with open(report_path_root, "w", encoding="utf-8") as f:
            f.write(report_text)

        self.logger.info(f"    > [Output] Validation report: {report_path}")

        # 콘솔에 간략 판정 출력
        self.logger.info("")
        self.logger.info(f"    === 판정: {verdict} "
                                  f"({n_pass}통과/{n_warn}주의/{n_fail}실패) ===")
        if n_fail > 0:
            for name, status, detail in checks:
                if status == 'FAIL':
                    self.logger.info(f"    [실패] {name}: {detail}")

    # ================================================================== #
    #  HTML Dashboard
    # ================================================================== #
    def generate_html_dashboard(self, ranking_df, final_csv_path: str) -> None:
        """Generate a single-file HTML dashboard with inline CSS/JS."""
        s = self.stats

        # Encode funnel plot images
        img_tags = {}
        for img_name, label in [('energy_funnel.png', 'Energy Funnel'),
                                ('energy_funnel_best.png', 'Energy Funnel (Best)')]:
            img_path = os.path.join(self.root_dir, img_name)
            if os.path.exists(img_path):
                with open(img_path, 'rb') as f:
                    b64 = base64.b64encode(f.read()).decode('ascii')
                img_tags[label] = f'<img src="data:image/png;base64,{b64}" style="max-width:100%;height:auto;" alt="{label}">'

        # Build ranking table data
        table_cols = ['Rank', 'Parent', 'dG_separated', 'dSASA', 'sc_value',
                      'packstat', 'dG_density', 'L_RMSD', 'I_RMSD']
        if 'Site_Group' in ranking_df.columns:
            table_cols.append('Site_Group')
        if 'key_contact_ratio' in ranking_df.columns:
            table_cols.append('key_contact_ratio')
        table_data = ranking_df[
            [c for c in table_cols if c in ranking_df.columns]
        ].to_dict('records')

        # Summary cards
        filter_ver = s.get('filter_version', 'v1.0')
        best_dG = ranking_df['dG_separated'].iloc[0] if len(ranking_df) > 0 else 0
        mean_sc = ranking_df['sc_value'].mean() if len(ranking_df) > 0 else 0
        n_clusters = s.get('num_clusters', 0)
        total_rt = s.get('total_runtime', 'N/A')

        # Boltzmann bar data
        boltz = s.get('boltzmann_probs', {})
        boltz_sorted = sorted(boltz.items(), key=lambda x: -x[1])[:15]

        # Site groups
        site_groups = s.get('site_groups')
        group_boltz = s.get('site_group_boltzmann', {})

        # Stage times
        stage_times = self.stage_times

        html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PPI Docking Dashboard - {self.filename}</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       background: #f5f7fa; color: #333; padding: 20px; }}
.header {{ background: linear-gradient(135deg, #1a365d, #2563eb); color: white;
           padding: 24px 32px; border-radius: 12px; margin-bottom: 20px; }}
.header h1 {{ font-size: 1.5em; margin-bottom: 8px; }}
.header .sub {{ opacity: 0.85; font-size: 0.9em; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
          gap: 12px; margin-bottom: 20px; }}
.card {{ background: white; border-radius: 10px; padding: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.card .label {{ font-size: 0.75em; color: #666; text-transform: uppercase; letter-spacing: 0.5px; }}
.card .value {{ font-size: 1.5em; font-weight: 700; margin-top: 4px; }}
.card .value.good {{ color: #059669; }}
.card .value.warn {{ color: #d97706; }}
.card .value.bad {{ color: #dc2626; }}
.section {{ background: white; border-radius: 10px; padding: 20px; margin-bottom: 16px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.section h2 {{ font-size: 1.1em; margin-bottom: 12px; color: #1a365d;
               border-bottom: 2px solid #e2e8f0; padding-bottom: 8px; }}
.img-row {{ display: flex; gap: 16px; flex-wrap: wrap; }}
.img-row img {{ border-radius: 8px; border: 1px solid #e2e8f0; }}
table {{ width: 100%; border-collapse: collapse; font-size: 0.85em; }}
th {{ background: #f1f5f9; padding: 8px 10px; text-align: left; cursor: pointer;
      user-select: none; position: sticky; top: 0; }}
th:hover {{ background: #e2e8f0; }}
td {{ padding: 6px 10px; border-bottom: 1px solid #f1f5f9; }}
tr:hover td {{ background: #f8fafc; }}
.bar-container {{ display: flex; align-items: center; gap: 8px; margin: 3px 0; }}
.bar-label {{ min-width: 40px; font-size: 0.85em; font-weight: 600; }}
.bar {{ height: 20px; background: linear-gradient(90deg, #3b82f6, #60a5fa);
        border-radius: 4px; transition: width 0.3s; }}
.bar-pct {{ font-size: 0.8em; color: #666; min-width: 50px; }}
.timing-row {{ display: flex; justify-content: space-between; padding: 4px 0;
               border-bottom: 1px solid #f1f5f9; font-size: 0.9em; }}
.group-badge {{ display: inline-block; padding: 2px 8px; border-radius: 10px;
                font-size: 0.75em; font-weight: 600; background: #e0e7ff; color: #3730a3; }}
</style>
</head>
<body>
<div class="header">
  <h1>PPI 결합 사이트 탐색 대시보드</h1>
  <div class="sub">{self.filename} | {self.total_global:,} models | Filter {filter_ver} | {datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
</div>

<div class="cards">
  <div class="card">
    <div class="label">Best dG</div>
    <div class="value {'good' if best_dG < -10 else 'warn' if best_dG < -3 else 'bad'}">{best_dG:.2f}</div>
  </div>
  <div class="card">
    <div class="label">Mean sc</div>
    <div class="value {'good' if mean_sc > 0.5 else 'warn' if mean_sc > 0.35 else 'bad'}">{mean_sc:.3f}</div>
  </div>
  <div class="card">
    <div class="label">사이트 수</div>
    <div class="value">{n_clusters}</div>
  </div>
  <div class="card">
    <div class="label">출력 모델</div>
    <div class="value">{len(ranking_df)}</div>
  </div>
  <div class="card">
    <div class="label">실행 시간</div>
    <div class="value" style="font-size:1.1em">{total_rt}</div>
  </div>
</div>
"""

        # Energy landscape section
        if img_tags:
            html += '<div class="section"><h2>에너지 랜드스케이프</h2><div class="img-row">\n'
            for label, tag in img_tags.items():
                html += f'<div style="flex:1;min-width:300px">{tag}<p style="text-align:center;font-size:0.8em;color:#666;margin-top:4px">{label}</p></div>\n'
            html += '</div></div>\n'

        # Boltzmann bar chart
        if boltz_sorted:
            html += '<div class="section"><h2>Boltzmann 사이트 확률</h2>\n'
            max_prob = boltz_sorted[0][1] if boltz_sorted else 1
            for cid, prob in boltz_sorted:
                bar_w = prob / max_prob * 100 if max_prob > 0 else 0
                html += f'<div class="bar-container"><span class="bar-label">{html_mod.escape(str(cid))}</span>'
                html += f'<div class="bar" style="width:{bar_w}%"></div>'
                html += f'<span class="bar-pct">{prob:.1%}</span></div>\n'
            html += '</div>\n'

        # Site groups
        if site_groups and site_groups.get('n_groups', 0) > 1:
            html += '<div class="section"><h2>사이트 그룹핑</h2>\n'
            html += f'<p style="font-size:0.85em;color:#666;margin-bottom:8px">Jaccard 임계값: {self.site_group_jaccard}</p>\n'
            for gid, members in enumerate(site_groups.get('groups', []), 1):
                member_str = html_mod.escape(', '.join(members))
                gb = group_boltz.get(gid, 0)
                html += f'<div style="margin:4px 0"><span class="group-badge">Group {gid}</span> '
                html += f'{member_str} (Boltzmann: {gb:.1%})</div>\n'
            html += '</div>\n'

        # Stage timing
        if stage_times:
            html += '<div class="section"><h2>단계별 소요 시간</h2>\n'
            stage_order = ['Relaxation', 'Global Docking', 'Scoring & Filtering',
                           'Full Scoring', 'Clustering', 'Selection & Save']
            for sname in stage_order:
                elapsed = stage_times.get(sname)
                if elapsed is not None:
                    html += f'<div class="timing-row"><span>{sname}</span><span>{timedelta(seconds=int(elapsed))}</span></div>\n'
            html += '</div>\n'

        # Metrics table
        html += '<div class="section"><h2>상세 메트릭 테이블</h2>\n'
        html += '<div style="overflow-x:auto"><table id="metrics-table"><thead><tr>\n'
        display_cols = [c for c in table_cols if c in ranking_df.columns]
        for col in display_cols:
            html += f'<th onclick="sortTable(this)">{html_mod.escape(str(col))}</th>\n'
        html += '</tr></thead><tbody>\n'
        for row in table_data:
            html += '<tr>'
            for col in display_cols:
                val = row.get(col, '')
                if isinstance(val, float):
                    val = f'{val:.3f}' if abs(val) < 10 else f'{val:.2f}'
                html += f'<td>{html_mod.escape(str(val))}</td>'
            html += '</tr>\n'
        html += '</tbody></table></div></div>\n'

        # Experimental data section
        if getattr(self, 'exp_data_configured', False):
            exp_corr = self._compute_exp_correlation(ranking_df)
            if exp_corr:
                html += '<div class="section"><h2>실험 데이터 상관</h2>\n'
                html += '<div class="cards" style="margin-bottom:12px">\n'
                for metric, label in [('sensitivity', '민감도'), ('specificity', '특이도'),
                                      ('f1', 'F1'), ('enrichment', '농축배율')]:
                    val = exp_corr.get(metric)
                    if val is not None:
                        fmt = f'{val:.2%}' if metric != 'enrichment' else f'{val:.2f}x'
                        html += f'<div class="card"><div class="label">{label}</div><div class="value">{fmt}</div></div>\n'
                html += '</div></div>\n'

        # Sort script
        html += """
<script>
function sortTable(th) {
  var table = th.closest('table');
  var tbody = table.querySelector('tbody');
  var rows = Array.from(tbody.querySelectorAll('tr'));
  var idx = Array.from(th.parentNode.children).indexOf(th);
  var asc = th.dataset.sort !== 'asc';
  th.parentNode.querySelectorAll('th').forEach(function(t){delete t.dataset.sort;});
  th.dataset.sort = asc ? 'asc' : 'desc';
  rows.sort(function(a, b) {
    var va = a.children[idx].textContent.trim();
    var vb = b.children[idx].textContent.trim();
    var na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
    return asc ? va.localeCompare(vb) : vb.localeCompare(va);
  });
  rows.forEach(function(r){tbody.appendChild(r);});
}
</script>
</body></html>"""

        dashboard_path = os.path.join(self.dir_final, "dashboard.html")
        with open(dashboard_path, 'w', encoding='utf-8') as f:
            f.write(html)
        self.logger.info(f"    > [Output] HTML Dashboard: {dashboard_path}")

    # ================================================================== #
    #  PIPELINE STEPS (private methods)
    # ================================================================== #

    def _step1_relax(self) -> Tuple[Optional[str], Optional[str]]:
        """Step 1: Check/run relaxed structure. Returns (relaxed_pdb, cache_path)."""
        self.logger.info("=" * 60)
        self.logger.info(">>> [Step 1/7] Check Relaxed Structure")
        self.logger.info("=" * 60)
        t_start = time.time()

        cache_filename = f"{os.path.splitext(self.filename)[0]}_relaxed.pdb"
        cache_path = os.path.join(self.dir_cache, cache_filename)

        if os.path.exists(cache_path):
            self.logger.info(f"    > [Cache Hit] Found: {cache_path}")
            with open(cache_path, 'r') as f:
                relaxed_pdb = f.read()
        else:
            self.logger.info("    > [Cache Miss] Running FastRelax...")
            relax_result = movers.run_relax_task(self.input_pdb)

            if relax_result['status'] == 'error':
                self.logger.critical(f"!!! Relax Error: {relax_result['error']}")
                if relax_result.get('trace'):
                    self.logger.critical(relax_result['trace'])
                return None, None

            relaxed_pdb = relax_result['pdb_data']
            with open(cache_path, 'w') as f:
                f.write(relaxed_pdb)
            self.logger.info(f"    > [Cache Saved] {cache_path}")

        elapsed = time.time() - t_start
        self.stage_times['Relaxation'] = elapsed
        self.logger.info(f"    [Time] Relax: {elapsed:.1f} sec")
        return relaxed_pdb, cache_path

    def _step2_global_docking(self, relaxed_pdb: str) -> Optional[List[Dict[str, Any]]]:
        """Step 2: Global docking. Returns list of docking result dicts or None."""
        # Log constraint configuration
        if self.excluded_residues_A:
            self.logger.info(
                f"    [Constraints] Excluded Chain A residues: "
                f"{len(self.excluded_residues_A)} residues, "
                f"early_rejection={self.enable_early_rejection}, "
                f"contact_dist={self.exclusion_contact_dist}A")
        if self.key_residues_B:
            self.logger.info(
                f"    [Constraints] Key Chain B residues: "
                f"{len(self.key_residues_B)} residues, "
                f"bonus_weight={self.key_residue_bonus_weight}")

        self.logger.info("=" * 60)
        self.logger.info(f">>> [Step 2/7] Global Docking ({self.total_global:,} models)")
        self.logger.info("=" * 60)
        t_start = time.time()

        # Build task args with constraint info for early rejection
        if self.enable_early_rejection and self.excluded_residues_A:
            tasks = [(i, relaxed_pdb, self.excluded_residues_A,
                      self.exclusion_contact_dist, self.max_excluded_contacts)
                     for i in range(1, self.total_global + 1)]
        else:
            tasks = [(i, relaxed_pdb) for i in range(1, self.total_global + 1)]

        docking_results = []
        count_done = 0
        count_rejected = 0
        count_errors = 0

        # Adaptive progress interval: ~10 updates total, at least every 500, at most every 5000
        _prog_interval = max(500, min(5000, self.total_global // 10))

        with multiprocessing.Pool(self.n_cpus) as pool:
            for res in pool.imap_unordered(movers.run_global_docking_task, tasks, chunksize=self.chunksize):
                count_done += 1
                if count_done % _prog_interval == 0 or count_done == self.total_global:
                    pct = count_done / self.total_global * 100
                    self.logger.info(
                        f"    > [Docking] {count_done:,}/{self.total_global:,} ({pct:.0f}%) "
                        f"accepted={len(docking_results):,}, rejected={count_rejected:,}, errors={count_errors}")
                if res['status'] == 'success':
                    docking_results.append(res)
                elif res['status'] == 'rejected':
                    count_rejected += 1
                else:
                    count_errors += 1
                    if count_done < 5 or count_done % 1000 == 0:
                        self.logger.warning(f"!!! Docking Task Error: {res.get('error')}")

        elapsed = time.time() - t_start
        self.stage_times['Global Docking'] = elapsed
        self.logger.info(f"    [Time] Global Docking: {elapsed:.1f} sec")
        if count_rejected > 0:
            rej_pct = count_rejected / self.total_global * 100
            self.logger.info(
                f"    > [Early Rejection] {count_rejected}/{self.total_global} "
                f"({rej_pct:.1f}%) rejected (excluded zone contact)")

        if not docking_results:
            self.logger.critical("!!! CRITICAL: 0 models generated.")
            return None

        self.stats['docking_success'] = len(docking_results)
        self.stats['docking_rejected'] = count_rejected
        self.stats['docking_errors'] = count_errors
        return docking_results

    def _step2_5_scoring_filtering(self, docking_results: List[Dict], cache_path: str) -> Tuple[Optional[List[Dict]], Dict]:
        """Step 2.5: Fast scoring & multi-criteria filtering. Returns (cluster_candidates, constraints_dict) or None."""
        self.logger.info("=" * 60)
        self.logger.info(f">>> [Step 3/7] Fast Scoring & Filtering ({len(docking_results):,} models)")
        self.logger.info("=" * 60)
        t_start = time.time()

        constraints_dict = {
            'excluded_residues_A': self.excluded_residues_A,
            'key_residues_B': self.key_residues_B,
            'exclusion_contact_dist': self.exclusion_contact_dist,
            'key_residue_weights': self.key_residue_weights,
        }
        # v2.0 + MiniRefine: sc는 리파인 후 re-score에서 계산되므로 Pass 1에서 스킵
        if self.filter_v2_enabled and self.stage1_enabled and self.mini_refine_enabled:
            constraints_dict['compute_sc'] = False

        if constraints_dict.get('compute_sc', True):
            metric_list = "dG + dSASA + sc_value + total_score"
        else:
            metric_list = "dG + dSASA + total_score (sc deferred to post-refinement)"
        self.logger.info(f"    > Fast scoring: {metric_list} + Center + constraint metrics")

        fast_tasks = ((d['pdb_data'], cache_path, self.contact_distance, constraints_dict)
                      for d in docking_results)

        combined_data = []
        n_dock = len(docking_results)
        _score_prog = max(100, min(2000, n_dock // 10))

        with multiprocessing.Pool(self.n_cpus) as pool:
            score_iter = pool.imap(scoring.run_fast_scoring_task,
                                   fast_tasks, chunksize=self.chunksize)

            for i, score_res in enumerate(score_iter):
                origin = docking_results[i]
                if score_res['status'] == 'success':
                    score_res['id'] = origin['id']
                    score_res['pdb_data'] = origin['pdb_data']
                    combined_data.append(score_res)
                else:
                    self.logger.warning(
                        f"Scoring Error (ID: {origin['id']}): {score_res.get('error')}")
                docking_results[i] = None  # Free memory incrementally
                done = i + 1
                if done % _score_prog == 0 or done == n_dock:
                    pct = done / n_dock * 100
                    self.logger.info(
                        f"    > [Fast Scoring] {done:,}/{n_dock:,} ({pct:.0f}%) "
                        f"passed={len(combined_data):,}")

        del docking_results
        gc.collect()

        if not combined_data:
            self.logger.critical("!!! CRITICAL: All scoring tasks failed.")
            return None, constraints_dict

        total_scored = len(combined_data)
        scores_np = np.array([d['dG_separated'] for d in combined_data])

        # --- Statistics ---
        mean_dG = float(np.mean(scores_np))
        median_dG = float(np.median(scores_np))
        std_dG = float(np.std(scores_np))
        min_dG = float(np.min(scores_np))
        max_dG = float(np.max(scores_np))

        self.logger.info(f"    > [Stats] Scored: {total_scored}/{self.total_global}")
        self.logger.info(f"    > [Stats] dG  Mean={mean_dG:.2f}  Std={std_dG:.2f}  Min={min_dG:.2f}  Max={max_dG:.2f}")

        # Save initial all-models CSV (filter_status filled later)
        self._save_scored_all_models_csv(combined_data)

        if self.filter_v2_enabled and self.stage1_enabled:
            cluster_candidates = self._filter_v2(
                combined_data, total_scored, mean_dG, median_dG, std_dG, min_dG, max_dG, scores_np, t_start)
        else:
            cluster_candidates = self._filter_v1(
                combined_data, total_scored, mean_dG, median_dG, std_dG, min_dG, max_dG, scores_np, t_start)

        if cluster_candidates is None:
            return None, constraints_dict

        # --- Save filtered models ---
        if self.save_filter_max > 0:
            count_save = min(len(cluster_candidates), self.save_filter_max)
            self.logger.info(f"    > [Output] Saving Top {count_save} Filtered Models...")
            for i, item in enumerate(cluster_candidates[:count_save]):
                fname = f"F{i+1:04d}_S{item['dG_separated']:.2f}.pdb"
                with open(os.path.join(self.dir_filter, fname), "w") as f:
                    f.write(item['pdb_data'])

        self.stage_times['Scoring & Filtering'] = time.time() - t_start
        return cluster_candidates, constraints_dict

    def _mini_refinement(self, stage1_pool: List[Dict]) -> List[Dict]:
        """Lightweight refinement of Stage 1 survivors to improve interface packing."""
        count = len(stage1_pool)
        n_rounds = self.mini_refine_n_rounds
        self.logger.info(
            f"    > [Mini Refinement] Refining {count} Stage 1 survivors "
            f"(mode={self.mini_refine_mode}, n_rounds={n_rounds})...")

        t_refine = time.time()

        # --- Pre-refinement metrics snapshot ---
        sc_deferred = not any(d.get('sc_value', 0) > 0 for d in stage1_pool)
        pre_dG = [d.get('dG_separated', 0) for d in stage1_pool]
        pre_dSASA = [d.get('dSASA', 0) for d in stage1_pool]
        pre_sc = [d.get('sc_value', 0) for d in stage1_pool]

        # --- Step 1: Repack/minimize interface ---
        refine_tasks = (
            (d['id'], d['pdb_data'], self.mini_refine_mode, n_rounds)
            for d in stage1_pool
        )

        refined_pool = []
        n_refined = 0
        with multiprocessing.Pool(self.n_cpus) as pool:
            refine_iter = pool.imap(
                movers.run_mini_refinement_task,
                refine_tasks, chunksize=self.chunksize)

            for i, result in enumerate(refine_iter):
                d = stage1_pool[i]
                if result['status'] == 'success':
                    d['pdb_data'] = result['pdb_data']
                    if result.get('refined', False):
                        n_refined += 1
                else:
                    self.logger.warning(
                        f"[MiniRefine] ID {d.get('id')} error: "
                        f"{result.get('error')}")
                refined_pool.append(d)

        self.logger.info(
            f"    > [Mini Refinement] {n_refined}/{count} successfully refined "
            f"({time.time() - t_refine:.1f} sec)")

        # --- Step 2: Re-score with fast scoring to update metrics ---
        self.logger.info(
            f"    > [Mini Refinement] Re-scoring {len(refined_pool)} refined models...")

        t_rescore = time.time()
        constraints_dict = {}
        if self.excluded_residues_A:
            constraints_dict['excluded_residues_A'] = self.excluded_residues_A
        if self.key_residues_B:
            constraints_dict['key_residues_B'] = self.key_residues_B
        if self.exclusion_contact_dist:
            constraints_dict['exclusion_contact_dist'] = self.exclusion_contact_dist
        if self.key_residue_weights:
            constraints_dict['key_residue_weights'] = self.key_residue_weights

        score_tasks = (
            (d['pdb_data'], None, self.contact_distance, constraints_dict)
            for d in refined_pool
        )

        with multiprocessing.Pool(self.n_cpus) as pool:
            score_iter = pool.imap(
                scoring.run_fast_scoring_task,
                score_tasks, chunksize=self.chunksize)

            for i, score_res in enumerate(score_iter):
                d = refined_pool[i]
                if score_res['status'] == 'success':
                    d['dG_separated'] = score_res['dG_separated']
                    d['dSASA'] = score_res['dSASA']
                    d['sc_value'] = score_res['sc_value']
                    d['total_score'] = score_res['total_score']
                    d['center_x'] = score_res['center_x']
                    d['center_y'] = score_res['center_y']
                    d['center_z'] = score_res['center_z']
                    d['dSASA_polar'] = score_res.get('dSASA_polar', 0.0)
                    d['dSASA_hydrophobic'] = score_res.get('dSASA_hydrophobic', 0.0)
                    d['excluded_contacts'] = score_res.get('excluded_contacts', 0)
                    d['key_contact_ratio'] = score_res.get('key_contact_ratio', 0.0)
                    # Recalculate adjusted_dG if key_residue_bonus is active
                    if self.key_residue_bonus_weight > 0:
                        d['adjusted_dG'] = (d['dG_separated']
                                            - self.key_residue_bonus_weight
                                            * d['key_contact_ratio'])
                    # Remove stale dG_density (will be recalculated in Stage 2)
                    d.pop('dG_density', None)
                else:
                    self.logger.warning(
                        f"[MiniRefine Re-score] ID {d.get('id')} error: "
                        f"{score_res.get('error')}")

        self.logger.info(
            f"    > [Mini Refinement] Re-scoring done ({time.time() - t_rescore:.1f} sec)")

        # --- Post-refinement metrics comparison ---
        t_total = time.time() - t_refine
        post_dG = [d.get('dG_separated', 0) for d in refined_pool]
        post_dSASA = [d.get('dSASA', 0) for d in refined_pool]
        post_sc = [d.get('sc_value', 0) for d in refined_pool]

        pre_dG_mean = float(np.mean(pre_dG)) if pre_dG else 0
        pre_dSASA_mean = float(np.mean(pre_dSASA)) if pre_dSASA else 0
        pre_sc_mean = float(np.mean(pre_sc)) if pre_sc else 0
        post_dG_mean = float(np.mean(post_dG)) if post_dG else 0
        post_dSASA_mean = float(np.mean(post_dSASA)) if post_dSASA else 0
        post_sc_mean = float(np.mean(post_sc)) if post_sc else 0

        delta_dG = post_dG_mean - pre_dG_mean
        delta_dSASA = post_dSASA_mean - pre_dSASA_mean
        delta_sc = post_sc_mean - pre_sc_mean

        self.logger.info(
            f"    > [Mini Refinement] Delta: "
            f"dG={delta_dG:+.2f}, dSASA={delta_dSASA:+.1f}, sc={delta_sc:+.3f}")

        self.stats.update({
            'mini_refine_enabled': True,
            'mini_refine_mode': self.mini_refine_mode,
            'mini_refine_n_rounds': n_rounds,
            'mini_refine_count': count,
            'mini_refine_success': n_refined,
            'mini_refine_time': t_total,
            'mini_refine_sc_deferred': sc_deferred,
            'mini_refine_pre_dG_mean': pre_dG_mean,
            'mini_refine_pre_dSASA_mean': pre_dSASA_mean,
            'mini_refine_pre_sc_mean': pre_sc_mean,
            'mini_refine_post_dG_mean': post_dG_mean,
            'mini_refine_post_dSASA_mean': post_dSASA_mean,
            'mini_refine_post_sc_mean': post_sc_mean,
            'mini_refine_delta_dG': delta_dG,
            'mini_refine_delta_dSASA': delta_dSASA,
            'mini_refine_delta_sc': delta_sc,
            'mini_refine_post_sc_median': float(np.median(post_sc)) if post_sc else 0,
            'mini_refine_post_sc_q25': float(np.percentile(post_sc, 25)) if post_sc else 0,
            'mini_refine_post_sc_q75': float(np.percentile(post_sc, 75)) if post_sc else 0,
        })

        return refined_pool

    def _save_filter_thresholds_csv(self, threshold_records: List[Dict]) -> None:
        """Save filter threshold records to a machine-readable CSV file.

        Each record has: filter_name, threshold_value, source, stage,
        n_input, n_passed, n_rejected.  Also stores the records in
        self.stats['filter_thresholds'] for the validation report.
        """
        self.stats['filter_thresholds'] = threshold_records

        csv_path = os.path.join(self.root_dir, "filter_thresholds.csv")
        fieldnames = [
            'filter_name', 'threshold_value', 'source', 'stage',
            'n_input', 'n_passed', 'n_rejected',
        ]
        try:
            with open(csv_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for rec in threshold_records:
                    writer.writerow(rec)
            self.logger.info(
                f"    > [Output] Filter thresholds saved to {csv_path} "
                f"({len(threshold_records)} rows)")
        except Exception as e:
            self.logger.warning(
                f"Failed to write filter_thresholds.csv: {e}")

    def _filter_v2(self, combined_data: List[Dict], total_scored: int, mean_dG: float, median_dG: float, std_dG: float, min_dG: float, max_dG: float, scores_np: Any, t_start: float) -> List[Dict]:
        """v2.0 multi-stage filtering pathway. Returns cluster_candidates list."""
        self.logger.info("    > [Filter] Using v2.0 multi-stage filtering")

        # Threshold tracking records
        _thresh_records = []

        # Filter status tracking: id -> rejection reason
        _filter_status_map: Dict[str, str] = {}

        # ---- Stage 1: Coarse Energy Filtering ---- #

        # 1a. Excluded residue hard filter (preserved from v1.0)
        if self.excluded_residues_A:
            pre_excl = len(combined_data)
            _excl_pass = []
            for d in combined_data:
                if d.get('excluded_contacts', 0) <= self.max_excluded_contacts:
                    _excl_pass.append(d)
                else:
                    _filter_status_map[str(d.get('id', ''))] = 'reject_exclusion'
            combined_data = _excl_pass
            excl_rejected = pre_excl - len(combined_data)
            self.logger.info(
                f"    > [S1-Excl] Excluded zone contacts <= {self.max_excluded_contacts}  "
                f"Passed={len(combined_data)}/{pre_excl}  "
                f"(rejected {excl_rejected})")
            self.stats['excl_filter_rejected'] = excl_rejected
            _thresh_records.append({
                'filter_name': 'stage1_excluded_zone_max_contacts',
                'threshold_value': self.max_excluded_contacts,
                'source': 'config',
                'stage': 'stage1',
                'n_input': pre_excl,
                'n_passed': len(combined_data),
                'n_rejected': excl_rejected,
            })
        else:
            self.stats['excl_filter_rejected'] = 0

        # 1b. dG > max_dG_separated removal (repulsive interfaces)
        pre_dG = len(combined_data)
        stage1_pool = []
        for d in combined_data:
            if d.get('dG_separated', 0) <= self.stage1_max_dG:
                stage1_pool.append(d)
            else:
                _filter_status_map[str(d.get('id', ''))] = 'reject_dG'
        self.logger.info(
            f"    > [S1-dG] dG <= {self.stage1_max_dG}  "
            f"Passed={len(stage1_pool)}/{pre_dG}")
        _thresh_records.append({
            'filter_name': 'stage1_max_dG_separated',
            'threshold_value': self.stage1_max_dG,
            'source': 'config',
            'stage': 'stage1',
            'n_input': pre_dG,
            'n_passed': len(stage1_pool),
            'n_rejected': pre_dG - len(stage1_pool),
        })

        # 1c. total_score percentile filter (top N%)
        ts_cutoff = None
        if stage1_pool and self.stage1_total_score_pct < 100:
            ts_arr = np.array([d.get('total_score', 0) for d in stage1_pool])
            ts_cutoff = float(np.percentile(ts_arr, self.stage1_total_score_pct))
            pre_ts = len(stage1_pool)
            _ts_pass = []
            for d in stage1_pool:
                if d.get('total_score', 0) <= ts_cutoff:
                    _ts_pass.append(d)
                else:
                    _filter_status_map[str(d.get('id', ''))] = 'reject_total_score'
            stage1_pool = _ts_pass
            self.logger.info(
                f"    > [S1-TS] total_score {self.stage1_total_score_pct}%  "
                f"Cutoff={ts_cutoff:.2f}  Passed={len(stage1_pool)}/{pre_ts}")
            _thresh_records.append({
                'filter_name': 'stage1_total_score_percentile',
                'threshold_value': ts_cutoff,
                'source': 'computed_percentile',
                'stage': 'stage1',
                'n_input': pre_ts,
                'n_passed': len(stage1_pool),
                'n_rejected': pre_ts - len(stage1_pool),
            })

        stage1_count = len(stage1_pool)
        self.logger.info(
            f"    > [Stage 1 Result] {stage1_count} survivors "
            f"({stage1_count/total_scored*100:.1f}%)")

        # ---- Mini Refinement (between Stage 1 and Stage 2) ---- #
        if self.mini_refine_enabled and stage1_count > 0:
            stage1_pool = self._mini_refinement(stage1_pool)
        else:
            self.stats['mini_refine_enabled'] = False

        # ---- Funnel Data (L_RMSD computation for energy landscape plot) ----
        if stage1_pool:
            self._compute_funnel_data(stage1_pool)

        # ---- Stage 2: Interface Quality Filtering ---- #

        # 2a. Cheap metrics (already computed in Pass 1)
        stage2_cheap = []
        _s2_rej_dSASA = 0
        _s2_rej_dG_density = 0
        _s2_rej_sc = 0
        for d in stage1_pool:
            dSASA = d.get('dSASA', 0)
            dG = d.get('dG_separated', 0)
            sc = d.get('sc_value', 0)
            _did = str(d.get('id', ''))

            if dSASA < self.stage2_min_dSASA:
                _s2_rej_dSASA += 1
                _filter_status_map[_did] = 'reject_stage2'
                continue
            if abs(dSASA) > 1e-6:
                dG_density = dG / dSASA * 100
            else:
                dG_density = 0.0
            d['dG_density'] = dG_density
            if dG_density > self.stage2_max_dG_density:
                _s2_rej_dG_density += 1
                _filter_status_map[_did] = 'reject_stage2'
                continue
            if sc < self.stage2_min_sc:
                _s2_rej_sc += 1
                _filter_status_map[_did] = 'reject_stage2'
                continue
            stage2_cheap.append(d)

        self.logger.info(
            f"    > [S2-Cheap] dSASA>={self.stage2_min_dSASA} & "
            f"dG_density<={self.stage2_max_dG_density} & "
            f"sc>={self.stage2_min_sc}  "
            f"Passed={len(stage2_cheap)}/{stage1_count}")

        # Track per-filter counts for cheap metrics (sequential cascade)
        _s2c_input = stage1_count
        _thresh_records.append({
            'filter_name': 'stage2_min_dSASA',
            'threshold_value': self.stage2_min_dSASA,
            'source': 'config',
            'stage': 'stage2_cheap',
            'n_input': _s2c_input,
            'n_passed': _s2c_input - _s2_rej_dSASA,
            'n_rejected': _s2_rej_dSASA,
        })
        _thresh_records.append({
            'filter_name': 'stage2_max_dG_density',
            'threshold_value': self.stage2_max_dG_density,
            'source': 'config',
            'stage': 'stage2_cheap',
            'n_input': _s2c_input - _s2_rej_dSASA,
            'n_passed': _s2c_input - _s2_rej_dSASA - _s2_rej_dG_density,
            'n_rejected': _s2_rej_dG_density,
        })
        _thresh_records.append({
            'filter_name': 'stage2_min_sc_value',
            'threshold_value': self.stage2_min_sc,
            'source': 'config',
            'stage': 'stage2_cheap',
            'n_input': _s2c_input - _s2_rej_dSASA - _s2_rej_dG_density,
            'n_passed': len(stage2_cheap),
            'n_rejected': _s2_rej_sc,
        })

        # 2b. Expensive metrics (Pass 2 - only on cheap survivors)
        stage2_full = stage2_cheap
        _s2e_rej_packstat = 0
        _s2e_rej_unsat = 0
        _s2e_rej_nres = 0
        _s2e_rej_hbonds = 0
        _s2e_errors = 0
        _stage2_all_scored: List[Dict] = []  # all models with expensive metrics (pass+fail)
        if self.stage2_expensive and stage2_cheap:
            self.logger.info(
                f"    > [S2-Pass2] Computing expensive metrics for "
                f"{len(stage2_cheap)} survivors...")

            inter_tasks = ((d['pdb_data'],) for d in stage2_cheap)

            with multiprocessing.Pool(self.n_cpus) as pool:
                inter_iter = pool.imap(
                    scoring.run_intermediate_scoring_task,
                    inter_tasks, chunksize=self.chunksize)

                stage2_full = []
                for i, inter_res in enumerate(inter_iter):
                    d = stage2_cheap[i]
                    _did = str(d.get('id', ''))
                    if inter_res['status'] != 'success':
                        self.logger.warning(
                            f"Intermediate scoring error (ID: {d.get('id')}): "
                            f"{inter_res.get('error')}")
                        _s2e_errors += 1
                        _filter_status_map[_did] = 'reject_stage2'
                        continue

                    d['packstat'] = inter_res.get('packstat', 0)
                    d['delta_unsatHbonds'] = inter_res.get('delta_unsatHbonds', 99)
                    d['nres_int'] = inter_res.get('nres_int', 0)
                    d['hbonds_int'] = inter_res.get('hbonds_int', 0)

                    _s2_passed = True
                    if (self.stage2_min_packstat > 0 and
                            d['packstat'] < self.stage2_min_packstat):
                        _s2e_rej_packstat += 1
                        _s2_passed = False
                    elif (self.stage2_max_unsat < 99 and
                            d['delta_unsatHbonds'] > self.stage2_max_unsat):
                        _s2e_rej_unsat += 1
                        _s2_passed = False
                    elif (self.stage2_min_nres > 0 and
                            d['nres_int'] < self.stage2_min_nres):
                        _s2e_rej_nres += 1
                        _s2_passed = False
                    elif (self.stage2_min_hbonds > 0 and
                            d['hbonds_int'] < self.stage2_min_hbonds):
                        _s2e_rej_hbonds += 1
                        _s2_passed = False

                    if _s2_passed:
                        d['_filter_status'] = 'pass'
                        stage2_full.append(d)
                    else:
                        d['_filter_status'] = 'reject_stage2'
                        _filter_status_map[_did] = 'reject_stage2'

                    _stage2_all_scored.append(d)

            # Save Stage 2 CSV (all models with expensive metrics)
            self._save_scored_stage2_models_csv(_stage2_all_scored)

            self.logger.info(
                f"    > [S2-Expensive] packstat>={self.stage2_min_packstat} & "
                f"unsatHb<={self.stage2_max_unsat} & "
                f"nres>={self.stage2_min_nres} & "
                f"hbonds>={self.stage2_min_hbonds}  "
                f"Passed={len(stage2_full)}/{len(stage2_cheap)}")

            # Track expensive filter thresholds (sequential cascade)
            _s2e_input = len(stage2_cheap)
            _s2e_scored = _s2e_input - _s2e_errors  # successfully scored
            _thresh_records.append({
                'filter_name': 'stage2_min_packstat',
                'threshold_value': self.stage2_min_packstat,
                'source': 'config',
                'stage': 'stage2_expensive',
                'n_input': _s2e_scored,
                'n_passed': _s2e_scored - _s2e_rej_packstat,
                'n_rejected': _s2e_rej_packstat,
            })
            _thresh_records.append({
                'filter_name': 'stage2_max_delta_unsatHbonds',
                'threshold_value': self.stage2_max_unsat,
                'source': 'config',
                'stage': 'stage2_expensive',
                'n_input': _s2e_scored - _s2e_rej_packstat,
                'n_passed': _s2e_scored - _s2e_rej_packstat - _s2e_rej_unsat,
                'n_rejected': _s2e_rej_unsat,
            })
            _thresh_records.append({
                'filter_name': 'stage2_min_nres_int',
                'threshold_value': self.stage2_min_nres,
                'source': 'config',
                'stage': 'stage2_expensive',
                'n_input': _s2e_scored - _s2e_rej_packstat - _s2e_rej_unsat,
                'n_passed': _s2e_scored - _s2e_rej_packstat - _s2e_rej_unsat - _s2e_rej_nres,
                'n_rejected': _s2e_rej_nres,
            })
            _thresh_records.append({
                'filter_name': 'stage2_min_hbonds_int',
                'threshold_value': self.stage2_min_hbonds,
                'source': 'config',
                'stage': 'stage2_expensive',
                'n_input': _s2e_scored - _s2e_rej_packstat - _s2e_rej_unsat - _s2e_rej_nres,
                'n_passed': len(stage2_full),
                'n_rejected': _s2e_rej_hbonds,
            })

        stage2_count = len(stage2_full)

        # ---- Graduated Fallback v2.0 ---- #
        fallback_level = 0
        min_surv = self.stage2_min_survivors

        if stage2_count >= min_surv:
            cluster_candidates = stage2_full
        elif len(stage2_cheap) >= min_surv:
            fallback_level = 1
            cluster_candidates = stage2_cheap
            self.logger.warning(
                f"    > [Fallback-1] Expensive metrics relaxed -> "
                f"{len(stage2_cheap)} models")
        elif stage1_count >= min_surv:
            fallback_level = 2
            softer_dSASA = self.stage2_min_dSASA * 0.5
            cluster_candidates = [d for d in stage1_pool
                                  if d.get('dSASA', 0) >= softer_dSASA]
            if len(cluster_candidates) < min_surv:
                cluster_candidates = stage1_pool
            self.logger.warning(
                f"    > [Fallback-2] Stage 2 relaxed (dSASA>={softer_dSASA:.0f}) -> "
                f"{len(cluster_candidates)} models")
        else:
            fallback_level = 3
            combined_data.sort(key=lambda x: x['dG_separated'])
            cluster_candidates = combined_data[:min_surv]
            self.logger.warning(
                f"    > [Fallback-3] All Stage 2 off -> top {min_surv} by dG")

        # Compute dG_density for candidates that don't have it yet
        for d in cluster_candidates:
            if 'dG_density' not in d:
                dSASA = d.get('dSASA', 0)
                if abs(dSASA) > 1e-6:
                    d['dG_density'] = d.get('dG_separated', 0) / dSASA * 100
                else:
                    d['dG_density'] = 0.0

        cluster_candidates.sort(key=lambda x: x['dG_separated'])

        pass_rate = len(cluster_candidates) / total_scored * 100
        self.logger.info(
            f"    > [Filter Result] {len(cluster_candidates)} candidates "
            f"({pass_rate:.1f}% of {total_scored})")
        self.logger.info(
            f"    [Time] Scoring & Filtering: {time.time() - t_start:.1f} sec")

        self.stats.update({
            'total_scored': total_scored,
            'mean_dG': mean_dG, 'median_dG': median_dG, 'std_dG': std_dG,
            'min_dG': min_dG, 'max_dG': max_dG,
            'dG_cutoff': None,
            'dG_passed': stage1_count,
            'dSASA_passed': len([d for d in stage1_pool
                                 if d.get('dSASA', 0) >= self.stage2_min_dSASA]),
            'fallback_used': fallback_level > 0,
            'fallback_level': fallback_level,
            'final_candidates': len(cluster_candidates),
            'pop_dSASA_mean': float(np.mean(
                [d.get('dSASA', 0) for d in combined_data])),
            'filter_version': 'v2.0',
            'stage1_survivors': stage1_count,
            'stage2_survivors': stage2_count,
        })

        # Record fallback information
        if fallback_level > 0:
            _thresh_records.append({
                'filter_name': f'graduated_fallback_level_{fallback_level}',
                'threshold_value': fallback_level,
                'source': 'fallback',
                'stage': 'fallback',
                'n_input': stage2_count if fallback_level == 1 else stage1_count,
                'n_passed': len(cluster_candidates),
                'n_rejected': 0,
            })

        # Save filter thresholds to CSV
        self._save_filter_thresholds_csv(_thresh_records)

        # Mark final candidates as "pass" in the filter status map
        # (overrides any prior rejection if rescued via fallback)
        for d in cluster_candidates:
            _did = str(d.get('id', ''))
            _filter_status_map[_did] = 'pass'

        # Update scored_all_models.csv with filter_status column
        self._update_scored_all_models_filter_status(_filter_status_map)

        del combined_data, stage1_pool, stage2_cheap, stage2_full, scores_np
        del _stage2_all_scored, _filter_status_map
        gc.collect()
        return cluster_candidates

    def _filter_v1(self, combined_data: List[Dict], total_scored: int, mean_dG: float, median_dG: float, std_dG: float, min_dG: float, max_dG: float, scores_np: Any, t_start: float) -> List[Dict]:
        """v1.0 filtering pathway (legacy). Returns cluster_candidates list."""
        if self.filter_v2_enabled and not self.stage1_enabled:
            self.logger.info(
                "    > [Filter] v2.0 sections present but Stage1 disabled, "
                "using v1.0 path")
        else:
            self.logger.info("    > [Filter] Using v1.0 filtering")

        # Filter status tracking: id -> rejection reason
        _filter_status_map: Dict[str, str] = {}

        # Threshold tracking records
        _thresh_records = []

        # --- Filter 1: dSASA (physical contact — applied first) ---
        dSASA_passed = []
        for d in combined_data:
            if d.get('dSASA', 0) >= self.min_dSASA:
                dSASA_passed.append(d)
            else:
                _filter_status_map[str(d.get('id', ''))] = 'reject_dSASA'
        self.logger.info(
            f"    > [Filter-1] dSASA >= {self.min_dSASA} A^2  "
            f"Passed={len(dSASA_passed)}/{total_scored}")
        _thresh_records.append({
            'filter_name': 'v1_min_dSASA',
            'threshold_value': self.min_dSASA,
            'source': 'config',
            'stage': 'filter1_dSASA',
            'n_input': total_scored,
            'n_passed': len(dSASA_passed),
            'n_rejected': total_scored - len(dSASA_passed),
        })

        # --- Filter 1.5: Excluded residue contact (hard rejection) ---
        if self.excluded_residues_A:
            pre_excl = len(dSASA_passed)
            excl_passed = []
            for d in dSASA_passed:
                if d.get('excluded_contacts', 0) <= self.max_excluded_contacts:
                    excl_passed.append(d)
                else:
                    _filter_status_map[str(d.get('id', ''))] = 'reject_exclusion'
            excl_rejected = pre_excl - len(excl_passed)
            self.logger.info(
                f"    > [Filter-1.5] Excluded zone contacts <= {self.max_excluded_contacts}  "
                f"Passed={len(excl_passed)}/{pre_excl}  "
                f"(rejected {excl_rejected} forbidden-zone poses)")
            dSASA_passed = excl_passed
            self.stats['excl_filter_rejected'] = excl_rejected
            _thresh_records.append({
                'filter_name': 'v1_excluded_zone_max_contacts',
                'threshold_value': self.max_excluded_contacts,
                'source': 'config',
                'stage': 'filter1.5_excluded',
                'n_input': pre_excl,
                'n_passed': len(excl_passed),
                'n_rejected': excl_rejected,
            })
        else:
            self.stats['excl_filter_rejected'] = 0

        # --- Filter 2: dG percentile (on contact-forming models only) ---
        dG_cutoff = None
        if not dSASA_passed:
            dG_passed = []
            _thresh_records.append({
                'filter_name': 'v1_dG_percentile',
                'threshold_value': '',
                'source': 'skipped_no_input',
                'stage': 'filter2_dG',
                'n_input': 0,
                'n_passed': 0,
                'n_rejected': 0,
            })
        elif std_dG < 1e-6:
            self.logger.warning("    > [Warning] dG std ~= 0. Skipping dG filter.")
            dG_passed = dSASA_passed
            _thresh_records.append({
                'filter_name': 'v1_dG_percentile',
                'threshold_value': '',
                'source': 'skipped_zero_std',
                'stage': 'filter2_dG',
                'n_input': len(dSASA_passed),
                'n_passed': len(dSASA_passed),
                'n_rejected': 0,
            })
        else:
            contact_dGs = np.array([d['dG_separated'] for d in dSASA_passed])
            dG_cutoff = float(np.percentile(contact_dGs, self.filter_percentile))
            dG_passed = []
            for d in dSASA_passed:
                if d['dG_separated'] <= dG_cutoff:
                    dG_passed.append(d)
                else:
                    _filter_status_map[str(d.get('id', ''))] = 'reject_dG'
            self.logger.info(
                f"    > [Filter-2] dG {self.filter_percentile}%  Cutoff={dG_cutoff:.2f}  "
                f"Passed={len(dG_passed)}/{len(dSASA_passed)}")
            _thresh_records.append({
                'filter_name': 'v1_dG_percentile',
                'threshold_value': dG_cutoff,
                'source': 'computed_percentile',
                'stage': 'filter2_dG',
                'n_input': len(dSASA_passed),
                'n_passed': len(dG_passed),
                'n_rejected': len(dSASA_passed) - len(dG_passed),
            })

        # ---- Funnel Data (L_RMSD computation for energy landscape plot) ----
        if dG_passed:
            self._compute_funnel_data(dG_passed)

        # --- Filter 3: sc_value (shape complementarity sanity check) ---
        sc_passed = []
        for d in dG_passed:
            if d.get('sc_value', 0) >= self.min_sc_value:
                sc_passed.append(d)
            else:
                _filter_status_map[str(d.get('id', ''))] = 'reject_sc'
        if len(dG_passed) - len(sc_passed) > 0:
            self.logger.info(
                f"    > [Filter-3] sc >= {self.min_sc_value}  "
                f"Passed={len(sc_passed)}/{len(dG_passed)}")
        _thresh_records.append({
            'filter_name': 'v1_min_sc_value',
            'threshold_value': self.min_sc_value,
            'source': 'config',
            'stage': 'filter3_sc',
            'n_input': len(dG_passed),
            'n_passed': len(sc_passed),
            'n_rejected': len(dG_passed) - len(sc_passed),
        })

        # --- Graduated Fallback ---
        fallback_level = 0
        if len(sc_passed) >= self.min_survivors:
            cluster_candidates = sc_passed
        elif len(dG_passed) >= self.min_survivors:
            fallback_level = 1
            cluster_candidates = dG_passed
            self.logger.warning(
                f"    > [Fallback-1] sc filter relaxed -> {len(dG_passed)} models")
        elif len(dSASA_passed) >= self.min_survivors:
            fallback_level = 2
            dSASA_passed.sort(key=lambda x: x['dG_separated'])
            cluster_candidates = dSASA_passed
            self.logger.warning(
                f"    > [Fallback-2] dG filter relaxed -> {len(dSASA_passed)} models")
        else:
            fallback_level = 3
            combined_data.sort(key=lambda x: x['dG_separated'])
            cluster_candidates = combined_data[:self.min_survivors]
            self.logger.warning(
                f"    > [Fallback-3] All filters relaxed -> top {self.min_survivors} by dG")

        cluster_candidates.sort(key=lambda x: x['dG_separated'])

        # Mark final candidates as "pass" in the filter status map
        for d in cluster_candidates:
            _did = str(d.get('id', ''))
            _filter_status_map[_did] = 'pass'

        pass_rate = len(cluster_candidates) / total_scored * 100
        self.logger.info(
            f"    > [Filter Result] {len(cluster_candidates)} candidates "
            f"({pass_rate:.1f}% of {total_scored})")
        self.logger.info(f"    [Time] Scoring & Filtering: {time.time() - t_start:.1f} sec")

        self.stats.update({
            'total_scored': total_scored,
            'mean_dG': mean_dG, 'median_dG': median_dG, 'std_dG': std_dG,
            'min_dG': min_dG, 'max_dG': max_dG,
            'dG_cutoff': dG_cutoff,
            'dG_passed': len(dG_passed),
            'dSASA_passed': len(dSASA_passed),
            'fallback_used': fallback_level > 0,
            'fallback_level': fallback_level,
            'final_candidates': len(cluster_candidates),
            'pop_dSASA_mean': float(np.mean([d.get('dSASA', 0) for d in combined_data])),
            'filter_version': 'v1.0',
        })

        # Record fallback information
        if fallback_level > 0:
            _thresh_records.append({
                'filter_name': f'graduated_fallback_level_{fallback_level}',
                'threshold_value': fallback_level,
                'source': 'fallback',
                'stage': 'fallback',
                'n_input': len(sc_passed) if fallback_level == 1 else (
                    len(dG_passed) if fallback_level == 2 else total_scored),
                'n_passed': len(cluster_candidates),
                'n_rejected': 0,
            })

        # Save filter thresholds to CSV
        self._save_filter_thresholds_csv(_thresh_records)

        # Update scored_all_models.csv with filter_status column
        self._update_scored_all_models_filter_status(_filter_status_map)

        del combined_data, dSASA_passed, dG_passed, sc_passed, scores_np
        del _filter_status_map
        gc.collect()
        return cluster_candidates

    def _step3_full_scoring(self, cluster_candidates: List[Dict], cache_path: str, constraints_dict: Dict) -> Optional[List[Dict]]:
        """Step 3: Full scoring of filter survivors. Returns fully_scored list or None."""
        n = len(cluster_candidates)
        self.logger.info("=" * 60)
        self.logger.info(f">>> [Step 4/7] Full Scoring ({n:,} filter survivors)")
        self.logger.info("=" * 60)
        t_start = time.time()

        score_tasks = ((d['pdb_data'], cache_path, self.contact_distance, constraints_dict)
                       for d in cluster_candidates)

        fully_scored = []
        n_errors = 0

        with multiprocessing.Pool(self.n_cpus) as pool:
            score_iter = pool.imap(scoring.run_scoring_task,
                                   score_tasks, chunksize=self.chunksize)

            for i, score_res in enumerate(score_iter):
                origin = cluster_candidates[i]
                if score_res['status'] == 'success':
                    score_res['id'] = origin.get('id', i + 1)
                    # Recalculate dG_density from full scoring values
                    dSASA = score_res.get('dSASA', 0)
                    dG = score_res.get('dG_separated', 0)
                    score_res['dG_density'] = (dG / dSASA * 100) if abs(dSASA) > 1e-6 else 0.0
                    fully_scored.append(score_res)
                else:
                    n_errors += 1
                    self.logger.warning(
                        f"Full scoring error (ID: {origin.get('id', i+1)}): "
                        f"{score_res.get('error')} [step={score_res.get('step')}]")

        del cluster_candidates
        gc.collect()

        self.logger.info(
            f"    > [Full Scoring] {len(fully_scored)}/{n} successful "
            f"({n_errors} errors)")
        elapsed = time.time() - t_start
        self.stage_times['Full Scoring'] = elapsed
        self.logger.info(f"    [Time] Full Scoring: {elapsed:.1f} sec")

        self.stats['total_full_scored'] = len(fully_scored)

        if not fully_scored:
            self.logger.critical("!!! CRITICAL: All full scoring tasks failed.")
            return None

        fully_scored.sort(key=lambda x: x['dG_separated'])
        return fully_scored

    def _step4_clustering(self, cluster_candidates: List[Dict]) -> Optional[List[Dict]]:
        """Step 4: L_RMSD greedy clustering. Returns final_representatives or None."""
        n_cands = len(cluster_candidates)
        self.logger.info("=" * 60)
        self.logger.info(f">>> [Step 5/7] L_RMSD Greedy Clustering ({n_cands:,} candidates)")
        self.logger.info("=" * 60)
        t_start = time.time()

        centers = np.array([[c.get('center_x', 0), c.get('center_y', 0),
                             c.get('center_z', 0)] for c in cluster_candidates])
        # Warn about models with unresolved center coordinates (all zeros)
        n_zero_center = int(np.sum(np.all(centers == 0, axis=1)))
        if n_zero_center > 0:
            self.logger.warning(
                f"    [Cluster WARN] {n_zero_center}/{n_cands} models have center (0,0,0) "
                f"— scoring may have failed for these; CoM pre-filter unreliable for them")
        com_threshold = max(self.cluster_threshold * 8, COM_THRESHOLD_MIN)

        leaders = []
        leader_centers = []
        cluster_members = []
        cluster_populations = []
        com_skipped = 0
        potential_extra_clusters = 0
        dropped_records = []
        all_membership_records = []  # Track every model→cluster assignment

        # Adaptive escalation state
        escalation_applied = False
        escalation_threshold = max(int(self.cluster_top_n * 0.2), 3)
        original_cluster_top_n = self.cluster_top_n

        for idx, candidate in enumerate(cluster_candidates):
            if (idx + 1) % 500 == 0:
                self.logger.info(
                    f"    > [Cluster] {idx+1}/{n_cands}  "
                    f"({len(leaders)} clusters, {com_skipped} CoM-skipped)")

            cand_center = centers[idx]

            potential_leaders = []
            for lid, lc in enumerate(leader_centers):
                if np.sqrt(np.sum((cand_center - lc) ** 2)) <= com_threshold:
                    potential_leaders.append(lid)

            if not potential_leaders and len(leaders) >= self.cluster_top_n:
                # Shadow archive: record metadata before dropping
                nearest_dist = float('inf')
                if leader_centers:
                    dists = [np.sqrt(np.sum((cand_center - lc) ** 2)) for lc in leader_centers]
                    nearest_dist = min(dists)
                dropped_records.append({
                    'id': candidate.get('id', idx + 1),
                    'center_x': round(float(cand_center[0]), 2),
                    'center_y': round(float(cand_center[1]), 2),
                    'center_z': round(float(cand_center[2]), 2),
                    'dG_separated': candidate.get('dG_separated', 0),
                    'dSASA': candidate.get('dSASA', 0),
                    'sc_value': candidate.get('sc_value', 0),
                    'L_RMSD_to_nearest_leader': '',
                    'dist_to_nearest_leader': round(nearest_dist, 2),
                    'drop_reason': 'com_skip_cap_full',
                })
                com_skipped += 1
                continue

            cand_pose = pyrosetta_init.string_to_pose(candidate['pdb_data'])
            if cand_pose is None or cand_pose.num_chains() < 2:
                continue

            assigned_cid = None
            min_lrmsd = float('inf')
            for lid in potential_leaders:
                _, leader_pose = leaders[lid]
                l_rmsd = self._compute_l_rmsd(cand_pose, leader_pose)
                if l_rmsd <= self.cluster_threshold and l_rmsd < min_lrmsd:
                    min_lrmsd = l_rmsd
                    assigned_cid = lid

            if assigned_cid is not None:
                cluster_populations[assigned_cid] += 1
                is_rep = False
                members = cluster_members[assigned_cid]
                if len(members) < self.members_per_cluster:
                    is_diverse = True
                    for _, existing_pose in members:
                        if self._compute_l_rmsd(cand_pose, existing_pose) < self.member_diversity_dist:
                            is_diverse = False
                            break
                    if is_diverse:
                        candidate['Cluster_ID'] = assigned_cid + 1
                        candidate['Member_ID'] = len(members) + 1
                        candidate['Parent'] = f"C{candidate['Cluster_ID']:02d}_M{candidate['Member_ID']:02d}"
                        members.append((candidate, cand_pose))
                        is_rep = True
                all_membership_records.append({
                    'model_id': candidate.get('id', idx + 1),
                    'cluster_id': f"C{assigned_cid + 1:02d}",
                    'is_representative': is_rep,
                    'dG_separated': candidate.get('dG_separated', 0),
                    'L_RMSD': candidate.get('L_RMSD', 0),
                    'center_x': round(float(cand_center[0]), 2),
                    'center_y': round(float(cand_center[1]), 2),
                    'center_z': round(float(cand_center[2]), 2),
                })
            elif len(leaders) < self.cluster_top_n:
                cid = len(leaders)
                candidate['Cluster_ID'] = cid + 1
                candidate['Member_ID'] = 1
                candidate['Parent'] = f"C{candidate['Cluster_ID']:02d}_M{candidate['Member_ID']:02d}"
                leaders.append((candidate, cand_pose))
                leader_centers.append(cand_center)
                cluster_members.append([(candidate, cand_pose)])
                cluster_populations.append(1)
                all_membership_records.append({
                    'model_id': candidate.get('id', idx + 1),
                    'cluster_id': f"C{cid + 1:02d}",
                    'is_representative': True,
                    'dG_separated': candidate.get('dG_separated', 0),
                    'L_RMSD': candidate.get('L_RMSD', 0),
                    'center_x': round(float(cand_center[0]), 2),
                    'center_y': round(float(cand_center[1]), 2),
                    'center_z': round(float(cand_center[2]), 2),
                })
            else:
                # Shadow archive: record with L_RMSD (pose already deserialized)
                nearest_dist = float('inf')
                min_lrmsd_any = float('inf')
                if leader_centers:
                    dists = [np.sqrt(np.sum((cand_center - lc) ** 2)) for lc in leader_centers]
                    nearest_dist = min(dists)
                for lid, (_, lp) in enumerate(leaders):
                    lrmsd = self._compute_l_rmsd(cand_pose, lp)
                    if lrmsd < min_lrmsd_any:
                        min_lrmsd_any = lrmsd
                dropped_records.append({
                    'id': candidate.get('id', idx + 1),
                    'center_x': round(float(cand_center[0]), 2),
                    'center_y': round(float(cand_center[1]), 2),
                    'center_z': round(float(cand_center[2]), 2),
                    'dG_separated': candidate.get('dG_separated', 0),
                    'dSASA': candidate.get('dSASA', 0),
                    'sc_value': candidate.get('sc_value', 0),
                    'L_RMSD_to_nearest_leader': round(min_lrmsd_any, 2),
                    'dist_to_nearest_leader': round(nearest_dist, 2),
                    'drop_reason': 'cap_full_no_match',
                })
                potential_extra_clusters += 1

                # Adaptive escalation: expand cluster_top_n once if too many extras
                if not escalation_applied and potential_extra_clusters >= escalation_threshold:
                    new_top_n = int(self.cluster_top_n * 1.5)
                    self.logger.info(
                        f"    > [Escalation] potential_extra={potential_extra_clusters} >= "
                        f"threshold={escalation_threshold}. "
                        f"Expanding cluster_top_n: {self.cluster_top_n} -> {new_top_n}")
                    self.cluster_top_n = new_top_n
                    escalation_applied = True

        final_representatives = []
        for members in cluster_members:
            for data, _ in members:
                final_representatives.append(data)

        extra_msg = ""
        if potential_extra_clusters > 0:
            total_possible = len(leaders) + potential_extra_clusters
            extra_msg = f", 실제 가능: {total_possible}개"
        self.logger.info(
            f"    > Found {len(leaders)} clusters{extra_msg}, "
            f"Total {len(final_representatives)} representatives. "
            f"(CoM skipped {com_skipped})")
        elapsed = time.time() - t_start
        self.stage_times['Clustering'] = elapsed
        self.logger.info(f"    [Time] Clustering: {elapsed:.1f} sec")

        # Save cluster PDBs + cluster summary CSV
        cluster_csv_rows = []
        for item in final_representatives:
            fname = f"C{item['Cluster_ID']:02d}_M{item['Member_ID']:02d}_S{item['dG_separated']:.2f}.pdb"
            with open(os.path.join(self.dir_cluster, fname), "w") as f:
                f.write(item['pdb_data'])
            row = {
                'Cluster_ID': item['Cluster_ID'],
                'Member_ID': item['Member_ID'],
                'dG_separated': item.get('dG_separated', 0),
                'Total_Score': item.get('Total_Score', 0),
                'dSASA': item.get('dSASA', 0),
                'dSASA_polar': item.get('dSASA_polar', 0.0),
                'dSASA_hydrophobic': item.get('dSASA_hydrophobic', 0.0),
                'dG_density': item.get('dG_density', 0),
                'sc_value': item.get('sc_value', 0),
                'packstat': item.get('packstat', 0),
                'delta_unsatHbonds': item.get('delta_unsatHbonds', ''),
                'nres_int': item.get('nres_int', ''),
                'hbonds_int': item.get('hbonds_int', ''),
                'L_RMSD': item.get('L_RMSD', 0),
                'I_RMSD': item.get('I_RMSD', ''),
                'L_RMSD_best': getattr(self, 'lrmsd_best_map', {}).get(
                    item.get('id'), ''),
                'Binding_Residues': item.get('Binding_Residues', ''),
                'Binding_Residues_A': item.get('Binding_Residues_A', ''),
                'Binding_Residues_B': item.get('Binding_Residues_B', ''),
                'Population': cluster_populations[item['Cluster_ID'] - 1],
                'center_x': item.get('center_x', 0),
                'center_y': item.get('center_y', 0),
                'center_z': item.get('center_z', 0),
                'File_PDB': fname,
            }
            if self.key_residues_B:
                row['key_contact_ratio'] = item.get('key_contact_ratio', 0)
            cluster_csv_rows.append(row)

        if cluster_csv_rows:
            cluster_csv_path = os.path.join(self.dir_cluster, "cluster_summary.csv")
            pd.DataFrame(cluster_csv_rows).to_csv(cluster_csv_path, index=False)
            self.logger.info(f"    > [Output] Cluster summary: {cluster_csv_path}")

        # Save full cluster membership CSV (all models → cluster mapping)
        if all_membership_records:
            membership_csv_path = os.path.join(self.dir_cluster, "cluster_membership.csv")
            membership_df = pd.DataFrame(all_membership_records, columns=[
                'model_id', 'cluster_id', 'is_representative',
                'dG_separated', 'L_RMSD',
                'center_x', 'center_y', 'center_z',
            ])
            membership_df.to_csv(membership_csv_path, index=False)
            self.logger.info(
                f"    > [Output] Cluster membership: {membership_csv_path} "
                f"({len(all_membership_records)} models mapped)")

        if not final_representatives:
            self.logger.critical("!!! CRITICAL: 0 cluster representatives.")
            return None

        # Shadow archive: save dropped candidates metadata
        if dropped_records:
            dropped_csv_path = os.path.join(self.dir_cluster, "dropped_candidates.csv")
            pd.DataFrame(dropped_records).to_csv(dropped_csv_path, index=False)
            self.logger.info(
                f"    > [Shadow Archive] {len(dropped_records)} dropped candidates "
                f"saved to {dropped_csv_path}")
            self.stats['dropped_candidates_path'] = dropped_csv_path
            self.stats['n_dropped_candidates'] = len(dropped_records)

        self.stats['num_clusters'] = len(leaders)
        self.stats['cluster_sizes'] = [len(m) for m in cluster_members]
        self.stats['cluster_populations'] = list(cluster_populations)
        self.stats['num_representatives'] = len(final_representatives)
        self.stats['potential_extra_clusters'] = potential_extra_clusters
        self.stats['total_possible_clusters'] = len(leaders) + potential_extra_clusters
        self.stats['cluster_top_n_escalated'] = escalation_applied
        if escalation_applied:
            self.stats['cluster_top_n_original'] = original_cluster_top_n
            self.stats['cluster_top_n_escalated_to'] = self.cluster_top_n

        del leaders, leader_centers, cluster_members, all_membership_records
        gc.collect()
        return final_representatives

    def _step5_selection_and_save(self, final_representatives: List[Dict]) -> Tuple[Any, str]:
        """Step 5: Diversity-aware selection, dedup, and save. Returns (ranking_df, final_csv_path)."""
        self.logger.info("=" * 60)
        self.logger.info(">>> [Step 6/7] Diversity-Aware Selection & Deduplication")
        self.logger.info("=" * 60)
        t_start_step6 = time.time()
        final_scores = final_representatives

        save_count = min(len(final_scores), self.save_top_n)

        # Group final scores by cluster
        cluster_groups = defaultdict(list)
        for item in final_scores:
            cid = item['Parent'].split('_')[0] if '_' in item.get('Parent', '') else item.get('Parent', 'UNK')
            item['_cluster_id'] = cid
            cluster_groups[cid].append(item)

        bonus_w = self.key_residue_bonus_weight
        for cid in cluster_groups:
            cluster_groups[cid].sort(
                key=lambda x: x['dG_separated'] - bonus_w * x.get('key_contact_ratio', 0))

        sorted_cids = sorted(cluster_groups.keys(),
                             key=lambda c: (cluster_groups[c][0]['dG_separated']
                                            - bonus_w * cluster_groups[c][0].get('key_contact_ratio', 0)))

        # Round-robin selection
        pool_target = min(len(final_scores), save_count * SELECTION_POOL_MULTIPLIER)
        selection_pool = []
        cluster_ptrs = {cid: 0 for cid in sorted_cids}

        while len(selection_pool) < pool_target:
            added_this_round = False
            for cid in sorted_cids:
                ptr = cluster_ptrs[cid]
                if ptr < len(cluster_groups[cid]):
                    selection_pool.append(cluster_groups[cid][ptr])
                    cluster_ptrs[cid] += 1
                    added_this_round = True
                    if len(selection_pool) >= pool_target:
                        break
            if not added_this_round:
                break

        n_clusters_in_pool = len(set(item['_cluster_id'] for item in selection_pool))
        self.logger.info(
            f"    > [Round-Robin] {len(selection_pool)} candidates from "
            f"{n_clusters_in_pool} clusters")

        # Post-selection L_RMSD deduplication
        dedup_threshold = self.member_diversity_dist
        deduped = []
        dedup_poses = []
        n_dedup_skipped = 0

        for item in selection_pool:
            pose = pyrosetta_init.string_to_pose(item['pdb_data'])
            if pose is None or pose.num_chains() < 2:
                continue

            is_duplicate = False
            for kept_pose in dedup_poses:
                if self._compute_l_rmsd(pose, kept_pose) < dedup_threshold:
                    is_duplicate = True
                    break

            if is_duplicate:
                n_dedup_skipped += 1
            else:
                deduped.append(item)
                dedup_poses.append(pose)

            if len(deduped) >= save_count:
                break

        # Compute L_RMSD_best for all final models (vs Rank 1 = best dG)
        # This guarantees every model gets a value, regardless of pipeline path
        final_lrmsd_best = []
        if dedup_poses:
            ref_pose = dedup_poses[0]
            for i, pose in enumerate(dedup_poses):
                if i == 0:
                    final_lrmsd_best.append(0.0)
                else:
                    final_lrmsd_best.append(self._compute_l_rmsd(pose, ref_pose))

        del dedup_poses
        gc.collect()

        if n_dedup_skipped > 0:
            self.logger.info(
                f"    > [Dedup] Removed {n_dedup_skipped} near-duplicates "
                f"(L_RMSD < {dedup_threshold}A)")

        if final_lrmsd_best:
            self.logger.info(
                f"    > [L_RMSD_best] Recomputed for {len(final_lrmsd_best)} models "
                f"(ref: Rank 1 = {deduped[0].get('Parent', 'N/A')})")
        self.logger.info(
            f"    > Saving {len(deduped)} diverse final models...")

        # Boltzmann-weighted cluster probabilities
        kT = BOLTZMANN_KT
        cluster_boltz = {}
        all_dGs = [item['dG_separated'] for item in final_scores]
        if all_dGs:
            min_dG_ref = min(all_dGs)
            for cid in sorted_cids:
                group = cluster_groups[cid]
                cluster_boltz[cid] = sum(
                    math.exp(-(item['dG_separated'] - min_dG_ref) / kT)
                    for item in group)
            z_total = sum(cluster_boltz.values())
            if z_total > 0:
                for cid in cluster_boltz:
                    cluster_boltz[cid] /= z_total
        self.stats['boltzmann_probs'] = cluster_boltz
        self.stats['n_dedup_removed'] = n_dedup_skipped

        # Site grouping
        site_groups = self._compute_site_groups(final_scores)
        self.stats['site_groups'] = site_groups
        if site_groups:
            c2g = site_groups['cluster_to_group']
            # Group Boltzmann probabilities
            group_boltz = defaultdict(float)
            for cid, prob in cluster_boltz.items():
                gid = c2g.get(cid, 0)
                group_boltz[gid] += prob
            self.stats['site_group_boltzmann'] = dict(group_boltz)

            # Add Site_Group to cluster_summary.csv
            cluster_csv_path = os.path.join(self.dir_cluster, "cluster_summary.csv")
            if os.path.exists(cluster_csv_path):
                cdf = pd.read_csv(cluster_csv_path)
                cdf['Site_Group'] = cdf['Cluster_ID'].apply(
                    lambda x: c2g.get(f"C{x:02d}", 0))
                cdf.to_csv(cluster_csv_path, index=False)

        # Save selected models
        csv_rows = []

        for rank, item in enumerate(deduped):
            fname_pdb = f"Rank{rank+1:02d}_{item['Parent']}_S{item['dG_separated']:.2f}.pdb"
            fpath_pdb = os.path.join(self.dir_final, fname_pdb)
            with open(fpath_pdb, "w") as f:
                f.write(item['pdb_data'])

            fname_csv = f"Rank{rank+1:02d}_{item['Parent']}_Energies.csv"
            fpath_csv = os.path.join(self.dir_final, fname_csv)
            with open(fpath_csv, "w") as f:
                f.write(item['csv_data'])

            interface_csv = item.get('interface_csv_data', '')
            if interface_csv and len(interface_csv.split('\n')) > 1:
                fname_iface = f"Rank{rank+1:02d}_{item['Parent']}_InterfaceEnergies.csv"
                with open(os.path.join(self.dir_final, fname_iface), "w") as f:
                    f.write(interface_csv)

            try:
                contact_pairs_csv = item.get('contact_pairs_csv_data', '')
                if contact_pairs_csv and len(contact_pairs_csv.split('\n')) > 1:
                    fname_contacts = f"Rank{rank+1:02d}_{item['Parent']}_ContactPairs.csv"
                    with open(os.path.join(self.dir_final, fname_contacts), "w") as f:
                        f.write(contact_pairs_csv)
            except Exception as e:
                self.logger.warning(
                    f"Error writing ContactPairs CSV for Rank {rank+1}: {e}")

            dSASA_val = item.get('dSASA', 0)
            dG_val = item.get('dG_separated', 0)
            dG_density = (dG_val / dSASA_val * 100) if abs(dSASA_val) > 1e-6 else 0.0

            row = {
                'Rank': rank + 1,
                'Parent': item.get('Parent', ''),
                'dG_separated': dG_val,
                'Total_Score': item.get('Total_Score', 0),
                'dSASA': dSASA_val,
                'dSASA_polar': item.get('dSASA_polar', 0.0),
                'dSASA_hydrophobic': item.get('dSASA_hydrophobic', 0.0),
                'dG_density': round(dG_density, 4),
                'sc_value': item.get('sc_value', 0),
                'packstat': item.get('packstat', 0),
                'delta_unsatHbonds': item.get('delta_unsatHbonds', ''),
                'nres_int': item.get('nres_int', ''),
                'hbonds_int': item.get('hbonds_int', ''),
                'L_RMSD': item.get('L_RMSD', 0),
                'I_RMSD': item.get('I_RMSD', ''),
                'L_RMSD_best': final_lrmsd_best[rank] if rank < len(final_lrmsd_best) else '',
                'Binding_Residues': item.get('Binding_Residues', ''),
                'Binding_Residues_A': item.get('Binding_Residues_A', ''),
                'Binding_Residues_B': item.get('Binding_Residues_B', ''),
                'center_x': item.get('center_x', 0),
                'center_y': item.get('center_y', 0),
                'center_z': item.get('center_z', 0),
                'File_PDB': fname_pdb,
                'File_CSV': fname_csv,
            }
            if self.key_residues_B:
                row['key_contact_ratio'] = item.get('key_contact_ratio', 0)
                row['key_contacted'] = item.get('key_contacted', 0)
                row['key_present'] = item.get('key_present', 0)
            if site_groups:
                cid = item.get('_cluster_id', '')
                row['Site_Group'] = site_groups['cluster_to_group'].get(cid, 0)
            csv_rows.append(row)

        ranking_df = pd.DataFrame(csv_rows)
        final_csv_path = os.path.join(self.dir_final, "final_ranking.csv")
        ranking_df.to_csv(final_csv_path, index=False)
        ranking_df.to_csv(os.path.join(self.root_dir, "final_ranking.csv"), index=False)
        self.logger.info(f"    > [Output] Final ranking: {final_csv_path}")

        # Archive: extended metadata CSV (no PDB files) for broader coverage
        if self.archive_top_n > self.save_top_n and len(selection_pool) > len(deduped):
            archive_count = min(len(selection_pool), self.archive_top_n)
            archive_rows = []
            for rank_a, item in enumerate(selection_pool[:archive_count]):
                dSASA_a = item.get('dSASA', 0)
                dG_a = item.get('dG_separated', 0)
                dG_density_a = (dG_a / dSASA_a * 100) if abs(dSASA_a) > 1e-6 else 0.0
                archive_rows.append({
                    'Archive_Rank': rank_a + 1,
                    'Parent': item.get('Parent', ''),
                    'Cluster_ID': item.get('Cluster_ID', ''),
                    'dG_separated': dG_a,
                    'dSASA': dSASA_a,
                    'dSASA_polar': item.get('dSASA_polar', 0.0),
                    'dSASA_hydrophobic': item.get('dSASA_hydrophobic', 0.0),
                    'sc_value': item.get('sc_value', 0),
                    'dG_density': round(dG_density_a, 4),
                    'packstat': item.get('packstat', 0),
                    'delta_unsatHbonds': item.get('delta_unsatHbonds', ''),
                    'nres_int': item.get('nres_int', ''),
                    'hbonds_int': item.get('hbonds_int', ''),
                    'L_RMSD': item.get('L_RMSD', 0),
                    'I_RMSD': item.get('I_RMSD', ''),
                    'center_x': item.get('center_x', 0),
                    'center_y': item.get('center_y', 0),
                    'center_z': item.get('center_z', 0),
                    'Binding_Residues_A': item.get('Binding_Residues_A', ''),
                })
            archive_csv_path = os.path.join(self.dir_final, "archive_ranking.csv")
            pd.DataFrame(archive_rows).to_csv(archive_csv_path, index=False)
            self.logger.info(
                f"    > [Archive] {len(archive_rows)} models in {archive_csv_path}")

        elapsed = time.time() - t_start_step6
        self.stage_times['Selection & Save'] = elapsed
        self.logger.info(f"    [Time] Selection & Save: {elapsed:.1f} sec")
        return ranking_df, final_csv_path

    def _step6_visualization(self, ranking_df: Any, final_csv_path: str, overall_start_time: float) -> None:
        """Step 6: Visualization & validation report."""
        self.logger.info("=" * 60)
        self.logger.info(">>> [Step 7/7] Generating Visualization & Validation Report")
        self.logger.info("=" * 60)
        self.plot_energy_funnel(final_csv_path)
        self.generate_cluster_overview()
        self.generate_per_cluster_views(final_csv_path)

        total_duration = time.time() - overall_start_time
        self.stats['total_runtime'] = str(timedelta(seconds=int(total_duration)))

        self.generate_validation_report(ranking_df)
        self.generate_html_dashboard(ranking_df, final_csv_path)

        self.logger.info("=" * 50)
        self.logger.info(f"    BENCHMARK REPORT")
        self.logger.info(f"    Total Time      : {self.stats['total_runtime']}")
        self.logger.info(">>> PIPELINE COMPLETED.")

    # ================================================================== #
    #  MAIN PIPELINE
    # ================================================================== #
    def execute(self) -> None:
        overall_start_time = time.time()
        try:
            if self.master_seed is not None:
                pyrosetta_init.set_master_seed(self.master_seed)
            pyrosetta_init.init_rosetta()
            self.logger.info(f"=== Pipeline Started (CPU: {self.n_cpus}) ===")
            self.logger.info(f"Target Input: {self.input_pdb}")
            self.logger.info(f"Root Directory: {self.root_dir}")
            if self.master_seed is not None:
                self.logger.info(f"Random Seed: {self.master_seed}")

            # Step 1: Relax
            relaxed_pdb, cache_path = self._step1_relax()
            if relaxed_pdb is None:
                return
            self.cache_path = cache_path

            self._resolve_chain_sizes(relaxed_pdb)

            # Preflight: validate chain labels, numbering, constraints
            if not self._preflight_check(relaxed_pdb):
                self.logger.critical(
                    "!!! Preflight check FAILED. Aborting pipeline. "
                    "Fix input PDB or config before re-running.")
                return

            # Log all filter thresholds for transparency
            self._log_filter_thresholds()

            self._save_input_validation_report()
            self._save_config_snapshot()
            self._save_run_metadata("started")

            # Step 2: Global Docking
            docking_results = self._step2_global_docking(relaxed_pdb)
            if docking_results is None:
                self._save_run_metadata("docking_failed")
                return

            # Step 2.5: Fast Scoring & Filtering
            result = self._step2_5_scoring_filtering(docking_results, cache_path)
            cluster_candidates, constraints_dict = result
            if cluster_candidates is None:
                self._save_run_metadata("filtering_failed")
                return

            # Step 3: Full Scoring of filter survivors
            fully_scored = self._step3_full_scoring(cluster_candidates, cache_path, constraints_dict)
            if fully_scored is None:
                self._save_run_metadata("full_scoring_failed")
                return
            self._export_decoy_scores(fully_scored)
            self._save_run_metadata("full_scoring_completed")

            # Step 4: Clustering
            final_representatives = self._step4_clustering(fully_scored)
            if final_representatives is None:
                self._save_run_metadata("clustering_failed")
                return

            # Step 5: Selection & Save
            ranking_df, final_csv_path = self._step5_selection_and_save(final_representatives)
            self._save_run_metadata("selection_completed")

            # Step 6: Visualization & Report
            self._step6_visualization(ranking_df, final_csv_path, overall_start_time)
            self._save_run_metadata("completed")

            gc.collect()

        except Exception as e:
            self.logger.critical(f"!!! CRITICAL UNHANDLED ERROR: {str(e)}")
            self.logger.critical(traceback.format_exc())
            try:
                self._save_run_metadata("unhandled_error")
            except Exception:
                pass


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python pipeline_manager.py <config_file> [input_pdb]")
        sys.exit(1)
    PipelineManager(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None).execute()
