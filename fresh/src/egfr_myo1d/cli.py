"""Minimal CLI for the fresh EGFR-MYO1D workflow skeleton."""

import argparse
import json
import os
import subprocess
import sys

from egfr_myo1d import WORKFLOW_STAGE, __version__


def _candidate_python_commands():
    commands = []
    env_python = os.environ.get("EGFR_MYO1D_PYTHON")
    if env_python:
        commands.append([env_python])

    path_value = os.environ.get("PATH", "")
    for directory in path_value.split(os.pathsep):
        pytest_exe = os.path.join(directory, "pytest.exe")
        if os.path.exists(pytest_exe):
            candidate = os.path.abspath(os.path.join(directory, os.pardir, "python.exe"))
            if os.path.exists(candidate):
                commands.append([candidate])

    commands.extend([["py", "-3"], ["python3"], ["python"]])
    return commands


def _candidate_supports_task2(command):
    if command[0] == sys.executable:
        return False
    try:
        proc = subprocess.Popen(
            command
            + [
                "-c",
                "import sys; import yaml, numpy, pandas; "
                "raise SystemExit(0 if sys.version_info[:2] >= (3, 9) else 1)",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        proc.communicate()
        return proc.returncode == 0
    except OSError:
        return False


def _bootstrap_python3_if_needed():
    if sys.version_info[:2] >= (3, 9):
        return
    for command in _candidate_python_commands():
        if _candidate_supports_task2(command):
            raise SystemExit(
                subprocess.call(command + ["-m", "egfr_myo1d.cli"] + sys.argv[1:])
            )
    sys.stderr.write(
        "Python 3.9+ with PyYAML, numpy, and pandas is required for Task 2 commands.\n"
    )
    raise SystemExit(1)


_bootstrap_python3_if_needed()


def build_parser():
    parser = argparse.ArgumentParser(
        prog="egfr-myo1d",
        description=(
            "Fresh EGFR-MYO1D workflow foundation. Task 2 provides run "
            "context, manifests, logging, status, and preflight only."
        ),
    )
    subparsers = parser.add_subparsers(dest="command")

    version_parser = subparsers.add_parser(
        "version", help="Print the fresh workflow skeleton version."
    )
    version_parser.set_defaults(func=_cmd_version)

    init_parser = subparsers.add_parser(
        "init-run", help="Create a fresh run directory with baseline manifests and logs."
    )
    init_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input", "mini", "prod", "production"])
    init_parser.add_argument("--run-id", help="Optional run identifier under fresh/runs/.")
    init_parser.set_defaults(func=_cmd_init_run)

    preflight_parser = subparsers.add_parser(
        "preflight", help="Run environment and input preflight checks."
    )
    preflight_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    preflight_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input"])
    preflight_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev warns on missing heavy tools; hpc_strict fails.",
    )
    preflight_parser.set_defaults(func=_cmd_preflight)

    tool_preflight_parser = subparsers.add_parser(
        "tool-preflight",
        help="Discover or smoke-test PPI-surface runtime tools without installing anything.",
    )
    tool_preflight_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    tool_preflight_parser.add_argument(
        "--mode",
        default="discover",
        choices=["discover", "smoke"],
        help="discover records availability; smoke also runs help/version/import probes.",
    )
    tool_preflight_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="hpc_strict turns missing required core tools into blockers.",
    )
    tool_preflight_parser.add_argument(
        "--registry",
        default=None,
        help="Optional path to a tool_registry.yaml override.",
    )
    tool_preflight_parser.add_argument(
        "--timeout",
        type=int,
        default=15,
        help="Per-tool smoke command timeout in seconds.",
    )
    tool_preflight_parser.set_defaults(func=_cmd_tool_preflight)

    status_parser = subparsers.add_parser("status", help="Summarize an existing fresh run.")
    status_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    status_parser.set_defaults(func=_cmd_status)

    structure_parser = subparsers.add_parser(
        "validate-structures",
        help="Validate structural input contracts and write Task 3 QC reports.",
    )
    structure_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    structure_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input"])
    structure_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev warns on missing future inputs; hpc_strict fails.",
    )
    structure_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing contract.json and referenced structure fixtures or inputs.",
    )
    structure_parser.set_defaults(func=_cmd_validate_structures)

    prepare_parser = subparsers.add_parser(
        "prepare-ppi-inputs",
        help="Prepare audited EGFR/MYO1D PPI input packs without docking.",
    )
    prepare_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    prepare_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input"])
    prepare_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev permits explicit fixture warnings; hpc_strict blocks them.",
    )
    prepare_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing ppi_input_contract.json and referenced structures.",
    )
    prepare_parser.add_argument(
        "--contract",
        help="Optional contract path relative to --input-root.",
    )
    prepare_parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat fixture warning classes as blockers unless explicitly handled by future policy.",
    )
    prepare_parser.set_defaults(func=_cmd_prepare_ppi_inputs)

    real_parser = subparsers.add_parser(
        "validate-real-inputs",
        help="Validate real EGFR/MYO1D/membrane-frame inputs for readiness without docking.",
    )
    real_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    real_parser.add_argument("--mode", default="smoke_input", choices=["smoke_env", "smoke_input"])
    real_parser.add_argument(
        "--profile",
        default="hpc_strict",
        choices=["codex_dev", "hpc_strict"],
        help="hpc_strict treats fixture-only warning classes as blockers/quarantine.",
    )
    real_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing real ppi_input_contract.json and referenced inputs.",
    )
    real_parser.add_argument(
        "--contract",
        help="Optional real contract path. Absolute paths and paths under --input-root are accepted.",
    )
    real_parser.add_argument(
        "--strict",
        action="store_true",
        help="Force production-like blocker policy even under codex_dev.",
    )
    real_parser.set_defaults(func=_cmd_validate_real_inputs)

    sampling_parser = subparsers.add_parser(
        "plan-ppi-sampling",
        help="Write spec-only EGFR-MYO1D PPI sampling jobs and future pose-QC policy.",
    )
    sampling_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    sampling_parser.add_argument("--mode", default="smoke_input", choices=["smoke_env", "smoke_input"])
    sampling_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev allows fixture quarantines; hpc_strict blocks production-unsafe inputs.",
    )
    sampling_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing ppi_input_contract.json and referenced inputs.",
    )
    sampling_parser.add_argument(
        "--contract",
        help="Optional PPI input contract path. Absolute paths and paths under --input-root are accepted.",
    )
    sampling_parser.add_argument(
        "--strict",
        action="store_true",
        help="Force production-like blocker policy even under codex_dev.",
    )
    sampling_parser.set_defaults(func=_cmd_plan_ppi_sampling)

    consensus_parser = subparsers.add_parser(
        "summarize-ppi-consensus",
        help="Summarize supplied PPI contact records into guarded EGFR-side consensus patches.",
    )
    consensus_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    consensus_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input"])
    consensus_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="Profile label recorded in consensus manifests; no external tools are run.",
    )
    consensus_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing accepted_ppi_contacts.csv or the selected contact table.",
    )
    consensus_parser.add_argument(
        "--contact-table",
        help="Optional contact table path relative to --input-root.",
    )
    consensus_parser.add_argument(
        "--input-kind",
        default="synthetic_fixture",
        choices=["synthetic_fixture", "real_input_derived", "future_accepted_pose_export"],
        help="Provenance label for supplied PPI contact records.",
    )
    consensus_parser.set_defaults(func=_cmd_summarize_ppi_consensus)

    pocket_parser = subparsers.add_parser(
        "plan-pocket-discovery",
        help="Plan PPI-guided EGFR pocket selection from Task 7 consensus evidence without running pocket tools.",
    )
    pocket_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    pocket_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input"])
    pocket_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="Profile label recorded in pocket-planning manifests; no external tools are run.",
    )
    pocket_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing ppi_consensus_patch.csv or the selected Task 7 consensus patch table.",
    )
    pocket_parser.add_argument(
        "--ppi-consensus",
        help="Optional Task 7 consensus patch CSV path relative to --input-root.",
    )
    pocket_parser.add_argument(
        "--input-kind",
        default="task7_consensus",
        choices=["task7_consensus", "synthetic_fixture", "real_input_derived"],
        help="Provenance label for supplied Task 7 PPI consensus evidence.",
    )
    pocket_parser.set_defaults(func=_cmd_plan_pocket_discovery)

    pocket_candidate_parser = subparsers.add_parser(
        "prioritize-pocket-candidates",
        help="Ingest provided pocket candidate records and prioritize PPI-guided non-ATP EGFR pockets without docking.",
    )
    pocket_candidate_parser.add_argument("--run-id", required=True, help="Run identifier under fresh/runs/.")
    pocket_candidate_parser.add_argument("--mode", default="smoke_env", choices=["smoke_env", "smoke_input"])
    pocket_candidate_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="Profile label recorded in pocket candidate manifests; no external tools are run.",
    )
    pocket_candidate_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="Directory containing pocket_discovery_plan.json and provided detector-style candidate CSV records.",
    )
    pocket_candidate_parser.add_argument(
        "--pocket-plan",
        help="Task 8 pocket_discovery_plan.json path relative to --input-root.",
    )
    pocket_candidate_parser.add_argument(
        "--candidate-pockets",
        help="Provided detector-style pocket candidate CSV path relative to --input-root.",
    )
    pocket_candidate_parser.add_argument(
        "--ppi-consensus",
        help="Optional Task 7 consensus CSV path relative to --input-root, recorded for provenance.",
    )
    pocket_candidate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Record strict mode policy; external runtimes still are not executed.",
    )
    pocket_candidate_parser.set_defaults(func=_cmd_prioritize_pocket_candidates)

    prepare_inputs_parser = subparsers.add_parser(
        "prepare-inputs",
        help="Orchestrator: preflight + prepare-receptor (per state) + compute-membrane-frame + prepare-myo1d + manifest-ligands.",
    )
    prepare_inputs_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/.",
    )
    prepare_inputs_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input"],
    )
    prepare_inputs_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="hpc_strict stops orchestration at first sub-step FAIL.",
    )
    prepare_inputs_parser.add_argument(
        "--input-root",
        default=None,
        help=(
            "Optional directory layout: "
            "<root>/{receptors,myo1d,ligands,private}/. "
            "Default: per-config paths (fresh/data/raw/...)."
        ),
    )
    prepare_inputs_parser.add_argument(
        "--states",
        default=None,
        help="Comma-separated state IDs. Default: all from receptor_states.yaml.",
    )
    prepare_inputs_parser.add_argument(
        "--skip-ligands",
        choices=["true", "false"],
        default="false",
        help="Skip the manifest-ligands sub-step.",
    )
    prepare_inputs_parser.add_argument(
        "--strict",
        action="store_true",
        help="Force production-like blocker policy even in codex_dev.",
    )
    prepare_inputs_parser.add_argument(
        "--compound-stage-enabled",
        choices=["true", "false"],
        default="false",
    )
    prepare_inputs_parser.set_defaults(func=_cmd_prepare_inputs)

    m2_ppi_parser = subparsers.add_parser(
        "generate-m2-ppi-inputs",
        help="Generate M2.1 PyRosetta/LightDock input packs/specs from M1 normalized outputs without execution.",
    )
    m2_ppi_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M1 prepare-inputs outputs should already exist in this run.",
    )
    m2_ppi_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.1 manifests.",
    )
    m2_ppi_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="hpc_strict treats skipped/non-primary states and missing M1 artifacts as blockers.",
    )
    m2_ppi_parser.add_argument(
        "--states",
        default=None,
        help="Comma-separated receptor states. Default: primary membrane-validated states from receptor_states.yaml.",
    )
    m2_ppi_parser.add_argument(
        "--strict",
        action="store_true",
        help="Force blocker behavior for non-primary/reference-control states.",
    )
    m2_ppi_parser.set_defaults(func=_cmd_generate_m2_ppi_inputs)

    m2_pyrosetta_parser = subparsers.add_parser(
        "prepare-m2-pyrosetta-harness",
        help="Prepare M2.2 PyRosetta adapter dry-run job manifests from M2.1 specs without docking.",
    )
    m2_pyrosetta_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M2.1 specs should already exist in this run.",
    )
    m2_pyrosetta_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Harness scale label. Smoke emits one state/seed dry-run job.",
    )
    m2_pyrosetta_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_pyrosetta_parser.add_argument(
        "--states",
        default=None,
        help="Optional comma-separated receptor states to include.",
    )
    m2_pyrosetta_parser.add_argument(
        "--m2-1-manifest",
        default=None,
        help="Optional M2.1 manifest path under the same run dir.",
    )
    m2_pyrosetta_parser.set_defaults(func=_cmd_prepare_m2_pyrosetta_harness)

    m2_collect_parser = subparsers.add_parser(
        "collect-m2-ppi-outputs",
        help="Collect M2.2 dry-run PPI job status and restore optional run-local raw contacts without docking.",
    )
    m2_collect_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M2.2 job manifests should already exist in this run.",
    )
    m2_collect_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.3 manifests.",
    )
    m2_collect_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_collect_parser.add_argument(
        "--job-manifest",
        default=None,
        help="Optional M2.2 job manifest JSONL path under the same run dir.",
    )
    m2_collect_parser.add_argument(
        "--raw-contact-table",
        default=None,
        help="Optional run-local raw contact CSV to restore through M1 receptor mapping.",
    )
    m2_collect_parser.add_argument(
        "--unmapped-fraction-fail-threshold",
        type=float,
        default=0.0,
        help="Fail if restored contact unmapped fraction is greater than this threshold.",
    )
    m2_collect_parser.set_defaults(func=_cmd_collect_m2_ppi_outputs)

    m2_consensus_parser = subparsers.add_parser(
        "build-m2-ppi-consensus-patch",
        help="Build M2.4 chain-resolved EGFR-side PPI consensus patch tables from M2.3 restored contacts.",
    )
    m2_consensus_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M2.3 contact tables should already exist in this run.",
    )
    m2_consensus_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.4 manifests.",
    )
    m2_consensus_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_consensus_parser.add_argument(
        "--contact-table",
        default=None,
        help="Optional contact table path inside the run dir. Default: phase1_ppi/tables/ppi_pose_contacts.csv.",
    )
    m2_consensus_parser.add_argument(
        "--job-manifest",
        default=None,
        help="Optional M2.2 job manifest JSONL path inside the run dir.",
    )
    m2_consensus_parser.set_defaults(func=_cmd_build_m2_ppi_consensus_patch)

    m2_atp_parser = subparsers.add_parser(
        "build-m2-atp-reference",
        help="Build M2.5 ATP-site reference tables for later non-ATP pocket gates without pocket discovery.",
    )
    m2_atp_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M1 normalized receptor outputs should already exist in this run.",
    )
    m2_atp_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.5 manifests.",
    )
    m2_atp_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_atp_parser.add_argument(
        "--states",
        default=None,
        help="Optional comma-separated receptor states. Default: primary states plus reference controls.",
    )
    m2_atp_parser.add_argument(
        "--include-reference-states",
        default="true",
        choices=["true", "false"],
        help="Include reference/control states such as 3GT8_raw in mapping outputs.",
    )
    m2_atp_parser.add_argument(
        "--config",
        default=None,
        help="Optional ATP reference config path inside fresh/ or the run dir.",
    )
    m2_atp_parser.add_argument(
        "--reference-pdb",
        default=None,
        help="Optional ligand-bearing ATP reference PDB inside fresh/ or the run dir.",
    )
    m2_atp_parser.add_argument(
        "--synthetic-fixture",
        default="false",
        choices=["true", "false"],
        help="Mark output reference_mode as synthetic_fixture for test fixtures.",
    )
    m2_atp_parser.set_defaults(func=_cmd_build_m2_atp_reference)

    m2_fpocket_parser = subparsers.add_parser(
        "run-m2-fpocket-discovery",
        help="Run or parse M2.6 fpocket pocket discovery and raw pocket normalization without M2.7 gates.",
    )
    m2_fpocket_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M2.4/M2.5 outputs should already exist for production interpretation.",
    )
    m2_fpocket_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.6 manifests.",
    )
    m2_fpocket_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_fpocket_parser.add_argument(
        "--states",
        default=None,
        help="Optional comma-separated receptor states. Default: primary states plus reference controls.",
    )
    m2_fpocket_parser.add_argument(
        "--include-reference-states",
        default="true",
        choices=["true", "false"],
        help="Include reference/control states such as 3GT8_raw in raw evidence outputs.",
    )
    m2_fpocket_parser.add_argument(
        "--execution-mode",
        default="parser_only",
        choices=["parser_only", "production"],
        help="parser_only parses existing fpocket output; production invokes fpocket on staged run-local receptors.",
    )
    m2_fpocket_parser.add_argument(
        "--fpocket-output-root",
        default=None,
        help="Parser-only fpocket output root inside fresh/ or the run dir.",
    )
    m2_fpocket_parser.add_argument(
        "--config",
        default=None,
        help="Optional pocket config path inside fresh/ or the run dir.",
    )
    m2_fpocket_parser.add_argument(
        "--fpocket-binary",
        default=None,
        help="Optional fpocket executable name/path for production mode.",
    )
    m2_fpocket_parser.add_argument(
        "--synthetic-fixture",
        default="false",
        choices=["true", "false"],
        help="Mark parser-only evidence as synthetic fixture data for tests.",
    )
    m2_fpocket_parser.set_defaults(func=_cmd_run_m2_fpocket_discovery)

    m2_gate_parser = subparsers.add_parser(
        "gate-m2-pockets",
        help="Apply M2.7 PPI/membrane/dimer hard gates to M2.6 pocket families without M2.8 export or docking.",
    )
    m2_gate_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M2.4, M2.5, and M2.6 outputs should already exist for production gating.",
    )
    m2_gate_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.7 manifests.",
    )
    m2_gate_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_gate_parser.add_argument("--config", default=None, help="Optional pocket config path inside fresh/ or the run dir.")
    m2_gate_parser.add_argument("--ppi-consensus", default=None, help="Optional M2.4 consensus patch CSV inside fresh/ or the run dir.")
    m2_gate_parser.add_argument("--atp-reference", default=None, help="Optional M2.5 ATP reference CSV inside fresh/ or the run dir.")
    m2_gate_parser.add_argument("--atp-centroids", default=None, help="Optional M2.5 ATP centroid CSV inside fresh/ or the run dir.")
    m2_gate_parser.add_argument("--pocket-families", default=None, help="Optional M2.6 merged pocket family CSV inside fresh/ or the run dir.")
    m2_gate_parser.add_argument("--raw-candidates", default=None, help="Optional M2.6 raw candidate CSV inside fresh/ or the run dir.")
    m2_gate_parser.add_argument("--membrane-frame", default=None, help="Optional membrane_frame.json inside fresh/ or the run dir.")
    m2_gate_parser.add_argument(
        "--synthetic-fixture",
        default="false",
        choices=["true", "false"],
        help="Allow synthetic/smoke fixtures with warnings instead of production blocker semantics.",
    )
    m2_gate_parser.set_defaults(func=_cmd_gate_m2_pockets)

    m2_export_parser = subparsers.add_parser(
        "export-m2-results",
        help="Build M2.8 accepted pocket export package and Milestone 2 report without starting M3.",
    )
    m2_export_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/. M2.4-M2.7 outputs should already exist for production export.",
    )
    m2_export_parser.add_argument(
        "--mode",
        default="smoke_input",
        choices=["smoke_env", "smoke_input", "mini", "production"],
        help="Mode label recorded in M2.8 manifests.",
    )
    m2_export_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    m2_export_parser.add_argument("--ppi-consensus", default=None, help="Optional M2.4 consensus patch CSV inside fresh/ or the run dir.")
    m2_export_parser.add_argument("--atp-reference", default=None, help="Optional M2.5 ATP reference CSV inside fresh/ or the run dir.")
    m2_export_parser.add_argument("--pocket-families", default=None, help="Optional M2.6 merged pocket family CSV inside fresh/ or the run dir.")
    m2_export_parser.add_argument("--pocket-gate-qc", default=None, help="Optional M2.7 pocket gate QC CSV inside fresh/ or the run dir.")
    m2_export_parser.add_argument("--accepted-families", default=None, help="Optional M2.7 accepted pocket families CSV inside fresh/ or the run dir.")
    m2_export_parser.add_argument("--rejected-families", default=None, help="Optional M2.7 rejected pocket families CSV inside fresh/ or the run dir.")
    m2_export_parser.add_argument(
        "--synthetic-fixture",
        default="false",
        choices=["true", "false"],
        help="Allow synthetic/smoke fixtures with warnings instead of production blocker semantics.",
    )
    m2_export_parser.set_defaults(func=_cmd_export_m2_results)

    m3_readiness_parser = subparsers.add_parser(
        "m3-readiness",
        help="Run M3-T0 readiness audit and create the safe phase3 directory skeleton without docking.",
    )
    m3_readiness_parser.add_argument(
        "--run-id",
        required=True,
        help="M3 run identifier under fresh/runs/. Outputs are written only to this run.",
    )
    m3_readiness_parser.add_argument(
        "--m2-run-id",
        required=True,
        help="M2 source run identifier under fresh/runs/.",
    )
    m3_readiness_parser.add_argument(
        "--mode",
        default="dry-run",
        choices=["dry-run", "docking"],
        help="dry-run warns on missing ligand files; docking requires them.",
    )
    m3_readiness_parser.add_argument(
        "--create-skeleton",
        dest="create_skeleton",
        action="store_true",
        default=True,
        help="Create the M3 phase3_compounds directory skeleton (default).",
    )
    m3_readiness_parser.add_argument(
        "--no-create-skeleton",
        dest="create_skeleton",
        action="store_false",
        help="Create only the output/log directories required for the readiness report.",
    )
    m3_readiness_parser.set_defaults(func=_cmd_m3_readiness)

    m3_ligand_qc_parser = subparsers.add_parser(
        "m3-ligand-qc",
        help="Run M3-T1 public-safe ligand manifest, confidentiality audit, and ligand file QC without preparation.",
    )
    m3_ligand_qc_parser.add_argument(
        "--run-id",
        required=True,
        help="M3 run identifier under fresh/runs/.",
    )
    m3_ligand_qc_parser.add_argument(
        "--mode",
        default="dry-run",
        choices=["dry-run", "docking"],
        help="dry-run warns on missing ligand files; docking requires them.",
    )
    m3_ligand_qc_parser.add_argument(
        "--private-map",
        default="fresh/data/private/compound_id_map.csv",
        help="Private public_id/internal_id map inside fresh/; internal IDs are never written to public outputs.",
    )
    m3_ligand_qc_parser.add_argument(
        "--update-gitignore",
        dest="update_gitignore",
        action="store_true",
        default=True,
        help="Append missing sensitive ligand/private-data ignore patterns (default).",
    )
    m3_ligand_qc_parser.add_argument(
        "--no-update-gitignore",
        dest="update_gitignore",
        action="store_false",
        help="Verify existing gitignore protection without modifying .gitignore.",
    )
    m3_ligand_qc_parser.set_defaults(func=_cmd_m3_ligand_qc)

    m3_prepare_ligands_parser = subparsers.add_parser(
        "m3-prepare-ligands",
        help="Run M3-T2 ligand PDBQT preparation/QC for public Cpd-A/B/C without receptor or Vina work.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--run-id",
        required=True,
        help="M3 run identifier under fresh/runs/.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--mode",
        default="prepare",
        choices=["dry-run", "prepare"],
        help="dry-run reports what would happen; prepare writes ligand PDBQT outputs when allowed.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--private-map",
        default="fresh/data/private/compound_id_map.csv",
        help="Private public_id/internal_id map inside fresh/; internal IDs are never written to public outputs.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Overwrite existing prepared ligand PDBQT files.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--no-force",
        dest="force",
        action="store_false",
        help="Do not overwrite existing prepared ligand PDBQT files (default).",
    )
    m3_prepare_ligands_parser.add_argument(
        "--keep-intermediates",
        dest="keep_intermediates",
        action="store_true",
        default=False,
        help="Keep conversion intermediates under run-local scratch when future routes create them.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--no-keep-intermediates",
        dest="keep_intermediates",
        action="store_false",
        help="Do not keep conversion intermediates (default).",
    )
    m3_prepare_ligands_parser.add_argument(
        "--update-gitignore",
        dest="update_gitignore",
        action="store_true",
        default=True,
        help="Append missing sensitive ligand/prepared-output ignore patterns (default).",
    )
    m3_prepare_ligands_parser.add_argument(
        "--no-update-gitignore",
        dest="update_gitignore",
        action="store_false",
        help="Verify existing gitignore protection without modifying .gitignore.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--allow-independent",
        action="store_true",
        help="Allow diagnostic/synthetic execution even when M3-T1 did not allow M3-T2.",
    )
    m3_prepare_ligands_parser.add_argument(
        "--require-m3-t1-pass",
        dest="require_m3_t1_pass",
        action="store_true",
        default=True,
        help="Require M3-T1 not to have failed before prepare mode proceeds (default).",
    )
    m3_prepare_ligands_parser.add_argument(
        "--no-require-m3-t1-pass",
        dest="require_m3_t1_pass",
        action="store_false",
        help="Do not require M3-T1 PASS/clean status; recorded for diagnostics only.",
    )
    m3_prepare_ligands_parser.set_defaults(func=_cmd_m3_prepare_ligands)

    for command_name, help_text, operation in [
        (
            "m3-prepare-receptors",
            "Run M3-T3 receptor PDBQT preparation and receptor provenance QC without ligand/Vina work.",
            "prepare-receptors",
        ),
        (
            "m3-validate-boxes",
            "Run M3-T3 M2 accepted docking-box validation without Vina/job generation.",
            "validate-boxes",
        ),
        (
            "m3-receptor-boxes",
            "Run combined M3-T3 receptor PDBQT preparation and docking-box validation.",
            "combined",
        ),
    ]:
        m3_receptor_box_parser = subparsers.add_parser(command_name, help=help_text)
        m3_receptor_box_parser.add_argument("--run-id", required=True, help="M3 run identifier under fresh/runs/.")
        m3_receptor_box_parser.add_argument("--m2-run-id", required=True, help="M2 source run identifier under fresh/runs/.")
        m3_receptor_box_parser.add_argument(
            "--mode",
            default="prepare",
            choices=["dry-run", "prepare"],
            help="dry-run reports what would happen; prepare writes receptor PDBQT and validation outputs when allowed.",
        )
        m3_receptor_box_parser.add_argument("--force", dest="force", action="store_true", default=False, help="Overwrite existing receptor PDBQT outputs.")
        m3_receptor_box_parser.add_argument("--no-force", dest="force", action="store_false", help="Do not overwrite existing receptor PDBQT outputs (default).")
        m3_receptor_box_parser.add_argument("--box-margin", type=float, default=None, help="Margin added to each valid box face after raw validation. Default reads fresh/configs/gates.yaml or 0.0.")
        m3_receptor_box_parser.add_argument("--states", default=None, help="Comma-separated state IDs. Default: all states represented in M2 accepted_pocket_boxes.csv.")
        m3_receptor_box_parser.add_argument("--update-gitignore", dest="update_gitignore", action="store_true", default=True, help="Append missing receptor PDBQT ignore patterns (default).")
        m3_receptor_box_parser.add_argument("--no-update-gitignore", dest="update_gitignore", action="store_false", help="Verify existing gitignore protection without modifying .gitignore.")
        m3_receptor_box_parser.add_argument("--allow-independent", action="store_true", help="Record diagnostic/synthetic execution outside normal M3-T2 pass gate.")
        m3_receptor_box_parser.add_argument("--require-m3-t2-pass", dest="require_m3_t2_pass", action="store_true", default=True, help="Require M3-T2 not to have failed before prepare mode proceeds (default).")
        m3_receptor_box_parser.add_argument("--no-require-m3-t2-pass", dest="require_m3_t2_pass", action="store_false", help="Do not require M3-T2 PASS/clean status; recorded for diagnostics only.")
        m3_receptor_box_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv", help="Private public_id/internal_id map inside fresh/; internal IDs are never written to public outputs.")
        m3_receptor_box_parser.set_defaults(func=_cmd_m3_receptor_boxes, operation=operation)

    m3_vina_smoke_parser = subparsers.add_parser(
        "m3-vina-smoke",
        help="Run M3-T4 focused Vina smoke planning or exactly one smoke command without production docking.",
    )
    m3_vina_smoke_parser.add_argument("--run-id", required=True, help="M3 run identifier under fresh/runs/.")
    m3_vina_smoke_parser.add_argument("--m2-run-id", required=True, help="M2 source run identifier under fresh/runs/.")
    m3_vina_smoke_parser.add_argument("--mode", default="run", choices=["dry-run", "run"], help="dry-run writes command manifest only; run invokes exactly one Vina command when allowed.")
    m3_vina_smoke_parser.add_argument("--force", dest="force", action="store_true", default=False, help="Overwrite existing smoke output files.")
    m3_vina_smoke_parser.add_argument("--no-force", dest="force", action="store_false", help="Do not overwrite existing smoke output files (default).")
    m3_vina_smoke_parser.add_argument("--compound-public-id", choices=["Cpd-A", "Cpd-B", "Cpd-C"], default=None, help="Optional explicit public compound ID.")
    m3_vina_smoke_parser.add_argument("--state-id", default=None, help="Optional explicit receptor state ID.")
    m3_vina_smoke_parser.add_argument("--pocket-family-id", default=None, help="Optional explicit pocket family ID.")
    m3_vina_smoke_parser.add_argument("--box-id", default=None, help="Optional explicit docking box ID.")
    m3_vina_smoke_parser.add_argument("--protomer-id", default=None, help="Optional explicit protomer ID.")
    m3_vina_smoke_parser.add_argument("--vina-executable", default="vina", help="Vina executable path or name.")
    m3_vina_smoke_parser.add_argument("--exhaustiveness", type=int, default=None, help="Vina exhaustiveness; default from gates.yaml or 1.")
    m3_vina_smoke_parser.add_argument("--num-modes", type=int, default=None, help="Vina num_modes; default from gates.yaml or 1.")
    m3_vina_smoke_parser.add_argument("--cpu", type=int, default=1, help="Vina CPU count for smoke test.")
    m3_vina_smoke_parser.add_argument("--seed", type=int, default=None, help="Vina seed; default from gates.yaml or 20260427.")
    m3_vina_smoke_parser.add_argument("--timeout-seconds", type=int, default=None, help="Vina timeout; default from gates.yaml or 300.")
    m3_vina_smoke_parser.add_argument("--require-m3-t2-pass", dest="require_m3_t2_pass", action="store_true", default=True, help="Require M3-T2 PASS before run mode proceeds (default).")
    m3_vina_smoke_parser.add_argument("--no-require-m3-t2-pass", dest="require_m3_t2_pass", action="store_false", help="Do not require M3-T2 PASS; diagnostic only.")
    m3_vina_smoke_parser.add_argument("--require-m3-t3-pass", dest="require_m3_t3_pass", action="store_true", default=True, help="Require M3-T3 PASS and m3_t4_allowed=true before run mode proceeds (default).")
    m3_vina_smoke_parser.add_argument("--no-require-m3-t3-pass", dest="require_m3_t3_pass", action="store_false", help="Do not require M3-T3 PASS; diagnostic only.")
    m3_vina_smoke_parser.add_argument("--allow-independent", action="store_true", help="Record diagnostic bypass and keep production gates closed.")
    m3_vina_smoke_parser.add_argument("--allow-reference-smoke", action="store_true", help="Allow 3GT8_raw/reference smoke diagnostically; M3-T5 remains closed.")
    m3_vina_smoke_parser.add_argument("--update-gitignore", dest="update_gitignore", action="store_true", default=True, help="Append missing Vina smoke output ignore patterns (default).")
    m3_vina_smoke_parser.add_argument("--no-update-gitignore", dest="update_gitignore", action="store_false", help="Verify gitignore protection without modifying .gitignore.")
    m3_vina_smoke_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv", help="Private public_id/internal_id map inside fresh/; internal IDs are never written to public outputs.")
    m3_vina_smoke_parser.set_defaults(func=_cmd_m3_vina_smoke)

    m3_vina_jobs_parser = subparsers.add_parser(
        "m3-vina-jobs",
        help="Generate M3-T5 focused Vina job manifests and PBS scripts without qsub or Vina execution.",
    )
    m3_vina_jobs_parser.add_argument("--run-id", required=True, help="M3 run identifier under fresh/runs/.")
    m3_vina_jobs_parser.add_argument("--m2-run-id", required=True, help="M2 source run identifier under fresh/runs/.")
    m3_vina_jobs_parser.add_argument("--mode", default="generate", choices=["dry-run", "generate"], help="dry-run writes planning outputs; generate writes submit-ready PBS files when gates pass.")
    m3_vina_jobs_parser.add_argument("--profile", default="mini", choices=["smoke", "mini", "scaling", "production"], help="Job-plan scope profile.")
    m3_vina_jobs_parser.add_argument("--force", dest="force", action="store_true", default=False, help="Regenerate existing planning outputs.")
    m3_vina_jobs_parser.add_argument("--no-force", dest="force", action="store_false", help="Do not overwrite planned outputs unless safe (default).")
    m3_vina_jobs_parser.add_argument("--compound-public-id", choices=["Cpd-A", "Cpd-B", "Cpd-C"], default=None)
    m3_vina_jobs_parser.add_argument("--state-id", default=None)
    m3_vina_jobs_parser.add_argument("--pocket-family-id", default=None)
    m3_vina_jobs_parser.add_argument("--box-id", default=None)
    m3_vina_jobs_parser.add_argument("--protomer-id", default=None)
    m3_vina_jobs_parser.add_argument("--include-reference", action="store_true")
    m3_vina_jobs_parser.add_argument("--include-3gt8", action="store_true")
    m3_vina_jobs_parser.add_argument("--allow-reference-jobs", action="store_true")
    m3_vina_jobs_parser.add_argument("--allow-independent", action="store_true", help="Record diagnostic bypass and keep submission gates closed.")
    m3_vina_jobs_parser.add_argument("--allow-production-plan-without-mini", action="store_true", help="Allow production planning estimate without mini PASS; submission remains closed.")
    m3_vina_jobs_parser.add_argument("--vina-executable", default="vina")
    m3_vina_jobs_parser.add_argument("--exhaustiveness", type=int, default=None)
    m3_vina_jobs_parser.add_argument("--num-modes", type=int, default=None)
    m3_vina_jobs_parser.add_argument("--cpu-per-vina", type=int, default=1)
    m3_vina_jobs_parser.add_argument("--repeats", type=int, default=None)
    m3_vina_jobs_parser.add_argument("--base-seed", type=int, default=None)
    m3_vina_jobs_parser.add_argument("--timeout-seconds", type=int, default=None)
    m3_vina_jobs_parser.add_argument("--nodes", default="node04,node05,node06")
    m3_vina_jobs_parser.add_argument("--ppn", type=int, default=None)
    m3_vina_jobs_parser.add_argument("--walltime", default=None)
    m3_vina_jobs_parser.add_argument("--queue", default=None)
    m3_vina_jobs_parser.add_argument("--conda-env", default=None)
    m3_vina_jobs_parser.add_argument("--conda-sh", default=None)
    m3_vina_jobs_parser.add_argument("--python-executable", default=None)
    m3_vina_jobs_parser.add_argument("--update-gitignore", dest="update_gitignore", action="store_true", default=True)
    m3_vina_jobs_parser.add_argument("--no-update-gitignore", dest="update_gitignore", action="store_false")
    m3_vina_jobs_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_vina_jobs_parser.set_defaults(func=_cmd_m3_vina_jobs)

    m3_run_vina_chunk_parser = subparsers.add_parser(
        "m3-run-vina-chunk",
        help="Run one pre-planned M3-T5 Vina manifest chunk; intended for generated PBS scripts.",
    )
    m3_run_vina_chunk_parser.add_argument("--run-id", required=True)
    m3_run_vina_chunk_parser.add_argument("--job-manifest", required=True)
    m3_run_vina_chunk_parser.add_argument("--chunk-id", required=True)
    m3_run_vina_chunk_parser.add_argument("--max-workers", type=int, default=1)
    m3_run_vina_chunk_parser.add_argument("--dry-run-runner", action="store_true", help="Validate chunk selection/path contract without invoking Vina.")
    m3_run_vina_chunk_parser.set_defaults(func=_cmd_m3_run_vina_chunk)

    m3_vina_collect_parser = subparsers.add_parser(
        "m3-vina-collect",
        help="Collect M3-T6 Vina outputs into raw pose and completion tables without running Vina or attribution.",
    )
    m3_vina_collect_parser.add_argument("--run-id", required=True)
    m3_vina_collect_parser.add_argument("--m2-run-id", required=True)
    m3_vina_collect_parser.add_argument("--job-manifest", default=None)
    m3_vina_collect_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_vina_collect_parser.add_argument("--mode", default="collect", choices=["dry-run", "collect"])
    m3_vina_collect_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_vina_collect_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_vina_collect_parser.add_argument("--allow-partial", action="store_true")
    m3_vina_collect_parser.add_argument("--require-all-jobs", dest="require_all_jobs", action="store_true", default=True)
    m3_vina_collect_parser.add_argument("--no-require-all-jobs", dest="require_all_jobs", action="store_false")
    m3_vina_collect_parser.add_argument("--allow-reference-only", action="store_true")
    m3_vina_collect_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_vina_collect_parser.add_argument("--update-gitignore", dest="update_gitignore", action="store_true", default=True)
    m3_vina_collect_parser.add_argument("--no-update-gitignore", dest="update_gitignore", action="store_false")
    m3_vina_collect_parser.add_argument("--max-log-bytes", type=int, default=200000)
    m3_vina_collect_parser.add_argument("--strict-log-parse", dest="strict_log_parse", action="store_true", default=True)
    m3_vina_collect_parser.add_argument("--no-strict-log-parse", dest="strict_log_parse", action="store_false")
    m3_vina_collect_parser.add_argument("--write-empty-raw-table", dest="write_empty_raw_table", action="store_true", default=False)
    m3_vina_collect_parser.add_argument("--no-write-empty-raw-table", dest="write_empty_raw_table", action="store_false")
    m3_vina_collect_parser.set_defaults(func=_cmd_m3_vina_collect)

    m3_pose_attribution_parser = subparsers.add_parser(
        "m3-pose-attribution",
        help="Run M3-T7 pose-level pocket/ATP/PPI/membrane geometry attribution without docking, clustering, or scoring.",
    )
    m3_pose_attribution_parser.add_argument("--run-id", required=True)
    m3_pose_attribution_parser.add_argument("--m2-run-id", required=True)
    m3_pose_attribution_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_pose_attribution_parser.add_argument("--mode", default="attribute", choices=["dry-run", "attribute"])
    m3_pose_attribution_parser.add_argument("--pose-raw-table", default=None)
    m3_pose_attribution_parser.add_argument("--job-completion-table", default=None)
    m3_pose_attribution_parser.add_argument("--accepted-pockets", default=None)
    m3_pose_attribution_parser.add_argument("--accepted-boxes", default=None)
    m3_pose_attribution_parser.add_argument("--atp-reference", default=None)
    m3_pose_attribution_parser.add_argument("--ppi-consensus-patch", default=None)
    m3_pose_attribution_parser.add_argument("--membrane-frame", default=None)
    m3_pose_attribution_parser.add_argument("--receptor-mapping", default=None)
    m3_pose_attribution_parser.add_argument("--receptor-manifest", default=None)
    m3_pose_attribution_parser.add_argument("--config", default=None)
    m3_pose_attribution_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_pose_attribution_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_pose_attribution_parser.add_argument("--allow-partial", action="store_true")
    m3_pose_attribution_parser.add_argument("--require-all-poses", dest="require_all_poses", action="store_true", default=True)
    m3_pose_attribution_parser.add_argument("--no-require-all-poses", dest="require_all_poses", action="store_false")
    m3_pose_attribution_parser.add_argument("--allow-reference-only", action="store_true")
    m3_pose_attribution_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_pose_attribution_parser.add_argument("--strict-geometry", dest="strict_geometry", action="store_true", default=True)
    m3_pose_attribution_parser.add_argument("--no-strict-geometry", dest="strict_geometry", action="store_false")
    m3_pose_attribution_parser.add_argument("--write-empty-attribution", dest="write_empty_attribution", action="store_true", default=False)
    m3_pose_attribution_parser.add_argument("--no-write-empty-attribution", dest="write_empty_attribution", action="store_false")
    m3_pose_attribution_parser.set_defaults(func=_cmd_m3_pose_attribution)

    m3_pose_clustering_parser = subparsers.add_parser(
        "m3-pose-clustering",
        help="Run M3-T8 pose clustering/convergence analysis without docking, attribution reruns, scoring, or candidate claims.",
    )
    m3_pose_clustering_parser.add_argument("--run-id", required=True)
    m3_pose_clustering_parser.add_argument("--m2-run-id", required=True)
    m3_pose_clustering_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_pose_clustering_parser.add_argument("--mode", default="cluster", choices=["dry-run", "cluster"])
    m3_pose_clustering_parser.add_argument("--pose-attribution-table", default=None)
    m3_pose_clustering_parser.add_argument("--pose-attribution-qc", default=None)
    m3_pose_clustering_parser.add_argument("--accepted-boxes", default=None)
    m3_pose_clustering_parser.add_argument("--config", default=None)
    m3_pose_clustering_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_pose_clustering_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_pose_clustering_parser.add_argument("--allow-partial", action="store_true")
    m3_pose_clustering_parser.add_argument("--require-attribution-ready", dest="require_attribution_ready", action="store_true", default=True)
    m3_pose_clustering_parser.add_argument("--no-require-attribution-ready", dest="require_attribution_ready", action="store_false")
    m3_pose_clustering_parser.add_argument("--allow-reference-only", action="store_true")
    m3_pose_clustering_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_pose_clustering_parser.add_argument("--cluster-method", default="hybrid", choices=["centroid", "rmsd", "hybrid"])
    m3_pose_clustering_parser.add_argument("--strict-convergence", dest="strict_convergence", action="store_true", default=True)
    m3_pose_clustering_parser.add_argument("--no-strict-convergence", dest="strict_convergence", action="store_false")
    m3_pose_clustering_parser.add_argument("--write-empty-clusters", dest="write_empty_clusters", action="store_true", default=False)
    m3_pose_clustering_parser.add_argument("--no-write-empty-clusters", dest="write_empty_clusters", action="store_false")
    m3_pose_clustering_parser.add_argument("--write-cluster-membership", dest="write_cluster_membership", action="store_true", default=True)
    m3_pose_clustering_parser.add_argument("--no-write-cluster-membership", dest="write_cluster_membership", action="store_false")
    m3_pose_clustering_parser.set_defaults(func=_cmd_m3_pose_clustering)

    m3_anchor_convergence_parser = subparsers.add_parser(
        "m3-anchor-convergence",
        help="Run M3-T9 compound-anchor convergence aggregation without docking, reruns, scoring, tiering, or candidate claims.",
    )
    m3_anchor_convergence_parser.add_argument("--run-id", required=True)
    m3_anchor_convergence_parser.add_argument("--m2-run-id", required=True)
    m3_anchor_convergence_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_anchor_convergence_parser.add_argument("--mode", default="converge", choices=["dry-run", "converge"])
    m3_anchor_convergence_parser.add_argument("--compound-pose-clusters", default=None)
    m3_anchor_convergence_parser.add_argument("--cluster-membership-table", default=None)
    m3_anchor_convergence_parser.add_argument("--pose-attribution-table", default=None)
    m3_anchor_convergence_parser.add_argument("--pose-clustering-qc", default=None)
    m3_anchor_convergence_parser.add_argument("--accepted-pockets", default=None)
    m3_anchor_convergence_parser.add_argument("--accepted-boxes", default=None)
    m3_anchor_convergence_parser.add_argument("--compound-manifest", default=None)
    m3_anchor_convergence_parser.add_argument("--config", default=None)
    m3_anchor_convergence_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_anchor_convergence_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_anchor_convergence_parser.add_argument("--allow-partial", action="store_true")
    m3_anchor_convergence_parser.add_argument("--require-clustering-ready", dest="require_clustering_ready", action="store_true", default=True)
    m3_anchor_convergence_parser.add_argument("--no-require-clustering-ready", dest="require_clustering_ready", action="store_false")
    m3_anchor_convergence_parser.add_argument("--allow-reference-only", action="store_true")
    m3_anchor_convergence_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_anchor_convergence_parser.add_argument("--strict-primary-state-support", dest="strict_primary_state_support", action="store_true", default=True)
    m3_anchor_convergence_parser.add_argument("--no-strict-primary-state-support", dest="strict_primary_state_support", action="store_false")
    m3_anchor_convergence_parser.add_argument("--write-empty-convergence", dest="write_empty_convergence", action="store_true", default=False)
    m3_anchor_convergence_parser.add_argument("--no-write-empty-convergence", dest="write_empty_convergence", action="store_false")
    m3_anchor_convergence_parser.set_defaults(func=_cmd_m3_anchor_convergence)

    m3_evidence_tiering_parser = subparsers.add_parser(
        "m3-evidence-tiering",
        help="Run M3-T10 integrated evidence table and cautious candidate-hypothesis tiering without upstream reruns.",
    )
    m3_evidence_tiering_parser.add_argument("--run-id", required=True)
    m3_evidence_tiering_parser.add_argument("--m2-run-id", required=True)
    m3_evidence_tiering_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_evidence_tiering_parser.add_argument("--mode", default="integrate", choices=["dry-run", "integrate"])
    m3_evidence_tiering_parser.add_argument("--accepted-pockets", default=None)
    m3_evidence_tiering_parser.add_argument("--accepted-boxes", default=None)
    m3_evidence_tiering_parser.add_argument("--pocket-gate-qc", default=None)
    m3_evidence_tiering_parser.add_argument("--ppi-consensus", default=None)
    m3_evidence_tiering_parser.add_argument("--ppi-to-pocket-evidence", default=None)
    m3_evidence_tiering_parser.add_argument("--pocket-ppi-relationship", default=None)
    m3_evidence_tiering_parser.add_argument("--pocket-membrane-geometry", default=None)
    m3_evidence_tiering_parser.add_argument("--atp-reference", default=None)
    m3_evidence_tiering_parser.add_argument("--compound-pose-attribution", default=None)
    m3_evidence_tiering_parser.add_argument("--compound-mechanism-classification", default=None)
    m3_evidence_tiering_parser.add_argument("--compound-pose-clusters", default=None)
    m3_evidence_tiering_parser.add_argument("--compound-pocket-support", default=None)
    m3_evidence_tiering_parser.add_argument("--compound-anchor-convergence", default=None)
    m3_evidence_tiering_parser.add_argument("--pose-attribution-qc", default=None)
    m3_evidence_tiering_parser.add_argument("--pose-clustering-qc", default=None)
    m3_evidence_tiering_parser.add_argument("--anchor-convergence-qc", default=None)
    m3_evidence_tiering_parser.add_argument("--compound-manifest", default=None)
    m3_evidence_tiering_parser.add_argument("--broad-anchor-scan-support", default=None)
    m3_evidence_tiering_parser.add_argument("--allow-broad-scan-support", dest="allow_broad_scan_support", action="store_true", default=False)
    m3_evidence_tiering_parser.add_argument("--no-allow-broad-scan-support", dest="allow_broad_scan_support", action="store_false")
    m3_evidence_tiering_parser.add_argument("--config", default=None)
    m3_evidence_tiering_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_evidence_tiering_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_evidence_tiering_parser.add_argument("--allow-partial", action="store_true")
    m3_evidence_tiering_parser.add_argument("--require-anchor-ready", dest="require_anchor_ready", action="store_true", default=True)
    m3_evidence_tiering_parser.add_argument("--no-require-anchor-ready", dest="require_anchor_ready", action="store_false")
    m3_evidence_tiering_parser.add_argument("--strict-tier1", dest="strict_tier1", action="store_true", default=True)
    m3_evidence_tiering_parser.add_argument("--no-strict-tier1", dest="strict_tier1", action="store_false")
    m3_evidence_tiering_parser.add_argument("--write-empty-evidence", dest="write_empty_evidence", action="store_true", default=False)
    m3_evidence_tiering_parser.add_argument("--no-write-empty-evidence", dest="write_empty_evidence", action="store_false")
    m3_evidence_tiering_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_evidence_tiering_parser.set_defaults(func=_cmd_m3_evidence_tiering)

    m3_report_parser = subparsers.add_parser(
        "m3-report",
        help="Run M3-T11 final milestone report and handoff packaging without rerunning upstream work.",
    )
    m3_report_parser.add_argument("--run-id", required=True)
    m3_report_parser.add_argument("--m2-run-id", required=True)
    m3_report_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_report_parser.add_argument("--mode", default="report", choices=["dry-run", "report"])
    m3_report_parser.add_argument("--allow-partial", dest="allow_partial", action="store_true", default=False)
    m3_report_parser.add_argument("--no-allow-partial", dest="allow_partial", action="store_false")
    m3_report_parser.add_argument("--require-t10-ready", dest="require_t10_ready", action="store_true", default=True)
    m3_report_parser.add_argument("--no-require-t10-ready", dest="require_t10_ready", action="store_false")
    m3_report_parser.add_argument("--strict-confidentiality", dest="strict_confidentiality", action="store_true", default=True)
    m3_report_parser.add_argument("--no-strict-confidentiality", dest="strict_confidentiality", action="store_false")
    m3_report_parser.add_argument("--strict-claims", dest="strict_claims", action="store_true", default=True)
    m3_report_parser.add_argument("--no-strict-claims", dest="strict_claims", action="store_false")
    m3_report_parser.add_argument("--config", default=None)
    m3_report_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_report_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_report_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_report_parser.set_defaults(func=_cmd_m3_report)

    m3_cleanup_parser = subparsers.add_parser(
        "m3-cleanup",
        help="Plan or execute M3-T11 cleanup with run-dir-only safety gates.",
    )
    m3_cleanup_parser.add_argument("--run-id", required=True)
    m3_cleanup_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_cleanup_parser.add_argument("--mode", default="plan", choices=["plan", "execute"])
    m3_cleanup_parser.add_argument("--cleanup-scope", default="conservative", choices=["test", "mini", "production", "conservative", "aggressive-test-only"])
    m3_cleanup_parser.add_argument("--config", default=None)
    m3_cleanup_parser.add_argument("--preserve-selected-poses", dest="preserve_selected_poses", action="store_true", default=True)
    m3_cleanup_parser.add_argument("--no-preserve-selected-poses", dest="preserve_selected_poses", action="store_false")
    m3_cleanup_parser.add_argument("--preserve-all-docking-outputs", dest="preserve_all_docking_outputs", action="store_true", default=True)
    m3_cleanup_parser.add_argument("--no-preserve-all-docking-outputs", dest="preserve_all_docking_outputs", action="store_false")
    m3_cleanup_parser.add_argument("--confirm-run-id", default=None)
    m3_cleanup_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_cleanup_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_cleanup_parser.set_defaults(func=_cmd_m3_cleanup)

    m3_report_cleanup_parser = subparsers.add_parser(
        "m3-report-cleanup",
        help="Run M3-T11 report generation plus cleanup plan or confirmed cleanup execution.",
    )
    m3_report_cleanup_parser.add_argument("--run-id", required=True)
    m3_report_cleanup_parser.add_argument("--m2-run-id", required=True)
    m3_report_cleanup_parser.add_argument("--profile", choices=["smoke", "mini", "scaling", "production"], default=None)
    m3_report_cleanup_parser.add_argument("--mode", default="report", choices=["dry-run", "report"])
    m3_report_cleanup_parser.add_argument("--cleanup-mode", default="plan", choices=["none", "plan", "execute"])
    m3_report_cleanup_parser.add_argument("--cleanup-scope", default="conservative", choices=["test", "mini", "production", "conservative", "aggressive-test-only"])
    m3_report_cleanup_parser.add_argument("--allow-partial", action="store_true")
    m3_report_cleanup_parser.add_argument("--require-t10-ready", dest="require_t10_ready", action="store_true", default=True)
    m3_report_cleanup_parser.add_argument("--no-require-t10-ready", dest="require_t10_ready", action="store_false")
    m3_report_cleanup_parser.add_argument("--strict-confidentiality", dest="strict_confidentiality", action="store_true", default=True)
    m3_report_cleanup_parser.add_argument("--no-strict-confidentiality", dest="strict_confidentiality", action="store_false")
    m3_report_cleanup_parser.add_argument("--strict-claims", dest="strict_claims", action="store_true", default=True)
    m3_report_cleanup_parser.add_argument("--no-strict-claims", dest="strict_claims", action="store_false")
    m3_report_cleanup_parser.add_argument("--config", default=None)
    m3_report_cleanup_parser.add_argument("--private-map", default="fresh/data/private/compound_id_map.csv")
    m3_report_cleanup_parser.add_argument("--confirm-run-id", default=None)
    m3_report_cleanup_parser.add_argument("--force", dest="force", action="store_true", default=False)
    m3_report_cleanup_parser.add_argument("--no-force", dest="force", action="store_false")
    m3_report_cleanup_parser.add_argument("--preserve-selected-poses", dest="preserve_selected_poses", action="store_true", default=True)
    m3_report_cleanup_parser.add_argument("--no-preserve-selected-poses", dest="preserve_selected_poses", action="store_false")
    m3_report_cleanup_parser.add_argument("--preserve-all-docking-outputs", dest="preserve_all_docking_outputs", action="store_true", default=True)
    m3_report_cleanup_parser.add_argument("--no-preserve-all-docking-outputs", dest="preserve_all_docking_outputs", action="store_false")
    m3_report_cleanup_parser.set_defaults(func=_cmd_m3_report_cleanup)

    ligand_parser = subparsers.add_parser(
        "manifest-ligands",
        help="Build ligand manifest shell (public IDs only; gitignored private mapping).",
    )
    ligand_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/.",
    )
    ligand_parser.add_argument(
        "--ligands-dir",
        default=None,
        help="Directory containing public-ID SDF files. Default: paths.yaml raw_ligands.",
    )
    ligand_parser.add_argument(
        "--private-mapping",
        default=None,
        help="Private mapping CSV path. Default: paths.yaml private_data/compound_id_map.csv.",
    )
    ligand_parser.add_argument(
        "--mode",
        default="smoke_env",
        choices=["smoke_env", "smoke_input"],
    )
    ligand_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    ligand_parser.add_argument(
        "--compound-stage-enabled",
        choices=["true", "false"],
        default="false",
        help="Whether compound docking stage is enabled (affects severity of missing files).",
    )
    ligand_parser.set_defaults(func=_cmd_manifest_ligands)

    pbs_parser = subparsers.add_parser(
        "prepare-pbs",
        help="Generate a concrete PBS job file under runs/<run_id>/scripts/. Does NOT call qsub.",
    )
    pbs_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/.",
    )
    pbs_parser.add_argument(
        "--job-name",
        required=True,
        help="PBS job name (alnum/underscore/dot/dash).",
    )
    pbs_parser.add_argument(
        "--mode",
        required=True,
        choices=["smoke_env", "smoke_input", "mini", "scaling", "production"],
        help="Mode label; resolves ppn/walltime defaults from hpc.yaml.",
    )
    pbs_parser.add_argument(
        "--node",
        default=None,
        choices=["node04", "node05", "node06"],
        help="Target node. Default: first entry in hpc.yaml nodes.",
    )
    pbs_parser.add_argument(
        "--ppn",
        type=int,
        default=None,
        help="Override ppn. Default: hpc.yaml ppn[mode].",
    )
    pbs_parser.add_argument(
        "--walltime",
        default=None,
        help="Override walltime as HH:MM:SS. Default: hpc.yaml walltime[mode].",
    )
    pbs_parser.add_argument(
        "--output-path",
        default=None,
        help="Optional output path. Default: runs/<run_id>/scripts/<job_name>.pbs.",
    )
    pbs_parser.add_argument(
        "--input-root",
        default="fresh/data/raw",
        help="--input-root used by smoke_input mode body.",
    )
    pbs_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
    )
    pbs_parser.set_defaults(func=_cmd_prepare_pbs)

    membrane_parser = subparsers.add_parser(
        "compute-membrane-frame",
        help="Compute state-aware membrane_frame.json from coordinates (no hardcoded vectors).",
    )
    membrane_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/.",
    )
    membrane_parser.add_argument(
        "--state",
        default="all",
        choices=["EGFR_160-185", "EGFR_170-200", "3GT8_raw", "all"],
        help="State to compute. Default 'all' iterates all three states.",
    )
    membrane_parser.add_argument(
        "--full-frame-source",
        default=None,
        help="Optional fallback PDB path; default uses plus10_full_frame.pdb from receptor_states.yaml.",
    )
    membrane_parser.add_argument(
        "--mode",
        default="smoke_env",
        choices=["smoke_env", "smoke_input"],
    )
    membrane_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev tolerates missing source (WARN); hpc_strict treats missing as FAIL.",
    )
    membrane_parser.set_defaults(func=_cmd_compute_membrane_frame)

    receptor_parser = subparsers.add_parser(
        "prepare-receptor",
        help="Normalize EGFR receptor input: explicit A/B split, 669-1014 dockable crop, +1000 runtime offset, mapping CSV.",
    )
    receptor_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/.",
    )
    receptor_parser.add_argument(
        "--state",
        required=True,
        choices=["EGFR_160-185", "EGFR_170-200", "3GT8_raw"],
        help="Receptor state ID. 3GT8_raw is a reference/control, not a primary membrane-validated state.",
    )
    receptor_parser.add_argument(
        "--source",
        required=True,
        help="Path to source receptor PDB.",
    )
    receptor_parser.add_argument(
        "--mode",
        default="smoke_env",
        choices=["smoke_env", "smoke_input"],
        help="Mode label recorded in manifests.",
    )
    receptor_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev tolerates fixture warnings; hpc_strict treats them as FAIL/quarantine.",
    )
    receptor_parser.add_argument(
        "--strict",
        action="store_true",
        help="Force production-like blocker policy even in codex_dev.",
    )
    receptor_parser.set_defaults(func=_cmd_prepare_receptor)

    myo1d_parser = subparsers.add_parser(
        "prepare-myo1d",
        help="Slice MYO1D source PDB to canonical M1 construct and emit QC.",
    )
    myo1d_parser.add_argument(
        "--run-id",
        required=True,
        help="Run identifier under fresh/runs/.",
    )
    myo1d_parser.add_argument(
        "--source",
        required=True,
        help="Path to MYO1D source PDB (real input or fixture).",
    )
    myo1d_parser.add_argument(
        "--construct",
        default=None,
        help="Residue range as 'start-end'. Defaults to gates.yaml myo1d.construct (955-1006).",
    )
    myo1d_parser.add_argument(
        "--mode",
        default="smoke_env",
        choices=["smoke_env", "smoke_input"],
        help="Mode label recorded in manifests.",
    )
    myo1d_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="codex_dev tolerates missing source (WARN); hpc_strict treats as FAIL.",
    )
    myo1d_parser.set_defaults(func=_cmd_prepare_myo1d)

    cleanup_parser = subparsers.add_parser(
        "cleanup",
        help="Delete intermediate/scratch files inside a run directory; preserves manifest/logs/qc/reports.",
    )
    cleanup_parser.add_argument(
        "--run-id",
        required=True,
        help="Existing run identifier under fresh/runs/.",
    )
    cleanup_parser.add_argument(
        "--mode",
        required=True,
        choices=["test", "production"],
        help="test deletes intermediates; production defaults to dry-run.",
    )
    cleanup_parser.add_argument(
        "--dry-run",
        dest="dry_run",
        choices=["true", "false"],
        default=None,
        help=(
            "Override default. Default is 'false' for --mode test and 'true' "
            "for --mode production. Pass 'true' or 'false' explicitly to override."
        ),
    )
    cleanup_parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="Profile label recorded in cleanup_report.json; no external tools are run.",
    )
    cleanup_parser.set_defaults(func=_cmd_cleanup)

    return parser


def _cmd_version(_args):
    print("egfr-myo1d-fresh {0} ({1})".format(__version__, WORKFLOW_STAGE))
    return 0


def _cmd_init_run(args):
    from egfr_myo1d.core.logging_utils import append_phase_status, initialize_logs
    from egfr_myo1d.core.manifest import initialize_manifests
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    initialize_manifests(ctx, args.mode)
    append_phase_status(
        ctx,
        phase="init_run",
        status="PASS",
        message="run context initialized",
        details={"mode": args.mode, "run_dir": str(ctx.run_dir)},
    )
    print("initialized run_id={0}".format(ctx.run_id))
    print("run_dir={0}".format(ctx.run_dir))
    return 0


def _cmd_preflight(args):
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.preflight import run_preflight

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = run_preflight(ctx, args.mode, args.profile)
    status = report["status"]
    counts = report.get("counts", {})
    print(
        "preflight {0}: PASS={1} WARN={2} FAIL={3}".format(
            status, counts.get("PASS", 0), counts.get("WARN", 0), counts.get("FAIL", 0)
        )
    )
    print("environment_report={0}".format(ctx.manifest_dir / "environment_report.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_tool_preflight(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.tools.tool_preflight import run_tool_preflight

    ctx = RunContext.create(args.run_id, "smoke_env")
    initialize_logs(ctx)
    report = run_tool_preflight(
        ctx,
        mode=args.mode,
        profile=args.profile,
        registry_path=Path(args.registry) if args.registry else None,
        timeout=args.timeout,
    )
    counts = report.get("counts", {})
    print(
        "tool-preflight {0}: tools={1} available_or_passed={2} not_installed={3} smoke_failed={4} blockers={5}".format(
            report["status"],
            counts.get("total", 0),
            counts.get("available_or_passed", 0),
            counts.get("not_installed", 0),
            counts.get("smoke_failed", 0),
            counts.get("blockers", 0),
        )
    )
    print("tool_status={0}".format(ctx.manifest_dir / "tool_status.json"))
    print("tool_report={0}".format(ctx.reports_dir / "tool_installation_report.md"))
    if report["status"] == "FAIL":
        return 1
    return 0


def _load_json(path):
    if not os.path.exists(str(path)):
        return None
    with open(str(path), "r") as handle:
        return json.load(handle)


def _read_phase_statuses(path):
    statuses = []
    if not os.path.exists(str(path)):
        return statuses
    with open(str(path), "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                statuses.append(json.loads(line))
            except ValueError:
                statuses.append({"status": "FAIL", "message": "invalid JSONL line"})
    return statuses


def _cmd_status(args):
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.for_existing(args.run_id)
    manifest_files = [
        ctx.manifest_dir / "run_manifest.json",
        ctx.manifest_dir / "input_manifest.json",
        ctx.manifest_dir / "environment_report.json",
        ctx.manifest_dir / "git_snapshot.json",
    ]
    statuses = _read_phase_statuses(ctx.logs_dir / "phase_status.jsonl")
    warn_count = sum(1 for record in statuses if record.get("status") == "WARN")
    fail_count = sum(1 for record in statuses if record.get("status") == "FAIL")
    last = statuses[-1] if statuses else None

    print("run_id={0}".format(ctx.run_id))
    print("run_dir={0}".format(ctx.run_dir))
    for path in manifest_files:
        label = "present" if os.path.exists(str(path)) else "missing"
        print("manifest {0}: {1}".format(path.name, label))
    if last:
        print(
            "last_phase={0} {1}: {2}".format(
                last.get("phase"), last.get("status"), last.get("message")
            )
        )
    else:
        print("last_phase=none")
    print("phase_WARN_count={0}".format(warn_count))
    print("phase_FAIL_count={0}".format(fail_count))
    print("master_log={0}".format(ctx.logs_dir / "master.log"))
    print("error_summary={0}".format(ctx.errors_dir / "error_summary.txt"))
    return 0


def _cmd_validate_structures(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.structure_inputs import validate_structure_inputs

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = validate_structure_inputs(ctx, args.mode, args.profile, Path(args.input_root))
    status = report["status"]
    print("structure validation {0}".format(status))
    print("structure_qc_report={0}".format(ctx.manifest_dir / "structure_qc_report.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_prepare_ppi_inputs(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.prepared_inputs import prepare_ppi_inputs

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = prepare_ppi_inputs(
        ctx,
        args.mode,
        args.profile,
        Path(args.input_root),
        Path(args.contract) if args.contract else None,
        strict=args.strict,
    )
    status = report["status"]
    print("prepare-ppi-inputs {0}".format(status))
    print("preparation_qc_report={0}".format(ctx.manifest_dir / "preparation_qc_report.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_validate_real_inputs(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.real_inputs import validate_real_inputs

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = validate_real_inputs(
        ctx,
        args.mode,
        args.profile,
        Path(args.input_root),
        Path(args.contract) if args.contract else None,
        strict=args.strict,
    )
    status = report["status"]
    print("validate-real-inputs {0}".format(status))
    print("real_input_readiness_report={0}".format(ctx.manifest_dir / "real_input_readiness_report.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_plan_ppi_sampling(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.ppi_sampling_plan import plan_ppi_sampling

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = plan_ppi_sampling(
        ctx,
        args.mode,
        args.profile,
        Path(args.input_root),
        Path(args.contract) if args.contract else None,
        strict=args.strict,
        command_line=" ".join(sys.argv),
    )
    status = report["status"]
    print("plan-ppi-sampling {0}".format(status))
    print("ppi_sampling_plan_report={0}".format(ctx.manifest_dir / "ppi_sampling_plan_report.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_summarize_ppi_consensus(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.ppi_consensus import summarize_ppi_consensus

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = summarize_ppi_consensus(
        ctx,
        args.mode,
        args.profile,
        Path(args.input_root),
        Path(args.contact_table) if args.contact_table else None,
        input_kind=args.input_kind,
    )
    status = report["status"]
    print("summarize-ppi-consensus {0}".format(status))
    print("ppi_consensus_qc_report={0}".format(ctx.manifest_dir / "ppi_consensus_qc_report.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_plan_pocket_discovery(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.pocket_discovery import plan_pocket_discovery

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = plan_pocket_discovery(
        ctx,
        args.mode,
        args.profile,
        Path(args.input_root),
        Path(args.ppi_consensus) if args.ppi_consensus else None,
        input_kind=args.input_kind,
    )
    status = report["status"]
    print("plan-pocket-discovery {0}".format(status))
    print("pocket_discovery_manifest={0}".format(ctx.manifest_dir / "pocket_discovery_manifest.json"))
    return 0


def _cmd_prioritize_pocket_candidates(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.validation.pocket_candidate_prioritization import prioritize_pocket_candidates

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    report = prioritize_pocket_candidates(
        ctx,
        args.mode,
        args.profile,
        Path(args.input_root),
        Path(args.pocket_plan) if args.pocket_plan else None,
        Path(args.candidate_pockets) if args.candidate_pockets else None,
        Path(args.ppi_consensus) if args.ppi_consensus else None,
        strict=args.strict,
    )
    status = report["status"]
    print("prioritize-pocket-candidates {0}".format(status))
    print("pocket_candidate_prioritization_manifest={0}".format(ctx.manifest_dir / "pocket_candidate_prioritization_manifest.json"))
    if status == "FAIL":
        return 1
    return 0


def _cmd_prepare_inputs(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.orchestrator.prepare_inputs import run_prepare_inputs

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)

    state_list = (
        [s.strip() for s in args.states.split(",") if s.strip()]
        if args.states
        else None
    )
    aggregate = run_prepare_inputs(
        ctx,
        mode=args.mode,
        profile=args.profile,
        input_root=Path(args.input_root) if args.input_root else None,
        states=state_list,
        skip_ligands=(args.skip_ligands == "true"),
        strict=args.strict,
        compound_stage_enabled=(args.compound_stage_enabled == "true"),
    )
    print(
        "prepare-inputs {0}: substeps={1} missing={2} blockers={3}".format(
            aggregate.status,
            len(aggregate.sub_steps),
            len(aggregate.missing_required_inputs),
            len(aggregate.blockers),
        )
    )
    print("aggregate_manifest={0}".format(aggregate.aggregate_manifest_path))
    print("summary_report={0}".format(aggregate.summary_report_path))
    if aggregate.status == "FAIL":
        return 1
    return 0


def _cmd_generate_m2_ppi_inputs(args):
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.m2.ppi_inputs import generate_m2_1_ppi_inputs

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    state_list = (
        [state.strip() for state in args.states.split(",") if state.strip()]
        if args.states
        else None
    )
    report = generate_m2_1_ppi_inputs(
        ctx,
        mode=args.mode,
        profile=args.profile,
        states=state_list,
        strict=args.strict,
    )
    print(
        "generate-m2-ppi-inputs {0}: packs={1} warnings={2} blockers={3}".format(
            report.status,
            len(report.packs),
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_1_manifest={0}".format(report.manifest_path))
    print("m2_1_qc_csv={0}".format(report.qc_csv_path))
    print("m2_1_summary={0}".format(report.summary_report_path))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_prepare_m2_pyrosetta_harness(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.ppi.pyrosetta_adapter import generate_pyrosetta_harness

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)
    state_list = (
        [state.strip() for state in args.states.split(",") if state.strip()]
        if args.states
        else None
    )
    report = generate_pyrosetta_harness(
        ctx,
        mode=args.mode,
        profile=args.profile,
        states=state_list,
        m2_1_manifest_path=Path(args.m2_1_manifest)
        if args.m2_1_manifest
        else None,
    )
    print(
        "prepare-m2-pyrosetta-harness {0}: jobs={1} warnings={2} blockers={3}".format(
            report.status,
            len(report.jobs),
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_2_manifest={0}".format(report.manifest_path))
    print("m2_2_job_manifest_csv={0}".format(report.job_manifest_csv))
    print("m2_2_launch_script={0}".format(report.launch_script_path))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_collect_m2_ppi_outputs(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.ppi.collect_ppi_outputs import collect_ppi_outputs

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = collect_ppi_outputs(
        ctx,
        mode=args.mode,
        profile=args.profile,
        job_manifest=Path(args.job_manifest) if args.job_manifest else None,
        raw_contact_table=Path(args.raw_contact_table)
        if args.raw_contact_table
        else None,
        unmapped_fraction_fail_threshold=args.unmapped_fraction_fail_threshold,
    )
    print(
        "collect-m2-ppi-outputs {0}: raw_pose_rows={1} restored_contacts={2} unmapped_contacts={3}".format(
            report.status,
            report.raw_pose_count,
            report.restored_contact_count,
            report.unmapped_contact_count,
        )
    )
    print("m2_3_manifest={0}".format(report.manifest_path))
    print("m2_3_raw_pose_table={0}".format(report.raw_pose_table))
    print("m2_3_pose_contacts={0}".format(report.pose_contacts_csv))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_build_m2_ppi_consensus_patch(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.ppi.consensus_patch import build_m2_4_consensus_patch

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = build_m2_4_consensus_patch(
        ctx,
        mode=args.mode,
        profile=args.profile,
        contact_table=Path(args.contact_table) if args.contact_table else None,
        job_manifest=Path(args.job_manifest) if args.job_manifest else None,
    )
    print(
        "build-m2-ppi-consensus-patch {0}: contacts={1} patches={2} warnings={3} blockers={4}".format(
            report.status,
            report.contact_count,
            report.patch_count,
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_4_manifest={0}".format(report.manifest_path))
    print("ppi_consensus_patch={0}".format(report.consensus_patch_csv))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_build_m2_atp_reference(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.pocket.atp_reference import build_atp_reference

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    states = (
        [state.strip() for state in args.states.split(",") if state.strip()]
        if args.states
        else None
    )
    report = build_atp_reference(
        ctx,
        mode=args.mode,
        profile=args.profile,
        states=states,
        include_reference_states=(args.include_reference_states == "true"),
        config_path=Path(args.config) if args.config else None,
        reference_pdb=Path(args.reference_pdb) if args.reference_pdb else None,
        synthetic_fixture=(args.synthetic_fixture == "true"),
    )
    print(
        "build-m2-atp-reference {0}: references={1} mappings={2} centroids={3} warnings={4} blockers={5}".format(
            report.status,
            report.reference_count,
            report.mapping_count,
            report.centroid_count,
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_5_status={0}".format(report.status_json))
    print("atp_site_reference={0}".format(report.reference_csv))
    print("atp_site_residue_mapping={0}".format(report.residue_mapping_csv))
    print("atp_site_centroid_by_state={0}".format(report.centroid_csv))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_run_m2_fpocket_discovery(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.pocket.fpocket_adapter import discover_and_normalize_pockets

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    states = (
        [state.strip() for state in args.states.split(",") if state.strip()]
        if args.states
        else None
    )
    report = discover_and_normalize_pockets(
        ctx,
        mode=args.mode,
        profile=args.profile,
        states=states,
        include_reference_states=(args.include_reference_states == "true"),
        execution_mode=args.execution_mode,
        fpocket_output_root=Path(args.fpocket_output_root)
        if args.fpocket_output_root
        else None,
        config_path=Path(args.config) if args.config else None,
        fpocket_binary=args.fpocket_binary,
        synthetic_fixture=(args.synthetic_fixture == "true"),
    )
    print(
        "run-m2-fpocket-discovery {0}: raw_pockets={1} raw_candidates={2} families={3} warnings={4} blockers={5}".format(
            report.status,
            report.raw_pocket_count,
            report.raw_candidate_count,
            report.merged_family_count,
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_6_raw_status={0}".format(report.raw_status_json))
    print("pocket_candidates_raw={0}".format(report.raw_candidate_csv))
    print("pocket_candidates_merged={0}".format(report.merged_candidate_csv))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_gate_m2_pockets(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.pocket.pocket_gate_apply import apply_pocket_gates

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = apply_pocket_gates(
        ctx,
        mode=args.mode,
        profile=args.profile,
        config_path=Path(args.config) if args.config else None,
        ppi_consensus=Path(args.ppi_consensus) if args.ppi_consensus else None,
        atp_reference=Path(args.atp_reference) if args.atp_reference else None,
        atp_centroids=Path(args.atp_centroids) if args.atp_centroids else None,
        pocket_families=Path(args.pocket_families) if args.pocket_families else None,
        raw_candidates=Path(args.raw_candidates) if args.raw_candidates else None,
        membrane_frame=Path(args.membrane_frame) if args.membrane_frame else None,
        synthetic_fixture=(args.synthetic_fixture == "true"),
    )
    print(
        "gate-m2-pockets {0}: families={1} accepted_primary={2} accepted_secondary={3} reference_only={4} rejected={5} warnings={6} blockers={7}".format(
            report.status,
            report.family_count,
            report.accepted_primary_count,
            report.accepted_secondary_count,
            report.reference_only_count,
            report.atp_reject_count
            + report.ppi_reject_count
            + report.membrane_reject_count
            + report.dimer_reject_count
            + report.mapping_or_origin_reject_count,
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_7_status={0}".format(report.status_json))
    print("pocket_gate_qc={0}".format(report.gate_qc_csv))
    print("accepted_pocket_families={0}".format(report.accepted_csv))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_export_m2_results(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.pocket.m2_export import export_m2_results

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = export_m2_results(
        ctx,
        mode=args.mode,
        profile=args.profile,
        ppi_consensus=Path(args.ppi_consensus) if args.ppi_consensus else None,
        atp_reference=Path(args.atp_reference) if args.atp_reference else None,
        pocket_families=Path(args.pocket_families) if args.pocket_families else None,
        pocket_gate_qc=Path(args.pocket_gate_qc) if args.pocket_gate_qc else None,
        accepted_families=Path(args.accepted_families) if args.accepted_families else None,
        rejected_families=Path(args.rejected_families) if args.rejected_families else None,
        synthetic_fixture=(args.synthetic_fixture == "true"),
    )
    print(
        "export-m2-results {0}: total_families={1} exported_primary={2} exported_secondary={3} reference_only={4} m3_allowed={5} warnings={6} blockers={7}".format(
            report.status,
            report.total_family_count,
            report.exported_primary_count,
            report.exported_secondary_count,
            report.reference_only_count,
            report.m3_docking_allowed,
            len(report.warnings),
            len(report.blockers),
        )
    )
    print("m2_8_status={0}".format(report.final_status_json))
    print("accepted_pockets_for_m3={0}".format(report.export_accepted_csv))
    print("m2_final_report={0}".format(report.export_report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_readiness(args):
    from egfr_myo1d.compound.m3_readiness import run_m3_readiness
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_readiness(
        ctx,
        m2_run_id=args.m2_run_id,
        mode=args.mode,
        create_skeleton=args.create_skeleton,
    )
    print(
        "m3-readiness {0}: m3_allowed={1} blockers={2} warnings={3} docking_ready_pockets={4} valid_boxes={5}".format(
            report.status,
            report.m3_docking_allowed,
            len(report.blockers),
            len(report.warnings),
            report.counts.get("docking_ready_pocket_rows", 0),
            report.counts.get("valid_box_rows", 0),
        )
    )
    print("m3_readiness_report={0}".format(report.report_csv))
    print("m3_readiness_summary={0}".format(report.summary_json))
    print("m3_task0_readiness={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_ligand_qc(args):
    from pathlib import Path

    from egfr_myo1d.compound.ligand_manifest import run_m3_ligand_qc
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = run_m3_ligand_qc(
        ctx,
        mode=args.mode,
        private_map=Path(args.private_map) if args.private_map else None,
        update_ignore=args.update_gitignore,
    )
    print(
        "m3-ligand-qc {0}: m3_t2_allowed={1} ligands_found={2} missing={3} failed={4} leaks={5} tracked_sensitive={6}".format(
            report.status,
            report.m3_t2_allowed,
            report.counts.get("ligands_found", 0),
            report.counts.get("ligands_missing", 0),
            report.counts.get("ligands_failed_qc", 0),
            report.counts.get("confidentiality_leaks", 0),
            report.counts.get("tracked_sensitive_files", 0),
        )
    )
    print("ligand_manifest_public={0}".format(report.manifest_csv))
    print("ligand_qc_summary={0}".format(report.qc_csv))
    print("confidentiality_audit={0}".format(report.audit_csv))
    print("m3_ligand_qc_summary={0}".format(report.summary_json))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_prepare_ligands(args):
    from pathlib import Path

    from egfr_myo1d.compound.ligand_prepare import run_m3_ligand_prepare
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = run_m3_ligand_prepare(
        ctx,
        mode=args.mode,
        private_map=Path(args.private_map) if args.private_map else None,
        force=args.force,
        keep_intermediates=args.keep_intermediates,
        update_ignore=args.update_gitignore,
        allow_independent=args.allow_independent,
        require_m3_t1_pass=args.require_m3_t1_pass,
    )
    print(
        "m3-prepare-ligands {0}: m3_t3_allowed={1} prepared={2} reused_supplied_pdbqt={3} failed={4} warnings={5}".format(
            report.status,
            report.m3_t3_allowed,
            report.counts.get("ligands_prepared", 0),
            report.counts.get("ligands_reused_supplied_pdbqt", 0),
            report.counts.get("ligands_failed", 0),
            len(report.warnings),
        )
    )
    print("ligand_preparation_manifest={0}".format(report.preparation_manifest_csv))
    print("prepared_ligand_hashes={0}".format(report.prepared_hashes_csv))
    print("ligand_prep_qc={0}".format(report.prep_qc_csv))
    print("m3_ligand_prep_summary={0}".format(report.summary_json))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_receptor_boxes(args):
    from pathlib import Path

    from egfr_myo1d.compound.receptor_prepare import run_m3_receptor_boxes
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    state_list = [state.strip() for state in args.states.split(",") if state.strip()] if args.states else None
    report = run_m3_receptor_boxes(
        ctx,
        m2_run_id=args.m2_run_id,
        mode=args.mode,
        force=args.force,
        box_margin=args.box_margin,
        states=state_list,
        update_ignore=args.update_gitignore,
        allow_independent=args.allow_independent,
        require_m3_t2_pass=args.require_m3_t2_pass,
        private_map=Path(args.private_map) if args.private_map else None,
        operation=args.operation,
    )
    print(
        "{0} {1}: m3_t4_allowed={2} states_prepared={3} boxes_pass={4} boxes_fail={5} warnings={6}".format(
            args.operation,
            report.status,
            report.m3_t4_allowed,
            report.counts.get("states_prepared", 0),
            report.counts.get("boxes_pass", 0),
            report.counts.get("boxes_fail", 0),
            len(report.warnings),
        )
    )
    print("receptor_preparation_manifest={0}".format(report.receptor_preparation_manifest_csv))
    print("receptor_pdbqt_hashes={0}".format(report.receptor_hashes_csv))
    print("docking_box_manifest={0}".format(report.docking_box_manifest_csv))
    print("m3_receptor_box_summary={0}".format(report.summary_json))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_vina_smoke(args):
    from pathlib import Path

    from egfr_myo1d.compound.vina_smoke import run_m3_vina_smoke
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_vina_smoke(
        ctx,
        m2_run_id=args.m2_run_id,
        mode=args.mode,
        force=args.force,
        compound_public_id=args.compound_public_id,
        state_id=args.state_id,
        pocket_family_id=args.pocket_family_id,
        box_id=args.box_id,
        protomer_id=args.protomer_id,
        vina_executable=args.vina_executable,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        cpu=args.cpu,
        seed=args.seed,
        timeout_seconds=args.timeout_seconds,
        require_m3_t2_pass=args.require_m3_t2_pass,
        require_m3_t3_pass=args.require_m3_t3_pass,
        allow_independent=args.allow_independent,
        allow_reference_smoke=args.allow_reference_smoke,
        update_ignore=args.update_gitignore,
        private_map=Path(args.private_map) if args.private_map else None,
    )
    print(
        "m3-vina-smoke {0}: m3_t5_allowed={1} compound={2} state={3} box={4} invoked={5} exit_code={6}".format(
            report.status,
            report.m3_t5_allowed,
            report.selection.get("compound_public_id"),
            report.selection.get("state_id"),
            report.selection.get("box_id"),
            report.vina.get("invoked"),
            report.vina.get("exit_code"),
        )
    )
    print("vina_smoke_manifest={0}".format(report.manifest_csv))
    print("vina_smoke_qc_json={0}".format(report.qc_json))
    print("vina_smoke_qc_csv={0}".format(report.qc_csv))
    print("m3_task4_vina_smoke={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_vina_jobs(args):
    from pathlib import Path

    from egfr_myo1d.compound.vina_jobs import run_m3_vina_jobs
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    node_list = [item.strip() for item in args.nodes.split(",") if item.strip()] if args.nodes else None
    report = run_m3_vina_jobs(
        ctx,
        m2_run_id=args.m2_run_id,
        mode=args.mode,
        profile=args.profile,
        force=args.force,
        compound_public_id=args.compound_public_id,
        state_id=args.state_id,
        pocket_family_id=args.pocket_family_id,
        box_id=args.box_id,
        protomer_id=args.protomer_id,
        include_reference=args.include_reference,
        include_3gt8=args.include_3gt8,
        allow_reference_jobs=args.allow_reference_jobs,
        allow_independent=args.allow_independent,
        allow_production_plan_without_mini=args.allow_production_plan_without_mini,
        vina_executable=args.vina_executable,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        cpu_per_vina=args.cpu_per_vina,
        repeats=args.repeats,
        base_seed=args.base_seed,
        timeout_seconds=args.timeout_seconds,
        nodes=node_list,
        ppn=args.ppn,
        walltime=args.walltime,
        queue=args.queue,
        conda_env=args.conda_env,
        conda_sh=args.conda_sh,
        python_executable=args.python_executable,
        update_ignore=args.update_gitignore,
        private_map=Path(args.private_map) if args.private_map else None,
    )
    print(
        "m3-vina-jobs {0}: profile={1} qsub_submission_allowed={2} production_submission_allowed={3} planned_vina_commands={4} chunks={5} pbs_files={6}".format(
            report.status,
            args.profile,
            report.qsub_submission_allowed,
            report.production_submission_allowed,
            report.counts.get("planned_vina_commands", 0),
            report.counts.get("planned_chunks", 0),
            report.counts.get("pbs_files", 0),
        )
    )
    print("vina_job_manifest={0}".format(report.job_manifest_csv))
    print("vina_profile_job_manifest={0}".format(report.profile_job_manifest_csv))
    print("vina_pbs_manifest={0}".format(report.pbs_manifest_csv))
    print("vina_job_plan_qc={0}".format(report.qc_json))
    print("m3_task5_vina_job_plan={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_run_vina_chunk(args):
    from pathlib import Path

    from egfr_myo1d.compound.vina_job_runner import run_vina_chunk
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)
    report = run_vina_chunk(
        ctx,
        job_manifest=Path(args.job_manifest),
        chunk_id=args.chunk_id,
        max_workers=args.max_workers,
        dry_run_runner=args.dry_run_runner,
    )
    print(
        "m3-run-vina-chunk {0}: chunk_id={1} rows_seen={2} rows_run={3} rows_failed={4} dry_run_runner={5}".format(
            report.status,
            report.chunk_id,
            report.rows_seen,
            report.rows_run,
            report.rows_failed,
            report.dry_run_runner,
        )
    )
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_vina_collect(args):
    from pathlib import Path

    from egfr_myo1d.compound.vina_collect import run_m3_vina_collect
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_vina_collect(
        ctx,
        m2_run_id=args.m2_run_id,
        job_manifest=Path(args.job_manifest) if args.job_manifest else None,
        profile=args.profile,
        mode=args.mode,
        force=args.force,
        allow_partial=args.allow_partial,
        require_all_jobs=args.require_all_jobs,
        allow_reference_only=args.allow_reference_only,
        private_map=Path(args.private_map) if args.private_map else None,
        update_ignore=args.update_gitignore,
        max_log_bytes=args.max_log_bytes,
        strict_log_parse=args.strict_log_parse,
        write_empty_raw_table=args.write_empty_raw_table,
    )
    print(
        "m3-vina-collect {0}: profile={1} expected_jobs={2} completed={3} parsed_pose_rows={4} m3_t7_attribution_ready={5}".format(
            report.status,
            report.profiles.get("selected_profile"),
            report.collection.get("expected_jobs", 0),
            report.collection.get("completed_jobs", 0),
            report.collection.get("parsed_pose_rows", 0),
            report.m3_t7_attribution_ready,
        )
    )
    print("compound_pose_raw={0}".format(report.raw_pose_csv))
    print("vina_job_completion={0}".format(report.completion_csv))
    print("docking_completion_qc_json={0}".format(report.qc_json))
    print("docking_completion_qc_csv={0}".format(report.qc_csv))
    print("m3_task6_vina_collection={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_pose_attribution(args):
    from pathlib import Path

    from egfr_myo1d.compound.pose_attribution import run_m3_pose_attribution
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_pose_attribution(
        ctx,
        m2_run_id=args.m2_run_id,
        profile=args.profile,
        mode=args.mode,
        pose_raw_table=Path(args.pose_raw_table) if args.pose_raw_table else None,
        job_completion_table=Path(args.job_completion_table) if args.job_completion_table else None,
        accepted_pockets=Path(args.accepted_pockets) if args.accepted_pockets else None,
        accepted_boxes=Path(args.accepted_boxes) if args.accepted_boxes else None,
        atp_reference=Path(args.atp_reference) if args.atp_reference else None,
        ppi_consensus_patch=Path(args.ppi_consensus_patch) if args.ppi_consensus_patch else None,
        membrane_frame=Path(args.membrane_frame) if args.membrane_frame else None,
        receptor_mapping=Path(args.receptor_mapping) if args.receptor_mapping else None,
        receptor_manifest=Path(args.receptor_manifest) if args.receptor_manifest else None,
        config=Path(args.config) if args.config else None,
        force=args.force,
        allow_partial=args.allow_partial,
        require_all_poses=args.require_all_poses,
        allow_reference_only=args.allow_reference_only,
        private_map=Path(args.private_map) if args.private_map else None,
        strict_geometry=args.strict_geometry,
        write_empty_attribution=args.write_empty_attribution,
    )
    print(
        "m3-pose-attribution {0}: profile={1} poses_classified={2} hard_gate_pass={3} rejected={4} m3_t8_clustering_ready={5}".format(
            report.status,
            report.collection_context.get("selected_profile"),
            report.attribution.get("poses_classified", 0),
            report.attribution.get("poses_pass", 0),
            report.attribution.get("poses_rejected", 0),
            report.m3_t8_clustering_ready,
        )
    )
    print("compound_pose_attribution={0}".format(report.attribution_csv))
    print("compound_pose_mechanism_classification={0}".format(report.mechanism_csv))
    print("pose_attribution_qc_json={0}".format(report.qc_json))
    print("pose_attribution_qc_csv={0}".format(report.qc_csv))
    print("m3_task7_pose_attribution={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_pose_clustering(args):
    from pathlib import Path

    from egfr_myo1d.compound.pose_clustering import run_m3_pose_clustering
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_pose_clustering(
        ctx,
        m2_run_id=args.m2_run_id,
        profile=args.profile,
        mode=args.mode,
        pose_attribution_table=Path(args.pose_attribution_table) if args.pose_attribution_table else None,
        pose_attribution_qc=Path(args.pose_attribution_qc) if args.pose_attribution_qc else None,
        accepted_boxes=Path(args.accepted_boxes) if args.accepted_boxes else None,
        config=Path(args.config) if args.config else None,
        force=args.force,
        allow_partial=args.allow_partial,
        require_attribution_ready=args.require_attribution_ready,
        allow_reference_only=args.allow_reference_only,
        private_map=Path(args.private_map) if args.private_map else None,
        cluster_method=args.cluster_method,
        strict_convergence=args.strict_convergence,
        write_empty_clusters=args.write_empty_clusters,
        write_cluster_membership=args.write_cluster_membership,
    )
    print(
        "m3-pose-clustering {0}: profile={1} eligible={2} clusters={3} converged={4} m3_t9_anchor_convergence_ready={5}".format(
            report.status,
            report.attribution_context.get("selected_profile"),
            report.attribution_context.get("eligible_pose_rows", 0),
            report.clustering.get("clusters_written", 0),
            report.clustering.get("converged_clusters", 0),
            report.m3_t9_anchor_convergence_ready,
        )
    )
    print("compound_pose_clusters={0}".format(report.clusters_csv))
    print("compound_pose_cluster_membership={0}".format(report.membership_csv))
    print("pose_clustering_qc_json={0}".format(report.qc_json))
    print("pose_clustering_qc_csv={0}".format(report.qc_csv))
    print("m3_task8_pose_clustering={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_anchor_convergence(args):
    from pathlib import Path

    from egfr_myo1d.compound.anchor_convergence import run_m3_anchor_convergence
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_anchor_convergence(
        ctx,
        m2_run_id=args.m2_run_id,
        profile=args.profile,
        mode=args.mode,
        compound_pose_clusters=Path(args.compound_pose_clusters) if args.compound_pose_clusters else None,
        cluster_membership_table=Path(args.cluster_membership_table) if args.cluster_membership_table else None,
        pose_attribution_table=Path(args.pose_attribution_table) if args.pose_attribution_table else None,
        pose_clustering_qc=Path(args.pose_clustering_qc) if args.pose_clustering_qc else None,
        accepted_pockets=Path(args.accepted_pockets) if args.accepted_pockets else None,
        accepted_boxes=Path(args.accepted_boxes) if args.accepted_boxes else None,
        compound_manifest=Path(args.compound_manifest) if args.compound_manifest else None,
        config=Path(args.config) if args.config else None,
        force=args.force,
        allow_partial=args.allow_partial,
        require_clustering_ready=args.require_clustering_ready,
        allow_reference_only=args.allow_reference_only,
        private_map=Path(args.private_map) if args.private_map else None,
        strict_primary_state_support=args.strict_primary_state_support,
        write_empty_convergence=args.write_empty_convergence,
    )
    print(
        "m3-anchor-convergence {0}: profile={1} support_rows={2} anchor_rows={3} evidence_ready={4}".format(
            report.status,
            report.cluster_context.get("selected_profile"),
            report.support.get("compound_pocket_support_rows", 0),
            report.anchor_convergence.get("anchor_rows", 0),
            report.m3_t10_evidence_integration_ready,
        )
    )
    print("compound_pocket_support={0}".format(report.support_csv))
    print("compound_anchor_convergence={0}".format(report.anchor_csv))
    print("anchor_convergence_qc_json={0}".format(report.qc_json))
    print("anchor_convergence_qc_csv={0}".format(report.qc_csv))
    print("m3_task9_anchor_convergence={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_evidence_tiering(args):
    from pathlib import Path

    from egfr_myo1d.compound.evidence_integration import run_m3_evidence_tiering
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    def opt_path(name):
        value = getattr(args, name, None)
        return Path(value) if value else None

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_evidence_tiering(
        ctx,
        m2_run_id=args.m2_run_id,
        profile=getattr(args, "profile", None),
        mode=getattr(args, "mode", "integrate"),
        accepted_pockets=opt_path("accepted_pockets"),
        accepted_boxes=opt_path("accepted_boxes"),
        pocket_gate_qc=opt_path("pocket_gate_qc"),
        ppi_consensus=opt_path("ppi_consensus"),
        ppi_to_pocket_evidence=opt_path("ppi_to_pocket_evidence"),
        pocket_ppi_relationship=opt_path("pocket_ppi_relationship"),
        pocket_membrane_geometry=opt_path("pocket_membrane_geometry"),
        atp_reference=opt_path("atp_reference"),
        compound_pose_attribution=opt_path("compound_pose_attribution"),
        compound_mechanism_classification=opt_path("compound_mechanism_classification"),
        compound_pose_clusters=opt_path("compound_pose_clusters"),
        compound_pocket_support=opt_path("compound_pocket_support"),
        compound_anchor_convergence=opt_path("compound_anchor_convergence"),
        pose_attribution_qc=opt_path("pose_attribution_qc"),
        pose_clustering_qc=opt_path("pose_clustering_qc"),
        anchor_convergence_qc=opt_path("anchor_convergence_qc"),
        compound_manifest=opt_path("compound_manifest"),
        broad_anchor_scan_support=opt_path("broad_anchor_scan_support"),
        allow_broad_scan_support=getattr(args, "allow_broad_scan_support", False),
        config=opt_path("config"),
        force=getattr(args, "force", False),
        allow_partial=getattr(args, "allow_partial", False),
        require_anchor_ready=getattr(args, "require_anchor_ready", True),
        strict_tier1=getattr(args, "strict_tier1", True),
        write_empty_evidence=getattr(args, "write_empty_evidence", False),
        private_map=Path(getattr(args, "private_map", "fresh/data/private/compound_id_map.csv")),
    )
    print(
        "m3-evidence-tiering {0}: profile={1} evidence_rows={2} candidate_rows={3} tier1={4} tier2={5} report_ready={6}".format(
            report.status,
            report.evidence_context.get("selected_profile"),
            report.evidence_integration.get("pocket_compound_evidence_rows", 0),
            report.candidate_tiering.get("candidate_rows", 0),
            report.candidate_tiering.get("tier1", 0),
            report.candidate_tiering.get("tier2", 0),
            report.m3_t11_report_ready,
        )
    )
    print("pocket_compound_evidence_table={0}".format(report.evidence_csv))
    print("final_m3_candidate_hypotheses={0}".format(report.candidates_csv))
    print("rejected_candidate_reasons={0}".format(report.reject_csv))
    print("final_candidate_gate_qc_json={0}".format(report.qc_json))
    print("final_candidate_gate_qc_csv={0}".format(report.qc_csv))
    print("m3_task10_evidence_tiering={0}".format(report.report_md))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_report(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.report.milestone3_report import run_m3_report

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_report(
        ctx,
        m2_run_id=args.m2_run_id,
        profile=args.profile,
        mode=args.mode,
        allow_partial=args.allow_partial,
        require_t10_ready=args.require_t10_ready,
        strict_confidentiality=args.strict_confidentiality,
        strict_claims=args.strict_claims,
        config=Path(args.config) if args.config else None,
        private_map=Path(args.private_map) if args.private_map else None,
        force=args.force,
        cleanup_mode="plan",
    )
    print(
        "m3-report {0}: profile={1} candidates={2} tier1={3} tier2={4} post_m3_ready={5} m4_ready={6}".format(
            report.status,
            report.profile,
            report.counts.get("candidate_rows", 0),
            report.counts.get("tier1", 0),
            report.counts.get("tier2", 0),
            report.post_m3_ready,
            report.m4_ready,
        )
    )
    print("milestone3_summary={0}".format(report.outputs.get("milestone3_summary")))
    print("final_candidate_summary={0}".format(report.outputs.get("final_candidate_summary")))
    print("m3_report_qc_json={0}".format(report.qc_json))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_cleanup(args):
    from pathlib import Path

    from egfr_myo1d.cleanup.milestone_cleanup import plan_m3_cleanup
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = plan_m3_cleanup(
        ctx,
        profile=args.profile,
        mode=args.mode,
        cleanup_scope=args.cleanup_scope,
        config=Path(args.config) if args.config else None,
        preserve_selected_poses=args.preserve_selected_poses,
        preserve_all_docking_outputs=args.preserve_all_docking_outputs,
        confirm_run_id=args.confirm_run_id,
        force=args.force,
    )
    print(
        "m3-cleanup {0}: mode={1} scope={2} planned={3} deleted={4}".format(
            report.status,
            report.report.get("mode"),
            report.report.get("cleanup_scope"),
            report.report.get("files_planned_for_deletion"),
            report.report.get("files_deleted"),
        )
    )
    print("cleanup_report={0}".format(report.report_root))
    print("cleanup_qc={0}".format(report.cleanup_qc))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_m3_report_cleanup(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.report.milestone3_report import run_m3_report

    ctx = RunContext.create(args.run_id, "smoke_input")
    initialize_logs(ctx)
    report = run_m3_report(
        ctx,
        m2_run_id=args.m2_run_id,
        profile=args.profile,
        mode=args.mode,
        allow_partial=args.allow_partial,
        require_t10_ready=args.require_t10_ready,
        strict_confidentiality=args.strict_confidentiality,
        strict_claims=args.strict_claims,
        config=Path(args.config) if args.config else None,
        private_map=Path(args.private_map) if args.private_map else None,
        force=args.force,
        cleanup_mode=args.cleanup_mode,
        cleanup_scope=args.cleanup_scope,
        confirm_run_id=args.confirm_run_id,
        preserve_selected_poses=args.preserve_selected_poses,
        preserve_all_docking_outputs=args.preserve_all_docking_outputs,
    )
    print(
        "m3-report-cleanup {0}: profile={1} reports={2} cleanup_mode={3} planned={4} deleted={5}".format(
            report.status,
            report.profile,
            report.counts.get("candidate_rows", 0),
            report.cleanup.report.get("mode"),
            report.cleanup.report.get("files_planned_for_deletion"),
            report.cleanup.report.get("files_deleted"),
        )
    )
    print("milestone3_summary={0}".format(report.outputs.get("milestone3_summary")))
    print("cleanup_report={0}".format(report.cleanup.report_root))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_manifest_ligands(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.ligand.manifest import build_ligand_manifest

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)

    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=Path(args.ligands_dir) if args.ligands_dir else None,
        private_mapping_path=Path(args.private_mapping) if args.private_mapping else None,
        profile=args.profile,
        compound_stage_enabled=(args.compound_stage_enabled == "true"),
    )
    print(
        "manifest-ligands {0}: present={1} missing={2} leak={3}".format(
            manifest.status,
            manifest.present_count,
            manifest.missing_count,
            manifest.internal_ids_leaked_into_outputs,
        )
    )
    print("ligand_manifest_qc_csv={0}".format(manifest.output_qc_csv))
    print("ligand_manifest_report_json={0}".format(manifest.output_manifest_json))
    if manifest.status == "FAIL":
        return 1
    return 0


def _cmd_prepare_pbs(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.hpc.pbs import generate_pbs

    # Use --mode 'smoke_env' as the run-context init mode for non-smoke modes
    init_mode = args.mode if args.mode in ("smoke_env", "smoke_input") else "smoke_env"
    try:
        ctx = RunContext.for_existing(args.run_id)
    except Exception:
        ctx = RunContext.create(args.run_id, init_mode)
    initialize_logs(ctx)

    output_path = Path(args.output_path) if args.output_path else None

    pbs = generate_pbs(
        ctx,
        job_name=args.job_name,
        mode=args.mode,
        node=args.node,
        ppn=args.ppn,
        walltime=args.walltime,
        output_path=output_path,
        profile=args.profile,
        input_root=args.input_root,
    )
    print(
        "prepare-pbs {0}: job={1} mode={2} node={3} ppn={4} walltime={5}".format(
            pbs.status, pbs.job_name, pbs.mode, pbs.node, pbs.ppn, pbs.walltime
        )
    )
    print("pbs_file={0}".format(pbs.output_path))
    if pbs.status == "FAIL":
        return 1
    return 0


def _cmd_compute_membrane_frame(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.model.membrane_frame import (
        ALL_STATES,
        run_membrane_frame_computation,
    )

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)

    state_ids = list(ALL_STATES) if args.state == "all" else [args.state]
    source = Path(args.full_frame_source) if args.full_frame_source else None

    frames, overall = run_membrane_frame_computation(
        ctx, state_ids=state_ids, full_frame_source=source, profile=args.profile
    )
    print(
        "compute-membrane-frame {0}: states={1} statuses={2}".format(
            overall,
            len(frames),
            ",".join("{0}={1}".format(f.state_id, f.status) for f in frames),
        )
    )
    print("membrane_frame_json={0}".format(ctx.manifest_dir / "membrane_frame.json"))
    if overall == "FAIL":
        return 1
    return 0


def _cmd_prepare_receptor(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.model.receptor_normalize import normalize_receptor

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)

    report = normalize_receptor(
        ctx,
        Path(args.source),
        state_id=args.state,
        profile=args.profile,
        strict=args.strict,
    )
    print(
        "prepare-receptor {0}: state={1} case={2} protomers={3} v924r_warn={4} warnings={5}".format(
            report.status,
            report.state_id,
            report.case,
            report.protomer_count,
            report.v924r_warn,
            len(report.warnings),
        )
    )
    print("receptor_manifest={0}".format(report.manifest_json))
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_prepare_myo1d(args):
    from pathlib import Path

    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext
    from egfr_myo1d.myo1d.qc import run_myo1d_qc

    ctx = RunContext.create(args.run_id, args.mode)
    initialize_logs(ctx)

    report = run_myo1d_qc(
        ctx,
        Path(args.source),
        construct_range=args.construct,
        profile=args.profile,
    )
    print(
        "prepare-myo1d {0}: construct={1} n_residues={2} caps={3} warnings={4}".format(
            report.status,
            report.construct_id,
            report.n_residues,
            report.ace_nme_caps_present,
            len(report.warnings),
        )
    )
    print(
        "myo1d_construct_manifest={0}".format(
            ctx.manifest_dir / "myo1d_construct_manifest.json"
        )
    )
    if report.status == "FAIL":
        return 1
    return 0


def _cmd_cleanup(args):
    from egfr_myo1d.core.cleanup import run_cleanup
    from egfr_myo1d.core.logging_utils import initialize_logs
    from egfr_myo1d.core.run_context import RunContext

    ctx = RunContext.for_existing(args.run_id)
    initialize_logs(ctx)

    if args.dry_run is None:
        dry_run = None
    else:
        dry_run = args.dry_run == "true"

    report = run_cleanup(ctx, mode=args.mode, dry_run=dry_run, profile=args.profile)

    print(
        "cleanup {0}: mode={1} dry_run={2} candidates={3} deleted={4} preserved={5} errors={6}".format(
            report.status,
            report.mode,
            report.dry_run,
            report.candidate_count,
            report.deleted_count,
            report.preserved_count,
            len(report.errors),
        )
    )
    print("cleanup_report={0}".format(ctx.manifest_dir / "cleanup_report.json"))
    if report.status == "FAIL":
        return 1
    return 0


def main(argv=None):
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        if exc.code == 0:
            return 0
        raise
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    try:
        return args.func(args)
    except Exception as exc:
        from egfr_myo1d.core.run_context import RunContextError

        if isinstance(exc, (RunContextError, RuntimeError, ValueError)):
            sys.stderr.write("ERROR: {0}\n".format(exc))
            return 2
        raise


if __name__ == "__main__":
    raise SystemExit(main())
