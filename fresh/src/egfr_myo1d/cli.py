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


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
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
