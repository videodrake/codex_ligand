from pathlib import Path


FRESH_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = FRESH_ROOT.parent


def test_required_skeleton_paths_exist():
    required_paths = [
        FRESH_ROOT,
        FRESH_ROOT / "README.md",
        FRESH_ROOT / "pyproject.toml",
        FRESH_ROOT / "configs" / "fresh_run.yaml",
        FRESH_ROOT / "configs" / "hpc.yaml",
        FRESH_ROOT / "configs" / "paths.yaml",
        FRESH_ROOT / "configs" / "gates.yaml",
        FRESH_ROOT / "configs" / "receptor_states.yaml",
        FRESH_ROOT / "docs" / "milestone1_foundation_plan.md",
        FRESH_ROOT / "docs" / "repo_audit.md",
        FRESH_ROOT / "docs" / "donor_module_list.csv",
        FRESH_ROOT / "docs" / "do_not_touch.md",
        FRESH_ROOT / "scripts" / "generate_pbs.py",
        FRESH_ROOT / "scripts" / "submit_smoke_env.sh",
        FRESH_ROOT / "scripts" / "submit_smoke_input.sh",
        FRESH_ROOT / "scripts" / "cleanup_run.py",
        FRESH_ROOT / "src" / "egfr_myo1d" / "__init__.py",
        FRESH_ROOT / "src" / "egfr_myo1d" / "cli.py",
        FRESH_ROOT / "runs" / ".gitkeep",
    ]

    missing = [path for path in required_paths if not path.exists()]
    assert missing == []


def test_required_fixture_files_exist():
    fixture_dir = FRESH_ROOT / "tests" / "fixtures"
    required_fixtures = [
        "mini_explicit_AB.pdb",
        "mini_duplicate_chain_X.pdb",
        "mini_atom_hetatm_cap.pdb",
        "mini_myo1d_955_1006.pdb",
    ]

    missing = [name for name in required_fixtures if not (fixture_dir / name).exists()]
    assert missing == []


def test_gitignore_mentions_fresh_guardrails():
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")
    expected_patterns = [
        "fresh/runs/*",
        "!fresh/runs/.gitkeep",
        "fresh/data/normalized/*",
        "fresh/data/raw/ligands/*",
        "fresh/data/private/*",
        "fresh/data/raw/receptors/*.xtc",
        "!fresh/data/raw/receptors/*.pdb",
        "!fresh/data/raw/myo1d/*.pdb",
    ]

    missing = [pattern for pattern in expected_patterns if pattern not in gitignore]
    assert missing == []
