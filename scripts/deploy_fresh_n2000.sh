#!/bin/bash
# ============================================================
# Fresh Deployment Script: n_poses=2000, pocket_cap=100
#
# Purpose:
#   Git clone → precheck → full Workflow A production
#   Vina n_poses=2000, postprocess pocket cap=100
#
# Usage (HPC login node):
#   bash scripts/deploy_fresh_n2000.sh
#
# Or from an existing repo:
#   bash /work4/hwang/onepack/codex_ligand/scripts/deploy_fresh_n2000.sh
#
# If git clone fails (network blocked), use local copy mode:
#   CLONE_MODE=copy bash scripts/deploy_fresh_n2000.sh
#   CLONE_MODE=copy SOURCE_DIR=/work4/hwang/onepack/codex_ligand bash scripts/deploy_fresh_n2000.sh
# ============================================================

set -eo pipefail

# ── Configuration ──
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BASE_DIR="${BASE_DIR:-/work4/hwang/onepack}"
DEPLOY_DIR="${DEPLOY_DIR:-${BASE_DIR}/codex_ligand_n2000_${TIMESTAMP}}"
REPO_URL="https://github.com/videodrake/codex_ligand.git"
CLONE_MODE="${CLONE_MODE:-git}"   # "git" or "copy"
SOURCE_DIR="${SOURCE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)}"

echo "============================================"
echo "  Fresh Deployment: Workflow A"
echo "    n_poses=2000, pocket_cap=100"
echo "  Target: ${DEPLOY_DIR}"
echo "  Clone mode: ${CLONE_MODE}"
echo "  Started: $(date)"
echo "============================================"

# ── Step 1: Clone or copy ──
echo ""
echo "[1/4] Setting up working directory..."

if [[ "${CLONE_MODE}" == "git" ]]; then
    echo "  Git cloning from ${REPO_URL}"
    git clone "${REPO_URL}" "${DEPLOY_DIR}"
else
    echo "  Copying from ${SOURCE_DIR}"
    if [[ ! -d "${SOURCE_DIR}" ]]; then
        echo "  [ERROR] Source directory not found: ${SOURCE_DIR}"
        exit 1
    fi
    mkdir -p "${DEPLOY_DIR}"
    # Copy code, config, scripts, input (exclude output/ and .git large objects)
    rsync -a \
        --exclude='output/' \
        --exclude='.git/objects/pack/' \
        --exclude='*.o' \
        --exclude='*.e' \
        --exclude='__pycache__/' \
        --exclude='.venv/' \
        "${SOURCE_DIR}/" "${DEPLOY_DIR}/"
    echo "  Copied. Initializing git..."
    cd "${DEPLOY_DIR}"
    if [[ -d .git ]]; then
        git checkout -- . 2>/dev/null || true
    fi
fi

cd "${DEPLOY_DIR}"
echo "  Working directory: $(pwd)"

# ── Step 2: Verify config ──
echo ""
echo "[2/4] Verifying configuration..."

# Check key parameters
N_POSES=$(python3 -c "
import yaml
with open('config/example-project.yaml') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('vina', {}).get('n_poses', '?'))
" 2>/dev/null || grep 'n_poses:' config/example-project.yaml | tail -1 | awk '{print $2}')

POCKET_CAP=$(python3 -c "
import yaml
with open('config/example-project.yaml') as f:
    cfg = yaml.safe_load(f)
print(cfg.get('postprocess', {}).get('max_poses_per_pocket', '?'))
" 2>/dev/null || grep 'max_poses_per_pocket:' config/example-project.yaml | awk '{print $2}')

echo "  vina.n_poses:                    ${N_POSES}"
echo "  postprocess.max_poses_per_pocket: ${POCKET_CAP}"
echo "  vina.exhaustiveness:             $(grep 'exhaustiveness:' config/example-project.yaml | awk '{print $2}')"
echo "  vina.energy_range:               $(grep 'energy_range:' config/example-project.yaml | awk '{print $2}')"

if [[ "${N_POSES}" != "2000" ]]; then
    echo "  [WARN] n_poses is ${N_POSES}, expected 2000"
fi
if [[ "${POCKET_CAP}" != "100" ]]; then
    echo "  [WARN] max_poses_per_pocket is ${POCKET_CAP}, expected 100"
fi

# Check input files
echo ""
echo "  Input files:"
for f in input/receptors/*.pdb; do
    [[ -f "$f" ]] && echo "    [OK] $f" || echo "    [MISSING] $f"
done
for f in input/ligands/*.sdf; do
    [[ -f "$f" ]] && echo "    [OK] $f" || echo "    [MISSING] $f"
done

# Check PPI partner files
for f in "input/PPI/beta_meander.pdb" "input/PPI/TH1 domain.pdb"; do
    [[ -f "$f" ]] && echo "    [OK] $f" || echo "    [MISSING] $f"
done

# Check Phase 1 INI configs
INI_COUNT=$(ls config/phase1/phase1_prod_*.ini 2>/dev/null | wc -l)
echo "  Phase 1 INI configs: ${INI_COUNT} production files"

# ── Step 3: Verify conda environment ──
echo ""
echo "[3/4] Checking conda environment..."

if command -v conda >/dev/null 2>&1; then
    set +u
    source ~/.bashrc 2>/dev/null || true
    set -u
    if conda activate pyrosetta 2>/dev/null; then
        echo "  [OK] conda env 'pyrosetta' activated"
        python -c "import yaml; print('  [OK] PyYAML available')" 2>/dev/null || echo "  [WARN] PyYAML not found"
        python -c "import numpy; print('  [OK] NumPy available')" 2>/dev/null || echo "  [WARN] NumPy not found"
        python -c "import pandas; print('  [OK] Pandas available')" 2>/dev/null || echo "  [WARN] Pandas not found"
        # Check Vina
        if command -v vina >/dev/null 2>&1; then
            echo "  [OK] Vina binary found: $(which vina)"
        else
            echo "  [WARN] Vina binary not found in PATH"
        fi
        # Check PDBQT conversion tools
        for tool in obabel prepare_receptor mk_prepare_ligand.py; do
            if command -v "$tool" >/dev/null 2>&1; then
                echo "  [OK] ${tool} available"
                break
            fi
        done
    else
        echo "  [WARN] Could not activate 'pyrosetta' conda env"
    fi
else
    echo "  [WARN] conda not found"
fi

# ── Step 4: Submit PBS jobs ──
echo ""
echo "[4/4] Submitting PBS jobs..."
echo ""
echo "  Submission chain:"
echo "    1. Pre-qsub checks (syntax, tests, environment)"
echo "    2. Production fresh (reset outputs + full Workflow A)"
echo ""

PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
echo "  [SUBMITTED] Precheck:   ${PRECHECK_JOB}"

PROD_JOB=$(qsub -W depend=afterok:${PRECHECK_JOB} config/run_production_fresh.pbs)
echo "  [SUBMITTED] Production: ${PROD_JOB} (depends on precheck)"

echo ""
echo "============================================"
echo "  Deployment Complete"
echo "============================================"
echo ""
echo "  Working directory: ${DEPLOY_DIR}"
echo "  Precheck job:      ${PRECHECK_JOB}"
echo "  Production job:    ${PROD_JOB}"
echo ""
echo "  Key parameters:"
echo "    Vina n_poses:    2000"
echo "    Pocket cap:      100 poses/pocket (round-robin by ligand)"
echo "    Exhaustiveness:  384"
echo "    Receptors:       3GT8_raw, EGFR_160-185, EGFR_170-200"
echo "    Ligands:         173940, 97806, VAX-C12_0"
echo "    PPI:             3 states x 5 seeds = 300K models"
echo ""
echo "  Monitor commands:"
echo "    qstat -a                                    # Job status"
echo "    tail -f ${DEPLOY_DIR}/production_fresh.o*   # PBS output"
echo "    # PPI seed progress:"
echo "    tail -f ${DEPLOY_DIR}/output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/PROGRESS.log"
echo ""
echo "  Pocket cap results (after Phase 4):"
echo "    cat ${DEPLOY_DIR}/output/workflow_a/phase4_vina_postprocess/egfr_myo1d_vina/vina_pocket_cap_report.csv"
echo "============================================"
