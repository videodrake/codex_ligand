# Optional AI Tool Runtime Status

This document records optional AI/interface tools tested on the HPC after the
core M2/M3 preflight. These tools are optional evidence layers. The fresh
workflow must not require them for core execution.

Verified on:

- Date: 2026-04-30
- Node: `node04`
- Branch: `codex/m2-ppi-input-generation-clean-v2`
- Repository root on HPC:
  `/work4/hwang/onepack/new/codex_ligand/codex_ligand/codex_ligand`

## Summary

| Tool | Status | Runtime | Use in M3 |
| --- | --- | --- | --- |
| PeSTo | Installed; dependency smoke OK; one-PDB prediction OK | `conda env: pesto` | Residue-level protein-interface likelihood |
| PocketMiner | Installed; source import OK; checkpoint restore OK; one-PDB prediction OK | `conda env: pocketminer` | Residue-level cryptic-pocket likelihood |
| MaSIF | Container pulled; basic container smoke OK; prediction smoke pending | rootless `podman` image | Surface/PPI-interface pattern evidence |
| InDeep | Not locally installed | external/server candidate | Optional ligandability/interface evidence |
| PASSer | External-server only | external server | Optional allosteric-site evidence |

## PeSTo

Purpose:

- Predict residue-level protein-binding interface likelihood from a PDB
  structure.
- Use as soft evidence that an EGFR surface patch is interface-like.
- Output predictions are written into PDB B-factor fields.

Installation location:

```text
/home/eunae/tools/PeSTo
```

Conda environment:

```text
/home/eunae/.conda/envs/pesto
```

Verified runtime:

```text
torch 2.5.1
torch_cuda_available False
gemmi 0.7.5
```

Verified model:

```text
/home/eunae/tools/PeSTo/model/save/i_v4_1_2021-09-07_11-21/model_ckpt.pt
```

Smoke result:

```text
PeSTo dependency import OK
PeSTo prediction smoke OK
```

Smoke output:

```text
fresh/runs/optional_ai_tool_discovery/pesto_prediction_smoke/output/
    pesto_smoke_input_i0.pdb
    pesto_smoke_input_i1.pdb
    pesto_smoke_input_i2.pdb
    pesto_smoke_input_i3.pdb
    pesto_smoke_input_i4.pdb
```

Run setup:

```bash
source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pesto
cd /home/eunae/tools/PeSTo
```

Adapter guidance:

- Use CPU by default.
- Use model `i_v4_1_2021-09-07_11-21`.
- For each input PDB, run PeSTo inference and capture the generated `_i*.pdb`
  files.
- Parse B-factor values from the relevant PeSTo output PDB and map them back to
  receptor residue identifiers.
- Treat PeSTo scores as soft evidence only; they must not override ATP,
  membrane-accessibility, dimer-accessibility, or PPI-consensus gates.

## PocketMiner

Purpose:

- Predict residue-level cryptic-pocket likelihood from a single protein
  structure.
- Use as soft evidence that a PPI-adjacent shallow groove may be openable or
  cryptic.

Installation location:

```text
/home/eunae/tools/gvp_pocketminer
```

Repository state:

```text
branch: pocket_pred
commit: 187062d
```

Conda environment:

```text
/home/eunae/.conda/envs/pocketminer
```

Verified model checkpoint:

```text
/home/eunae/tools/gvp_pocketminer/models/pocketminer.index
/home/eunae/tools/gvp_pocketminer/models/pocketminer.data-00000-of-00001
```

Required runtime setup:

```bash
source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pocketminer
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=""
export TF_CPP_MIN_LOG_LEVEL=2
cd /home/eunae/tools/gvp_pocketminer
```

The `LD_LIBRARY_PATH` line is required because the HPC shell otherwise loads the
system `/lib64/libstdc++.so.6` first, which lacks `GLIBCXX_3.4.29` needed by
conda-installed pandas/scipy extensions. This is a session-local environment
variable and does not modify the shared server system.

Verified source imports:

```text
import OK gvp
import OK models
import OK util
```

Verified prediction entrypoint:

```text
/home/eunae/tools/gvp_pocketminer/src/xtal_predict.py
```

Smoke result:

```text
CHECKPOINT RESTORED FROM /home/eunae/tools/gvp_pocketminer/models/pocketminer
pred_shape (1, 292)
pred_min 0.025487031787633896
pred_max 0.930513858795166
PocketMiner prediction smoke OK
```

Smoke output:

```text
fresh/runs/optional_ai_tool_discovery/pocketminer_prediction_smoke/output/
    pocketminer_smoke_preds.npy
    pocketminer_smoke_predictions.txt
```

Adapter guidance:

- Use CPU by default.
- Load checkpoint basename `/home/eunae/tools/gvp_pocketminer/models/pocketminer`.
- Use `src/xtal_predict.py` / `make_predictions` for one or more PDB inputs.
- Save raw `.npy`, a flat text score file, and a manifest mapping score index to
  residue identity.
- Treat scores as soft evidence only. PocketMiner must not nominate a pocket
  without PPI adjacency and non-ATP/membrane/dimer gate support.

## MaSIF

Purpose:

- Predict and characterize surface/interface patterns from molecular surface
  features.
- Use as soft PPI-surface evidence for EGFR PPI-adjacent patches.

Runtime:

```text
rootless podman 4.2.0
image: docker.io/pablogainza/masif:latest
image id: 6b3c808b7bf7fabdadfee4c6dc2a48c4761b4a118d94983131f34e5a76754a12
container repo: /masif
container commit: dde04c8
container python: 3.6.6
```

Verified inside the container:

```text
/masif/data/masif_site/data_prepare_one.sh
/masif/data/masif_site/predict_site.sh
masif smoke OK
```

Podman/GPU status:

- `podman run hello-world`: OK.
- CUDA image pull: OK.
- CDI GPU passthrough: not configured.
- GPU container execution: not verified.

Adapter guidance:

- Treat MaSIF as container-backed optional evidence.
- Run CPU/basic prediction first. Do not assume GPU passthrough.
- Keep all generated files under `fresh/runs/<run_id>/`.
- Prediction smoke with a small PDB is still pending before production adapter
  integration.

## Not Local

InDeep:

- No local command or conda env found.
- Treat as external/server or future isolated-install candidate.

PASSer:

- External-server only.
- Use only if needed for independent allosteric-site evidence.

## Safety Rules

- Do not use `sudo`.
- Do not modify `/usr`, `/bin`, `/lib64`, `/opt`, system Java, CUDA drivers, or
  lab-wide modules.
- Keep external tool clones under `$HOME/tools` or another user-owned external
  tool directory.
- Keep optional tool outputs under `fresh/runs/<run_id>/`.
- Optional AI evidence cannot override hard workflow gates.

## Preflight Note

`python -m egfr_myo1d.cli tool-preflight` currently checks tools from the active
environment and `PATH`. It does not yet activate `pesto`, `pocketminer`, or run
the MaSIF podman image by itself. Until env-specific optional-tool adapters are
implemented, use this document and `fresh/configs/tool_envs.yaml` as the
authoritative runtime record for these optional tools.
