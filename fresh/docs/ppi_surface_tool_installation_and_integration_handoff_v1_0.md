# EGFR-MYO1D Fresh Workflow — PPI-Surface Tool Installation, Testing, and Integration Handoff v1.0

**Project:** EGFR-MYO1D membrane-compatible PPI/pocket/fresh-VS discovery
**Audience:** Codex or another implementation agent
**Repository:** `https://github.com/videodrake/codex_ligand`
**Execution target:** HPC server, PBS/qsub, nodes `node04`, `node05`, `node06`, 32 cores per node target
**Status:** Agent handoff document for installing, testing, and integrating additional PPI-surface pocket prioritization tools through Milestone 3
**Important correction:** Professor-provided private compounds are **not** part of the public Milestone 3 workflow. They may only be used as a separate private internal confidence check after public pocket zones are frozen.

---

## 0. One-sentence purpose

This document explains how to install/test additional pocket, surface, hotspot, and AI-interface tools and how to attach them to the EGFR-MYO1D workflow so that PyRosetta-derived PPI-adjacent sites are filtered down to one or two defensible **PPI-compatible EGFR surface pocket zones** for public fresh virtual screening.

The goal is not to find a generic deep ligand pocket. The goal is to find a **surface-accessible pocket/groove/zone capable of perturbing MYO1D binding to EGFR**.

---

## 1. Scientific scope that every agent must preserve

### 1.1 What the project is trying to find

The project is based on the following evidence chain:

```text
Membrane-compatible EGFR inactive dimer model
    -> MYO1D TH1 beta-meander PPI consensus patch
    -> PPI-compatible non-ATP EGFR dimer-surface pocket zone
    -> public fresh virtual screening
    -> public EGFR-MYO1D PPI-disruptive candidate hypotheses
```

The target is a **non-ATP, C-lobe / dimer-surface / lower-lateral / membrane-compatible pocket zone** near the MYO1D PPI surface. The workflow must not drift toward ordinary ATP-site kinase inhibitor discovery.

### 1.2 What the added tools are for

PyRosetta PPI docking will likely generate many possible EGFR-MYO1D poses and many nearby pockets. The added tools are not meant to create more final candidates. They are **filters** that reduce many possible PPI-adjacent sites to a few defensible pocket zones.

The filtering logic is:

```text
Many PyRosetta PPI sites and adjacent surface pockets
    -> keep repeated MYO1D PPI consensus patches
    -> keep nearby real surface cavities/grooves
    -> keep dynamically persistent or repeatedly observed cavities
    -> keep chemically ligandable probe hotspots
    -> keep PPI-interface-like surfaces
    -> keep non-ATP, lower/lateral, dimer-accessible sites
    -> export one or two public PPI-compatible pocket zones to Milestone 3
```

### 1.3 Hard constraints

The following rules override any tool score:

```text
1. ATP pocket = hard reject for PPI-modulator candidate.
2. PPI-unrelated pocket = hard reject for this project.
3. Dimer-buried or membrane-inaccessible pocket = reject or quarantine.
4. 3GT8_raw-only support = reference/control only, not primary final evidence.
5. AI/ML scores are soft prioritization only, never hard proof.
6. Private professor compounds are not public workflow evidence.
7. All outputs must stay under fresh/runs/<run_id>/.
8. Every install/test/run must produce centralized logs and cleanup reports.
```

---

## 2. Current known environment

From the user-provided HPC environment:

```text
conda env: pyrosetta
python: /home/eunae/.conda/envs/pyrosetta/bin/python
python version: 3.9.25
pyrosetta: import OK
rdkit: import OK
BioPython: import OK
numpy: import OK
pandas: import OK
vina: /usr/local/anaconda/3/2023.09/bin/vina
vina version: AutoDock Vina 52ec525-mod
fpocket: /home/eunae/.conda/envs/pyrosetta/bin/fpocket
fpocket version: 4.0
obabel: /usr/local/anaconda/3/2023.09/bin/obabel
Open Babel version: 3.1.0
PBS/qsub: available
nodes: node04, node05, node06
cores target: 32 cores/node
```

Do not assume that `mdpocket`, `gmx`, `P2Rank`, `InDeep`, `PeSTo`, `MaSIF`, `PocketMiner`, or PASSer local tools are installed until preflight confirms them.

---

## 3. Global implementation strategy

### 3.1 Do not install everything blindly

The agent must implement and run a **tool installation/testing ladder**:

```text
Stage A: discover existing tools without installing anything
Stage B: test already-installed core tools
Stage C: install only low-risk core additions
Stage D: add optional AI/surface tools one by one
Stage E: add heavy validation tools only after M2 core passes
```

Every tool gets one of these statuses:

```text
available
installed_and_smoke_passed
installed_but_smoke_failed
not_installed
optional_disabled
requires_user_review
external_server_only
excluded_from_core
```

### 3.2 Keep optional tools isolated

Do not damage the working `pyrosetta` environment. Prefer separate environments for optional tools.

Recommended environment policy:

```text
pyrosetta env:
    PyRosetta PPI docking, Vina, RDKit, fpocket, OpenBabel, core fresh workflow.

ppi_surface env:
    pyKVFinder, simple Python analysis tools, optional adapters.

pesto env:
    PeSTo only.

masif env or container:
    MaSIF only; optional.

indeep env:
    InDeep only; optional.

msmd/md_env:
    GROMACS/MSMD validation only; optional.
```

### 3.3 Do not vendor external tool repositories into the public repo

External tool clones should go outside the public repo or into a gitignored folder.

Recommended:

```text
/work4/eunae/external_tools/
    fpocket_src_if_needed/
    pyKVFinder_src_if_needed/
    PeSTo/
    masif/
    InDeep/
    exprorer_msmd/
    HOTPocket_or_PocketMiner/
```

If a tool path must be tracked, store only the path in:

```text
fresh/configs/tool_registry.yaml
```

Do not commit external tool source, large models, ligand files, trajectory files, or private probe outputs.

---

## 4. Tool registry and preflight output

### 4.1 Create `fresh/configs/tool_registry.yaml`

The agent should add a registry like this:

```yaml
tools:
  fpocket:
    required_level: core
    env: pyrosetta
    binary: /home/eunae/.conda/envs/pyrosetta/bin/fpocket
    preflight: "fpocket --help or fpocket 2>&1 | head"

  mdpocket:
    required_level: core_if_available
    env: pyrosetta
    binary: null
    discover_from: fpocket_directory
    preflight: "mdpocket --help"

  pyKVFinder:
    required_level: core_addition
    env: ppi_surface_or_pyrosetta
    python_import: pyKVFinder
    install_candidate: "pip install pyKVFinder"

  mini_ftmap:
    required_level: core_custom
    env: pyrosetta
    depends_on: [vina, rdkit, obabel]
    install_candidate: "no external package; implement in fresh/src"

  indeep:
    required_level: optional_ai
    env: indeep
    external_path: /work4/eunae/external_tools/InDeep
    license_note: "GNU AGPLv3; do not vendor into public repo"

  pesto:
    required_level: optional_ai
    env: pesto
    external_path: /work4/eunae/external_tools/PeSTo
    license_note: "non-commercial/share-alike license; use as optional external adapter"

  masif_site:
    required_level: optional_heavy_ai
    env: masif_or_container
    external_path: /work4/eunae/external_tools/masif
    install_mode: "docker/container preferred if available"

  msmd:
    required_level: optional_heavy_validation
    env: md_env
    external_path: /work4/eunae/external_tools/exprorer_msmd
    depends_on: [gromacs]

  pocketminer:
    required_level: optional_cryptic
    env: pocketminer
    external_path: /work4/eunae/external_tools/HOTPocket_or_PocketMiner

  passer:
    required_level: optional_allosteric
    mode: external_server_or_api
    privacy_note: "Do not upload private receptor models without user approval"
```

### 4.2 Create tool preflight outputs

For each run:

```text
fresh/runs/<run_id>/manifest/tool_status.json
fresh/runs/<run_id>/logs/phase_tool_preflight.log
fresh/runs/<run_id>/reports/tool_installation_report.md
```

Minimum `tool_status.json` schema:

```json
{
  "run_id": "test_tool_preflight_YYYYMMDD_HHMMSS",
  "node": "node04",
  "tools": {
    "fpocket": {
      "status": "available",
      "version": "4.0",
      "path": "/home/eunae/.conda/envs/pyrosetta/bin/fpocket",
      "smoke_test": "passed"
    },
    "mdpocket": {
      "status": "not_tested",
      "path": null,
      "smoke_test": "pending"
    }
  }
}
```

---

## 5. Tool-by-tool purpose, install/test plan, and workflow integration

## 5.1 MDpocket / fpocket dynamic pocket tracking

### Why it matters

`fpocket` finds pockets from a single structure. `mdpocket` tracks pocket occurrence across MD frames, trajectory snapshots, or aligned conformational ensembles. This is useful because the project has two primary MD-derived receptor states, `EGFR_160-185` and `EGFR_170-200`, and possibly full trajectory frames later.

### What information it gives

```text
- transient or persistent cavity regions
- frequency grid: how often a pocket is open
- density grid: where alpha-sphere-like pocket points accumulate
- selected pocket descriptors per frame
- pocket volume variation over conformations
- dynamic cavity persistence score
```

### Install/discovery test

Do not reinstall `fpocket` unless necessary. First discover whether `mdpocket` was installed with `fpocket`.

```bash
source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

which fpocket
fpocket 2>&1 | head

which mdpocket || true
ls -l $(dirname $(which fpocket))/mdpocket 2>/dev/null || true
mdpocket --help 2>&1 | head || true
```

If `mdpocket` is missing but needed, build from fpocket source in an external tool directory, not in the repo:

```bash
mkdir -p /work4/eunae/external_tools
cd /work4/eunae/external_tools
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
# Do not sudo install unless user explicitly approves.
./bin/fpocket --help 2>&1 | head || true
./bin/mdpocket --help 2>&1 | head || true
```

If the build requires static standard C++ libraries or other system packages, stop and record `requires_user_review` rather than using `sudo`.

### Real smoke test

`mdpocket` needs aligned frames. When only the two MD-derived state PDBs are available, make a **cross-state PDB ensemble** after aligning them to the same receptor frame.

Preferred real test when trajectory is available:

```bash
mdpocket --trajectory_file aligned_egfr.xtc --trajectory_format xtc -f reference.pdb
```

Expected outputs include DX grid files and descriptor files. The adapter should parse at least:

```text
mdpout_freq_grid.dx
mdpout_dens_grid.dx
mdpout_descritpors.txt or mdpout_descriptors.txt, depending on version spelling
```

### Integration into Milestone 2

Add:

```text
M2-T8.6 mdpocket cross-state/dynamic pocket tracking
```

Outputs:

```text
fresh/runs/<run_id>/phase2_pockets/mdpocket/
fresh/runs/<run_id>/phase2_pockets/tables/mdpocket_persistence.csv
```

Suggested columns:

```text
state_group
pocket_family_id
freq_grid_support
mean_volume
volume_std
persistence_score
near_ppi_patch
near_ppi_hotspot
non_atp_region
lower_lateral_class
```

Use as a **dynamic persistence filter**, not as a final pocket selector.

---

## 5.2 pyKVFinder local PPI-zone cavity characterization

### Why it matters

`pyKVFinder` is useful for measuring whether a PPI-adjacent surface groove is actually cavity-like. It gives physical descriptors such as volume, area, depth, and hydropathy. Unlike generic global pocket ranking, this workflow should run pyKVFinder in **local PPI-zone mode**.

### What information it gives

```text
- local cavity volume
- local cavity area
- cavity depth
- hydropathy or polarity descriptors
- residue list around the cavity
- grid points defining the cavity
```

### Install test

First check existing import:

```bash
source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta
python - <<'PY'
try:
    import pyKVFinder
    print('pyKVFinder OK', getattr(pyKVFinder, '__version__', 'unknown'))
except Exception as e:
    print('pyKVFinder FAIL', repr(e))
PY
```

If not installed, install in a separate environment if possible:

```bash
conda create -n ppi_surface python=3.9 -y
conda activate ppi_surface
pip install pyKVFinder
python -c "import pyKVFinder; print('pyKVFinder OK')"
```

If the agent chooses to install into `pyrosetta`, it must first create a rollback record and log the package list before and after.

### Smoke test

Use a tiny receptor fixture or a real normalized EGFR dockable receptor if available.

The agent should implement a Python adapter rather than relying on a one-line command. The adapter should:

```text
1. read a PDB
2. restrict analysis to a local box around PPI patch residues
3. run pyKVFinder cavity detection
4. write cavity descriptors to CSV
```

### Integration into Milestone 2

Add:

```text
M2-T8.5 local PPI-zone cavity detection with pyKVFinder
```

Outputs:

```text
fresh/runs/<run_id>/phase2_pockets/pykvfinder/
fresh/runs/<run_id>/phase2_pockets/tables/pykvfinder_local_cavities.csv
```

Suggested columns:

```text
pocket_family_id
state_id
protomer_id
local_box_id
volume
area
depth
hydropathy
surface_openness
cavity_residues
near_ppi_hotspot
near_myo1d_approach_corridor
non_atp_pass
```

Use as a **surface-cavity geometry filter**.

---

## 5.3 mini-FTMap-like local probe hotspot mapping

### Why it matters

This is a custom local probe mapping method using tools already present: Vina, RDKit, and OpenBabel. It is more appropriate for PPI-surface pockets than generic deep ligand docking because it asks: “Do diverse small chemical probes repeatedly cluster in the MYO1D PPI-adjacent surface groove?”

### What information it gives

```text
- probe hotspot clusters
- probe chemical diversity at each hotspot
- hydrophobic vs polar probe preference
- overlap with EGFR PPI hotspot residues
- overlap with MYO1D approach corridor
- subpocket assignment within a larger composite zone
```

### Install test

No external install. Confirm dependencies:

```bash
source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta
which vina
vina --version 2>&1 | head
which obabel
obabel -V
python - <<'PY'
from rdkit import Chem
print('RDKit OK')
PY
```

### Probe library

Create a small public probe set in code/config, not as private data.

Recommended initial probe set:

```text
methanol
ethanol
isopropanol
acetamide
acetonitrile
urea
benzene
toluene
phenol
aniline
pyridine
imidazole
thiophene
dimethyl ether
methylamine-like neutral/charged variants if preparation is robust
acetate-like variants only if charge handling is robust
```

Start with neutral probes for smoke tests. Add charged probes only after ligand preparation and charge sanity checks are stable.

### Smoke test

Use a single accepted PPI patch local box or fixture box.

```text
1. Generate 3 probes: methanol, benzene, acetamide.
2. Convert to PDBQT with OpenBabel/RDKit pipeline.
3. Run Vina in a small local box around fixture pocket.
4. Cluster probe poses by coordinate proximity.
5. Write probe_hotspot_clusters.csv.
```

### Integration into Milestone 2

Add:

```text
M2-T8.7 local probe hotspot mapping / mini-FTMap
```

Outputs:

```text
fresh/runs/<run_id>/phase2_pockets/mini_ftmap/
fresh/runs/<run_id>/phase2_pockets/tables/mini_ftmap_probe_hotspots.csv
```

Suggested columns:

```text
probe_cluster_id
pocket_family_id
subpocket_id
state_id
protomer_id
probe_pose_count
probe_type_count
probe_chemical_diversity
mean_probe_score
best_probe_score
overlaps_ppi_hotspot
overlaps_approach_corridor
non_atp_region
surface_zone_rank
```

Use as a **chemical ligandability filter**. This is one of the most important additions.

---

## 5.4 InDeep PPI ligandability / epitope-site prediction

### Why it matters

InDeep is relevant because it was built for PPI drug discovery and predicts interactibility/ligandability over protein structures. It is more suitable than generic ligand docking AI for this project because the target is a PPI surface pocket zone.

### What information it gives

```text
- PPI ligandability-like site score
- epitope-binding-like site prediction
- 3D score maps or predicted binding-site regions
- possible small-molecule-targetable surface areas within PPI regions
```

### Install/test policy

InDeep is optional and must not be required for the core workflow. Its GitLab repository is AGPLv3; do not vendor it into the public repo.

Install in isolated external path:

```bash
mkdir -p /work4/eunae/external_tools
cd /work4/eunae/external_tools
git clone https://gitlab.pasteur.fr/InDeep/InDeep.git
cd InDeep
# The agent must inspect README and use the repository-provided installation instructions.
# Do not invent install commands beyond the README.
```

Preflight status can be one of:

```text
optional_disabled
requires_user_review
installed_and_smoke_passed
```

### Smoke test

Run only the repository-provided example first. Then run a small EGFR dockable receptor test only if the example passes.

Expected adapter output:

```text
indeep_ligandability.csv
```

Suggested columns:

```text
pocket_family_id
state_id
protomer_id
indeep_ligandability_score
indeep_epitope_score
overlap_with_ppi_patch
overlap_with_probe_hotspot
priority_note
```

### Integration into Milestone 2

Add optional:

```text
M2-OPT2 InDeep PPI-surface ligandability scoring
```

Use as a **PPI-ligandability AI filter**, but never as a hard gate.

---

## 5.5 PeSTo protein-interface support

### Why it matters

PeSTo predicts protein interaction interfaces from protein structures. It does not find pockets. Its role is to answer whether the PyRosetta-derived EGFR-side patch looks like a protein-protein interaction surface.

### What information it gives

```text
- residue-level protein-interface probability
- output PDB with predicted interface confidence encoded in B-factor field
- optional interface type probabilities, depending on model
```

### Install test

Install in isolated environment:

```bash
mkdir -p /work4/eunae/external_tools
cd /work4/eunae/external_tools
git clone https://github.com/LBM-EPFL/PeSTo.git
cd PeSTo
conda env create -f pesto.yml
conda activate pesto
python -c "import torch, gemmi; print('PeSTo env basic imports OK')"
```

If GPU dependencies fail, try CPU mode or mark optional. Do not block the core workflow.

### Smoke test

Use PeSTo’s own example or a small PDB from its `pdbs_test` folder first. Then run on one normalized EGFR dockable receptor.

Expected output interpretation:

```text
PDB B-factor values 0 to 1 = predicted interface confidence.
```

### Integration into Milestone 2

Add optional:

```text
M2-T5.6 PeSTo PPI surface plausibility scoring
```

Outputs:

```text
fresh/runs/<run_id>/phase1_ppi/ai_interface/pesto/
fresh/runs/<run_id>/phase1_ppi/tables/pesto_interface_scores.csv
```

Suggested columns:

```text
egfr_residue
state_id
protomer_id
pesto_ppi_score
pyrosetta_contact_frequency
combined_ppi_surface_score
```

Use as a **PPI patch confidence filter**, not as pocket selection.

---

## 5.6 MaSIF-site molecular surface fingerprinting

### Why it matters

MaSIF-site predicts surface patches with propensity for protein-protein interactions using geometric and chemical surface features. It is philosophically well-matched to PPI surfaces, but is heavier to install and run than PeSTo.

### What information it gives

```text
- surface patch interaction score
- molecular surface fingerprint
- PPI-site-like surface patches
```

### Install/test policy

MaSIF is optional heavy AI. Do not require it for the core pipeline. Install only if the HPC allows the needed environment or container approach.

Recommended first option:

```text
Use the repository's Docker/container tutorial if container execution is available.
```

If container execution is not available, install manually following the MaSIF README in an isolated env.

```bash
mkdir -p /work4/eunae/external_tools
cd /work4/eunae/external_tools
git clone https://github.com/lpdi-epfl/masif.git
cd masif
# Follow README / docker_tutorial.md. Do not modify system packages.
```

### Smoke test

Run only MaSIF’s own demo first. If it passes, run on one small normalized EGFR state.

### Integration into Milestone 2

Add optional:

```text
M2-OPT3 MaSIF-site surface-pattern validation
```

Outputs:

```text
masif_site_scores.csv
```

Suggested columns:

```text
surface_patch_id
mapped_egfr_residues
masif_site_score
overlap_with_pyrosetta_ppi_patch
overlap_with_probe_hotspot
```

Use as a **surface-pattern filter**. Disable if installation is unstable.

---

## 5.7 Mixed-solvent MD / MSMD

### Why it matters

MSMD can map dynamic chemical hotspots by simulating protein in water plus small probe molecules. This is useful for top PPI-compatible surface zones because it includes protein motion and competition between water and probes.

### What information it gives

```text
- dynamic probe density maps
- probe-specific hotspot regions
- transient or cryptic surface hotspot support
- dynamic ligandability under motion
```

### Install/test policy

MSMD is **not** core. It is heavy validation for the top 1-3 PPI-compatible zones after Milestone 2. It should not block M2 or M3.

Candidate tool:

```text
EXPRORER_MSMD, a GROMACS-based MSMD automation repository.
```

Install in external path only:

```bash
mkdir -p /work4/eunae/external_tools
cd /work4/eunae/external_tools
git clone https://github.com/keisuke-yanagisawa/exprorer_msmd.git
cd exprorer_msmd
# Follow repository README. Confirm gromacs/gmx path first.
```

Preflight:

```bash
which gmx || which gmx_mpi || true
gmx --version 2>&1 | head || true
```

### Smoke test

Do not run a full MSMD simulation during tool installation. First run the repository’s minimal test if available. Then define a tiny 10-100 ps test on a small fixture or one top pocket zone after the receptor is available.

### Integration point

Add optional after M2 pocket freeze:

```text
M2.5-OPT MSMD validation of top PPI-compatible zones
```

Outputs:

```text
msmd_hotspots.csv
msmd_probe_density_summary.csv
```

Suggested columns:

```text
msmd_hotspot_id
pocket_family_id
probe_type
probe_density
dynamic_support
overlap_with_ppi_patch
overlap_with_mini_ftmap
```

Use as a **dynamic chemical hotspot validation filter** for the final 1-3 zones.

---

## 5.8 PocketMiner cryptic pocket prediction

### Why it matters

PocketMiner predicts residues likely to participate in cryptic pocket opening from a single protein structure. It can help when a PPI-adjacent surface groove is shallow but may open under dynamics.

### What information it gives

```text
- residue-level cryptic pocket opening probability
- cryptic region support near existing surface pockets
```

### Install/test policy

PocketMiner is optional and should not be treated as a robust core executable unless the agent confirms an installable implementation. Some repositories expose PocketMiner through broader benchmark or ensemble pipelines, not a simple official CLI.

Possible path:

```text
HOTPocket / PocketMiner-related environment
```

The agent must verify whether a direct PocketMiner model can be run locally. If not, status should be:

```text
optional_disabled: no stable local CLI confirmed
```

### Integration point

Add optional:

```text
M2-OPT4 PocketMiner cryptic-support scoring
```

Outputs:

```text
pocketminer_cryptic_scores.csv
```

Suggested columns:

```text
egfr_residue
state_id
protomer_id
pocketminer_score
mapped_pocket_family_id
near_ppi_patch
near_local_cavity
priority_note
```

Use only as a **cryptic possibility filter**. It cannot rescue a PPI-unrelated or ATP pocket.

---

## 5.9 PASSer allosteric site prediction

### Why it matters

PASSer predicts allosteric sites using machine-learning models. It can help interpret allosteric-near pockets, but it is not PPI-specific.

### What information it gives

```text
- top allosteric pocket probability or rank score
- residue list for top predicted pockets
- model-specific allosteric score
```

### Install/test policy

PASSer is primarily a web/API service. Use with caution because project receptor models may be private or not appropriate for external upload.

Safe test using public PDB ID:

```bash
curl -X POST -d pdb=5dkk -d chain=A https://passer.smu.edu/api
```

Do not upload custom EGFR-MYO1D receptor models unless the user explicitly approves.

If using the API, record:

```text
external_server_used: true
server_url: https://passer.smu.edu
input_type: public PDB ID or uploaded PDB
privacy_approved: yes/no
```

### Integration point

Add optional:

```text
M2-OPT5 PASSer allosteric-priority scoring
```

Outputs:

```text
passer_allosteric_scores.csv
```

Suggested columns:

```text
pocket_family_id
state_id
protomer_id
passer_rank
passer_probability
ppi_relationship_class
near_ppi_patch
non_atp_pass
allosteric_priority
privacy_status
```

Use as an **allosteric possibility filter** only. It cannot select final pockets alone.

---

## 6. Recommended priority order

### 6.1 Core to implement first

```text
1. Existing fpocket + mdpocket discovery/preflight
2. pyKVFinder local PPI-zone cavity analysis
3. mini-FTMap local probe hotspot mapping
4. PyRosetta PPI hotspot/interface-energy analysis
```

These are most aligned with the PPI-surface pocket objective and most likely to be installable.

### 6.2 Strong optional tools

```text
5. PeSTo
6. InDeep
```

PeSTo supports PPI interface plausibility. InDeep supports PPI ligandability. Both are useful but optional.

### 6.3 Heavy optional validation

```text
7. MaSIF-site
8. MSMD / mixed-solvent MD
9. PocketMiner
10. PASSer
```

These should be used only after M2 core succeeds, mainly for tie-breaking or validation of top zones.

### 6.4 Excluded from core workflow

```text
DiffDock
DynamicBind
AlphaFold3 ligand complex prediction
generic blind AI ligand docking
```

Reason: these methods tend to favor conventional protein-ligand pockets and are not optimized for shallow PPI-surface pocket zones in this workflow.

---

## 7. Updated Milestone 2 integration map

The current Milestone 2 document should be extended as follows.

### 7.1 Existing M2 core remains

```text
M2-T0  Readiness check and M2 config extension
M2-T1  PPI input and job manifest generation
M2-T2  PyRosetta PPI adapter and smoke execution harness
M2-T3  PPI pose collection and residue mapping restoration
M2-T4  MYO1D artifact and orientation QC
M2-T5  PPI consensus patch builder
M2-T6  PPI report and manual review checkpoint
M2-T7  ATP reference builder
M2-T8  fpocket pocket discovery adapter
M2-T9  Pocket normalization and pocket family merge
M2-T10 PPI-pocket relationship classifier
M2-T11 Membrane/lateral/dimer accessibility gates
M2-T12 Accepted pocket export and Milestone 2 report
```

### 7.2 Add PPI-surface filtering tasks

```text
M2-T5.5 PyRosetta PPI hotspot residue analysis
M2-T5.6 Optional PeSTo PPI surface plausibility scoring
M2-T5.7 Optional InDeep PPI ligandability scoring
M2-T8.5 pyKVFinder local PPI-zone cavity characterization
M2-T8.6 mdpocket cross-state/dynamic pocket persistence
M2-T8.7 mini-FTMap local probe hotspot mapping
M2-T8.8 Optional MaSIF-site surface-pattern validation
M2-T8.9 Optional PocketMiner cryptic pocket support
M2-T8.10 Optional PASSer allosteric pocket support
M2-T11.5 PPI-compatible pocket-zone evidence integration
```

### 7.3 New M2 primary output

Replace a simple accepted pocket list with a pocket-zone evidence table:

```text
fresh/runs/<run_id>/phase2_pockets/tables/ppi_surface_zone_evidence.csv
fresh/runs/<run_id>/phase2_pockets/tables/accepted_ppi_surface_zones_for_m3.csv
fresh/runs/<run_id>/phase2_pockets/export_for_m3/zone_composite_boxes.csv
```

Minimum columns:

```text
zone_id
pocket_family_id
state_support
symmetry_support
ppi_consensus_support
ppi_hotspot_support
surface_interface_support
cavity_geometry_support
dynamic_persistence_support
probe_hotspot_support
cryptic_support
allosteric_support
non_atp_pass
lower_lateral_pass
dimer_accessibility_pass
primary_state_support_pass
final_zone_tier
decision_reason
```

### 7.4 Hard gates remain unchanged

Hard gates:

```text
G1 non-ATP pass
G2 PPI-related pass
G3 lower/lateral membrane-compatible pass
G4 dimer-accessible pass
G5 mapping/QC pass
G6 primary-state support or explicit review flag
```

AI, probe, or cryptic scores cannot override hard gate failure.

---

## 8. Updated Milestone 3 integration map

### 8.1 Corrected official Milestone 3 scope

Milestone 3 should be renamed or interpreted as:

```text
Milestone 3: Fresh Virtual Screening and Public Candidate Nomination
```

The professor-provided private compounds are not part of the public workflow and must not appear in public reports or public evidence tables.

Official M3 flow:

```text
accepted_ppi_surface_zones_for_m3.csv
    -> public fresh library preparation
    -> focused docking to composite PPI-compatible zones
    -> pose clustering and subpocket attribution
    -> ATP migration filtering
    -> PPI-disruption geometry classification
    -> public hit ranking
    -> final public candidate shortlist
```

Private internal layer:

```text
private professor compounds
    -> run only after public M2 zone list is frozen
    -> output only in private_runs/ or external private directory
    -> not used in public scoring
```

### 8.2 M3 uses M2 zones, not generic pockets

Milestone 3 should dock to:

```text
zone_composite_boxes.csv
```

not to arbitrary fpocket boxes.

The docking box must represent:

```text
PPI-compatible pocket zone = surface groove/rim/subpocket cluster near MYO1D PPI patch
```

### 8.3 M3 outputs

```text
fresh/runs/<run_id>/phase3_fresh_vs/tables/public_vs_pose_raw.csv
fresh/runs/<run_id>/phase3_fresh_vs/tables/public_vs_pose_qc.csv
fresh/runs/<run_id>/phase3_fresh_vs/tables/zone_subpocket_attribution.csv
fresh/runs/<run_id>/phase3_fresh_vs/tables/public_hit_ranking.csv
fresh/runs/<run_id>/phase3_fresh_vs/tables/final_public_candidate_shortlist.csv
fresh/runs/<run_id>/reports/milestone3_public_vs_summary.md
```

### 8.4 M3 selection logic

A public compound hit is not ranked by Vina score alone. It must satisfy:

```text
1. docks inside an M2 accepted PPI-compatible zone
2. does not migrate to ATP pocket
3. remains near PPI hotspot residues or approach corridor
4. clusters reproducibly in a subpocket
5. has acceptable chemistry/QC
6. preferably appears in both primary receptor states or maps to the same equivalent zone
```

---

## 9. Installation/testing task sequence for another agent

### Task TOOL-0 — Add tool registry and preflight skeleton

Create:

```text
fresh/configs/tool_registry.yaml
fresh/src/egfr_myo1d/tools/tool_preflight.py
fresh/scripts/preflight_tools.sh
fresh/scripts/submit_tool_preflight.sh
fresh/docs/tool_installation_report_template.md
```

Acceptance:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.tools.tool_preflight --run-id test_tool_preflight --mode discover
```

Expected:

```text
fresh/runs/test_tool_preflight/manifest/tool_status.json
fresh/runs/test_tool_preflight/logs/phase_tool_preflight.log
```

### Task TOOL-1 — Test existing core tools

Test:

```text
fpocket
mdpocket discovery
vina
obabel
rdkit
pyrosetta
gromacs discovery only
```

Do not install yet.

### Task TOOL-2 — Install/test pyKVFinder

Use separate env if possible.

Acceptance:

```bash
conda activate ppi_surface
python -c "import pyKVFinder; print('pyKVFinder OK')"
```

Then run adapter smoke on fixture.

### Task TOOL-3 — Implement mini-FTMap

No install. Use Vina/RDKit/OpenBabel.

Acceptance:

```text
mini_ftmap_probe_hotspots.csv exists for fixture/local box
probe clusters are assigned
ATP escape is not observed in fixture smoke
```

### Task TOOL-4 — Test PeSTo optional

Install isolated `pesto` env. Run repository example. If it fails, mark optional disabled. If it passes, create adapter to read B-factor output.

### Task TOOL-5 — Test InDeep optional

Clone in external path. Follow README. Run repository example only. If stable, create adapter.

### Task TOOL-6 — Evaluate MaSIF/MSMD/PocketMiner/PASSer

Do not force these into core. Create status rows and only implement adapters if installation is clean.

### Task TOOL-7 — Integrate evidence table

Implement:

```text
ppi_surface_zone_evidence.csv builder
accepted_ppi_surface_zones_for_m3.csv exporter
zone_composite_boxes.csv exporter
```

---

## 10. Logging, cleanup, and safety requirements

Every tool installation/test must produce:

```text
install log
smoke test log
status JSON
error summary
cleanup report
```

Suggested layout:

```text
fresh/runs/<run_id>/logs/tools/<tool>.install.log
fresh/runs/<run_id>/logs/tools/<tool>.smoke.log
fresh/runs/<run_id>/manifest/tool_status.json
fresh/runs/<run_id>/reports/tool_installation_report.md
fresh/runs/<run_id>/cleanup_report.json
```

Cleanup policy:

```text
Delete:
    temporary tool scratch
    temporary probe PDBQT files from smoke tests
    intermediate test docking files
    failed partial outputs outside summary tables

Keep:
    logs
    manifests
    status JSON
    summary CSVs
    representative fixture outputs if needed
```

Never delete outside:

```text
fresh/runs/<run_id>/
```

unless the user explicitly approves.

---

## 11. How this narrows many sites to one or two zones

The final narrowing strategy is:

```text
1. PyRosetta PPI consensus identifies where MYO1D repeatedly contacts EGFR.
2. PPI hotspot analysis identifies which EGFR residues matter energetically for MYO1D binding.
3. PeSTo/MaSIF/InDeep optionally ask whether that surface looks like a PPI/ligandable PPI site.
4. fpocket/pyKVFinder identify whether nearby surface grooves/cavities exist.
5. mdpocket asks whether these grooves persist across MD-derived states or trajectory frames.
6. mini-FTMap asks whether small chemical probes repeatedly cluster in the same surface subpocket.
7. MSMD optionally validates top zones dynamically with explicit probe solvent.
8. PocketMiner/PASSer optionally support cryptic/allosteric interpretation.
9. Hard gates remove ATP, membrane-inaccessible, dimer-buried, and PPI-unrelated sites.
10. The remaining top one or two zones become official M3 fresh virtual screening targets.
```

This is a filter stack, not a voting contest. A tool can raise confidence, but it cannot override project hard gates.

---

## 12. Deliverables expected from the agent

The agent must deliver:

```text
1. Updated configs:
   fresh/configs/tool_registry.yaml
   fresh/configs/pocket_zone_scoring.yaml

2. Preflight code:
   fresh/src/egfr_myo1d/tools/tool_preflight.py
   fresh/scripts/preflight_tools.sh
   fresh/scripts/submit_tool_preflight.sh

3. Tool adapters, only where installation passes:
   mdpocket_adapter.py
   pykvfinder_adapter.py
   mini_ftmap_adapter.py
   optional_pesto_adapter.py
   optional_indeep_adapter.py
   optional_masif_adapter.py
   optional_msmd_adapter.py
   optional_pocketminer_adapter.py
   optional_passer_adapter.py

4. Output tables:
   tool_status.json
   ppi_surface_zone_evidence.csv
   accepted_ppi_surface_zones_for_m3.csv
   zone_composite_boxes.csv

5. Reports:
   tool_installation_report.md
   milestone2_surface_zone_filtering_summary.md
   milestone3_public_vs_target_definition.md
```

---

## 13. Acceptance criteria

### 13.1 Tool preflight acceptance

```text
- Existing fpocket, Vina, OpenBabel, RDKit, PyRosetta status recorded.
- mdpocket status recorded, whether available or not.
- pyKVFinder installed or marked with clear failure reason.
- mini-FTMap dependencies pass.
- Optional tools are marked available/disabled with reasons.
- tool_status.json exists.
- No external tool source committed to repo.
```

### 13.2 M2 integration acceptance

```text
- PPI consensus patch table exists.
- PPI hotspot table exists.
- local cavity table exists from pyKVFinder or fallback fpocket.
- mini-FTMap probe hotspot table exists.
- mdpocket persistence table exists or is explicitly skipped due no trajectory/mdpocket.
- ppi_surface_zone_evidence.csv exists.
- accepted_ppi_surface_zones_for_m3.csv contains only hard-gate-passing zones.
```

### 13.3 M3 integration acceptance

```text
- M3 targets accepted PPI-compatible zones, not generic pockets.
- Fresh VS docking boxes come from zone_composite_boxes.csv.
- ATP migration filter exists.
- final_public_candidate_shortlist.csv contains no private compound IDs.
- private probe outputs, if any, are outside public workflow and not mixed into public scoring.
```

---

## 14. References and source notes for the agent

Project-specific basis:

- `PROJECT_KNOWLEDGE_FINAL_CLEAN.md` defines the EGFR-MYO1D objective, receptor model, MYO1D construct, and non-ATP/C-lobe PPI strategy.
- `milestone1_foundation_codex_handoff_v0_5.md` defines the fresh workflow, HPC, logging, cleanup, and input normalization foundation.
- `milestone2_ppi_to_pocket_detailed_plan_v0_2.md` defines the original PPI-to-pocket flow.
- `milestone3_compound_anchor_integration_detailed_plan_v0_3.md` should be superseded for public workflow by fresh virtual screening. Private compounds are internal only.

Tool source notes:

- fpocket/mdpocket: fpocket suite; mdpocket tracks pocket occurrence over MD trajectories and can output frequency/density grids and descriptors.
- pyKVFinder: Python package for biomolecular cavity detection and characterization.
- PeSTo: geometric deep learning for protein binding interface prediction; output interface scores can be stored in B-factor field.
- MaSIF: molecular surface interaction fingerprints; includes MaSIF-site for PPI-site surface prediction.
- InDeep: PPI drug discovery toolbox for interactibility/ligandability site prediction.
- MSMD: mixed-solvent MD for dynamic probe hotspot mapping; EXPRORER_MSMD is one GROMACS-based automation option.
- PocketMiner: cryptic pocket opening predictor using graph neural networks; optional until stable local execution is confirmed.
- PASSer: allosteric site predictor; mainly web/API and must respect privacy constraints.

---

## 15. Final instruction to the implementation agent

Do not treat this as a request to install every tool at once. Treat it as a controlled, logged, qsub-tested, cleanup-safe integration plan. The core workflow must remain usable even if all optional AI tools fail to install.

The minimum scientifically useful extension is:

```text
PyRosetta PPI consensus
+ PyRosetta PPI hotspot analysis
+ fpocket/mdpocket status
+ pyKVFinder local cavity descriptors
+ mini-FTMap local probe hotspots
+ hard gates
= accepted PPI-compatible pocket zones for M3 fresh virtual screening
```

Everything else is optional evidence.
