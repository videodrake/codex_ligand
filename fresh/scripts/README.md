# Fresh Scripts

These scripts generate safe workflow helpers. They do not submit jobs or run docking unless an explicit submit flag and confirmation flag are both provided by the user.

Future PBS snippets must include:

```bash
#PBS -q workq
#PBS -l nodes=node04:ppn=4
#PBS -l walltime=02:00:00
```

Future PBS jobs and shell smoke wrappers must also include:

```bash
export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta
```

`submit_smoke_env.sh` and `submit_smoke_input.sh` document smoke behavior but do not call `qsub`.

`submit_m2_pyrosetta_real_jobs.py` writes the M2.2 real PyRosetta PBS plan for node-level chunked execution. The production default is intended for `node04,node05,node06` with `ppn=32`; it keeps the useful `PBS_NP` worker-cap contract but runs each chunk as an isolated PyRosetta subprocess. Chunk size is automatic by default, roughly one worker wave per state/seed. It does not call `qsub` unless both `--submit` and `--i-understand-this-submits-hpc-jobs` are supplied.
