# Fresh Scripts

These scripts are safe placeholders for Milestone 1 Task 1. They do not submit jobs or run docking.

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

`submit_smoke_env.sh` and `submit_smoke_input.sh` document later smoke behavior but do not call `qsub` in Task 1.
