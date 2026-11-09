# HPC site settings

Everything site-specific lives in this file and in the environment of the shell you launch from,
never in a script: commit b531c14 removed the last site paths on purpose. Fill in the blanks once
per cluster, then export the variables before launching (for example from a small file you
`source`, kept outside the repo).

## Site settings (fill in)

| setting | variable | value at this site | used by |
|---|---|---|---|
| partition | `SLURM_PARTITION` | ________ | every sbatch; the launcher refuses to run without it |
| account | `SLURM_ACCOUNT` | ________ | every sbatch; the launcher refuses to run without it |
| log directory | `LOG_DIR` | default `<repo>/logs` | every job's `--output` / `--error` (always an absolute path) |
| largest array | (none) | ________ (`scontrol show config \| grep MaxArraySize`) | the 22-task chromosome arrays |
| QOS limits | (none) | ________ (`sacctmgr show qos`) | how many arrays can run at once |

If the site has no Slurm accounts, set `SLURM_ACCOUNT` to your default account
(`sacctmgr show user $USER`), so the run record still names one.

## Launching

```bash
export SLURM_PARTITION=<partition> SLURM_ACCOUNT=<account>
slurm/submit_full_pipeline.sh --dry-run   # prints every sbatch line, submits nothing
slurm/submit_full_pipeline.sh
```

The launcher can be started from any directory: it changes to the repo root first. It adds
`--partition`, `--account`, `--export=ALL` and `--kill-on-invalid-dep=yes` to every sbatch it
issues, and gives every job an absolute `--output` (and `--error`) under `LOG_DIR`. An `#SBATCH`
line cannot read a shell variable, which is why the partition and account are on the command line
and not in the scripts.

## Submitting one script by hand

The scripts under `slurm/` carry no partition, account or log path. Submit them from the repo root
and give those three yourself:

```bash
cd <repo>
mkdir -p logs
sbatch --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
       --output="$PWD/logs/gwas_%A_%a.out" --error="$PWD/logs/gwas_%A_%a.err" \
       slurm/submit_gwas_array.sh EUR OS survival
```

Or let sbatch pick the first two up by itself: it reads `SBATCH_PARTITION` and `SBATCH_ACCOUNT`
from the environment.

```bash
export SBATCH_PARTITION="$SLURM_PARTITION" SBATCH_ACCOUNT="$SLURM_ACCOUNT"
sbatch --output="$PWD/logs/pcair_%j.out" slurm/submit_pcair.sh
```

Without `--output`, Slurm writes `slurm-<jobid>.out` in the directory you submitted from.

Why the repo root: Slurm runs a job script from a copy in its own spool directory, so a script
cannot find the repo from its own path. Each script looks for the repo in `GWAS_REPO` (exported by
the launcher), then in `SLURM_SUBMIT_DIR` (the directory sbatch was run in), and stops with
"cannot find the repo" if neither holds it.
