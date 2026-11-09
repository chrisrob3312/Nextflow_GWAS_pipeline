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
| scratch | `GWAS_SCRATCH` | ________ | `TMPDIR` is `$GWAS_SCRATCH/tmp` unless the site sets `SLURM_TMPDIR` (see TMPDIR below) |
| does the site set `SLURM_TMPDIR` in jobs? | `SLURM_TMPDIR` | yes / no, filesystem: ________ | first choice for `TMPDIR` |
| what `/tmp` is on compute nodes | (none) | ________ (`stat -f -c %T /tmp` in a job) | if tmpfs, `GWAS_SCRATCH` is required |
| conda env directory | `GWAS_ENVS_DIR` | ________ | holds the four env prefixes, `$GWAS_ENVS_DIR/<env>` (envs.tsv) |
| how conda gets into a job | `GWAS_CONDA_MODULE` | ________ (module name, or blank if `conda` is on PATH) | `bin/activate_env.sh` runs `module load` with it when `conda` is missing |
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

## TMPDIR (M-04)

Every job sources `bin/activate_env.sh`, which exports `TMPDIR` before anything else runs. The
rule: `TMPDIR` is on disk, never on a tmpfs. `data.table::fread()` decompresses a `.gz` into
`TMPDIR` before reading it, so a chromosome-sized file in a tmpfs `/tmp` fills RAM and fails with
"No space left on device" (ENOSPC). The first of these that is a writable directory and not
tmpfs/ramfs becomes `TMPDIR`:

1. `$SLURM_TMPDIR`, if the site sets it (a per-job directory the site cleans up);
2. `$GWAS_SCRATCH/tmp`, created if missing;
3. the `TMPDIR` the job already had, or `/tmp`, with a note in the log.

If none qualifies the job stops before running anything, and says what it rejected. R
(`tempdir()`), Python (`tempfile`) and `sort` all honour `TMPDIR`. Java does not, so a step that
runs FLARE should pass `-Djava.io.tmpdir="$TMPDIR"` itself.

## Building the envs

The four specs live in `envs/` (`gwas-r`, `gwas-py`, `lai`, `mrmega`); `envs.tsv` lists each
env's prefix and the sha256 of the yml it was built from. Build them once per site, into
`GWAS_ENVS_DIR`:

```bash
export GWAS_ENVS_DIR=<site path>
for e in gwas-r gwas-py lai mrmega; do
    CONDA_CHANNEL_PRIORITY=strict conda env create -p "$GWAS_ENVS_DIR/$e" -f "envs/$e.yml"
done
bin/envs_tsv.sh        # records the sha256 of each yml; no change if the ymls are the committed ones
```

`mamba env create` takes the same arguments. Strict channel priority is what bioconda asks for.
After changing a yml: rebuild that env (`conda env remove -p ...`, then create it again), run
`bin/envs_tsv.sh`, and commit the yml and `envs.tsv` together.

Two things are not conda packages and go in by hand:

- **flare.jar** for the `lai` env: download it as the FLARE README says
  (`https://faculty.washington.edu/browning/flare.jar`) to `$GWAS_ENVS_DIR/lai/share/flare/`,
  and Tractor's `scripts/extract_tracts_flare.py` (github.com/Atkinson-Lab/Tractor, commit
  2860728) to `$GWAS_ENVS_DIR/lai/bin/`. Record the versions you fetched in the site notes below.
- **env-MR-MEGA** for the `mrmega` env (CP-29): build it from the authors' release (Wang et al.,
  Communications Biology 2024, PMID 39543362) and copy the binary to `$GWAS_ENVS_DIR/mrmega/bin/`.
  MR-MEGA itself comes from bioconda (`mr-mega=0.2`).

Inside a job, an env is used with `source bin/activate_env.sh <env>`; it stops with the env name
and the missing path if the prefix is not there.

## Checking the envs and a job's environment

From the repo root, with the settings above exported:

```bash
mkdir -p logs
sbatch --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
       --output="$PWD/logs/check_envs_%j.out" bin/check_envs.sh
sbatch --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
       --output="$PWD/logs/test_job_env_%A_%a.out" bin/test_job_env.sh
```

`check_envs.sh` loads every R package of the SLURM path in `gwas-r` (the 14 of the plan plus
`R.utils`, which `fread` needs for any `.gz`), the Python imports in `gwas-py` and `lai`, `java`
in `lai` and `MR-MEGA` in `mrmega`, and checks `envs.tsv` against the ymls. Exit 0 and a last
line `check_envs: PASS` mean all of it loaded. A missing `coxme` or `metafor` is called out on its
own line: without them the adapter and the meta-analysis would still run, on a weaker model.

`test_job_env.sh` prints `SLURM_JOB_ID`, `SLURM_ARRAY_TASK_ID` and `SLURM_CPUS_PER_TASK` after
`conda activate`, from bash and from R, checks that R's `tempdir()` is inside `TMPDIR`, and reads a
~2 GB `.gz` with `fread()`. Set `TEST_GZ=/path/to/file.gz` to read a real chromosome file instead
of the generated one. Its last line is `test_job_env: PASS` or `FAIL`.

## Site notes

(What was fetched by hand and when: flare.jar version, Tractor commit, env-MR-MEGA version.)
