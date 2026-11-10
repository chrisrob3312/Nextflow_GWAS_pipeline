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

## Steps

The steps (the site variables, building the four conda envs, `check_envs`, `test_job_env`, the dry
run, the run-record self-test) are in the combined runbook: the harmonization repo's
`genotyping_pipeline/benchmarking/SMOKE_RUNBOOK.md`, part "GWAS pipeline".

## What each setting does

- **`SLURM_PARTITION`, `SLURM_ACCOUNT`.** The launcher puts both on every sbatch line and refuses
  to start without them; an `#SBATCH` line cannot read a variable, so the scripts carry neither.
- **`LOG_DIR`.** Every job's `--output` / `--error` is an absolute path under it.
- **`GWAS_SCRATCH`, `SLURM_TMPDIR`.** See TMPDIR below.
- **`GWAS_ENVS_DIR`.** `envs.tsv` gives each env's prefix as `${GWAS_ENVS_DIR}/<env>`;
  `bin/activate_env.sh` stops, naming the env and the path, if the prefix is not there.
- **`GWAS_CONDA_MODULE`.** If a job's shell has no `conda`, `bin/activate_env.sh` runs
  `module load` with it first.

## Submitting one script by hand

The scripts under `slurm/` carry no partition, account or log path. Submit them from the repo root
and give those yourself:

```bash
sbatch --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
       --output="$PWD/logs/gwas_%A_%a.out" --error="$PWD/logs/gwas_%A_%a.err" \
       slurm/submit_gwas_array.sh EUR OS survival
```

(or export `SBATCH_PARTITION` and `SBATCH_ACCOUNT`, which sbatch reads by itself). From the repo
root, because Slurm runs a job script from a copy in its spool directory: each script finds the
repo through `GWAS_REPO` (exported by the launcher) or `SLURM_SUBMIT_DIR`, and stops with "cannot
find the repo" otherwise.

## TMPDIR

Every job sources `bin/activate_env.sh`, which exports `TMPDIR` on disk, never on a tmpfs:
`data.table::fread()` decompresses a `.gz` into `TMPDIR` first, so a chromosome-sized file in a
tmpfs `/tmp` fills RAM and fails with "No space left on device". `TMPDIR` is the first of these that
is a writable directory and not tmpfs/ramfs: `$SLURM_TMPDIR`; `$GWAS_SCRATCH/tmp` (created if
missing); the job's existing `TMPDIR`, or `/tmp`. If none qualifies the job stops before running
anything. R, Python and `sort` honour `TMPDIR`; Java does not, so a FLARE step passes
`-Djava.io.tmpdir="$TMPDIR"` itself.

## Building the envs

The build steps are in the runbook (see Steps). Installed by hand, once per site, after the envs
are built:

- **hyprcoloc** into `gwas-r` (needed only with `RUN_COLOC=true`): the command is in the header of
  `envs/gwas-r.yml`;
- **flare.jar** and Tractor's **extract_tracts_flare.py** into `lai`: where to get them is in the
  header of `envs/lai.yml`;
- **env-MR-MEGA** into `mrmega` (CP-29): it is in no conda channel; build it from the authors'
  release (Wang et al., Communications Biology 2024, PMID 39543362) and copy the binary to
  `$GWAS_ENVS_DIR/mrmega/bin/`. MR-MEGA itself comes with the env (bioconda `mr-mega=0.2`).

## Site notes

(What was fetched by hand and when: flare.jar version, Tractor commit, hyprcoloc commit,
env-MR-MEGA version.)
