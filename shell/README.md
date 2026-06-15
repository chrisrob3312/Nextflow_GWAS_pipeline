# Shell Wrappers for Ancestry-Aware GWAS Pipeline

## Dual Shell + Nextflow Structure

This directory contains shell wrapper scripts that provide command-line access to the same core analysis scripts used by the Nextflow modules. This dual structure prevents code drift between interactive debugging and production pipeline runs.

### Architecture

```
pipeline/
├── bin/                          # Core analysis scripts (R, Python)
│   ├── tractor_genesis_adapter.R # Core Tractor→GENESIS logic
│   ├── run_coloc_analysis.R      # Core colocalization
│   └── ...
│
├── shell/                        # Shell wrappers for interactive use
│   ├── tractor_genesis_adapter.sh
│   ├── run_coloc.sh
│   └── ...
│
└── modules/local/                # Nextflow thin wrappers
    ├── tractor_genesis_adapter/
    │   └── main.nf               # Calls bin/tractor_genesis_adapter.R
    └── ...
```

### Key Principles

1. **Single Source of Truth**: Core analysis logic lives in `bin/` scripts only
2. **Thin Wrappers**: Both shell scripts and Nextflow modules call the same `bin/` script
3. **Equivalent I/O**: Shell and Nextflow produce identical outputs for identical inputs
4. **Debug-Friendly**: Use shell scripts for interactive debugging, then run via Nextflow
5. **No Drift**: Changes to analysis logic only need to be made in one place

### Usage

#### Interactive Debugging
```bash
# Debug a specific sample with shell wrapper
./shell/tractor_genesis_adapter.sh \
    -g data/sample.gds \
    -p data/phenotypes.tsv \
    -t relapse \
    -m survival \
    -a EUR,AFR,AMR \
    --time_col os_time \
    --event_col os_event \
    -v

# Examine outputs, iterate on parameters
less results/tractor_genesis.log
```

#### Production Pipeline
```bash
# Run full pipeline via Nextflow
nextflow run main.nf \
    --input samplesheet.csv \
    --phenotype_file phenotypes.tsv \
    --run_tractor true \
    --tractor_groups 'LATINO' \
    -profile slurm
```

### Available Shell Wrappers

| Script | Description | Corresponding Module |
|--------|-------------|---------------------|
| `tractor_genesis_adapter.sh` | Tractor→GENESIS survival/binary/quantitative | `tractor_genesis_adapter` |
| `run_coloc.sh` | Colocalization analysis (planned) | `coloc` |
| `run_env_mr_mega.sh` | env-MR-MEGA meta-analysis (planned) | `mr_mega_env` |
| `calculate_ld.sh` | Cohort-specific LD calculation (planned) | `ld_calculation` |

### Adding New Wrappers

When adding a new analysis:

1. **Create core script in `bin/`**
   - Self-contained with all analysis logic
   - Uses argparse/optparse for CLI arguments
   - Produces consistent output formats
   - Includes version reporting

2. **Create shell wrapper in `shell/`**
   - Handles argument parsing and validation
   - Calls the core `bin/` script
   - Provides helpful usage messages
   - Logs output for debugging

3. **Create Nextflow module in `modules/local/`**
   - Thin wrapper that calls `bin/` script
   - Handles Nextflow-specific I/O
   - Uses same arguments as shell wrapper

### Example: Debugging a Failed Run

```bash
# 1. Identify the failed process from Nextflow logs
# Process: TRACTOR_GENESIS_SURVIVAL (sample_123.LATINO.os)

# 2. Find the work directory
ls -la work/ab/cd1234.../

# 3. Recreate the run with shell wrapper
./shell/tractor_genesis_adapter.sh \
    -g work/ab/cd1234.../sample_123.tractor_gds \
    -p work/ab/cd1234.../sample_123.phenotype.tsv \
    -t os \
    -m survival \
    -a EUR,AFR,AMR \
    -v

# 4. Debug interactively, fix issues in bin/ script
# 5. Re-run Nextflow - fix applies to both shell and Nextflow
```

### Testing Equivalence

```bash
# Run shell wrapper
./shell/tractor_genesis_adapter.sh -g test.gds -p test.tsv -t trait \
    -o results_shell/test

# Run Nextflow (single process)
nextflow run main.nf -entry test_tractor_genesis \
    --gds test.gds --phenotype test.tsv --trait trait \
    --outdir results_nextflow

# Compare outputs
diff <(zcat results_shell/test.joint.tsv.gz | sort) \
     <(zcat results_nextflow/test.joint.tsv.gz | sort)
```
