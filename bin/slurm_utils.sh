#!/bin/bash
# ============================================================================
# SLURM Utility Functions for Multi-Ancestry GWAS Pipeline
# ============================================================================
# Source this file in SLURM job scripts:
#   source /path/to/slurm_utils.sh
# ============================================================================

# The run record (rr_begin, ...; bin/run_record.sh), so every job that sources this file has it
source "$(dirname "${BASH_SOURCE[0]}")/run_record.sh"

# Check if running under SLURM
is_slurm_job() {
    [[ -n "${SLURM_JOB_ID:-}" ]]
}

# Get array task info
get_array_info() {
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
        echo "SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT:-1}"
        echo "SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN:-0}"
        echo "SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX:-0}"
    fi
}

# Get item from array based on SLURM task ID
# Usage: item=$(get_array_item "EUR,AAC,LAT1,LAT2,OTHER")
get_array_item() {
    local items="$1"
    local delimiter="${2:-,}"

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "ERROR: Not running as SLURM array job" >&2
        return 1
    fi

    # Convert to array
    IFS="$delimiter" read -ra arr <<< "$items"

    # Get item (0-indexed adjustment if needed)
    local idx=$((SLURM_ARRAY_TASK_ID))
    if [[ $idx -ge ${#arr[@]} ]]; then
        echo "ERROR: Array index $idx out of range (max: $((${#arr[@]}-1)))" >&2
        return 1
    fi

    echo "${arr[$idx]}"
}

# Print SLURM job info
print_slurm_info() {
    echo "=============================================="
    echo "SLURM Job Information"
    echo "=============================================="
    echo "Job ID:        ${SLURM_JOB_ID:-N/A}"
    echo "Job Name:      ${SLURM_JOB_NAME:-N/A}"
    echo "Partition:     ${SLURM_JOB_PARTITION:-N/A}"
    echo "Node:          ${SLURMD_NODENAME:-N/A}"
    echo "CPUs:          ${SLURM_CPUS_PER_TASK:-N/A}"
    echo "Memory:        ${SLURM_MEM_PER_NODE:-N/A} MB"
    echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
    echo "Working Dir:   $(pwd)"
    echo "=============================================="
}

# Set thread counts for common tools
set_thread_env() {
    local threads="${SLURM_CPUS_PER_TASK:-4}"

    export OMP_NUM_THREADS=$threads
    export MKL_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=$threads
    export NUMEXPR_NUM_THREADS=$threads
    export VECLIB_MAXIMUM_THREADS=$threads

    # R
    export R_THREADS=$threads
    export MC_CORES=$threads

    # PLINK
    export PLINK_THREADS=$threads

    echo "Set thread count to $threads for parallel operations"
}

# Memory-aware chunk processing
# Usage: process_with_memory_limit <total_items> <mem_per_item_mb>
calculate_chunk_size() {
    local total_items=$1
    local mem_per_item_mb=${2:-100}

    local available_mem_mb=${SLURM_MEM_PER_NODE:-16000}
    local max_items=$((available_mem_mb / mem_per_item_mb))

    # Leave some headroom
    max_items=$((max_items * 8 / 10))

    if [[ $max_items -gt $total_items ]]; then
        echo $total_items
    else
        echo $max_items
    fi
}

# Check if we should run this stratum based on sample size
# Usage: should_run_stratum "EUR" 150 30
should_run_stratum() {
    local stratum="$1"
    local n_samples="$2"
    local min_n="${3:-30}"

    if [[ $n_samples -ge $min_n ]]; then
        return 0  # Yes, run it
    else
        echo "Stratum $stratum has N=$n_samples < min_n=$min_n, skipping" >&2
        return 1  # No, skip it
    fi
}

# Print a command so it can be read, and pasted back into a shell: plain words as they are,
# anything with spaces or special characters quoted. Used by the DRY_RUN modes.
print_cmd() {
    local a out=""
    for a in "$@"; do
        if [[ "$a" =~ ^[A-Za-z0-9_./:=,%@+-]+$ ]]; then out+=" $a"; else out+=" $(printf '%q' "$a")"; fi
    done
    echo "${out# }"
}

# Ancestry labels for one stratum's Tractor decomposition (the adapter's --ancestries), read
# from data/tractor/<stratum>/ancestries.txt: one label per line, as in the handoff bundle's
# local_ancestry/ancestry_labels.txt (commas also accepted; blank lines and # comments ignored).
# Prints the labels comma-separated. Fails, naming the file, if it is missing or lists nothing,
# or if a label is malformed or repeated. There is no default list on purpose.
# Usage: anc="$(read_stratum_ancestries data/tractor/EUR/ancestries.txt)" || exit 1
read_stratum_ancestries() {
    local file="$1" labels bad dup
    if [[ ! -f "$file" ]]; then
        echo "ERROR: ${file} does not exist; every stratum needs its ancestries.txt (docs/HANDOFF_CONTRACT.md)" >&2
        return 1
    fi
    labels="$(awk '{ sub(/#.*/, ""); gsub(/\r/, ""); n = split($0, a, /[ \t,]+/)
                     for (i = 1; i <= n; i++) if (a[i] != "") print a[i] }' "$file")"
    if [[ -z "$labels" ]]; then
        echo "ERROR: ${file} lists no ancestry" >&2
        return 1
    fi
    bad="$(awk '!/^[A-Za-z][A-Za-z0-9_]*$/ { print; exit }' <<< "$labels")"
    if [[ -n "$bad" ]]; then
        echo "ERROR: ${file}: '${bad}' is not an ancestry label" >&2
        return 1
    fi
    dup="$(awk 'seen[$0]++ { print; exit }' <<< "$labels")"
    if [[ -n "$dup" ]]; then
        echo "ERROR: ${file} lists ${dup} twice" >&2
        return 1
    fi
    paste -sd, - <<< "$labels"
}

# GWAS types to run
GWAS_TYPES="tractor_la,standard,ancestry_specific"

# PRS methods
PRS_METHODS="prs_csx,gaudi,disco_divas,sdpr_admix,mussel,prosper"

# QTL types
QTL_TYPES="eqtl,sqtl,pqtl,mqtl,caqtl,hqtl"

# Export for use in other scripts. There is no cohort-wide ancestry list here (each stratum's comes
# from its own ancestries.txt, read_stratum_ancestries above) and no trait list: traits come only
# from the launcher's switch table (bin/run_record.sh); the stale list that was here is gone.
export GWAS_TYPES PRS_METHODS QTL_TYPES
