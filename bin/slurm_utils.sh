#!/bin/bash
# ============================================================================
# SLURM Utility Functions for Multi-Ancestry GWAS Pipeline
# ============================================================================
# Source this file in SLURM job scripts:
#   source /path/to/slurm_utils.sh
# ============================================================================

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

# Create SLURM array job script
# Usage: create_array_job "my_job" "EUR,AAC,LAT1,LAT2" "run_analysis.sh"
create_array_job() {
    local job_name="$1"
    local items="$2"
    local script="$3"
    local partition="${4:-normal}"
    local time="${5:-24:00:00}"
    local mem="${6:-32G}"
    local cpus="${7:-8}"

    # Count items
    IFS=',' read -ra arr <<< "$items"
    local n_items=${#arr[@]}
    local max_idx=$((n_items - 1))

    cat << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --array=0-${max_idx}
#SBATCH --partition=${partition}
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --output=logs/${job_name}_%A_%a.out
#SBATCH --error=logs/${job_name}_%A_%a.err

# Source utilities
source \$(dirname \$0)/slurm_utils.sh

# Print job info
print_slurm_info

# Set thread environment
set_thread_env

# Get item for this array task
ITEM=\$(get_array_item "${items}")
echo "Processing: \$ITEM"

# Run the analysis script
${script} \$ITEM

echo "Completed: \$ITEM"
EOF
}

# Ancestry group definitions for this cohort
ANCESTRIES="EUR,AAC,LAT1,LAT2,OTHER"
PRIMARY_ANCESTRIES="EUR,AAC,LAT1,LAT2"
MIN_STRATUM_N=30

# GWAS types to run
GWAS_TYPES="tractor_la,standard,ancestry_specific"

# PRS methods
PRS_METHODS="prs_csx,gaudi,disco_divas,sdpr_admix,mussel,prosper"

# QTL types
QTL_TYPES="eqtl,sqtl,pqtl,mqtl,caqtl,hqtl"

# Traits
TRAITS="OS,relapse,MRD"

# Export for use in other scripts
export ANCESTRIES PRIMARY_ANCESTRIES MIN_STRATUM_N
export GWAS_TYPES PRS_METHODS QTL_TYPES TRAITS
