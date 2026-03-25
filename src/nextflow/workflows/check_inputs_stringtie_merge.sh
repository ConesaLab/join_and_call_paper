#!/bin/bash

# Checks that all metadata_ind files needed by the StringTie merge workflow
# exist and that the per-sample GTFs referenced inside them are reachable.
# Run interactively on garnatxa (no SLURM needed).

PB_BASE="/home/fabianje/repos/documenting_NIH/fabian/data/output/isoseq"
ONT_BASE="/home/fabianje/repos/documenting_NIH/fabian/data/output/ont"
SCRIPT_DIR=$(dirname "$(realpath "$0")")
CORRECTED_META="${SCRIPT_DIR}/../scripts/stringtie_merge/metadata"

declare -A PATHS

# PacBio (corrected metadata for flair/isoseq/mandalorion)
PATHS[pb_flair]="${CORRECTED_META}/pb_flair_metadata_ind.tsv"
PATHS[pb_isoquant]="${PB_BASE}/isoquant/run3_data/isoquant/isoquant_out_metadata_ind.tsv"
PATHS[pb_bambu]="${PB_BASE}/bambu/run3_data/bambu/bambu_out_metadata_ind.tsv"
PATHS[pb_talon]="${PB_BASE}/talon/run3_data/talon/talon_out_metadata_ind.tsv"
PATHS[pb_isoseq]="${CORRECTED_META}/pb_isoseq_metadata_ind.tsv"
PATHS[pb_mandalorion]="${CORRECTED_META}/pb_mandalorion_metadata_ind.tsv"

# ONT (corrected metadata for flair)
PATHS[ont_flair]="${CORRECTED_META}/ont_flair_metadata_ind.tsv"
PATHS[ont_isoquant]="${ONT_BASE}/isoquant/run2_data/isoquant/isoquant_out_metadata_ind.tsv"
PATHS[ont_bambu]="${ONT_BASE}/bambu/run2_data/bambu/bambu_out_metadata_ind.tsv"
PATHS[ont_talon]="${ONT_BASE}/talon/run1_data/talon/talon_out_metadata_ind.tsv"

any_failures=0

echo "=============================================="
echo " Checking metadata_ind files"
echo "=============================================="
for key in $(echo "${!PATHS[@]}" | tr ' ' '\n' | sort); do
    path="${PATHS[$key]}"
    if [ -f "$path" ]; then
        n_samples=$(grep -vc '^#' "$path")
        echo "[OK]         $key ($n_samples data lines)"
    elif [ -L "$path" ]; then
        echo "[DEAD LINK]  $key: $path"
        any_failures=1
    else
        echo "[MISSING]    $key: $path"
        any_failures=1
    fi
done

echo ""
echo "=============================================="
echo " Checking per-sample .expr.gtf files (preferred input)"
echo "=============================================="
for key in $(echo "${!PATHS[@]}" | tr ' ' '\n' | sort); do
    path="${PATHS[$key]}"
    if [ ! -f "$path" ]; then
        echo "[SKIP]       $key (metadata file missing)"
        continue
    fi

    expr_ok=0
    expr_miss=0
    gtf_fallback=0
    both_miss=0
    while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq_file sample_gtf rest; do
        expr_gtf="${sample_gtf%.gtf}.expr.gtf"
        if [ -f "$expr_gtf" ]; then
            ((expr_ok++))
        elif [ -f "$sample_gtf" ]; then
            ((gtf_fallback++))
            if [ "$gtf_fallback" -le 3 ]; then
                echo "  [NO .expr] $key $sample_id: $expr_gtf (will fall back to .gtf)"
            fi
        else
            ((both_miss++))
            if [ "$both_miss" -le 3 ]; then
                echo "  [MISSING]  $key $sample_id: neither $expr_gtf nor $sample_gtf"
            fi
        fi
    done < <(grep -v '^#' "$path")

    total=$((expr_ok + gtf_fallback + both_miss))
    if [ "$both_miss" -eq 0 ] && [ "$gtf_fallback" -eq 0 ]; then
        echo "[OK]         $key: all $expr_ok .expr.gtf found"
    elif [ "$both_miss" -eq 0 ]; then
        echo "[WARN]       $key: $expr_ok .expr.gtf, $gtf_fallback .gtf fallback (usable but not filtered)"
    else
        echo "[PROBLEM]    $key: $expr_ok .expr.gtf, $gtf_fallback .gtf fallback, $both_miss MISSING"
        any_failures=1
    fi
done

echo ""
echo "=============================================="
if [ "$any_failures" -eq 0 ]; then
    echo " All checks passed."
else
    echo " Some checks FAILED. Review the output above."
fi
echo "=============================================="

exit $any_failures
