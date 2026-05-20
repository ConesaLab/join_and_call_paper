# ONT R10 SY5Y — Pychopper alternative branch

Parallel preprocessing to [`../ont_r10/`](../ont_r10/) (Dorado-trimmed FASTQs + minimap2 `-ub`). This branch runs **Pychopper** (`-k LSK114`), merges full-length + rescued reads, then maps with minimap2 **`-uf`**. The default `ont_r10` scripts are **not modified**.

Cluster data root: `/storage/gge/Fabian/ont_r10_sy5y/`

| Path | Role |
|------|------|
| `fastq/` | Input (Dorado-trimmed; shared, read-only) |
| `fastq_pychopper/` | Pychopper outputs + `*_for_map.fastq` |
| `bam_pychopper/` | Alignments |
| `gff_pychopper/` | Transcript GFF / SIRV subsets |
| `analysis/logs_pychopper/` | Slurm logs |
| `analysis/*_pychopper/` | QC, SQANTI, SQANTI-reads, read_qc |

Paths are in [`config.sh`](config.sh). Job scripts use relative `list_fastqs*.fof` and `source config.sh`, same pattern as [`../ont_r10/`](../ont_r10/).

## Kit flag

`pychopper -k LSK114` (SQK-LSK114). Do **not** pass `"ONT SQK-LSK114"` as `-k`.

## Submit order

Run `sbatch` from **this directory** (Slurm sets the job cwd to `SLURM_SUBMIT_DIR`):

```bash
cd src/preprocessing/ont_r10_pychopper

sbatch 0_pychopper.sh
sbatch 1_merge_pychopper_fastqs.sh   # after 0 completes
sbatch 2_map.sh                      # after 1 completes

# Optional downstream (after 2_map)
sbatch 3_fastq_qc.sh
sbatch 4_bam_qc.sh
sbatch 5_run_sqanti.sh
sbatch 6_run_sqanti_SIRVS.sh
sbatch 7_sqanti_reads.sh
sbatch 8_sqanti_reads_merge.sh
sbatch 9_concat_samples.sh
sbatch 10_count_reads_joint.sh
sbatch 11_get_read_lengths_bam.sh
```

Conda: `pychopper` (step 0), `SQANTI3.env` (mapping and SQANTI steps).

## Validation checklist

1. **Pychopper QC:** `fastq_pychopper/*_pychopper_report.pdf`
2. **Read counts:** `wc -l` on `*_for_map.fastq` vs `fastq/*.fastq`
3. **Mapping:** `samtools flagstat bam_pychopper/*` vs `bam/*`
4. **Strand:** `samtools view -f 16` spot-check on pychopper BAMs
5. **Kit:** Confirm ENA library prep; re-run step 0 with correct `-k` in `config.sh` if needed

## References

- Mouse Pychopper + `-uf`: [`../ont/3_pychopper.sh`](../ont/3_pychopper.sh), [`../ont/4_map.sh`](../ont/4_map.sh)
- Default R10 Dorado + `-ub`: [`../ont_r10/2_map.sh`](../ont_r10/2_map.sh)
