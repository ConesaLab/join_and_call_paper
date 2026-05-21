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

All jobs send Slurm mail on BEGIN, END, and FAIL to `fjetzinger@biobam.com`.

## Slurm resources (QoS)

| Step | cpus | mem | qos | time | Notes |
|------|------|-----|-----|------|-------|
| 0 pychopper | 4 | 100gb | long-mem | 7-00:00:00 | `-B 25000` batch size |
| 1 merge | 4 | 50gb | short | 24:00:00 | Matches mouse `ont/3c_merge_fastqs.sh` |
| 2 map | 8 | 50gb | medium | 2-00:00:00 | Matches `ont_r10/2_map.sh` |
| 3 fastqc | 8 | 10gb | short | 5:00:00 | Matches `ont_r10/1_fastq_qc.sh` |
| 4 bam_qc | 2 | 10gb | short | 5:00:00 | Matches `ont_r10/3_bam_qc.sh` |
| 5 sqanti | 1 | 100gb | medium | 2-00:00:00 | Matches `ont_r10/4_run_sqanti.sh` |
| 6 sqanti SIRV | 1 | 20gb | short | 10:00:00 | Matches `ont_r10/5_run_sqanti_SIRVS.sh` |
| 7 sqanti_reads | 2 | 50gb | medium | 2-00:00:00 | Matches `ont_r10/6_sqanti_reads.sh` |
| 8 sqanti merge | 12 | 100gb | medium | 7-00:00:00 | Matches `ont_r10/7_sqanti_reads_merge.sh` |
| 9 concat | 1 | 40gb | short | 24:00:00 | Matches `ont_r10/8_concat_samples.sh` |
| 10 count reads | 8 | 16gb | short | 24:00:00 | Matches `ont_r10/9_count_reads_joint.sh` |
| 11 read lengths | 8 | 16gb | medium | 2-00:00:00 | Matches `ont_r10/10_get_read_lengths_bam.sh` |

All jobs request **at most 7 days** wall time (cluster maintenance ~May 30 / Jun 1; longer reservations are rejected). Scripts otherwise match [`ont_r10`](../ont_r10/) mem/cpu/qos. If step 0 still OOMs, lower `PYCHOPPER_BATCH_SIZE` in `config.sh` or raise `--mem` in `0_pychopper.sh`.

## Validation checklist

1. **Pychopper QC:** `fastq_pychopper/*_pychopper_report.pdf`
2. **Read counts:** `wc -l` on `*_for_map.fastq` vs `fastq/*.fastq`
3. **Mapping:** `samtools flagstat bam_pychopper/*` vs `bam/*`
4. **Strand:** `samtools view -f 16` spot-check on pychopper BAMs
5. **Kit:** Confirm ENA library prep; re-run step 0 with correct `-k` in `config.sh` if needed

## References

- Mouse Pychopper + `-uf`: [`../ont/3_pychopper.sh`](../ont/3_pychopper.sh), [`../ont/4_map.sh`](../ont/4_map.sh)
- Default R10 Dorado + `-ub`: [`../ont_r10/2_map.sh`](../ont_r10/2_map.sh)
