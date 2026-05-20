# Shared paths for the ONT R10 Pychopper branch (source from Slurm scripts).
# Does not modify the default ont_r10 Dorado + minimap2 -ub pipeline.

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
fastq_in="${base_dir}/fastq"
fastq_out="${base_dir}/fastq_pychopper"
bam_dir="${base_dir}/bam_pychopper"
gff_dir="${base_dir}/gff_pychopper"
ref_dir="${base_dir}/ref"
analysis_dir="${base_dir}/analysis"
logs_dir="${analysis_dir}/logs_pychopper"

assembly="${ref_dir}/GRCh38_SIRV.fa"
ref_annotation="${ref_dir}/gencode.v49_SIRV.gtf"

PYCHOPPER_KIT="LSK114"

# Analysis subdirs (isolated from default ont_r10 branch)
fastqc_dir="${analysis_dir}/fastqc_pychopper"
flagstat_dir="${analysis_dir}/flagstat_pychopper"
sqanti_dir="${analysis_dir}/run_SQANTI_pychopper"
sqanti_reads_dir="${analysis_dir}/sqanti_reads_pychopper"
read_qc_dir="${analysis_dir}/read_qc_pychopper"
