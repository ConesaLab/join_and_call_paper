#!/bin/bash
#SBATCH --job-name=isoquant_ind
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=80gb # Per processor memory
#SBATCH -t 4-00:00:00     # Walltime
#SBATCH -o log/isoquant_ind_%A_%a.out

module load anaconda
source activate isoquant

# Input
n_cores="$SLURM_CPUS_PER_TASK"

fastq_files=()

# chr8 only
# fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/lrgasp_h1_h1de_chr8_1_60000000/endodermal_cdna_pacbio_chr8/ENCFF235QXW/ENCFF235QXW_chr8.fastq.gz")
# fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/lrgasp_h1_h1de_chr8_1_60000000/endodermal_cdna_pacbio_chr8/ENCFF712CBL/ENCFF712CBL_chr8.fastq.gz")
# fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/lrgasp_h1_h1de_chr8_1_60000000/endodermal_cdna_pacbio_chr8/ENCFF731HST/ENCFF731HST_chr8.fastq.gz")
# fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/lrgasp_h1_h1de_chr8_1_60000000/h1_cdna_pacbio_chr8/ENCFF400BQQ/ENCFF400BQQ_chr8.fastq.gz")
# fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/lrgasp_h1_h1de_chr8_1_60000000/h1_cdna_pacbio_chr8/ENCFF436GKZ/ENCFF436GKZ_chr8.fastq.gz")
# fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/lrgasp_h1_h1de_chr8_1_60000000/h1_cdna_pacbio_chr8/ENCFF853OFP/ENCFF853OFP_chr8.fastq.gz")

# genome="/storage/gge/Fabian/lrgasp/homo_sapiens/reference/chr8/GRCh38.primary_assembly.genome_chr8.fa"
# annotation="/storage/gge/Fabian/lrgasp/homo_sapiens/reference/chr8/gencode.v39.annotation_chr8.gtf"
# prefix="lrgasp_h1_h1de_chr8_1_60000000"
# out_iso_detect="./isoquant_lrgasp_h1_h1de_chr8_1_60000000/"

# full data
fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/endodermal/pacbio_cdna/ENCFF235QXW.fastq.gz")
fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/endodermal/pacbio_cdna/ENCFF712CBL.fastq.gz")
fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/endodermal/pacbio_cdna/ENCFF731HST.fastq.gz")
fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/h1/pacbio_cdna/ENCFF400BQQ.fastq.gz")
fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/h1/pacbio_cdna/ENCFF436GKZ.fastq.gz")
fastq_files+=("/storage/gge/Fabian/lrgasp/homo_sapiens/h1/pacbio_cdna/ENCFF853OFP.fastq.gz")

genome="/storage/gge/Fabian/lrgasp/homo_sapiens/reference/GRCh38.primary_assembly.genome.fa"
annotation="/storage/gge/Fabian/lrgasp/homo_sapiens/reference/gencode.v39.annotation.gtf"
prefix="lrgasp_h1_h1de"
out_iso_detect="./isoquant_lrgasp_h1_h1de/"

# Print Input
echo "[isoquant] FASTQ reads: ${fastq_files[*]}"
echo "[isoquant] Genome: $genome"
echo "[isoquant] Annotation: $annotation"
echo "[isoquant] N cores: $n_cores"
echo "[isoquant] Out dir: $out_iso_detect"
echo "[isoquant] Prefix: $prefix"

if [ ! -d "$out_iso_detect" ]; then
	mkdir -p "$out_iso_detect"
	echo "[isoquant] Directory created: $out_iso_detect"
else
	echo "[isoquant] Directory already exists: $out_iso_detect"
fi

# Convert array to space-separated string
fastq_files_str="${fastq_files[*]}"

isoquant.py -t $n_cores --reference $genome --genedb $annotation --fastq $fastq_files_str --data_type nanopore --fl_data -p $prefix -o $out_iso_detect
