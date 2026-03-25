#!/bin/bash
#SBATCH --job-name=create_subset
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/create_subset_%A.out

module load anaconda
conda activate flair

# Input
metadata_samples="$1"
target="$2"

reference_genome="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

# Targets
metadata_samples_target="${metadata_samples%.tsv}_subset.tsv"
metadata_concat_target="${metadata_samples%samples.tsv}concat_samples_subset.tsv"
target_bams="${target}/bams"
target_fastqs="${target}/fastqs"

mkdir -p $target_bams
mkdir -p $target_fastqs

> $metadata_samples_target
> $metadata_concat_target

declare -A unique_cond_fastqs

while IFS=$'\t' read -r cond sample_id pool bam fastq; do
    # Skip lines that start with '#'
    if [[ $cond == \#* ]]; then
        echo -e "$cond\t$sample_id\t$pool\t$bam\t$fastq" >> "$metadata_samples_target"
        echo -e "$cond\t$sample_id\t$pool\t$bam\t$fastq" >> "$metadata_concat_target"
        continue
    fi

    echo -e "\n---------------------------"
    echo "Condition: $cond"
    echo "Sample ID: $sample_id"
    echo "Pool: $pool"
    echo "BAM: $bam"
    echo "FASTQ: $fastq"
    echo -e "---------------------------\n"

    target_sam="${target_bams}/$(basename "${bam%.bam}")_subset.sam"
    target_bam="${target_bams}/$(basename "${bam%.bam}")_subset.bam"
    target_fastq="${target_fastqs}/$(basename "${fastq%.fastq}")_subset.fastq"

    target_concat_fastq="${target_fastqs}/${cond}_concat_subset.fastq"

    if [[ ! -v unique_cond_fastqs["$cond"] ]]; then
        unique_cond_fastqs["$cond"]="$target_concat_fastq"
        > $target_concat_fastq
    fi

    # take the first 1000 reads of $fastq
    head -n 4000 "$fastq" > "$target_fastq"
    head -n 4000 "$fastq" >> "$target_concat_fastq"

    # align
    minimap2 -ax splice -uf --MD -t 4 ${reference_genome} $target_fastq > $target_sam
    samtools sort -@ 4 $target_sam -o $target_bam --write-index

    # write metadata
    echo -e "${cond}\t${sample_id}\t${pool}\t${target_bam}\t${target_fastq}" >> $metadata_samples_target

done < "$metadata_samples"

for cond in "${!unique_cond_fastqs[@]}"; do
    target_concat_fastq="${unique_cond_fastqs[$cond]}"
    
    target_concat_sam="${target_bams}/${cond}_concat_subset.sam"
    target_concat_bam="${target_concat_sam%.sam}.bam"

    # align
    minimap2 -ax splice -uf --MD -t 4 ${reference_genome} $target_concat_fastq > $target_concat_sam 
    samtools sort -@ 4 $target_concat_sam -o $target_concat_bam --write-index

    case "$cond" in
        B0K100)
            sample_id="Kconcat"
            ;;
        B100K0)
            sample_id="Bconcat"
            ;;
        *)
            echo "Error: Unknown condition '$cond'" >&2
            continue  # Skip this record due to error
            ;;
    esac

    # write metadata
    echo -e "${cond}\t${sample_id}\t${pool}\t${target_concat_bam}\t${target_concat_fastq}" >> $metadata_concat_target
done
