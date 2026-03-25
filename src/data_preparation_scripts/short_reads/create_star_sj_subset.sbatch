#!/bin/bash
#SBATCH --job-name=create_star_sj_subset
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/create_star_sj_subset_%A.out

# create subsets of STAR SJ output files with the first 1000 lines, but guaranteeing to include at least the first 10 unannotated lines

# Directory source and target variables
source=$1
target=$2

# Pattern for files
pattern="$source/[BK]3[1-5]/*SJ.out.tab"

# Iterate over each file matching the pattern
for file in $pattern; do
    # Check if it's an actual file
    if [ -f "$file" ]; then
        echo "Processing file: $file"

        # Read the first 1000 lines
        head_lines=$(head -n 1000 "$file")

        # Filter lines where the 6th column is 0 and count them
        unannotated_count=$(echo "$head_lines" | awk '$6 == 0 {count++} END {print count+0}')

        # Determine the number of lines needed if unannotated_count is less than 10
        if [ "$unannotated_count" -lt 10 ]; then
            missing_count=$((10 - unannotated_count))
            
            # Extract lines beyond the first 1000 where the 6th column is 0
            replacement_lines=$(tail -n +1001 "$file" | awk '$6 == 0' | head -n "$missing_count")
            
            # Count the actual number of replacement lines found
            actual_replacement_count=$(echo "$replacement_lines" | wc -l)
            
            # Remove the last (actual_replacement_count) lines from head_lines
            trimmed_head_lines=$(echo "$head_lines" | head -n -$actual_replacement_count)
    
            # Combine the trimmed lines with the replacement lines
            final_lines=$(echo -e "${trimmed_head_lines}\n${replacement_lines}")
        else
            final_lines="$head_lines"
        fi

        # Construct the new filename
        base_filename=$(basename "$file" SJ.out.tab)
        new_filename="${target}/${base_filename}_subset_SJ.out.tab"

        # Write to the new file
        echo "$final_lines" > "$new_filename"

        echo "Written to: $new_filename"
    else
        echo "No valid files found matching the pattern."
    fi
done
