#!/bin/bash
#SBATCH --job-name=tusco_config
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=2gb # Per processor memory
#SBATCH -t 1:00:00     # Walltime
#SBATCH -o log/tusco_config_%A_%a.out

module load anaconda
source activate tusco_env

# input
base_dirs_file=$1

# Check if input file exists
if [ ! -f "$base_dirs_file" ]; then
    echo "Error: Base directories file not found: $base_dirs_file"
    exit 1
fi

# Create config file header (only once)
echo "classification_file,data_type,tool_name,tissue_name,sample_name" > tusco_config.csv

# Process each base directory
while IFS= read -r base_data_dir || [ -n "$base_data_dir" ]; do
    # Skip empty lines
    [ -z "$base_data_dir" ] && continue
    
    echo "Processing directory: $base_data_dir"
    
    # Determine data_type and tool_name from base_data_dir
    # Check paper_stuff first (more specific case)
    if [[ "$base_data_dir" == *"/paper_stuff/"* ]]; then
        data_type="isoseq"
    elif [[ "$base_data_dir" == *"/ont/"* ]]; then
        data_type="ont"
    elif [[ "$base_data_dir" == *"/isoseq/"* ]]; then
        data_type="isoseq"
    else
        echo "Warning: Skipping directory (must contain '/ont/', '/isoseq/', or '/paper_stuff/'): $base_data_dir"
        continue
    fi

    # Extract tool_name from path
    # Check paper_stuff first (more specific case)
    if [[ "$base_data_dir" =~ /paper_stuff/[^/]+/([^/]+)/isoform_detection ]]; then
        tool_name="${BASH_REMATCH[1]}"
    elif [[ "$base_data_dir" =~ /(ont|isoseq)/([^/]+)/ ]]; then
        tool_name="${BASH_REMATCH[2]}"
    else
        echo "Warning: Skipping directory (could not extract tool name): $base_data_dir"
        continue
    fi

    # Replace isoform_detection and everything after it with tama_condition
    if [[ "$base_data_dir" =~ (.*)/[^/]+/isoform_detection ]]; then
        tama_dir="${BASH_REMATCH[1]}/tama_condition"
    else
        echo "Warning: Directory must contain 'isoform_detection': $base_data_dir"
        continue
    fi

    # Process isoform detection files
    echo "Processing isoform detection files..."
    
    # Set patterns based on special case detection
    if [[ "$base_data_dir" == *"/paper_stuff/"* && "$tool_name" == "isoseq" ]]; then
        echo "Special case: paper_stuff/isoseq with collapse directory structure"
        # Special case patterns for paper_stuff/isoseq
        brain_find_pattern="*/collapse/*_sq3_orth/*_classification.txt"
        kidney_find_pattern="*/collapse/*_sq3_orth/*_classification.txt"
        brain_regex="/B100K0/([^/]+)/collapse/[^/]+_sq3_orth/[^/]+_classification\.txt$"
        kidney_regex="/B0K100/([^/]+)/collapse/[^/]+_sq3_orth/[^/]+_classification\.txt$"
    else
        # Standard case patterns
        brain_find_pattern="*/[^/]*_sq3/*_classification.txt -o -path */[^/]*_sq3_orth/*_classification.txt"
        kidney_find_pattern="*/[^/]*_sq3/*_classification.txt -o -path */[^/]*_sq3_orth/*_classification.txt"
        brain_regex="/B100K0/([^/]+)/[^/]+_sq3(_orth)?/[^/]+_classification\.txt$"
        kidney_regex="/B0K100/([^/]+)/[^/]+_sq3(_orth)?/[^/]+_classification\.txt$"
    fi
    
    # Find classification files using the determined patterns
    class_files_brain_raw=($(find "$base_data_dir/B100K0" -path $brain_find_pattern))
    class_files_kidney_raw=($(find "$base_data_dir/B0K100" -path $kidney_find_pattern))
    
    # Filter out excluded files for paper_stuff/mandalorion case
    class_files_brain=()
    class_files_kidney=()
    
    # Filter brain files
    for file in "${class_files_brain_raw[@]}"; do
        if [[ "$file" =~ $brain_regex ]]; then
            sample_name="${BASH_REMATCH[1]}"
            if [[ "$base_data_dir" == *"/paper_stuff/"* && "$tool_name" == "mandalorion" && "$sample_name" == "B100K0" ]]; then
                echo "Excluding file where sample_name equals tissue_name: $file"
            else
                class_files_brain+=("$file")
            fi
        fi
    done
    
    # Filter kidney files
    for file in "${class_files_kidney_raw[@]}"; do
        if [[ "$file" =~ $kidney_regex ]]; then
            sample_name="${BASH_REMATCH[1]}"
            if [[ "$base_data_dir" == *"/paper_stuff/"* && "$tool_name" == "mandalorion" && "$sample_name" == "B0K100" ]]; then
                echo "Excluding file where sample_name equals tissue_name: $file"
            else
                class_files_kidney+=("$file")
            fi
        fi
    done

    # Process brain files
    for classification_file in "${class_files_brain[@]}"; do
        # Extract sample_name from path using the appropriate regex
        if [[ "$classification_file" =~ $brain_regex ]]; then
            sample_name="${BASH_REMATCH[1]}"
            tissue_name="B100K0"
            
            # Add to config file
            echo "$classification_file,$data_type,$tool_name,$tissue_name,$sample_name" >> tusco_config.csv
        fi
    done

    # Process kidney files
    for classification_file in "${class_files_kidney[@]}"; do
        # Extract sample_name from path using the appropriate regex
        if [[ "$classification_file" =~ $kidney_regex ]]; then
            sample_name="${BASH_REMATCH[1]}"
            tissue_name="B0K100"
            
            # Add to config file
            echo "$classification_file,$data_type,$tool_name,$tissue_name,$sample_name" >> tusco_config.csv
        fi
    done

    # Process TAMA files
    echo "Processing TAMA files..."
    # Find TAMA classification files for brain and kidney
    tama_files_brain=($(find "$tama_dir/B100K0_TAMA" -path "*/B100K0_TAMA_sq3/*_classification.txt"))
    tama_files_kidney=($(find "$tama_dir/B0K100_TAMA" -path "*/B0K100_TAMA_sq3/*_classification.txt"))

    # Process brain files (TAMA)
    for classification_file in "${tama_files_brain[@]}"; do
        if [[ "$classification_file" =~ /B100K0_TAMA/B100K0_TAMA_sq3/[^/]+_classification\.txt$ ]]; then
            tissue_name="B100K0"
            sample_name="TAMA"
            
            # Add to config file
            echo "$classification_file,$data_type,$tool_name,$tissue_name,$sample_name" >> tusco_config.csv
        fi
    done

    # Process kidney files (TAMA)
    for classification_file in "${tama_files_kidney[@]}"; do
        if [[ "$classification_file" =~ /B0K100_TAMA/B0K100_TAMA_sq3/[^/]+_classification\.txt$ ]]; then
            tissue_name="B0K100"
            sample_name="TAMA"
            
            # Add to config file
            echo "$classification_file,$data_type,$tool_name,$tissue_name,$sample_name" >> tusco_config.csv
        fi
    done

    # Sanity checking - count files per tissue
    isoform_brain_count=${#class_files_brain[@]}
    isoform_kidney_count=${#class_files_kidney[@]}
    tama_brain_count=${#tama_files_brain[@]}
    tama_kidney_count=${#tama_files_kidney[@]}
    
    total_isoform=$((isoform_brain_count + isoform_kidney_count))
    total_tama=$((tama_brain_count + tama_kidney_count))
    
    echo "Summary for $data_type/$tool_name:"
    echo "  B100K0 - isoform_detection: $isoform_brain_count files, tama_condition: $tama_brain_count files"
    echo "  B0K100 - isoform_detection: $isoform_kidney_count files, tama_condition: $tama_kidney_count files"
    
    # Check expected counts and warn if different
    if [ $isoform_brain_count -ne 6 ]; then
        echo "WARNING: Expected 6 isoform_detection files for B100K0, found $isoform_brain_count"
    fi
    if [ $isoform_kidney_count -ne 6 ]; then
        echo "WARNING: Expected 6 isoform_detection files for B0K100, found $isoform_kidney_count"
    fi
    if [ $tama_brain_count -ne 1 ]; then
        echo "WARNING: Expected 1 tama_condition file for B100K0, found $tama_brain_count"
    fi
    if [ $tama_kidney_count -ne 1 ]; then
        echo "WARNING: Expected 1 tama_condition file for B0K100, found $tama_kidney_count"
    fi
    
    echo "Added $total_isoform isoform detection files and $total_tama TAMA files from $base_data_dir"
done < "$base_dirs_file"

total_entries=$(($(wc -l < tusco_config.csv) - 1))
echo "Created tusco_config.csv with $total_entries total entries"