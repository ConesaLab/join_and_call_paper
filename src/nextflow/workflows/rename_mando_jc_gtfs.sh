#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 metadata_concat"
    exit 1
fi

# Assign the first argument to a variable for clarity
metadata_file="$1"
output_file="${metadata_file%.*}_renamed.csv"

# Ensure the input file exists
if [ ! -f "$metadata_file" ]; then
    echo "Error: File '$metadata_file' not found!"
    exit 1
fi

# Get the header and locate the "gtf" column index
header=$(head -n 1 "$metadata_file")
IFS=$'\t' read -r -a columns <<< "$header"
gtf_index=-1

for i in "${!columns[@]}"; do
    if [ "${columns[i]}" == "gtf" ]; then
        gtf_index=$i
        break
    fi
done

if [ "$gtf_index" -eq -1 ]; then
    echo "Error: 'gtf' column not found in the header!"
    exit 1
fi

# Process each line of the file starting from the second line
{
    echo "$header" # Write the header to the output file
    tail -n +2 "$metadata_file" | while IFS=$'\t' read -r -a fields; do
        original_filename="${fields[$gtf_index]}"
        
        # Replace "Bconcat" with "B100K0" and "Kconcat" with "B0K100"
        new_filename="${original_filename/Bconcat/B100K0}"
        new_filename="${new_filename/Kconcat/B0K100}"

        # If the filename was changed, rename the actual file
        if [[ "$original_filename" != "$new_filename" ]]; then
            if [ -f "$original_filename" ]; then
                mv "$original_filename" "$new_filename"
            else
                echo "Warning: File '$original_filename' not found for renaming!"
            fi
        fi

        # Update the gtf field in the current line
        fields[$gtf_index]="$new_filename"
        
        # Print the updated line
        (IFS=$'\t'; echo "${fields[*]}")
    done
} > "$output_file"

echo "Processed metadata saved to '$output_file'."
