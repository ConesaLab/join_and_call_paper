#!/usr/bin/env bash

# Name: find_conda_envs.sh
# Usage: ./find_conda_envs.sh <directory_path>

# Exit on error
set -e

SEARCH_PATH="$1"

# Make sure a directory was provided
if [ -z "$SEARCH_PATH" ] || [ ! -d "$SEARCH_PATH" ]; then
  echo "Usage: $0 <path_to_directory>"
  exit 1
fi

# Create an associative array to store environment -> list_of_files
declare -A ENV_TO_FILES

# Find all .sh and .sbatch files, then examine their contents
while IFS= read -r file
do
  # Read each file line by line to find lines with "source activate X" or "conda activate X"
  while IFS= read -r line
  do
    # Match lines that have "source activate X" or "conda activate X" 
    # where X is a sequence of non-whitespace characters
    if [[ "$line" =~ (source|conda)[[:space:]]activate[[:space:]]([^[:space:]]+) ]]; then
      env_name="${BASH_REMATCH[2]}"
      # Accumulate file paths in the associative array
      if [[ -n "${ENV_TO_FILES[$env_name]}" ]]; then
        # Append the new file with a comma separator
        ENV_TO_FILES[$env_name]="${ENV_TO_FILES[$env_name]},$file"
      else
        # Initialize with this file path
        ENV_TO_FILES[$env_name]="$file"
      fi
    fi
  done < "$file"

done < <(find "$SEARCH_PATH" -type f \( -name "*.sh" -o -name "*.sbatch" \))

# Write results to conda_occurrences.tsv
{
  for env_name in "${!ENV_TO_FILES[@]}"; do
    # Print the environment, a tab, then the comma-separated list of files
    echo -e "${env_name}\t${ENV_TO_FILES[$env_name]}"
  done
} > conda_occurrences.tsv

echo "Occurrences saved to conda_occurrences.tsv"
