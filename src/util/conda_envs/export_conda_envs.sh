#!/usr/bin/env bash

# Name: export_conda_envs.sh
# Usage: ./export_conda_envs.sh <conda_occurrences_file>

set -e

CONDA_OCCURRENCES_FILE="$1"

# Check that we received a valid file
if [ -z "$CONDA_OCCURRENCES_FILE" ] || [ ! -f "$CONDA_OCCURRENCES_FILE" ]; then
  echo "Usage: $0 <path_to_conda_occurrences.tsv>"
  exit 1
fi

# For each line in the file, extract the environment name and export it
while IFS= read -r line; do
  # Get the environment name from the first tab-separated column
  env_name="$(echo "$line" | cut -f1)"

  # If for some reason the line is empty, skip
  [ -z "$env_name" ] && continue

  # Optionally, check if that environment actually exists
  if conda env list | awk '{print $1}' | grep -q "^${env_name}$"; then
    echo "Exporting environment: $env_name"
    conda env export --name "$env_name" > "${env_name}.yaml"
  else
    echo "Environment '$env_name' not found. Skipping..."
  fi

done < "$CONDA_OCCURRENCES_FILE"

echo "Done. YAML files for found environments have been created."
