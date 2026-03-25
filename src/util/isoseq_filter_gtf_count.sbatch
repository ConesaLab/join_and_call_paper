#!/bin/bash
#SBATCH --job-name=check_count_gtf
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=NONE      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=2gb     # Per processor memory
#SBATCH -t 1-00:00:00         # Walltime
#SBATCH -o log/check_count_gtf_%A.out

#############################
# Function Definitions
#############################
# Trim leading and trailing whitespace
trim() {
  local var="$*"
  # Remove leading whitespace
  var="${var#"${var%%[![:space:]]*}"}"
  # Remove trailing whitespace
  var="${var%"${var##*[![:space:]]}"}"
  echo "$var"
}

#############################
# Input Validation
#############################
count_file=$1
gtf_file=$2
out_gtf=$3

if [ -z "$count_file" ] || [ -z "$gtf_file" ]; then
  echo "Usage: $0 <count_file> <gtf_file>"
  exit 1
fi

if [ ! -f "$count_file" ]; then
  echo "Count file not found: $count_file"
  exit 1
fi

if [ ! -f "$gtf_file" ]; then
  echo "GTF file not found: $gtf_file"
  exit 1
fi

echo "Processing count file: $count_file"
echo "Using GTF file: $gtf_file"

#############################
# Process Count File
#############################
# Use an associative array for transcript IDs seen in the count file
declare -A count_ids_map
# And an array to keep track of transcript IDs that have all zero counts
zero_count_ids=()

# Use process substitution so the while loop runs in the current shell
while IFS=$'\t' read -r transcript_id counts; do
  # Skip header lines starting with '#' (extra safety)
  if [[ "$transcript_id" == \#* ]]; then
    continue
  fi

  # Record this transcript_id in the count file map
  count_ids_map["$transcript_id"]=1

  # Split counts into an array (Bash splits on whitespace by default)
  read -ra count_values <<< "$counts"

  all_zeros=true
  for count in "${count_values[@]}"; do
    cleaned_count=$(trim "$count")
    # If the cleaned count is not "0" and does not start with "0." then mark as nonzero.
    if [[ "$cleaned_count" != "0" && "$cleaned_count" != 0.* ]]; then
      all_zeros=false
      break
    fi
  done

  if $all_zeros; then
    zero_count_ids+=("$transcript_id")
  fi
done < <(tail -n +2 "$count_file")

#############################
# Process GTF File
#############################
# Extract transcript IDs from the GTF file using awk.
# We use mapfile to load them into an array.
mapfile -t gtf_transcript_ids < <(awk -F'\t' '$3 == "transcript" {
  if (match($9, /transcript_id "([^"]+)"/, arr)) {
    print arr[1]
  }
}' "$gtf_file" | sort -u)

# Build an associative array from the GTF transcript IDs
declare -A gtf_ids_map
for id in "${gtf_transcript_ids[@]}"; do
  gtf_ids_map["$id"]=1
done

echo "Number of transcript IDs in count file: ${#count_ids_map[@]}"
echo "Number of transcript IDs in GTF file: ${#gtf_ids_map[@]}"

#############################
# Determine Missing IDs
#############################
# IDs present in GTF but not in count file
missing_in_count=()
for id in "${!gtf_ids_map[@]}"; do
  if [ -z "${count_ids_map[$id]}" ]; then
    missing_in_count+=("$id")
  fi
done

# IDs present in count file but not in GTF
missing_in_gtf=()
for id in "${!count_ids_map[@]}"; do
  if [ -z "${gtf_ids_map[$id]}" ]; then
    missing_in_gtf+=("$id")
  fi
done

#############################
# Report Results
#############################
echo ""
echo "Transcript IDs with all zero counts:"
for id in "${zero_count_ids[@]}"; do
  echo "$id"
done

echo ""
echo "Transcript IDs present in GTF but not in count file:"
for id in "${missing_in_count[@]}"; do
  echo "$id"
done

echo ""
echo "Transcript IDs present in count file but not in GTF:"
for id in "${missing_in_gtf[@]}"; do
  echo "$id"
done

#############################
# Filter the GTF File
#############################
# Combine IDs to be filtered: those with all zero counts and those missing in count.
filter_ids=("${zero_count_ids[@]}" "${missing_in_count[@]}")


if [ ${#filter_ids[@]} -gt 0 ]; then
  # Write the filter IDs to a temporary file, one per line.
  tmp_filter=$(mktemp)
  if [ ! -f "$tmp_filter" ]; then
    echo "Failed to create a temporary file for filter IDs."
    exit 1
  fi
  for id in "${filter_ids[@]}"; do
    echo "$id" >> "$tmp_filter"
  done
  
  echo ""
  echo "Filtering GTF file to remove lines with transcript IDs in the filter list."
  
  # AWK script:
  # 1. In the BEGIN block, read the temporary filter file and build an associative array "filter".
  # 2. For each line, if a transcript_id is found in the filter array, skip the line.
  # 3. Otherwise, print the line.
  awk -v filter_file="$tmp_filter" 'BEGIN {
    while ((getline line < filter_file) > 0) {
      filter[line] = 1
    }
    close(filter_file)
  }
  {
    if (match($0, /transcript_id "([^"]+)"/, arr)) {
      if (arr[1] in filter)
        next
    }
    print
  }' "$gtf_file" > "$out_gtf"
  
  rm "$tmp_filter"
  echo "Filtered GTF saved to $out_gtf"
else
  echo ""
  echo "No transcript IDs to filter from GTF. Copying original GTF to filtered file."
  cp "$gtf_file" "$out_gtf"
  echo "Copied GTF saved as $out_gtf"
fi
