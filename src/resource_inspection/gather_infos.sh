#!/bin/bash

# Check if the search path is provided as the first argument
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_search> [output_directory]"
    exit 1
fi

search_path="$1"

# If a second parameter is provided, use it as the base output directory.
if [ -n "$2" ]; then
    output_base="$2"
    echo "Using output directory: $output_base"
else
    output_base=""
fi

# Function to convert memory strings with units to kilobytes
convert_to_kb() {
    local mem_str="$1"
    local num unit
    num=$(echo "$mem_str" | grep -oE '^[0-9]+')
    unit=$(echo "$mem_str" | grep -oE '[A-Za-z]+$')

    case "${unit^^}" in
        K) echo "$num" ;;
        M) echo $((num * 1024)) ;;
        G) echo $((num * 1024 * 1024)) ;;
        T) echo $((num * 1024 * 1024 * 1024)) ;;
        *) echo "$num" ;;  # Default to kilobytes if unit is missing or unrecognized
    esac
}

# Find all joblog.tsv files and process them
find "$search_path" -type f -name 'joblog.tsv' | while read -r joblog_file; do

    echo "found $joblog_file; processing ..."

    joblog_dir=$(dirname "$joblog_file")

    # If an output directory is provided, replicate the subdirectory structure
    if [ -n "$output_base" ]; then
        # Compute the relative path of joblog_dir with respect to search_path
        rel_path=${joblog_dir#"$search_path"}
        rel_path=${rel_path#/}  # remove any leading slash
        new_dir="$output_base/$rel_path"
        mkdir -p "$new_dir"
        output_file="$new_dir/job_details.tsv"
    else
        output_file="$joblog_dir/job_details.tsv"
    fi

    echo "writing to $output_file ..."

    # Write header to the output file
    printf "job name\tjob id\tstart\tend\telapsed\treqcpus\treqmem\tmaxrss\ttotalcpu\tstatus\texitcode\tnodelist\tallocnodes\n" > "$output_file"

    # Process each line in the joblog.tsv file
    while IFS=$'\t' read -r jobname jobid; do
        # Fetch job details using sacct
        sacct_output=$(sacct -j "$jobid" --format=JobID,JobName,Start,End,Elapsed,ReqCPUs,ReqMem,MaxRSS,TotalCPU,State,ExitCode,NodeList,AllocNodes -P -n)

        # Process sacct output line by line
        while IFS='|' read -r sacct_jobid sacct_jobname sacct_start sacct_end sacct_elapsed \
                              sacct_reqcpus sacct_reqmem sacct_maxrss sacct_totalcpu sacct_state sacct_exitcode \
                              sacct_nodelist sacct_allocnodes; do
            # Check if JobName is "batch" to capture MaxRSS
            if [[ "$sacct_jobname" == "batch" ]]; then
                # Convert MaxRSS to kilobytes if available
                if [[ -n "$sacct_maxrss" ]]; then
                    maxrss_kb=$(convert_to_kb "$sacct_maxrss")
                else
                    maxrss_kb=""
                fi

                # Write the collected information to the output file
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                    "$jobname" "$sacct_jobid" "$sacct_start" "$sacct_end" "$sacct_elapsed" \
                    "$sacct_reqcpus" "$sacct_reqmem" "$maxrss_kb" "$sacct_totalcpu" "$sacct_state" "$sacct_exitcode" \
                    "$sacct_nodelist" "$sacct_allocnodes" >> "$output_file"
            fi
        done <<< "$sacct_output"

    done < "$joblog_file"

done
