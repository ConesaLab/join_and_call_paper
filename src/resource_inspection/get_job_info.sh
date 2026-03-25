#!/bin/bash
# This script accepts a single job ID, a comma-separated list,
# or a file piped via STDIN (one job ID per line), and outputs
# job details in TSV format to stdout.

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

# If no parameters and no STDIN, print usage.
if [ -t 0 ] && [ "$#" -eq 0 ]; then
    echo "Usage: $0 <jobid(s)>"
    echo "Provide a single job ID, a comma-separated list, or pipe a file with one job ID per line."
    exit 1
fi

# Read job IDs from command-line or STDIN.
if [ "$#" -gt 0 ]; then
    # Combine all command-line arguments (can be comma- or space-separated)
    jobids="$*"
else
    jobids=$(cat)
fi

# Replace commas with spaces so that we can iterate over each job ID
jobids=$(echo "$jobids" | tr ',' ' ')

# Print header for TSV output to stdout.
printf "job name\tjob id\tstart\tend\telapsed\treqcpus\treqmem\tmaxrss\ttotalcpu\tstatus\texitcode\tnodelist\tallocnodes\n"

# Process each job ID
for jobid in $jobids; do
    # Trim any extra whitespace
    jobid=$(echo "$jobid" | xargs)
    [ -z "$jobid" ] && continue

    # Get job details from sacct in pipe-separated format (-P) and no header (-n)
    sacct_output=$(sacct -j "$jobid" --format=JobID,JobName,Start,End,Elapsed,ReqCPUs,ReqMem,MaxRSS,TotalCPU,State,ExitCode,NodeList,AllocNodes -P -n)

    # Process each line of sacct output
    while IFS='|' read -r sacct_jobid sacct_jobname sacct_start sacct_end sacct_elapsed \
                              sacct_reqcpus sacct_reqmem sacct_maxrss sacct_totalcpu sacct_state \
                              sacct_exitcode sacct_nodelist sacct_allocnodes; do
        # The "batch" line usually holds the aggregated resource usage info.
        if [[ "$sacct_jobname" == "batch" ]]; then
            if [[ -n "$sacct_maxrss" ]]; then
                maxrss_kb=$(convert_to_kb "$sacct_maxrss")
            else
                maxrss_kb=""
            fi
            # Output a TSV line with the collected details.
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$sacct_jobname" "$sacct_jobid" "$sacct_start" "$sacct_end" "$sacct_elapsed" \
                "$sacct_reqcpus" "$sacct_reqmem" "$maxrss_kb" "$sacct_totalcpu" "$sacct_state" \
                "$sacct_exitcode" "$sacct_nodelist" "$sacct_allocnodes"
        fi
    done <<< "$sacct_output"
done
