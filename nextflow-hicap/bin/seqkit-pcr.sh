#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <primer_info> <sample_IDs_file>"
    echo "++++++++++"
    echo "Make sure to: "
    echo "install seqkit (conda install -c bioconda seqkit)" 
    echo "before running the script."
    echo "++++++++++"
    exit 1
fi

# Get the start time
start_time=$(date +%s)

# Read the primer_info file and store the data in an array
declare -A primers
echo "Reading primer_info file..."
while IFS= read -r line || [[ -n "$line" ]]; do
    if $first_line; then
        first_line=false
        continue
    fi
    serotype=$(echo "$line" | awk '{print $1}')
    forward=$(echo "$line" | awk '{print $2}')
    reverse=$(echo "$line" | awk '{print $3}')
    primers[$serotype]="$forward $reverse"

    echo "Serotype: $serotype"
    echo "Forward: $forward"
    echo "Reverse: $reverse"
    echo "Primers: ${primers[$serotype]}"
done < "$1"

echo "Number of serotypes: ${#primers[@]}"

# Create a new file to store the detection results
detection_results="detection_results.csv"

# Write a header line to the detection results file
echo "SampleFile, Serotype, Detection_Status" > "$detection_results"

# Ensure the output directory exists
mkdir -p detailed_results

# Read the sample IDs from the input file
while IFS= read -r sample_id || [[ -n "$sample_id" ]]; do
    # Get the path to the FASTA file
    file=$(mdu contigs --sample_id "$sample_id" | awk '{print $2}')

    # Print a message indicating which file is being tested
    echo "Testing file: $file"

    # Loop through all the serotypes and check for their presence in the sample
    for serotype in "${!primers[@]}"; do

        # Create a new output file for the current sample and serotype
        output_file="detailed_results/${sample_id}_${serotype}_pcr.txt"

        # Extract the forward and reverse primers for the current serotype
        forward=$(echo "${primers[$serotype]}" | awk '{print $1}')
        reverse=$(echo "${primers[$serotype]}" | awk '{print $2}')

        # Print a message indicating which serotype is being tested with which primers
        echo "Testing serotype: $serotype"
        echo "Forward primer: $forward"
        echo "Reverse primer: $reverse"

        # Store the command in a variable
        command="seqkit amplicon -m 0 -F \"$forward\" -R \"$reverse\" < \"$file\" > \"$output_file\""

        # Print the command
        echo "$command"

        # Execute the command
        eval "$command"

        # Check if the output file contains the string '>NODE'
        if grep -q '>NODE' "$output_file"; then
        # The serotype is present in the sample
        echo "$sample_id, $serotype, +ve" >> "$detection_results"
        
        else
        # The serotype is not present in the sample
        echo "$sample_id, $serotype, NA" >> "$detection_results"
        fi

    done
done < "$2"

# Get the end time
end_time=$(date +%s)

# Calculate the time elapsed
time_elapsed=$((end_time - start_time))

echo "Time elapsed: $time_elapsed seconds"

# Print a success message
echo "Run ended"