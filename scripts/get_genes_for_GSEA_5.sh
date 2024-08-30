#!/bin/bash
#this script takes a specific SV type and lists affected genes for GSEA
# Function to print usage
usage() {
  echo "Usage: $0 -type [DEL|DUP|all] -i input_file -o output_file"
  exit 1
}

# Parse user input
if [ "$#" -ne 6 ]; then
  usage
fi

# Initialize variables
sv_type=""
input_file=""
output_file=""

# Parse the command line arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -type)
      sv_type="$2"
      shift 2
      ;;
    -i)
      input_file="$2"
      shift 2
      ;;
    -o)
      output_file="$2"
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

# Validate the sv_type
if [[ "$sv_type" != "DEL" && "$sv_type" != "DUP" && "$sv_type" != "all" ]]; then
  usage
fi

# Validate input and output files
if [ -z "$input_file" ] || [ -z "$output_file" ]; then
  usage
fi

# Extract unique gene IDs based on SV type and remove .x suffix
if [ "$sv_type" == "all" ]; then
  awk -F'\t' '{for (i=1; i<=NF; i++) if ($i ~ /gene_name /) {match($i, /gene_name "([^"]+)"/, arr); print arr[1]}}' "$input_file" | sed 's/\.[0-9]\+$//' | sort | uniq > "$output_file"
else
  grep "SVTYPE=${sv_type}" "$input_file" | awk -F'\t' '{for (i=1; i<=NF; i++) if ($i ~ /gene_name /) {match($i, /gene_name "([^"]+)"/, arr); print arr[1]}}' | sed 's/\.[0-9]\+$//' | sort | uniq > "$output_file"
fi

# Notify the user
echo "Unique gene names for SV type $sv_type have been written to $output_file"