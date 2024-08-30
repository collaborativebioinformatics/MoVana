#!/bin/bash
#thi script filters the vcf file to obtain events with below specific VAF (mosaic)
# Default values
AF="0.4"
INPUT_FILE="icgc_with_af_with_SVLEN.vcf"
OUTPUT_FILE="filtered_icgc_with_SVLEN.vcf"

# Function to display usage information
usage() {
  echo "Usage: $0 [-a allele_frequency] [-i input_file] [-o output_file]"
  exit 1
}

# Parse command-line options
while getopts "a:i:o:" opt; do
  case $opt in
    a)
      AF=$OPTARG
      ;;
    i)
      INPUT_FILE=$OPTARG
      ;;
    o)
      OUTPUT_FILE=$OPTARG
      ;;
    *)
      usage
      ;;
  esac
done

# Check if input and output files are provided
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
  echo "Error: Input and output files must be specified."
  usage
fi

# Run bcftools with the user-specified or default options
bcftools view -i "INFO/AF<$AF" "$INPUT_FILE" -o "$OUTPUT_FILE"