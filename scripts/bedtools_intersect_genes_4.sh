#!/bin/bash
#this script intersects the SV positions with genes

# Default values
INPUT_FILE="filtered_icgc_with_SVLEN2_sampled.vcf"
OUTPUT_FILE="icgc_with_genes.txt"
BED_FILE="genes_cds_only.bed"

# Function to display usage information
usage() {
  echo "Usage: $0 [-i input_file] [-o output_file]"
  exit 1
}

# Parse command-line options
while getopts "i:o:" opt; do
  case $opt in
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

# Run bedtools intersect with the user-specified or default options
bedtools intersect -a "$INPUT_FILE" -b "$BED_FILE" -wo > "$OUTPUT_FILE"
