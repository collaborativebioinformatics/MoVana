#!/bin/bash
# Input parameters
input_vcf=$1
output_vcf=$2

# Define the allele frequency threshold
af_threshold=0.4

# Run bcftools command to filter VCF
bcftools view -i 'INFO/SVTYPE="INV" || INFO/SVTYPE="DUP" || INFO/SVTYPE="DEL" && INFO/AF<${af_threshold}' "${input_vcf}" -o "${output_vcf}"
