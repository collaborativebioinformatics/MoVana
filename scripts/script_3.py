import random
import pandas as pd

# File paths
input_vcf = '/g/korbel/olisov/hackathon/filtered_icgc_with_SVLEN.vcf'
output_vcf_sampled = '/g/korbel/olisov/hackathon/filtered_icgc_with_SVLEN2_sampled.vcf'

# Initialize list to store data lines
data_lines = []

# Read the VCF file and separate header and data lines
with open(input_vcf, 'r') as infile:
    with open(output_vcf_sampled, 'w') as outfile:
        for line in infile:
            # Collect header lines
            if line.startswith('#'):
                outfile.write(line)
            else:
                data_lines.append(line.strip())

# Perform random sampling of 1000 lines
sampled_lines = random.sample(data_lines, min(1000, len(data_lines)))

# Write the sampled data lines
with open(output_vcf_sampled, 'a') as outfile:
    for line in sampled_lines:
        outfile.write(line + '\n')

print(f"VCF file '{output_vcf_sampled}' containing 1000 randomly sampled rows has been created.")