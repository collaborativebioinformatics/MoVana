import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Generate a Simulated Distribution of allele frequencies 
# Define parameters for the mixture of Gaussians
np.random.seed(1)  
n_samples = 10000  # Number of samples to generate

# Mixture components (means and standard deviations)
means = [0.05, 0.25, 0.5]
std_devs = [0.02, 0.05, 0.04]
weights = [0.3, 0.2, 0.5]  # Weights for each Gaussian

# Generate the samples
samples = []
for mean, std_dev, weight in zip(means, std_devs, weights):
    samples.append(np.random.normal(loc=mean, scale=std_dev, size=int(n_samples * weight)))

# Combine all samples
af_distribution = np.concatenate(samples)

# Clip the values to be within [0, 1]
af_distribution = np.clip(af_distribution, 0, 1)

# Step 2: Plotting with the specified colors and legend
# Define the bins for the histogram
bins = np.linspace(0, 1, 100)

# Plot the first two peaks in light orange ('mosaic')
plt.hist(af_distribution[(af_distribution <= 0.1) | ((af_distribution > 0.1) & (af_distribution <= 0.4))], 
         bins=bins, density=True, alpha=0.6, color='coral', label='Mosaic')

# Plot the third peak in light blue ('clonal')
plt.hist(af_distribution[(af_distribution > 0.4)], 
         bins=bins, density=True, alpha=0.6, color='blue', label='Clonal')

# Adding titles and labels
plt.title('Simulated VAF Distribution of SVs')
plt.xlabel('VAF')
plt.ylabel('Density')

# Adding the legend
plt.legend()

plt.savefig(f'/g/korbel/olisov/hackathon/ICGC_VAF_distribution.png', bbox_inches='tight', dpi=300)

# Show the plot
plt.show()

#Read the VCF file and add AF values
vcf_file_path = '/g/korbel/olisov/hackathon/icgc_filtered_sorted.vcf'
vcf_data = pd.read_csv(vcf_file_path, sep='\t', comment='#', header=None)

# Find the header line and set the column names
with open(vcf_file_path, 'r') as f:
    for line in f:
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            break

vcf_data.columns = headers

# Sample AF values from the generated distribution
af_values = np.random.choice(af_distribution, size=len(vcf_data))

# Step 3: Add the AF values to the VCF
vcf_data['INFO'] = vcf_data['INFO'] + ';AF=' + af_values.astype(str)

# Define the output file path
output_vcf_path = '/g/korbel/olisov/hackathon/icgc_with_af.vcf'

# Write the updated VCF to a new file, preserving the header
with open(output_vcf_path, 'w') as out_file:
    with open(vcf_file_path, 'r') as in_file:
        for line in in_file:
            if line.startswith('##'):
                out_file.write(line)
            elif line.startswith('#CHROM'):
                out_file.write(line)
                break

    vcf_data.to_csv(out_file, sep='\t', index=False, header=False)

print(f"VCF file with AF values saved to: {output_vcf_path}")

#add SVLEN for clinical_SVs
# File paths
input_vcf = '/g/korbel/olisov/hackathon/icgc_with_af.vcf'
output_vcf = '/g/korbel/olisov/hackathon/icgc_with_af_with_SVLEN.vcf'

def add_svlen_to_info(info, pos):
    # Parse the existing INFO field
    info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
    
    # Calculate SVLEN
    end = int(info_dict['END'])
    svlen = end - int(pos)
    
    # Add SVLEN to the INFO field
    info_dict['SVLEN'] = str(svlen)
    
    # Reconstruct the INFO field
    new_info = ";".join([f"{key}={value}" for key, value in info_dict.items()])
    
    return new_info

# Read the VCF file
with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
    for line in infile:
        # Write the header lines as-is
        if line.startswith('#'):
            outfile.write(line)
            continue
        
        # Split the line into columns
        columns = line.strip().split('\t')
        
        # Extract the necessary fields
        chrom = columns[0]
        pos = columns[1]
        info = columns[7]
        format_col = columns[8]
        sample = columns[9]
        
        # Update the INFO field with SVLEN
        updated_info = add_svlen_to_info(info, pos)
        
        # Replace the old INFO field with the updated one
        columns[7] = updated_info
        
        # Change GT in the FORMAT and SAMPLE columns
        if format_col == "GT":
            columns[9] = "0/1"
        
        # Write the updated line to the output VCF
        outfile.write('\t'.join(columns) + '\n')

print(f"VCF file '{output_vcf}' has been created with SVLEN added to the INFO field.")

##! next i use a bash script to filter the file based on AF to have only mosaic events

