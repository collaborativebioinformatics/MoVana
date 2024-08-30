version 1.0

# Task 1: Generate a Simulated Distribution of allele frequencies and modify VCF file
task GenerateSimulatedDistribution {
  input {
    File vcf_file
  }

  command <<<
    python3 <<CODE
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# Plotting with the specified colors and legend
bins = np.linspace(0, 1, 100)
plt.hist(af_distribution[(af_distribution <= 0.1) | ((af_distribution > 0.1) & (af_distribution <= 0.4))], 
         bins=bins, density=True, alpha=0.6, color='coral', label='Mosaic')
plt.hist(af_distribution[(af_distribution > 0.4)], 
         bins=bins, density=True, alpha=0.6, color='blue', label='Clonal')
plt.title('Simulated VAF Distribution of SVs')
plt.xlabel('VAF')
plt.ylabel('Density')
plt.legend()
plt.savefig('ICGC_VAF_distribution.png', bbox_inches='tight', dpi=300)
plt.show()

# Read the VCF file and add AF values
vcf_data = pd.read_csv('{vcf_file}', sep='\t', comment='#', header=None)

# Find the header line and set the column names
with open('{vcf_file}', 'r') as f:
    for line in f:
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            break

vcf_data.columns = headers

# Sample AF values from the generated distribution
af_values = np.random.choice(af_distribution, size=len(vcf_data))

# Add the AF values to the VCF
vcf_data['INFO'] = vcf_data['INFO'] + ';AF=' + af_values.astype(str)

# Write the updated VCF to a new file
output_vcf_path = 'icgc_with_af.vcf'
with open(output_vcf_path, 'w') as out_file:
    with open('{vcf_file}', 'r') as in_file:
        for line in in_file:
            if line.startswith('##'):
                out_file.write(line)
            elif line.startswith('#CHROM'):
                out_file.write(line)
                break
    vcf_data.to_csv(out_file, sep='\t', index=False, header=False)

print(f"VCF file with AF values saved to: {output_vcf_path}")

# Add SVLEN for clinical_SVs
def add_svlen_to_info(info, pos):
    info_dict = dict(item.split("=") for item in info.split(";") if "=" in item)
    end = int(info_dict['END'])
    svlen = end - int(pos)
    info_dict['SVLEN'] = str(svlen)
    new_info = ";".join([f"{key}={value}" for key, value in info_dict.items()])
    return new_info

input_vcf = 'icgc_with_af.vcf'
output_vcf = 'icgc_with_af_with_SVLEN.vcf'

with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
            continue
        columns = line.strip().split('\t')
        chrom = columns[0]
        pos = columns[1]
        info = columns[7]
        format_col = columns[8]
        sample = columns[9]
        updated_info = add_svlen_to_info(info, pos)
        columns[7] = updated_info
        if format_col == "GT":
            columns[9] = "0/1"
        outfile.write('\t'.join(columns) + '\n')

print(f"VCF file '{output_vcf}' has been created with SVLEN added to the INFO field.")
CODE
  >>>

  output {
    File updated_vcf = "icgc_with_af_with_SVLEN.vcf"
    File plot = "ICGC_VAF_distribution.png"
  }

  runtime {
    docker: "python:3.8"
    cpu : 1
    memory : "4 GiB"
    maxRetries : 1
  }
}

workflow MoVana_Workflow {
  input {
    File vcf_file
    String af
    File bed_file
    String sv_type
  }

  call GenerateSimulatedDistribution {
    input:
      vcf_file = vcf_file
  }

  output {
    File updated_vcf = GenerateSimulatedDistribution.updated_vcf
    File plot = GenerateSimulatedDistribution.plot
  }
}

