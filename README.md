# MoVana 
MOsaic structural Variants ANonation in cAncer is a pipeline for identifying mosaic DNA rearrangements in the call set and linking them to functional impact.

## Background


## Workflow 
![17B195E0-10F4-45B8-9F95-8EB99F8BBC56_1_201_a](https://github.com/user-attachments/assets/0a23cc97-28a1-470e-b572-0211579c35d4)

## Roles

*Dmitrii - get the dataset

*Farha - identify mosaic SVs and affected genes

*Gobi - integrate gene function DB into the pipeline for annotation

*Gil - write the paper 

## Filtration

- we will filter the VCF data file, based on AF (Allele frequency) < 0.4, VCF contains three types of structural variations (inversions, duplications and deletions)
- create shell script for filtering, in this step we used "bcftools" tool for filtration that takes three arguments (input.vcf, output.vcf, AF threshold)
- save shell script in the same directory of the input.vcf data to run command line for filtration './filter_vcf.sh input.vcf filtered_output.vcf' (this includes: shell script, input data file, resulted file)
