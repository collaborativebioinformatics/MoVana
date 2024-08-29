# MoVana 
MOsaic structural Variants ANonation in cAncer is a pipeline for identifying mosaic DNA rearrangements in the call set and linking them to functional impact.

## Background


## Workflow 

Mosaic SV calls in VCF (based on VAF) --> identify affected genes --> annotate the functional effect, e.g. deletion - downregulation, amplification - upregulation --> outline selective advantages of mosaic SVs knowing their implications in cancer 

- ''SnpEff'', ''ANNOVAR'', or ''VEP''  to annotate structural variants
- ''KEGG'', ''Reactome'', or ''Gene Ontology (GO)'' to find associations between genes and known biological pathways?
- Integrating epigenomic data?
- Compare the integrated epigenomic data with normal data (human genome) to see the variations 
- Visualise our results in heatmaps/ plots 

## Roles (future Hackaton team):

*Dmitrii - get the dataset

*Farha - identify mosaic SVs and affected genes

*Gobi - integrate gene function DB into the pipeline for annotation

*Gil - write the paper 

## Filtration

- we will filter the VCF data file, based on AF (Allele frequency) < 0.4, VCF contains three types of structural variations (inversions, duplications and deletions)
- create shell script for filtering, in this step we used "bcftools" tool for filtration that takes three arguments (input.vcf, output.vcf, AF threshold)
- save shell script in the same directory of the input.vcf data to run command line for filtration './filter_vcf.sh input.vcf filtered_output.vcf' (this includes: shell script, input data file, resulted file)

