# cancer_mosaic_SVs
## Our idea:

Linking mosaic SVs in cancer to affected genes and pathways can contribute to the projection of tumor evolution 

## Data use

SRP186687: CCLE_cancer_cell_line
SRP300100: stage_II/III_colorectal_cancer
SRP100797: mosaicisms_in_clinically_unremarkable_individuals
SRP411687: BRCA1/2_Tumor
SRP426587: synchronous_bilateral_Wilms_tumor
SRP115159: simulated_tumor_genomes
SRP297028: Colorectal_Cancer
SRP198194: hepatocellular_carcinoma_(HCC)

## Pipeline structure: 

Mosaic SV calls in VCF (based on VAF) --> identify affected genes --> annotate the functional effect, e.g. deletion - downregulation, amplification - upregulation --> outline selective advantage of mosaic SVs knowing their implication in cancer 

- ''SnpEff'', ''ANNOVAR'', or ''VEP''  to annotate structural variants
- ''KEGG'', ''Reactome'', or ''Gene Ontology (GO)'' to find associations between genes and known biological pathways?
- Integrating epigenomic data?
- Compare the integrated epigenomic data with normal data (human genome) to see the variations 
- Visualise our results in heatmaps/ plots 

## Roles:

*Dmitrii - get the dataset and identify mosaic SVs

*Farha - identify genes? 

*Franceso, Gobi - integrate VEP/other gene function DB into the pipeline for annotation? 

*Gil - write the paper 
