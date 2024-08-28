# cancer_mosaic_SVs

## Our idea:

Linking mosaic SVs in cancer to affected genes and pathways can contribute to the projection of tumor evolution 

## Data use:

i) SRP186687: CCLE_cancer_cell_line
ii) SRP300100: stage_II/III_colorectal_cancer
iii) SRP100797: mosaicisms_in_clinically_unremarkable_individuals
iv) SRP411687: BRCA1/2_Tumor
v) SRP426587: synchronous_bilateral_Wilms_tumor
vi) SRP115159: simulated_tumor_genomes
vii) SRP297028: Colorectal_Cancer
viii) SRP198194: hepatocellular_carcinoma_(HCC)

## Pipeline structure: 

Mosaic SV calls in VCF (based on VAF) --> identify affected genes --> annotate the functional effect, e.g. deletion - downregulation, amplification - upregulation --> outline selective advantages of mosaic SVs knowing their implications in cancer 

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
