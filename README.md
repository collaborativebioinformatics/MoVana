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

*Farha - identify mosaic SVs and affected genes? 

*Franceso, Gobi - integrate VEP/other gene function DB into the pipeline for annotation? 

*Gil - write the paper 
