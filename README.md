# MoVana 
**MO**saic structural **V**ariants **AN**nonation in c**A**ncer is a pipeline for identifying mosaic DNA rearrangements in the call set and linking them to functional impact.

## Background
Cancer is a highly heterogeneous system that arises from healthy cells through a series of point mutations and large DNA rearrangements. Sporadic mutagenesis gives rise to tumor subclones that have distinct sets of genomic alterations, which can promote tumor growth, metastasis, or treatment resistance. In comparison to single nucleotide variants, the characterization of more complex events that contribute to intratumoral genetic heterogeneity was lacking until recent efforts in deep whole genome sequencing of tumors and the development of mutation callers. In this project, we specifically focus on mosaic structural variants in cancer and have designed a tool for their functional annotation by identifying which genes and biological pathways are affected by these mosaic structural variants. While simple, this tool should serve as a stepping stone for further studies on the contribution of rare variants to the emergence of treatment-resistant subclones and the recurrence of disease.

## Workflow 
![17B195E0-10F4-45B8-9F95-8EB99F8BBC56_1_201_a](https://github.com/user-attachments/assets/0a23cc97-28a1-470e-b572-0211579c35d4)

Overall, the pipeline is designed to select mosaic events based on their allele frequency (AF), annotate them with overlapping genes, and perform gene set enrichment analysis to infer functional impact. MoVana takes as input a VCF file of SV calls, where each entry must have a reported AF (allele frequency) value. Before running the tool, we advise plotting the distribution of the variant allele frequencies (as in the example below). Putative mosaic events should form one or multiple separate peaks at frequencies below 0.5, and the border between clusters can be chosen as the VAF threshold. In this example, we chose 0.4. The tool uses bcftools to filter out the entries above the user-specified threshold. 

![ICGC_VAF_distribution](https://github.com/user-attachments/assets/7154eb1d-d17c-4665-a278-4f3fdca548c3)

Next, known SV breakpoints and gene coordinates are used to find overlaps with bedtools and annotate each event with the respective affected genes. The disruptive effect of various rearrangements on gene function depends on the type of SV; therefore, an extra filtering step is required before the enrichment analysis. For this, a filtering step based on the SV type is included, and its output can be submitted for gene set enrichment analysis or a search among known genes implicated in cancer. The example output of the GSEA analysis shows genes overlapping duplications in mock data.

![GSEA](https://github.com/user-attachments/assets/cae8f075-046b-4e60-9ab4-f1874b3a3469)

## Installation

## Usage

## Future directions
*Comparing the composition and frequency of mosaic SVs in primary tumors vs metastasis and relapse

*Identifying recurrent mosaic events implicated in treatment resistance 

## Test data
Test data is located in the files and was downloaded from the ICGC bucket available at AWS 
s3://icgc25k-open/PCAWG/consensus_sv/final_consensus_sv_bedpe_passonly.icgc.public.tgz. It includes deletions, duplications and inversions. The calls were assigned simulated AF values and filtered with the threshold of 0.4 VAF.
