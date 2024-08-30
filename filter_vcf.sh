#!/bin/sh

## Allele frequency : 0.4
input=SVclone_SVs_2.vcf
output=filtered_SVclone_SVs_2.vcf

bcftools view -i 'INFO/SVTYPE="INV" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" && INFO/AF<0.4' "$input" -o "$output"
