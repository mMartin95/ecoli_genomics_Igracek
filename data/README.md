# Data Directory

This directory stores the essential reference genome, raw variant calls, structural variant datasets, and functional annotation outputs required for the *E. coli* K-12 genomic analysis pipeline.

## Contents

- `ATCC25922.fasta`  
  **Reference genome**: The standard *Escherichia coli* ATCC 25922 reference sequence in FASTA format, serving as the coordinate backbone for read mapping (`bwa`) or (`bowtie2`), pileup generation, and visual validation in IGV.

- `ATCC25922.gff3`  
**Structural Annotation**: Genomic features, gene coordinates, and locus tags (`DR76_RS...`) utilized by annotation tools (`snpEff`).

  - `ATCC25922.gb`  
**GenBank Record**: Complete GenBank flatfile formatted specifically as the reference input for `breseq` structural variant and junction analysis.

- `all_samples_joint_ raw_variants.vcf` 
  **Raw Variants**: - Raw variants for all samples C40, ZnO40, ZnONPs40, ZnO20+20, ZnONPs20+20.

  - `all_samples_joint_ filtered_annotated.vcf` 
  **Filtered and annotated variants**: Filtered and annotated variants for all samples C40, ZnO40, ZnONPs40, ZnO20+20, ZnONPs20+20.

- `ZnO40_variants.vcf`
  **Variants**: Fixed variants for ZnO40. 

- `ZnONPs40_variants.vcf`
  **Variants**: Fixed variants for ZnONPs40.

  - `ZnO20+20_variants.vcf`
  **Variants**: Fixed variants for ZnO20+20.

  - `ZnONPs20+20_variants.vcf`
  **Variants**: Fixed variants for ZnONPs20+20. 
