# Xchr_inactivation_mouse

RNA-seq processing and differential expression analysis contributed to a study of imprinted X-chromosome inactivation in mouse.

Published as Maclary E, Hinten M, Harris C, **Sethuraman S**, Gayen S, Kalantry S. *PRC2 represses transcribed genes on the imprinted inactive X chromosome in mice.* Genome Biology (2017);18:82. [doi:10.1186/s13059-017-1211-5](https://doi.org/10.1186/s13059-017-1211-5)

## Background

In imprinted X-chromosome inactivation, one X chromosome is silenced according to parental origin rather than at random. Studying it requires distinguishing transcripts from the two X chromosomes within the same cell, which is done with an allele-specific approach: reads are assigned to the maternal or paternal chromosome using sequence variants that differ between the strains in a cross, rather than pooled as in a standard expression analysis.

The study examined how genes that remain transcribed on the inactive X are repressed, and identified a role for PRC2.

## Scope of this repository

This repository holds the RNA-seq processing and differential expression analysis contributed to the study. The remaining analyses in the paper were performed by co-authors. The published methods section is the authoritative description of the full study.

| File | Purpose |
|---|---|
| `mapping.sh` | Read alignment |
| `cufflinks.sh` | Transcript assembly and expression quantification |
| `HTcommand.sh` | Read counting |
| `DESeq_Lab.R` | Differential expression analysis in R |
| `totrc.sh`, `paste command.sh` | Intermediate file handling |
| `pipeline2.txt` | Pipeline notes and run order |

## Note

This is early work from a laboratory rotation, published as a record of the analysis contributed rather than as a maintained pipeline. Scripts retain absolute paths from the compute environment where they were originally run.
