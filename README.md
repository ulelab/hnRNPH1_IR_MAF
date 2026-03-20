## hnRNPH1 Intron 4 Decoy Analysis

#### This repository contains data, scripts, and outputs for the hnRNPH1 intron 4 retention decoy analysis workflow.

## Decoy Dataset Generation

#### The database of predicted decoy loci was generated as a BED file using the bedtools intersect workflow described in `decoy_splice_site_flowchart.pdf`

#### Briefly, all human intron coordinates were collected from Vast-DB PSI-TABLE-hg38.tab.gz as any EVENT with an ID beginning with 'HsaIN'. These coordinates were fed in a strandwise fashion to bedtools getfasta, then SpliceAI to predict splice donor scores for every intronic nucleotide, using ***nt flanks for internal normalization by the canonical 5' splice site. Scores below 0.01 were filtered out, then all coordinates were intersected against PRPF8-binding datasets of RBPnet predictions generated with the same scaled workflow of all introns, and eCLIP and iCLIP signal from experimental data. Predictions within 100nt of a canonical splice site were removed from the dataset. Finally this dataset `DecoySpliceSites_f50.bed` was filtered to keep only 'protein-coding' introns by overlapping with `Gencode.v49.annotation.gtf` using GenomicRanges in the `decoy_exon_overlaps.Rmd` script to output 

### bedtools intersect commands

#### `$ bedtools intersect -a Splice_All.filtered.01min.bed -b rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8_iCLIP_HepG2_xlinks.bed PRPF8_eCLIP_HepG2_rollmean10_minHeightAdjust1. 0_minPromAdjust1.0_minGeneCount5_Peaks.bed -c -s | awk '$NF > 0' > SupportedSpliceSites_f50.bed`

#### `$ bedtools intersect -a SupportedSpliceSites_f50.bed -b Wide_Canonical_splice_sites.bed -v -s > DecoySpliceSites_f50.bed`



## Calculate peak SpliceAI score around decoy site

#### The final `decoy.bed` dataset input to `SpliceAI_Inference.py` which slops the coordinates 24 nt wider on either side of the proposed decoy site, then strand-aware extracted fasta sequence from bedtools getfasta is passed to SpliceAI. Max Donor score for the 49 nt segment is kept for the locus and an updated BED file is output. 

## Figure 1 R Markdown (`scripts/hnRNPH1_figure1.rmd`)

#### This R Markdown document generates Figure 1: a comparative multiple-sequence alignment view of the hnRNPH1 intron 4 decoy region (`chr5:179620560-179620600`, hg38). It starts from an extracted UCSC multiz MAF alignment (converted to FASTA), plots the raw alignment, then creates a manuscript-ready alignment by renaming taxa, converting DNA bases from `T` to `U`, removing selected outlier species, and dropping columns that are gaps/missing across all taxa.

The final plot highlights the proposed decoy site (positions 27-33 in the processed alignment; labeled as genomic interval `179,620,576-179,620,582`) and is intended for direct use in manuscript figure generation.

### Inputs used by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.fa`
- `data/tree_to_clade_mapping.tsv`

### Output written by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.processed.fa`

