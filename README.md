# RNAseq_analysis
<p>This repository contains custom programs and scripts used in molecular evolution analyses of RNA-seq datasets
The RNA-seq pipeline consists of many different programs and these scripts link programs together.
Prior steps created a <i>de novo</i> transcriptome assembly and identified putative open reading frames.
This pipeline starts with those open reading frame sequences.</p>
1. Translate to protein sequences with EMBOSS
2. Use OrthoMCL to identify putative orthologs
3. Use custom program om_groups_to_fasta to take the OrthoMCL output file groups.txt and convert it to fasta format.
4. Align grouped sequences to each other using MAFFT
5. Convert the aligned sequences in the .fasta output files from MAFFT to PAML format
6. Perform maximum-likelihood phylogenetic analysis with PAML

<p>Three custom programs are contained in this directory:</p>
1. om_groups_to_fasta
<p>This program takes the groups.txt output and creates fasta files for each group using the original sequence fasta file </p>
2. fasta_to_paml
<p>This program converts fasta files to input paml format. </p>
3. find_single_copy_genes
<p>This program searches through either fasta or paml files to identify files (intended to be ortholog groups) containing single-copy genes only </p>

