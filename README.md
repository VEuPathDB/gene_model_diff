# Geneset merge script for Apollo annotations
============================================

## Project description

### Features

The main perl script gene_model_diff.pl compares two genesets in GFF3 and classifies each gene_model change. 
It outputs a GFF3 file with the genes that should go on the merged dataset and a file with the list of former gene ids to delete from the original geneset.

### Usage

The software expects a directory for each species with 3 files cap.gff, cap.fasta, core.gff. Each species directory name should be listed on a file (speciesList parameter).
Bioperl and the gene_model_diff/modules must be included as perl libraries.

To run the software: perl run_gene_model_diff.pl --config gene_model_diff.conf --speciesFile list_of_species.txt
This is an early release and the software is still under development. 

### Known Issues

This Software only works with Bioperl 1.6 (it does not work with 1.7).