# Geneset merge script for Apollo annotations

## Project description

### Features

The main perl script gene_model_diff.pl compares two genesets in GFF3 from the same assembly and classifies each gene model change. The original gene models are in the file "core.gff" while the updated gene models are located on the "cap.gff" file. Both GFF3 files reference to the sequences on the file "cap.fasta".
It outputs a new GFF3 file with the genes that should go on the merged dataset and a file with the list of former gene ids to delete from the original geneset.

### Usage

The software expects a directory for each species with 3 files "cap.gff", "cap.fasta"", "core.gff"". Each species directory name should be listed on a file (speciesList parameter).
Bioperl and the gene_model_diff/modules must be included as perl libraries.

MySQL server connection details need to be provided in the config file. The modules directory have all the necessary code to store, compare, validate and retrieve gene models from a custom made database.  

To run the software:  

    perl run_gene_model_diff.pl --config gene_model_diff.conf --speciesFile list_of_species.txt

This is an early release and the software is still under development. 

### Known Issues

This software only works with Bioperl 1.6 (it does not work with 1.7).