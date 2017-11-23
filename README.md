gene_model_diff compares two genesets in GFF3 and classifies each gene_model change
The software expects a directory for each species with 3 files cap.gff, cap.fasta, core.gff.
	
To run the software: perl run_gene_model_diff.pl --config gene_model_diff.conf --speciesFile list_of_species.txt
This is an early release and the software is still under development. 
Known Bugs:
This Software only works with Bioperl 1.6 (it does not work with 1.7)