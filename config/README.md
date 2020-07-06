# Gene model diff config file
The config file is in the format of Config::IniFiles see CPAN for more information.

## Database
The database connection information (*host*, *port*, *user*, *pass*) is used to create a MySQL database per species. 
The script uses *prefix* and *database* to create a database with the name "prefix_species_database" which is used to process the gene models
and store the final annotation events. The load_database needs to be created manually and is used to load GFF.

## Data
The Data section contains paths to the *datadir* where the input and output files are stored. The *scriptdir* is where run_gene_model_diff.pl is located.
*log_file* is the path to the Log::Log4perl config file.

## Validate
The error values (*GFF_mRNA*, *GFF_exon*, *CDS_start*, *CDS_internal_stop*, *runtime_error*) are used in Validate.pm. 
The error codes are added together, so a transcript without a start and stop codon will have a error value of CDS_start + CDS_stop. Do not change them without consulting the code.
*approved_email* can be used to override the validation by a superuser.  

# Log file configiration
The Log::Log4perl config file gene_model_logFile.conf can be used to change the behaviour of the logging see CPAN more information.
*log4perl.appender.LOGFILE.filename=* must be set or run_gene_model_diff.pl will fail.
