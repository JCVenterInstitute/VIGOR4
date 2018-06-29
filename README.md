# VIGOR - Viral Genome ORF Reader
VIGOR4 (Viral Genome ORF Reader) is a Java application to predict protein sequences encoded in viral genomes.
VIGOR4 determines the protein coding sequences by sequence similarity searching against curated viral protein databases.
This project is funded by The National Institute of Allergy and Infectious Diseases (NIH / DHHS) under Contract No.HHSN272201400028C and is a collaboration 
between Northrop Grumman Health IT, J. Craig Venter Institute, and Vecna Technologies.Research reported in this publication was supported by the National Institute Of Allergy And Infectious Diseases of the 
National Institutes of Health under Award Number U19AI110819. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

### VIGOR4 currently annotates the following viruses:
* Influenza (A & B for human, avian, and swine)
* West Nile Virus
* Zika Virus
* Chikungunya Virus
* Eastern Equine Encephalitis Virus
* Respiratory Syncytial Virus
* Rotavirus
* Enterovirus

## Installing VIGOR4
#### Build Dependencies
Vigor4 uses 'Maven' to build & package the program <br>
#### Runtime Dependencies
VIGOR4 runs on Unix environment <br>
JAVA8 or above<br>
Exonerate2.2.0<br>
Vigor Viral databse (Available at [VIGOR_DB](https://github.com/JCVenterInstitute/VIGOR_DB))

Refer to [INSTALL.md](https://github.com/JCVenterInstitute/VIGOR4/blob/master/INSTALL.md) for detailed instructions

#### Running VIGOR4

*usage*: vigor4 -i inputfasta -o outputprefix  -d refdb

 Named Arguments:
```
   -h, --help             show this help message and exit
   -i <input fasta>, --input-fasta <input fasta>
                            path to fasta  file  of  genomic  sequences  to be
                            annotated, (-I is a synonym for this option)
   -I <input fasta>       synonym for -i/--input-fasta)
   -o <output prefix>, --output-prefix <output prefix>
                            prefix for outputfile  files,  e.g.  if  the ouput
                            prefix is  /mydir/anno  VIGOR  will  create output
                            files /mydir/anno.tbl, /mydir/anno.stats, etc.,
                            (-O is a synonym for this option)
   -O <output prefix>     synonym for -o/--output-prefix
   -c MIN_COVERAGE, --min-coverage MIN_COVERAGE
                            minimum  coverage  of  reference  product  (0-100)
                            required to report a gene,  by default coverage is
                            ignored
   -P <parameter=value~~...~~parameter=value>, --parameter <parameter=value~~...~~parameter=value>
                            ~~ separated list of  VIGOR parameters to override
                            default values.  Use  --list-config-parameters  to
                            see settable parameters.
   -s <gene size>, --min-gene-size <gene size>
                            minimum size (aa) of product  required to report a
                            gene, by default size is ignored
   -v, --verbose           verbose logging (default=terse)
   --list-config-parameters
                            list available configuration parameters and exit
   --config-file CONFIG_FILE
                            config file to use
   --reference-database-path REFERENCE_DATABASE_PATH
                            reference database path
   --virus-config VIRUSSPECIFICCONFIG
                            Path to virus specific configuration
   --virus-config-path VIRUSSPECIFICCONFIGPATH
                            Path  to  directory   containing   virus  specific
                            config files.
   --overwrite-output     overwrite existing output files if they exist
   --temporary-directory  TEMPORARYDIRECTORY
                            Root directory to use for temporary directories
```
#### Outputs:
```
 outputprefix.rpt   -  summary of program results
 outputprefix.cds   -  fasta file of predicted CDSs
 outputprefix.pep  -  fasta file of predicted proteins
 outputprefix.tbl    -  predicted features in GenBank tbl format
 outputprefix.aln   -  alignment  of  predicted  protein  to  reference, and
                                    reference protein to genome
```
#### Reference Databases:

 | Name  | Description |
 | :-----: | :----------:|
 | flua |  Flu A |
 | flub |  Flu B|
 |fluc  |   Flu C |
 |rsv |        Respiratory syntactical virus (RSV)   |
 | rtva   |     Rotavirus A   |                                
  |rtvb    |    Rotavirus B     |                              
  |rtvc    |    Rotavirus C     |                              
 | rtvf   |     Rotavirus F  |
 | rtvg   |     Rotavirus G  |
 | sapo    |    Sapovirus  |
 | veev     |   Alphaviruses (VEEV/EEEV)|
 | wnvI     |   West Nile Virus - Lineage I |
 | wnvII    |   West Nile Virus - Lineage II |
 | zikv    |    Zika virus |

