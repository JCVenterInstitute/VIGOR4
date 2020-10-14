# VIGOR - Viral Genome ORF Reader
VIGOR4 (Viral Genome ORF Reader) is a Java application to predict protein sequences encoded in viral genomes.<br>
VIGOR4 determines the protein coding sequences by sequence similarity searching against curated viral protein databases.<br>

This project is funded by the National Institute Of Allergy And Infectious Diseases of the National Institutes of Health under Award Number U19AI110819 and Contract Number HHSN272201400028C and is a collaboration between Northrop Grumman Health IT, J. Craig Venter Institute, and Vecna Technologies.

### Annotatable Viruses

Vigor4 uses the [VIGOR_DB](https://github.com/JCVenterInstitute/VIGOR_DB) project which currently has databases for the following viruses:

* Influenza (A & B for human, avian, and swine, and C for human)
* SARS-CoV-2
* West Nile Virus
* Zika Virus
* Chikungunya Virus
* Eastern Equine Encephalitis Virus
* Respiratory Syncytial Virus
* Rotavirus
* Enterovirus
* Lassa Mammarenavirus

## Installing VIGOR4

#### Quick Start

```
mvn -DskipTests clean package
unzip target/vigor-4.0.VERSION.zip -d INSTALL_DIRECTORY
# set reference database path, exonerate path and temporary directory
vi INSTALL_DIRECTORY/vigor-4.0.VERSION/config/vigor.ini
INSTALL_DIRECTORY/vigor-4.0.VERSION/bin/vigor4 -i VIRUS.fasta -o OUTPUT_DIRECTORY -d VIRUS_DB
```

#### Build Dependencies
Vigor4 uses 'Maven' to build & package the program <br>
#### Runtime Dependencies
VIGOR4 runs on Linux environment <br>
JAVA8 or above<br>
[Exonerate2.2.0](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)<br>
Vigor Viral database (Available at [VIGOR_DB](https://github.com/JCVenterInstitute/VIGOR_DB))

Refer to [INSTALL.md](https://github.com/JCVenterInstitute/VIGOR4/blob/master/INSTALL.md) for detailed instructions

## Running VIGOR4

*usage*: vigor4 -i inputfasta -o outputprefix  -d refdb

 Named Arguments:
```
  -h, --help             show this help message and exit
  -i <input fasta>, --input-fasta <input fasta>
                         path to fasta  file  of  genomic  sequences  to be
                         annotated
  -o <output prefix>, --output-prefix <output prefix>
                         prefix for outputfile  files,  e.g.  if  the output
                         prefix is  /mydir/anno  VIGOR  will  create output
                         files /mydir/anno.tbl, /mydir/anno.stats, etc.
  -c MIN_COVERAGE, --min-coverage MIN_COVERAGE
                         minimum  coverage  of  reference  product  (0-100)
                         required to report a gene,  by default coverage is
                         ignored
  -P <parameter=value~~...~~parameter=value>, --parameter <parameter=value~~...~~parameter=value>
                         ~~ separated list of  VIGOR parameters to override
                         default values.  Use  --list-config-parameters  to
                         see settable parameters.
  -v, --verbose          verbose logging (default=terse)
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
  --temporary-directory TEMPORARYDIRECTORY
                         Root directory to use for temporary directories

reference database:
  -d <ref db>, --reference-database <ref db>
                         specify the reference database to  be used, (-D is
                         a synonym for this option)
  
locus tag usage:
  -l, --no-locus-tags    do  NOT  use   locus_tags   in   TBL  file  output
                         (incompatible with -L)
  -L [<locus_tag_prefix>], --locus-tags [<locus_tag_prefix>]
                         USE locus_tags in  TBL  file  output (incompatible
                         with -l). If  no  prefix  is  provided, the prefix
                         "vigor_" will be used.
```
#### Outputs:
```
 outputprefix.rpt   -  summary of program results
 outputprefix.cds   -  fasta file of predicted CDSs
 outputprefix.pep   -  fasta file of predicted proteins
 outputprefix.tbl   -  predicted features in GenBank tbl format
 outputprefix.aln   -  alignment  of  predicted  protein  to  reference, and
                       reference protein to genome
 outputprefix.gff3  -  predicted features in GFF3 format
```
#### Currently unimplemented VIGOR3 Command Line Options:
```
  -0, --circular         complete circular  genome  (allows  gene  to  span
                         origin). This feature is currently unimplemented
  -f {0,1,2}, --frameshift-sensitivity {0,1,2}
                         frameshift  sensitivity,   0=ignore   frameshifts,
                         1=normal (default), 2=sensitive. 
  -m, --ignore-reference-requirements
                         ignore      reference      match      requirements
                         (coverage/identity/similarity),  sometimes  useful
                         when running VIGOR  to  evaluate  raw  contigs and
                         rough draft sequences
  -x <ref_id,...,ref_id>, --ignore-refID <ref_id,...,ref_id>
                         comma separated list of  reference sequence IDs to
                         ignore  (useful   when   debugging   a   reference
                         database). Not currently implemented
```


## Reference Databases:

 | Name  | Description |
 | :------------ | :----------|
 | flua |  Influenza A |
 | flub |  Influenza B|
 |fluc  |   Influenza C |
 | lassa  |     Lassa Mammarenavirus   |                                
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

