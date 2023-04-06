# COMETS
COMETS (COmpare METagenomeS)

COMETS is an automated tool to compare multiple metagenomes. It can automatically process quality filters for raw reads, perform taxonomic classification, format taxonomic group tables (either Prokaryotes, Bacteria, Archaea, Eukaryota, Viruses or All) and output rarefaction curves, relative abundance, alpha diversity and NMDS Bray-Curtis dissimilarity plots. It can be run in a modular fashion, starting the whole process at any point or running one process at a time.

## Dependencies:
* fastp (https://github.com/OpenGene/fastp)
* kaiju-multi (https://github.com/bioinformatics-centre/kaiju)
* kaiju2table (https://github.com/bioinformatics-centre/kaiju)
* R (>= 4.1). Required libraries:
  + tibble, ggplot2 and dplyr (tidyverse: https://www.tidyverse.org/packages/)
  + phyloseq (https://bioconductor.org/packages/release/bioc/html/phyloseq.html)
  + mirlyn (https://github.com/escamero/mirlyn)

## Installation and setup:
comets is a bash script. Just download comets and comets_plots.R and place them in your local bin directory. If you don't have a bin directory just type within your HOME directory `mkdir bin`. If your bin directory is not in $PATH add the following to the end of your .profile file:

```{bash, eval=FALSE, echo=TRUE}
if [ -d "$HOME/bin" ] ; then
    PATH="$HOME/bin:$PATH"
fi
```

Then reload .profile:
```{bash, eval=FALSE, echo=TRUE}
source .profile
```

Make sure that the first line in comets_plots.R has the correct path to Rscript. To check your Rscript path type:
```{bash, eval=FALSE, echo=TRUE}
which Rscript
```

If necessary change the first line of comets_plots.R for the actual Rscript path

## Usage:
IMPORTANT NOTE: Before starting any process including -r option please provide a file named phyloseq_tables/SAM.table. This is a tab separated metadata table to be loaded by phyloseq. It must contain at least three columns: Files, Sample and SampleType  
Samples must be listed in alphabetical order (as in a simple ls and sort listing of files). Example: If you want to compare sample1 (paired end, seawater) and sample2 (single end freshwater) you will have three fastq.gz files (sample1_1.fastq.gz, sample1_2.fastq.gz and sample2.fastq.gz). The SAM.table within the phyloseq_tables/ directory should look like this:

Files|Sample|SampleType
--|--|--
sample1.kaiju.tsv|sample1|seawater
sample2.kaiju.tsv|sample2|freshwater

### Whole process starting from quality filter of raw reads:
```{bash, eval=FALSE, echo=TRUE}
comets -p pe|se|pe_se -f auto|<path_to_fastp> -k auto|<path_to_kaiju-multi> [-z <int>] [-d <path_to_kaiju_db>] [-K <path_to_kaiju2table>] [-n <path_to_nodes.dmp>] [-N <path_to_names.dmp>] -t -g All|Prokaryotes|Bacteria|Archaea|Eukaryota|Viruses -r
```
Example paired end files: sample1_1.fastq.gz, sample1_2.fastq.gz, sample2_1.fastq.gz, sample2_2.fastq.gz  
Example single end files: sample1.fastq.gz, sample2.fastq.gz  
NOTE: Single end files name must not have _1 or _2 suffix before .fastq.gz extension. Otherwise the script will wet very confused  

### Whole process starting from taxonomic classification of quality filtered reads:
```{bash, eval=FALSE, echo=TRUE}
comets -p pe|se|pe_se -k auto|<path_to_kaiju-multi> [-z <int>] [-d <path_to_kaiju_db>] [-K <path_to_kaiju2table>] [-n <path_to_nodes.dmp>] [-N <path_to_names.dmp>] -t -g All|Prokaryotes|Bacteria|Archaea|Eukaryota|Viruses -r
```
Example paired end files: paired_end/fastp/sample1_1.fastp.fastq.gz, paired_end/fastp/sample1_2.fastp.fastq.gz, paired_end/fastp/sample2_1.fastp.fastq.gz, paired_end/fastp/sample2_2.fastp.fastq.gz  
Example single end files: single_end/fastp/sample1.fastp.fastq.gz, single_end/fastp/sample2.fastp.fastq.gz  
NOTE: Single end files name must not have _1 or _2 suffix before .fastp.fastq.gz extension. Otherwise the script will wet very confused 
NOTE: Input files must be within a paired_end/fastp/ or single_end/fastp/ directory, respectively

### Whole process starting from creation of phyloseq OTU and TAX tables from taxonomic classification tables:
```{bash, eval=FALSE, echo=TRUE}
comets -p pe|se|pe_se -t -g All|Prokaryotes|Bacteria|Archaea|Eukaryota|Viruses -r
```
Example input files: kaiju2table/sample1.kaiju.tsv, kaiju2table/sample2.kaiju.tsv  
NOTE: Input files must be within a kaiju2table/ directory  

### Whole process starting from creation of R plots from phyloseq objects:
```{bash, eval=FALSE, echo=TRUE}
comets -p pe|se|pe_se -g All|Prokaryotes|Bacteria|Archaea|Eukaryota|Viruses -r
```
NOTE: comets_plots.R will save final dataframes as .RData files to phyloseq_tables directory. These can be loaded in RStudio with load("filename.RData") for custom processing

## Options:
* -h: help message
* -p: paired-end (pe), single-end (se), or both (pe_se). Always needed
* -f: path to fastp (if "auto", will automatically look for fastp in ~/, which may take a while. If option is not provided, will not run fastp quality filter)
* -k: path to kaiju-multi (if "auto", will automatically look for kaiju-multi, kaiju2table, nodes.dmp and names.dmp in ~/, which may take a while. If path provided, will also need to provide -K, -n, and -N. If option is not provided, will not run kaiju taxonomic assignment)
* -z: number of threads to use for kaiju taxonomic assignment (recomendation is 25 if server has enough cores. If option is not provided will use default: 1)
* -d: path to kaiju database (if option is not provided will automatically look for kaiju_db_nr_euk.fmi in ~/, which may take a while)
* -K: path to kaiju2table (needed if -k path is provided)
* -n: path to nodes.dmp (needed if -k path is provided)
* -N: path to names.dmp (needed if -k path is provided)
* -t: create phyloseq tables (no arguments needed. If provided, will also need to provide -g. If option is not provided, will not create phyloseq tables)
* -g: taxonomic group to create phyloseq tables for (Options are All, Prokaryotes, Bacteria, Archaea, Eukaryota, Viruses. Needed if -t or -r is provided)
* -r: create R plots (no arguments needed. If provided, will also need to provide -g. If option is not provided, will not create R plots)

Example output tree structure of a full process:
```{bash, eval=FALSE, echo=TRUE}
.
├── paired_end/
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   ├── sample2_1.fastq.gz
│   ├── sample2_2.fastq.gz
│   └── fastp/
│       ├── sample1_1.fastp.fastq.gz
│       ├── sample1_2.fastp.fastq.gz
│       ├── sample2_1.fastp.fastq.gz
│       ├── sample2_2.fastp.fastq.gz
│       └── kaiju/
│           ├── sample1.kaiju.out
│           └── sample2.kaiju.out
├── single_end/
│    ├── sample1.fastq.gz
│    ├── sample2.fastq.gz
│    └── fastp/
│        ├── sample1.fastp.fastq.gz
│        ├── sample2.fastp.fastq.gz
│        └── kaiju/
│            ├── sample1.kaiju.out
│            └── sample2.kaiju.out
├── kaiju2table/
│  ├── sample1.kaiju.tsv
│  ├── sample2.kaiju.tsv
│  ├── unclassified.txt
│  └── Viruses/
│      ├── sample1_Viruses.tsv
│      └── sample2_Viruses.tsv
└── phyloseq_tables/
    ├── Viruses_OTU.table
    ├── Viruses_TAX.table
    └── SAM.table
```
