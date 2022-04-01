# COMETS
COMETS (COmpare METagenomeS)

COMETS is an automated tool to compare multiple metagenomes. It can automatically process quality filters for raw reads, perform taxonomic classification, format taxonomic group tables (either Bacteria or Archaea or Eukaryota or Viruses) and output rarefaction curves, relative abundance, alpha diversity and NMDS Bray-Curtis dissimilarity plots.

## Dependencies:
* fastp (https://github.com/OpenGene/fastp)
* kaiju-multi (https://github.com/bioinformatics-centre/kaiju)
* kaiju2table (https://github.com/bioinformatics-centre/kaiju)
* R (>= 4.1). Required libraries:
  + tibble, ggplot2 and dplyr (tidyverse: https://www.tidyverse.org/packages/)
  + phyloseq (https://bioconductor.org/packages/release/bioc/html/phyloseq.html)
  + mirlyn (https://github.com/escamero/mirlyn)

## Installation and setup:
comets is a bash script. Just download comets and comets_plots.R and place them in your local bin directory. If you don't have a bin directory just type mkdir bin. If your bin directory is not in $PATH add the following to the end of your .bashrc file:

```{bash, eval=FALSE, echo=TRUE}
export PATH="$HOME/bin:$PATH"
```

Then reload .bashrc:
```{bash, eval=FALSE, echo=TRUE}
source ./bashrc
```

Make sure that the first line in comets_plots.R has the correct path to Rscript. To check your Rscript path type:
```{bash, eval=FALSE, echo=TRUE}
which Rscript
```

If necessary change the first line of comets_plots.R for the actual Rscript path

## Usage:
IMPORTANT NOTE: Before start please provide a file named phyloseq_tables/SAM.table. This is a tab separated metadata table to be loaded by phyloseq. It must contain at least three columns: Files, Sample and SampleType  
Samples must be listed in alphabetical order. Example: If you want to compare sample1 (paired end, seawater) and sample2 (single end freshwater) you will have three fastq.gz files (sample1_1.fastq.gz, sample1_2.fastq.gz and sample2.fastq.gz). The SAM.table within the phyloseq_tables/ directory should look like this:

Files|Sample|SampleType
--|--|--
sample1.kaiju.out|sample1|seawater
sample2.kaiju.out|sample2|freshwater

### Starting from raw reads:
```{bash, eval=FALSE, echo=TRUE}
comets [-z INTEGER] [-d path/to/kaiju_db.fmi] -f pe|se|pe_se -g Bacteria|Archaea|Eukaryota|Viruses
```
This mode will perform the whole process starting from fastp processing of raw reads.

* -z - Optional argument to specify the number of parallel threads for the taxonomic classification. Default=1
* -d - Optional argument to specify the path to the kaiju databse for the taxonomic classification. Default=kaiju_db_nr_euk.fmi
* -f - Mandatory argument to specify if files are paired end (pe), single end (se) or both (pe_se). Such files must be in the working directory
* -g - Mandatory argument to specify the taxonomic group to analyze (either Bacteria or Archaea or Eukaryota or Viruses). Names must be provided with the first letter in uppercase

Example paired end files: sample1_1.fastq.gz, sample1_2.fastq.gz, sample2_1.fastq.gz, sample2_2.fastq.gz  
Example single end files: sample1.fastq.gz, sample2.fastq.gz  
NOTE: Single end files name must not have _1 or _2 suffix before .fastq.gz extension. Otherwise the script will wet very confused  

### Starting from fastp quality filtered reads:
```{bash, eval=FALSE, echo=TRUE}
comets [-z INTEGER] [-d path/to/kaiju_db.fmi] -k pe|se|pe_se -g Bacteria|Archaea|Eukaryota|Viruses
```
This mode will assume that reads are already quality filtered and will start the process from the kaiju taxonomic classification step.

* -z - Optional argument to specify the number of parallel threads for the taxonomic classification. Default=1
* -d - Optional argument to specify the path to the kaiju databse for the taxonomic classification. Default=kaiju_db_nr_euk.fmi
* -k - Mandatory argument to specify if files are paired end (pe), single end (se) or both (pe_se)
* -g - Mandatory argument to specify the taxonomic group to analyze (either Bacteria or Archaea or Eukaryota or Viruses). Names must be provided with the first letter in uppercase

Example paired end files: paired_end/fastp/sample1_1.fastp.fastq.gz, paired_end/fastp/sample1_2.fastp.fastq.gz, paired_end/fastp/sample2_1.fastp.fastq.gz, paired_end/fastp/sample2_2.fastp.fastq.gz  
Example single end files: single_end/fastp/sample1.fastp.fastq.gz, single_end/fastp/sample2.fastp.fastq.gz  
NOTE: Single end files name must not have _1 or _2 suffix before .fastp.fastq.gz extension. Otherwise the script will wet very confused 
NOTE: Input files must be within a paired_end/fastp/ or single_end/fastp/ directory, respectively

### Starting from kaiju.out files:
```{bash, eval=FALSE, echo=TRUE}
comets -t pe|se|pe_se -g Bacteria|Archaea|Eukaryota|Viruses
```
This mode will assume that taxonomic classification is already done and will start the process from the kaiju2table classification summary.

* -t - Mandatory argument to specify if files are paired end (pe), single end (se) or both (pe_se)
* -g - Mandatory argument to specify the taxonomic group to analyze (either Bacteria or Archaea or Eukaryota or Viruses). Names must be provided with the first letter in uppercase

Example paired end files: paired_end/fastp/kaiju/sample1.kaiju.out, paired_end/fastp/kaiju/sample2.kaiju.out  
Example single end files:  single_end/fastp/kaiju/sample1.kaiju.out,  single_end/fastp/kaiju/sample2.kaiju.out  
NOTE: Input files must be within a paired_end/fastp/kaiju/ or single_end/fastp/kaiju/ directory, respectively

### Starting from classification summary tables.kaiju.tsv:
```{bash, eval=FALSE, echo=TRUE}
comets -g Bacteria|Archaea|Eukaryota|Viruses
```
This mode will assume that classification summary tables exist and will start the process from the taxonomic group table formatting to produce taxonomic group specific phyloseq tables. The user can loop this step with any of the valid options to get the phyloseq tables and plots for all taxonomic groups. If all phyloseq tables are already made and want to proceed only to plot generation the user must specify which taxonomic group to plot.

* -g - Mandatory argument to specify the taxonomic group to analyze (either Bacteria or Archaea or Eukaryota or Viruses). Names must be provided with the first letter in uppercase

Example input files: kaiju2table/sample1.kaiju.tsv, kaiju2table/sample2.kaiju.tsv  
NOTE: Input files must be within a kaiju2table/ directory  
NOTE: comets_plots.R will save final dataframes as .RData files to phyloseq_tables directory. These can be loaded in RStudio with load("filename.RData") for custom processing

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
