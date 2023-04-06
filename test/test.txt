This folder contains test data to run COMETS

It includes four publicly available metagenomes from Hypersaline environments:
2010Bt1.fastq.gz
2010Bt3.fastq.gz

and Marine environments:
Kingman_atoll.fastq.gz
Sargasso_sea.fastq.gz

Note that the phyloseq_tables directory is also included along with the SAM.table file

Since all files are single-end raw reads, a full test command to analyze viral reads from scratch should look like:

comets -p se -f auto -k auto -z 4 -t -g Viruses -r
