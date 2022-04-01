#!/usr/local/bin/Rscript
# Load libraries
print("Loading libraries")
library(tibble) # Needed for converting column to row names
library(phyloseq) # Needed to import, store, analyze, and graphically display complex phylogenetic sequencing data
library(ggplot2) # For high quality graphics
library(mirlyn) # For rarefaction curves
library(dplyr) # To manipulate dataframes

# Set arguments
args <- commandArgs(TRUE)

# Read tables
print("Reading tables")
otu_mat <- read.delim(args[1], check.names = F) # check.names is set to FALSE to allow sample names starting with numbers
tax_mat <- read.delim(args[2])
sam_tab <- read.delim(args[3])

# Define row names
print("Creating phyloseq object")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("otu")
sam_tab <- sam_tab %>%
  tibble::column_to_rownames("Files")

# Transform into matrixes (except for sample table)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(sam_tab)
PhyloData <- phyloseq(OTU, TAX, samples)

# Rarefaction curves
print("Making rarefaction curves")
rarefy_data <- rarefy_whole_rep(PhyloData, rep = 100)
colnames(rarefy_data)[1] <- "Files"
colnames(rarefy_data)[2] <- "OTUcount"
save(rarefy_data, file = "rarefy_data.RData")
pdf("rarefaction.pdf")
ggplot(rarefy_data, aes(x=LibSize, y=OTUcount, color=Sample)) +
  geom_line()
dev.off()

# Normalization by median sequencing depth
# 1. Get median depth from all samples
# 2. For each sample, multiply the proportion of every OTU by the median depth
print("Normalizing read counts")
total = median(sample_sums(PhyloData))
standf = function(x, t=total) round(t * (x / sum(x)))
PhyloData = transform_sample_counts(PhyloData, standf)

# Bar plots
print("Drawing Phylum bar plot")
PhyloData_abund <- filter_taxa(PhyloData, function(x) sum(x > total*0.01) > 0, TRUE)
PhyloData.frame <- psmelt(PhyloData_abund)
phylumPhyloData.frame <- PhyloData.frame[!is.na(PhyloData.frame$Phylum), ]
save(phylumPhyloData.frame, file = "phylumPhyloData.RData")
pdf("phylumBarPlot.pdf")
ggplot(phylumPhyloData.frame, aes(x=sample_Sample, y=Abundance, fill=Phylum)) +
  geom_bar(position = "fill", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative abundance") +
  xlab("")
dev.off()

print("Drawing Family bar plot")
familyPhyloData.frame <- PhyloData.frame[!is.na(PhyloData.frame$Family), ]
save(familyPhyloData.frame, file = "familyPhyloData.RData")
pdf("familyBarPlot.pdf")
ggplot(familyPhyloData.frame, aes(x=sample_Sample, y=Abundance, fill=Family)) +
  geom_bar(position = "fill", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative abundance") +
  xlab("")
dev.off()

# Alpha diversity plots
print("Drawing Shannon alpha diversity plot")
save(PhyloData, file = "PhyloData_norm.RData")
pdf("Shannon.pdf")
plot_richness(PhyloData, measures = "Shannon", x = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# NMDS Bray-Curtis plots
print("Drawing Bray-Curtis NMDS plot")
PhyloData.ord <- ordinate(PhyloData, "NMDS", "bray")
save(PhyloData.ord, file = "PhyloData_NMDS.RData")
pdf("NMDS.pdf")
plot_ordination(PhyloData, PhyloData.ord, type = "samples", color = "Sample", label = "Sample") + 
  theme_bw() + 
  theme(legend.position = "none")
dev.off()