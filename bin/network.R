#This script represents an example of alternative processing for the bray-curtis dissimilarity matrix resulting from bray.sh script. 

otu_df <- read.delim("/path/to/otu_matrix.txt", check.names = F, row.names = 1)
#otu_df <- otu_df[1:(length(otu_df)-1)]
otu_df_mat <- as.matrix(otu_df)
otu_df_mat_sym <- as.matrix(Matrix::forceSymmetric(otu_df_mat,uplo="L"))

library(pheatmap)
pheatmap(otu_df_mat_sym, cellheight = 4, cellwidth = 4, fontsize_row = 4, fontsize_col = 4, legend_breaks = c(0.25,0.5,0.75))

library(phangorn)
nnet <- neighborNet(otu_df_mat_sym)
write.nexus.networx(nnet, file = "otu_nnet.nexus") # This file can be visualized with SplitsTree software

library(igraph)
library(RColorBrewer)
otu_df_r <- 1-otu_df
otu_df_r_mat <- as.matrix(otu_df_r)
otu_df_r_mat_sym <- as.matrix(Matrix::forceSymmetric(otu_df_r_mat,uplo="L"))
# write.csv(otu_df_r_mat_sym, file = "otu_df_r_mat_sym.csv")
metadata <- read.delim("/path/to/metadata.txt") # This file can be derived from SAM.table. It must contain two columns: Sample and SampleType. Make sure
# the samples are in the same order as in otu_matrix.txt
coul <- brewer.pal(nlevels(as.factor(metadata$SampleType)), "Set2")
my_color <- coul[as.numeric(as.factor(metadata$SampleType))]

network_r_100 <- graph_from_adjacency_matrix(otu_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_100, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

library(reshape2)
otu_df_r_mat.melt <- melt(otu_df_r_mat)
otu_df_r_mat.melt <- subset(otu_df_r_mat.melt, Var1 != Var2)
otu_df_r_mat.melt <- otu_df_r_mat.melt %>% filter(!is.na(value))

Q1 <- unname(quantile(otu_df_r_mat.melt$value, prob=c(.25))) # 75%
MED <- median(otu_df_r_mat.melt$value) # 50%
Q3 <- unname(quantile(otu_df_r_mat.melt$value, prob=c(.75))) # 25%
SD1 <- sd(otu_df_r_mat.melt$value)+mean(otu_df_r_mat.melt$value) # 15.9%

otu_df_r_mat_sym[otu_df_r_mat_sym<Q1] <- 0
network_r_75 <- graph_from_adjacency_matrix(otu_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_75, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

otu_df_r_mat_sym[otu_df_r_mat_sym<MED] <- 0
network_r_50 <- graph_from_adjacency_matrix(otu_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_50, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

otu_df_r_mat_sym[otu_df_r_mat_sym<Q3] <- 0
network_r_25 <- graph_from_adjacency_matrix(otu_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_25, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

otu_df_r_mat_sym[otu_df_r_mat_sym<SD1] <- 0
network_r_16 <- graph_from_adjacency_matrix(otu_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_16, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")


