#This script represents an example of alternative processing for the bray-curtis dissimilarity matrix resulting from bray.sh script. 

species_df <- read.delim("/Volumes/TOSHIBA EXT/DomosArqueano/figuras/bray_new/network/species_matrix.txt", check.names = F, row.names = 1)
species_df <- species_df[1:(length(species_df)-1)]
species_df_mat <- as.matrix(species_df)
species_df_mat_sym <- as.matrix(Matrix::forceSymmetric(species_df_mat,uplo="L"))
# library(pheatmap)
pheatmap(species_df_mat_sym, cellheight = 4, cellwidth = 4, fontsize_row = 4, fontsize_col = 4, legend_breaks = c(0.25,0.5,0.75))
# library(phangorn)
nnet <- neighborNet(species_df_mat_sym)
write.nexus.networx(nnet, file = "species_nnet.nexus")

# library(igraph)
# library(RColorBrewer)
species_df_r <- 1-species_df
species_df_r_mat <- as.matrix(species_df_r)
species_df_r_mat_sym <- as.matrix(Matrix::forceSymmetric(species_df_r_mat,uplo="L"))
# write.csv(species_df_r_mat_sym, file = "species_df_r_mat_sym.csv")
metadata <- read.delim("/Volumes/TOSHIBA EXT/DomosArqueano/Otros_metagenomas/Bray-Curtis_UPGMA/netweork/metadata.txt")
coul <- brewer.pal(nlevels(as.factor(metadata$SampleType)), "Set2")
my_color <- coul[as.numeric(as.factor(metadata$SampleType))]

network_r_100 <- graph_from_adjacency_matrix(species_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_100, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

# library(reshape2)
species_df_r_mat.melt <- melt(species_df_r_mat)
species_df_r_mat.melt <- subset(species_df_r_mat.melt, Var1 != Var2)
species_df_r_mat.melt <- species_df_r_mat.melt %>% filter(!is.na(value))

Q1 <- unname(quantile(species_df_r_mat.melt$value, prob=c(.25))) # 75%
MED <- median(species_df_r_mat.melt$value) # 50%
Q3 <- unname(quantile(species_df_r_mat.melt$value, prob=c(.75))) # 25%
SD1 <- sd(species_df_r_mat.melt$value)+mean(species_df_r_mat.melt$value) # 15.9%

species_df_r_mat_sym[species_df_r_mat_sym<Q1] <- 0
network_r_75 <- graph_from_adjacency_matrix(species_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_75, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

species_df_r_mat_sym[species_df_r_mat_sym<MED] <- 0
network_r_50 <- graph_from_adjacency_matrix(species_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_50, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

species_df_r_mat_sym[species_df_r_mat_sym<Q3] <- 0
network_r_25 <- graph_from_adjacency_matrix(species_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_25, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")

species_df_r_mat_sym[species_df_r_mat_sym<SD1] <- 0
network_r_16 <- graph_from_adjacency_matrix(species_df_r_mat_sym, mode='undirected', weighted = TRUE, diag = F)
par(mar=c(0,0,0,0))
set.seed(4)
plot(network_r_16, layout=layout.fruchterman.reingold, main="", vertex.color=my_color, vertex.size=6, vertex.label.cex=0.5, vertex.label.color="black", vertex.frame.color="black")


