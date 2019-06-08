library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(reshape2)
library(LEA)
library(netview)
library(networkD3)

##
##Data Handling
##

#set working directory
setwd(dir = "C:/Users/Brandon/Google Drive/WVU/RRISSR/ISSRseq/C_bentleyi/")

#import GATK filtered VCF file
filtered_vcf <- read.vcfR("C_bentleyi_filtered_SNPs.vcf")

#create optional matrix of subsetted random SNPs from VCF
set.seed(123)
subset.1 <- sample(size = 100, x = c(1:nrow(filtered_vcf)), replace = FALSE)

subset_filtered_vcf <- filtered_vcf[subset.1,]

#import sample ID and population map
pop_data <- read.table("C_bentleyi_populations.txt", sep = "\t", header = FALSE)

#convert VCF object into genlight object
gl_vcf_subset <- vcfR2genlight(subset_filtered_vcf)

gl_vcf_full <- vcfR2genlight(filtered_vcf)

#try to extract a single random SNP from each locus
#locus_split <- seploc(gl_vcf_full, block.size = 1)
#set.seed(123)
#subset.1 <- sample(size = 100, x = c(1:nrow(locus_split)), replace = FALSE)

#set ploidy
ploidy(gl_vcf_subset) <- 2
ploidy(gl_vcf_full) <- 2

#assign populations and color by pop
pop(gl_vcf_subset) <- pop_data$V2
pop(gl_vcf_full) <- pop_data$V2

cols_subset <- brewer.pal(n = nPop(gl_vcf_subset), name = "Dark2")
cols_full <- brewer.pal(n = nPop(gl_vcf_full), name = "Dark2")

write.geno(as.matrix(gl_vcf_subset), "C_bentleyi_subset_genofile.geno")
write.geno(as.matrix(gl_vcf_full), "C_bentleyi_full_genofile.geno")

##
##Data Analysis
##

##create UPGMA tree

upgma_tree <- aboot(gl_vcf_full, tree = "upgma", distance = bitwise.dist, missing = "loci", mcutoff = 0.1, sample = 10, showtree = FALSE, cutoff = 49, quiet = T)

plot.phylo(upgma_tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols_full[pop(gl_vcf_full)])

nodelabels(upgma_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)

#create minimum spanning network

# dist_subset_filtered_vcf <- bitwise.dist(gl_vcf_full, missing_match = FALSE)
# 
# msn_subset_filtered_vcf <- poppr.msn(gl_vcf_full, dist_subset_filtered_vcf, showplot = FALSE, include.ties = T, clustering.algorithm = "nearest_neighbor")
# 
# node.size <- rep(2, times = nInd(gl_vcf_full))
# names(node.size) <- indNames(gl_vcf_full)
# vertex.attributes(msn_subset_filtered_vcf$graph)$size <- node.size
# plot_poppr_msn(gl_vcf_full, msn_subset_filtered_vcf , palette = brewer.pal(n = nPop(gl_vcf_full), name = "Dark2"), gadj = 70)

#conduct PCA analysis

pca_vcf <- glPca(gl_vcf_full, nf =3)
vcf_pca_scores <- as.data.frame(pca_vcf$scores)
vcf_pca_scores$pop <- pop(gl_vcf_full)

#plot PCA results

# p <- ggplot(vcf_pca_scores, aes(x=PC1, y=PC2, colour=pop)) 
# p <- p + geom_point(size=2)
# p <- p + stat_ellipse(level = 0.95, size = 1)
# p <- p + scale_color_manual(values = cols) 
# p <- p + geom_hline(yintercept = 0) 
# p <- p + geom_vline(xintercept = 0) 
# p <- p + theme_bw()
# p

#conduct Discriminant Analysis of Principle Components

vcf_dapc <- dapc(gl_vcf_full, n.pca = 3, n.da = 2)

#plot DAPC results
scatter(vcf_dapc, col = cols_full, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = FALSE,
        posi.pca = "topleft", cleg = 0.75)

#plot DAPC estimated population memberships
compoplot(vcf_dapc, col = cols_full, posi = 'top')

dapc.results <- as.data.frame(vcf_dapc$posterior)
dapc.results$pop <- pop(gl_vcf_full)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- melt(dapc.results)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

dapc_plot <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
dapc_plot <- dapc_plot + geom_bar(stat='identity') 
dapc_plot <- dapc_plot + scale_fill_manual(values = cols_full) 
dapc_plot <- dapc_plot + facet_grid(~Original_Pop, scales = "free")
dapc_plot <- dapc_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dapc_plot

######conduct structure-like analysis http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf

#run for a range of K values -- alpha is a non-negative regularization parameter, 10 is a good choice
#tolerance is the difference between two runs to measure convergence
obj.snmf = snmf("C_bentleyi_full_genofile.geno", K = 1:10, ploidy = 2, entropy = T,
                percentage = 0.25, alpha = 100, tolerance = 0.0000000001, project = "new")

#plot the cross entropy 
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)

#plot the appropriate K value
qmatrix = Q(obj.snmf, K = 2)

barplot(t(qmatrix), col = c("orange","violet","lightgreen", "blue", "red", "black", "gold","tan"), border = NA, space = 0,
        ylab = "Admixture coefficients", names.arg = gl_vcf_full$ind.names, las=2)

######conduct NETVIEW analysis

#create distance matrix
dist_matrix <- bitwise.dist(gl_vcf_full, mat = TRUE, missing_match = FALSE)

#create metadata file
net_metadata <- as.data.frame(cbind(pop_data, cols_full[pop(gl_vcf_full)]))
names(net_metadata)[1] <- "ID"
names(net_metadata)[2] <- "Population"
names(net_metadata)[3] <- "Colour"

#conduct analysis
netview_options <- netviewOptions(mknnWeights=TRUE, nodeGroup="Population", 
                                  nodeID="ID", nodeColour="Colour", communityAlgorithms = "Infomap")

netview_output <- netview(distMatrix = dist_matrix, metaData = net_metadata, 
                          k=1:20, cluster = TRUE, options = netview_options, mst = FALSE)

#plot netview results
netview_K_plot <- plotSelection(netview_output, options = netview_options)

netview_K_plot

netview_selected_K <- netview_output$k12
plot(netview_selected_K, vertex.size=10, vertex.label=net_metadata$ID)

######conduct AMOVA

