#Package Install
#Package for data processing from VCF file
library(vcfR)

#Packages for population genetic analyses (AMOVA, bootstrapping ect.)
library(poppr); library(adegenet); library(ape); library(mmod);
library(magrittr); library(hierfstat); library(pegas); library(qvalue)

#Package for local adaptation analysis
library(pcadapt)

#Package for admixture analysis
library(LEA)

#Packages for figures and plots
library(ggplot2); library(plot3D); library(plotly); library(igraph);
library(RColorBrewer); library(grid); library(R.devices)

#Package for Wes Anderson themed color palettes (Range 4-5 colors). 
#Choices are: BottleRocket1, BottleRocket2, Rushmore1, Royal1, Royal2, Zissou1, Darjeeling1, Darjeeling2, Chevalier1,
#FantasticFox1, Moonrise1, Moonrise2, Moonrise3, Cavalcanti1, GrandBudapest1, GrandBudapest2, IsleofDogs1, IsleofDogs2
library(wesanderson)

#Importing and formatting datasets for analysis
# 1. Import GATK filtered VCF file.
filtered_vcf <- read.vcfR("SNP-Matrix.vcf")

# 2. VCF can be stored as different object types including genlight (binary SNPs), genid (allelic counts), and genclone (multilocus genotype for clonal populations).
my_genlight <- vcfR2genlight(filtered_vcf)
my_genind <- vcfR2genind(filtered_vcf)

# 3. Import sample ID and population map.
pop_data <- read.table("Sample-Populations.txt", sep = "\t", header = FALSE)

# 5. Make sure to double check that the order of individuals between your population sample file and SNP matrix match.
my_genlight@ind.names

# 6. Set the ploidy of your individuals and assign the population groups you want to use for PCA. 
ploidy(my_genlight) <- 2
ploidy(my_genind) <- 2
pop(my_genlight) <- pop_data$V3
pop(my_genind) <- pop_data$V3

# 7. Next assign strata to your SNP matrix for use in building hierarchies in AMOVA. 
samplingstrata <- as.data.frame(pop_data$V4)
strata(my_genlight) <- samplingstrata
strata(my_genind) <- samplingstrata
splitStrata(my_genlight) <- ~Individual/Locality/Region
splitStrata(my_genind) <- ~Individual/Locality/Region

# 8. Finally, you can create datasets from your genlight and genind objects that are helpful in analyzing clonal populations.
my_snpclone <- poppr::as.snpclone(my_genlight)
my_genclone <- poppr::as.genclone(my_genind)


#Analysis of molecular variance (AMOVA) with poppr
#AMOVA with hierarchy of Sampling locality within Region
species_AMOVA <- poppr.amova(my_genclone, ~Region/Locality, within = TRUE)
species_AMOVA

#Permutational significance test for AMOVA
set.seed(1999)
speciessignif <- randtest(species_AMOVA, nrepet = 999)

#Plot of AMOVA results
#eps(file="AMOVARegion.eps", width = 10, height = 8)
plot(speciessignif)
#dev.off()

#Convert genind objects from adegenet into a hierfstat data.
my_genindFst <- genind2hierfstat(my_genind)
basicstat <- basic.stats(my_genindFst, diploid = TRUE, digits = 4) 
basicstat$overall

#Calculate pairwise Fst with method options: "Dch","Da","Ds","Fst","Dm","Dr","Cp", "WC84, "Nei87", or "X2"
genet.dist(my_genind, method = "WC84")

#Principle component analysis (PCA) with adegenet package
aPCA <- glPca(my_genlight)

#Percent variance explained by PC axes
eig.perc <- 100*aPCA$eig/sum(aPCA$eig)

#Plot PCA results
colors <- wes_palette(n = 4, name = "Moonrise2")
colors <- colors[as.numeric(my_genlight@pop)]
scatter3D(aPCApoints$PC1, aPCApoints$PC2, aPCApoints$PC3, theta = 290, phi = 40, bty = "g",pch = 20, cex = 2, ticktype = "detailed", colvar = NULL, col = colors)

#Discriminant Analysis of Principle Components with the adegenet package
#Determine an appropriate number of axes and clusters with the find.clusters() command.
grp <- find.clusters(my_snpclone, max.n.clust = 40)

#Run DAPC to determine membership probabilities
vcf_dapc <- dapc(my_snpclone, n.pca = 4, n.da = 5)

#Plot DAPC results
pdf(file="DAPCdotplotAxis13.pdf", width = 10, height = 8)
scatter(vcf_dapc, col = wes_palette(n = 4, name = "Moonrise2"), xax = 1, yax = 3, cex = 2, scree.da=FALSE, legend = TRUE, grp = my_snpclone@pop)
dev.off()

#Component plot of membership probabilities created with ggplot2
#Creating dataset and ordering individuals
DPACmatrix <- as.data.frame(vcf_dapc$posterior)
DPACmatrix$pop <- vcf_dapc$grp
DPACmatrix$individual<-my_snpclone$ind.names
DPACmatrix <- arrange(DPACmatrix, individual)
DPACmatrix <- arrange(DPACmatrix, pop)

#eps(file="species_PCAdapt_p-values.eps", width = 10, height = 8)
hist(pcadapt_result$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#dev.off()
