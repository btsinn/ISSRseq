##Sandra Simon & Brandon Sinn -- 2020
##Generalizable R script to conduct 
##common population genetic analyses
##on ISSRseq pipeline output.

##Must be executed in the output directory of an ISSRseq analysis

#Packages for data processing from VCF file
library(vcfR)

#Package for filtering SNPs in VCF files
library(radiator)

#Package for population genetic analyses (AMOVA, bootstrapping ect.)
library(poppr)
library(adegenet)
library(ape)
library(mmod)
library(magrittr)
library(hierfstat)
library(pegas)

#Package for local adaptation analysis
library(pcadapt)

#Package for admixture analysis
library(LEA)

#Package for dendrogram graphing
library(networkD3)

#Packages for figures and plots
library(ggplot2)
library(plot3D)
library(plotly)
library(igraph)
library(RColorBrewer)
library(grid);
library(R.devices);

#Package for Wes Anderson themed color palettes (Range 4-5 colors). 
#Choices are: BottleRocket1, BottleRocket2, Rushmore1, Royal1, Royal2, Zissou1, Darjeeling1, Darjeeling2, Chevalier1,
#FantasticFox1, Moonrise1, Moonrise2, Moonrise3, Cavalcanti1, GrandBudapest1, GrandBudapest2, IsleofDogs1, IsleofDogs2
library(wesanderson)

#create output directory
dir.create(path = population_genetics_outputs)

#import GATK filtered VCF file
filtered_vcf <- read.vcfR("filtered_SNPs.vcf")
#Variants = 7378

#create optional matrix of subsetted random SNPs from VCF
set.seed(123)
subset.1 <- sample(size = 100, x = c(1:nrow(filtered_vcf)), replace = FALSE)

subset_filtered_vcf <- filtered_vcf[subset.1,]

#import sample ID and population map
pop_data <- read.table("samples.txt", sep = "\t", header = FALSE)

#VCF can be stored as different object types including genlight (binary SNPs), genid (allelic counts), and genclone (multilocus genotype for clonal populations)
#example convert VCF object into genlight object
my_genlightsubset <- vcfR2genlight(subset_filtered_vcf)
my_genlight <- vcfR2genlight(filtered_vcf)

my_genindsubset <- vcfR2genind(subset_filtered_vcf)
my_genind <- vcfR2genind(filtered_vcf)

#Sorghum_16S.t<-as.data.frame(t(my_genlight$gen))
write.csv(file="Individuals.csv", Individuals, row.names=TRUE)

#set ploidy
ploidy(my_genlightsubset) <- 2
ploidy(my_genlight) <- 2

ploidy(my_genindsubset) <- 2
ploidy(my_genind) <- 2

#assign populations and color by pop
pop(my_genlightsubset) <- pop_data$V2
pop(my_genlight) <- pop_data$V2

pop(my_genindsubset) <- pop_data$V2
pop(my_genind) <- pop_data$V2

#assign strata for analysis
samplingstrata <- as.data.frame(pop_data$V4)

strata(my_genlight) <- samplingstrata
strata(my_genind) <- samplingstrata

splitStrata(my_genlight) <- ~Individual/Pop/Region
splitStrata(my_genind) <- ~Individual/Pop/Region

#strata(my_genindsubset) <- samplingstrata
#strata(my_genind) <- samplingstrata

#splitStrata(my_genindsubset) <- ~Local/Pop/Individual
#splitStrata(my_genind) <- ~Local/Pop/Individual

#create snpclone and genclone datasets
my_snpclonesubset <- poppr::as.snpclone(my_genlightsubset)
my_snpclone <- poppr::as.snpclone(my_genlight)

my_genclonesubset <- poppr::as.genclone(my_genindsubset)
my_genclone <- poppr::as.genclone(my_genind)

#assign population colors for figure generation
#cols_subset <- brewer.pal(n = nPop(my_genlightsubset), name = "Dark2")
#cols_full <- brewer.pal(n = nPop(my_genlight), name = "Dark2")

striata_AMOVA <- poppr.amova(my_genclone, ~Pop, within = TRUE)
#results for Striata
#Df    Sum Sq  Mean Sq
#Between samples 10  534.8965 53.48965
#Within samples  63  857.8029 13.61592
#Total           73 1392.6994 19.07807
#componentsofcovariance
#                            Sigma         %
#Variations  Between samples  6.910202  33.66541
#Variations  Within samples  13.615919  66.33459
#Total variations            20.526121 100.00000
#statphi
#Phi
#Phi-samples-total 0.3366541

set.seed(1999)
striatasignif   <- randtest(striata_AMOVA, nrepet = 999)

#eps(file="IndividualAMOVA.eps", width = 10, height = 8)
plot(striatasignif)
#dev.off()

#Monte-Carlo test for Bentleyi
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 15.29718 
#Based on 999 replicates
#Simulated p-value: 0.001 
#Alternative hypothesis: greater 
#Std.Obs Expectation    Variance 
#9.18337132 -0.01812962  2.78129718 

#Monte-Carlo test for Striata
#Call: as.randtest(sim = res, obs = sigma[1])
#Observation: 6.910202 
#Based on 999 replicates
#Simulated p-value: 0.001 
#Alternative hypothesis: greater 
#Std.Obs Expectation    Variance 
#17.59408942 -0.02002478  0.15515342

striata_AMOVA <- poppr.amova(my_genclone, ~Pop/Region, within = TRUE)
#results
#Df    Sum Sq  Mean Sq
#Between Pop                10  534.8965 53.48965
#Between samples Within Pop  3  140.7333 46.91110
#Within samples             60  717.0696 11.95116
#Total                      73 1392.6994 19.07807
#componentsofcovariance
#Sigma         %
#Variations  Between Pop                 3.702303  18.19888
#Variations  Between samples Within Pop  4.690111  23.05451
#Variations  Within samples             11.951160  58.74661
#Total variations                       20.343574 100.00000
#statphi
#Phi
#Phi-samples-total 0.4125339
#Phi-samples-Pop   0.2818361
#Phi-Pop-total     0.1819888

set.seed(1999)
striatasignif   <- randtest(striata_AMOVA, nrepet = 999)

#eps(file="PopAMOVA.eps", width = 10, height = 8)
plot(striatasignif)
#dev.off()

#Monte-Carlo tests
#Call: randtest.amova(xtest = striata_AMOVA, nrepet = 999)
#Number of tests:   3 
#Adjustment method for multiple comparisons:   none 
#Permutation number:   999 
#Test       Obs    Std.Obs   Alter Pvalue
#1  Variations within samples 11.951160 -18.715494    less  0.001
#2 Variations between samples  4.690111   9.722416 greater  0.001
#3     Variations between Pop  3.702303   1.161586 greater  0.141


pegas.fstat <- Fst(as.loci(my_genind), pop = my_genind$pop)
#apply(pegas.fstat, 2, mean)
#Cannot average across loci if NaNs produced.
#Creating excel file to average manually

Fstatist <- cbind(rownames(pegas.fstat), data.frame(pegas.fstat,row.names=NULL))
write.csv(file="StriataFst.csv", Fstatist, row.names=TRUE)
#Fit = overall inbreeding coefficient of an individual relative to total population (Ind/total),
#Fst = population substructure/differentiation among populations (pop/total) 
#Fis = inbreeding coefficient (ind/pop)
#Striata F-statistics
#Fit          Fst                  Fis 
#0.250276603  0.008313517          0.240801479

#PCAdapt
pcadapt_input <- read.pcadapt(input = "output_filtered061520.vcf", type = c("vcf"))

#Determine number of dimensions (K = 3)
x <- pcadapt(input = pcadapt_input, K = 20)
plot(x, option = "screeplot")

pcadapt_result <- pcadapt(pcadapt_input, K = 3)

pcadapt_scores <- as.data.frame(pcadapt_result$scores) 
Population <- pop(my_genlight)

PCAplot <- ggplot(data = pcadapt_scores, aes(V1, V2)) + 
  geom_point(size = 4, aes(color = Population, text =      paste("Individual:", my_genlight$ind.names)))+
  scale_color_manual(values = wes_palette(11, name = "FantasticFox1", type = "continuous")) + 
  #scale_shape_manual(values=c(0:2,15:17))+
  #xlim(-0.15,0.1) + ylim(-0.4,0.25) +
  theme_minimal()

#eps(file="PCAplot.eps", width = 10, height = 8)
PCAplot
#dev.off()

InteractivePCAplot <- ggplotly(PCAplot)

plot(pcadapt_result, option = "manhattan")

Manhattan <- plot(pcadapt_result, option = "manhattan")
Histpvalue <- hist(pcadapt_result$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
qqplot <- plot(pcadapt_result, option = "qqplot")
StatDist <- plot(pcadapt_result, option = "stat.distribution")

#eps(file="Stats.eps", width = 10, height = 8)
#plot(pcadapt_result, option = "stat.distribution")
#dev.off()


#detecting outlier loci analysis using q-values
BiocManager::install("qvalue")
library(qvalue)
qval <- qvalue(pcadapt_result$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
SNPoutliers <- as.data.frame(outliers)

library(dplyr)
SNP_pvalues <- as.data.frame(pcadapt_result$pvalues)
SNP_pvalues <- as.data.frame(SNP_pvalues[!is.na(SNP_pvalues$`pcadapt_result$pvalues`), ])
SNP_pvalues <- SNP_pvalues %>% rename(pvalues = "SNP_pvalues[!is.na(SNP_pvalues$`pcadapt_result$pvalues`), ]")

SNPs <- as.data.frame(pcadapt_result$pass)
SNP_pvalues$outliers <- SNPs$`pcadapt_result$pass`

FinalOutliers <- merge(SNP_pvalues, SNPoutliers, by = "outliers")

write.csv(file="SNPoutlierscheck.csv", FinalOutliers, row.names=TRUE)

#conduct Discriminant Analysis of Principle Components
grp <- find.clusters(my_snpclone, max.n.clust = 40)
vcf_dapc <- dapc(my_snpclone, n.pca = 4, n.da = 5)

#plot DAPC results

pdf(file="DPACplot.pdf", width = 10, height = 8)
scatter(vcf_dapc, col = cols_full, cex = 2, scree.da=FALSE, legend = TRUE)
dev.off()

scatter(vcf_dapc, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)

pdf(file="DensityPlot.pdf", width = 10, height = 8)
scatter(vcf_dapc,1,1, col=cols_full, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
dev.off()

#component plot
compoplot(vcf_dapc, posi = 'top')

DPACmatrix <- as.data.frame(vcf_dapc$posterior)
DPACmatrix$pop <- vcf_dapc$grp
DPACmatrix$individual<-my_snpclone$ind.names
DPACmatrix <- arrange(DPACmatrix, CA)
DPACmatrix <- arrange(DPACmatrix, pop)

#eps(file="DPAC.eps", width = 8, height = 5)
barplot(t(DPACmatrix[,1:11]), col = wes_palette(n = 11, name = "FantasticFox1", type = "continuous"), border = NA, space = 0,
        ylab = "Membership probability", name = DPACmatrix$individual, las=2, legend = TRUE, cex.axis = 0.5, cex.names = 0.3)
#dev.off()

