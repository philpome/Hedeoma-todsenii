library(poppr)
library(mmod)
library(magrittr)
library(treemap)
setwd("C:/Users/Megan/Desktop")
hed2018<-read.genalex("hed-cleanup-final.csv")
#create a new dataset with only informative loci
ihed2018<-informloci(hed2018)
(ihedlt <- locus_table(ihed2018))
poppr(ihed2018)

basic_stats <- poppr(ihed)
write.table(basic_stats, file = "hed_stats.csv", sep = ",", quote = FALSE, row.names = F)
#generate a genotype accumulation curve
gac<-genotype_curve(ihed2018,sample=1000,quiet=TRUE)
#number of alleles, simpson index, and gene diversity
(hedlt <- locus_table(hed2018))
#missing data
info_table(hed2018, type = "missing", plot = TRUE)
#diversity measures
#MLG=#multilocus genotypes, eMLG=expected number of multilocus genotypes, SE=standard error for MLG
#H=shannon-wiener, G=stoddart and taylor, lambda=simpson's index, E.5=evenness, 
#Hexp=nei's, Ia=association index, rbard=standardized Ia
poppr(hed2018)
locus_table(hed2018)
#split the strata
splitStrata(ihed2018) <- ~Loc/Year/Pop
#look at the strata distributions to check
library("dplyr")
hedstrata <- strata(ihed2018) %>%     
  group_by(Year, Loc, Pop) %>%
  summarize(Count = n())

hedstrata
#check HW equilibrium
library(pegas)
(hedhwe.full <- hw.test(ihed2018, B = 1000)) 
#visualize population level HW
(hedhwe.pop <- seppop(ihed2018) %>% lapply(hw.test, B = 0))
(hedhwe.mat <- sapply(hedhwe.pop, "[", i = TRUE, j = 3))
#visualize the HW by pop
alpha  <- 0.05
newmat <- hedhwe.mat
newmat[newmat > alpha] <- 1
library("lattice")
levelplot(t(newmat))
#to rotate the x-axis labels (but be careful running this code for some reason)
levelplot(t(newmat), scales=list(x=list(rot=45)))
#genotypic richness
library(vegan)
#check genotypic richness by population
setPop(ihed2018) <- ~Pop
(ihed2018_diversity <- poppr(ihed2018))
#plot MLG's (multilocus genotypes) to check richness
hed.tab <- mlg.table(ihed2018, plot=FALSE)
min_sample <- min(rowSums(hed.tab))
rarecurve(hed.tab, sample = min_sample, xlab = "Sample Size", ylab = "Expected MLGs")
title("Rarefaction of Hedeoma Populations")
#genotypic diversity
#first, calculate simpson's index for all pop's
N      <- ihed2018_diversity$N      # number of samples
lambda <- ihed2018_diversity$lambda # Simpson's index
(N/(N - 1)) * lambda              # Corrected Simpson's index
#calculate genotypic evenness
hed.tab <- mlg.table(ihed2018)
#check for linkage disequilibrium (clonal or sexual pops?)
WS1 <- popsub(ihed2018, "WS_17_1")
WS1
#WS1 shows us clones - if we have 100 original MLGs, we should have 100 individuals
ia(WS1, sample = 999)
WS2 <- popsub(ihed2018, "WS_17_2")
WS2
ia(WS2, sample = 999)
#but there were clones! so let's correct that:
WS2 %>% clonecorrect %>% ia(sample = 999)
#it's still in linkage disequilibrium!
LNF3 <- popsub(ihed2018, "LNF_17_3")
LNF3
ia(LNF3, sample = 999)
WS4 <- popsub(ihed2018, "WS_16_4")
WS4
ia(WS4, sample = 999)
#disequilibrium! but let's correct for clones
WS4 %>% clonecorrect %>% ia(sample = 999)
WS5 <- popsub(ihed2018, "WS_16_5")
WS5
ia(WS5, sample = 999)
WS6 <- popsub(ihed2018, "WS_15_6")
WS6
ia(WS6, sample = 999)
WS6 %>% clonecorrect %>% ia(sample = 999)
WS7 <- popsub(ihed2018, "WS_15_7")
WS7
ia(WS7, sample = 999)
#calculate pairwise r-bar for all loci
WSall <- popsub(ihed2018loc, "WS")
wspair <- WSall %>% clonecorrect(strata = ~Loc) %>% pair.ia
LNFall <- popsub(ihed2018loc, "LNF")
lnfpair <- LNFall %>% clonecorrect(strata = ~Loc) %>% pair.ia
plotrange <- range(c(wspair, lnfpair), na.rm = TRUE)
plot(wspair, limits = plotrange)
plot(lnfpair, limits = plotrange)
#let's visualize by WS vs LNF
ihed2018loc<-ihed2018
splitStrata(ihed2018loc) <- ~Loc/Year/Pop
setPop(ihed2018loc) <- ~Loc
ihed2018loc
(ihed2018loc_diversity <- poppr(ihed2018loc))
ihed2018loc.tab <- mlg.table(ihed2018loc, plot = FALSE)
min_sample <- min(rowSums(ihed2018loc.tab))
rarecurve(ihed2018loc.tab, sample = min_sample, xlab = "Sample Size", ylab = "Expected MLGs")
title("Rarefaction of White Sands & Lincoln National Forest populations")
#correct Simpson's Index for sample size
N      <- ihed2018loc_diversity$N      # number of samples
lambda <- ihed2018loc_diversity$lambda # Simpson's index
(N/(N - 1)) * lambda              # Corrected Simpson's index
#calculate the clonal fraction
hedMLG <- ihed2018loc_diversity$MLG
1-(hedMLG/N)
#population differentiation with Gst
library("mmod")
Gst_Hedrick(ihed2018)
#look at genetic distance for 10 sample individuals
library("ape")
library("magrittr")
ten_samples <- sample(nInd(ihed2018), 10)
mic10       <- ihed2018[ten_samples]
(micdist    <- provesti.dist(mic10))
#make a tree
theTree <- micdist %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.
#now do it by stratified data
ihed2018strat <- ihed2018
strata(ihed2018strat) <- data.frame(other(ihed2018strat))
ihed2018strat
nameStrata(ihed2018strat) <- ~Loc/Year/Pop

# Analysis
set.seed(999)
ihed2018strat %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist)
#K-means hierarchical clustering to check clonality
#first, perform a cluster analysis
#choose 50 PCs to retain at first, then choose the number of clusters with the lowest BIC
WS1clust <- find.clusters(WS1)
WS1clust
WS2clust <- find.clusters(WS2)
WS2clust
LNF3clust <- find.clusters(LNF3)
LNF3clust
WS4clust <-find.clusters(WS4)
WS4clust
WS5clust <- find.clusters(WS5)
WS5clust
WS6clust <- find.clusters(WS6)
WS6clust
WS7clust <- find.clusters(WS7)
WS7clust
#now make a tree
WS <- popsub(ihed2018loc, "WS")
WSclust <- find.clusters(WS)
LNF <- popsub(ihed2018loc, "LNF")
LNFclust <- find.clusters(LNF)
#length of microsatellite repeats per locus
hedreps <- c(2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
LNFtree <- bruvo.boot(LNF, replen = hedreps, cutoff = 50, quiet = TRUE, tree = "njs")
cols <- rainbow(4)
plot.phylo(LNFtree, cex = 0.8, font = 2, adj = 0, tip.color = cols[LNFclust$grp],
           label.offset = 0.0125)
nodelabels(LNFtree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)
#amova
ihed2018split<-ihed2018
table(strata(ihed2018split, ~Loc))
table(strata(ihed2018split, ~Loc/Pop, combine = FALSE))
#AMOVA without clone correction
Hedamova <- poppr.amova(ihed2018, ~Pop)
#AMOVA with clone correction
Aeutamovacc <- poppr.amova(Aeut, ~Pop/Subpop, clonecorrect = TRUE)
#visualize pop diff
library("adegenet")
library("pegas")
library("mmod")
library("reshape2")
library("ggplot2")
diff_hed <- diff_stats(ihed2018)
per.locus <- melt(diff_hed$per.locus, varnames = c("Locus", "Statistic"))
stats     <- c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
glob      <- data.frame(Statistic = stats, value = diff_hed$global)
head(per.locus)
head(glob)
ggplot(per.locus, aes(x = Statistic, y = value)) +
  geom_boxplot() +
  geom_point() +
  geom_point(size = rel(3), color = "red", data = glob) +
  ggtitle("Estimates of population differentiation")
#bootstrap it for confidence intervals
set.seed(20151219) # Be sure to set a seed for any random analysis!
bs_reps <- chao_bootstrap(ihed2018, nreps = 100)
summarise_bootstrap(bs_reps, D_Jost) # Using the D_Jost function to summarize.
summarise_bootstrap(bs_reps, Gst_Hedrick)
#pairwise g'st
pairwise_Gst_Hedrick(ihed2018, linearized = FALSE)
#pairwise D
pairwise_D(ihed2018, linearized = FALSE)
#look at observed vs expected heterozygosity
div <- summary(ihed2018)
div
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
plot(div$Hobs, div$Hexp, xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
#is there a significant difference b/w observed & expected?
bartlett.test(list(div$Hexp, div$Hobs))
#check within pops
div1 <- summary(WS1)
div1
bartlett.test(list(div1$Hexp, div1$Hobs))
div2 <- summary(WS2)
div2
bartlett.test(list(div2$Hexp, div2$Hobs))
div3 <- summary(LNF3)
div3
bartlett.test(list(div3$Hexp, div3$Hobs))
