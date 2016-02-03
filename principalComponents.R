#!/usr/bin/Rscript

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

# All files will be named "*.pcaRes.rdd" and "*.pcaRes.pdf" regardless of args
# So RENAME files after running to distinguish different runs using different subsets of fish.

# This code assumes that we are working only with ONE chromosome at a time
# Note: At most 3 ALT alleles permitted per snp in this version

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn", "paxl", "paxb", "marine-pac")
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb", "marine-pac")
# args <- c("BenlimAllMarine", "chrXXI", "paxl","paxb","pril","prib","qryl","qryb","ensl","ensb")
# args <- c("BenlimAllMarine", "chrXXI", "paxl","paxb","pril","prib","qryl","qryb","ensl","ensb","marine-pac","marine-atl","marine-jap","solitary")

project 		<- args[1]
chrname 		<- args[2]
groupnames 	<- args[3:length(args)]

pcaResFile 	<- paste(project, chrname, "pcaRes.rdd", sep = ".")

Glazerize 	<- TRUE # Requires file "glazerFileS4 NewScaffoldOrder.csv" in current directory

cat("\nProject, chrname, groupnames:\n")
cat(project, chrname, groupnames,"\n")

if(Glazerize){
	vcfresultsfile <- paste(project, ".", chrname, ".vcfresultsNew.rdd", sep = "")
	} else {
	vcfresultsfile <- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")}

# vcfresultsfile
# [1] "BenlimAllMarine.chrXXI.vcfresultsNew.rdd"

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcfresults"

# object.size(vcfresults)
# 2,392,863,888 bytes

dropRareAlleles	<- vcfresults$control$snpOptions["dropRareAlleles"]
saveBiAllelic 	<- vcfresults$control$snpOptions["saveBiAllelic"]
nMaxAlt			<- vcfresults$control$nMaxAlt

gc()
            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  14201705  758.5   17756025  948.3  16864700  900.7
# Vcells 146561218 1118.2  211583004 1614.3 159676022 1218.3

library(VariantAnnotation, quietly = TRUE)

if(Glazerize){  # grab fish names
	fishnames <- samples(header(vcfresults$vcf[[1]]))
	} else {
	fishnames <- samples(header(vcfresults$vcf))
	}

groupcodes <- rep(0, length(fishnames))		 # initialize
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
cat("\ngroupcodes:\n")
print(groupcodes)
 # [1]  8  8  8  8  8  8  7  7  7  7  7  7 10 10 10 11 11 11 11  9  9  9  9  9  9
# [26]  9  4  4  4  4  4  3  3  3  3  3  2  2  2  2  2  1  1  1  1  1  4  4  4  4
# [51]  3  3  3  3  6  6  6  6  6  6  5  5  5  5  5  5 12 12 12 12 12 12 12 12  2
# [76]  2  2  2  2  2  1  1  1  1  1  1

groups <- groupcodes[groupcodes > 0]
fish <- fishnames[groupcodes > 0]
# groups
# fish

# -----
cat("\nExtracting genotypes and deleting vcfresults\n")
# Pull out genotypes 
if(Glazerize){
	genotypes <- lapply(vcfresults$vcf, function(x){geno(x)$GT[ , groupcodes > 0]})
	genotypes <- do.call("rbind", genotypes)
	pos <- vcfresults$newPos
	} else {
	genotypes <- geno(vcfresults$vcf)$GT[ , groupcodes > 0]
	pos <- start(rowData(vcfresults$vcf))
	}

# Don't need this any more
rm(vcfresults)

gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  3916383 209.2   18306168  977.7  18683826  997.9
# Vcells 86723924 661.7  319657529 2438.8 362764926 2767.7


cat("\nNumber of snp\n")
print(nrow(genotypes))
# [1] 915040


# -----
# Keep only snp with all non-missing genotypes

z <- apply(genotypes, 1, function(x){sum(is.na(x))})
genotypes <- genotypes[ z == 0, ]
pos <- pos [z == 0]

cat("\nDropped all snp with any genotypes missing\n")
cat("\nNumber of snp remaining\n")
print(nrow(genotypes))
# [1] 405478

gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  3918636 209.3   14644934  782.2  18683826  997.9
# Vcells 43513899 332.0  255726023 1951.1 362764926 2767.7


# -----
# Tabulate allele frequencies and keep only bi-allelic snps

cat("\nTabulating allele frequencies\n")
alleleFreq <- apply(genotypes, 1, function(x){table(unlist(strsplit(x, split = "/")))})
# Result is a list of one-way tables

cat("\nKeeping only snp with exactly two alleles\n")
# count number of unique alleles
# "dim" works here to count columns (ncol does not, maybe because only 1 row)

z <- sapply(alleleFreq, dim) 
# table(z) # number of alleles per snp
     # 1      2      3      4 
# 105190 279159  20209    920 

genotypes <- genotypes[ z == 2, ]
alleleFreq <- alleleFreq[ z == 2 ]
pos <- pos[ z == 2 ]

cat("\nNumber of snp remaining\n")
nBiallelicSnpComplete <- nrow(genotypes)
print(nBiallelicSnpComplete)
# [1] 279159

# Grab the name of the second allele -- then count how many such alleles each genotype has as a score

z <- sapply(alleleFreq, function(x){dimnames(x)[[1]][2]})

# head(z)
      # chrXXI:58916_C/A       chrXXI:58947_G/A       chrXXI:58965_G/T 
                   # "1"                    "1"                    "1" 
      # chrXXI:58971_T/C chrXXI:58974_CACAGAA/C       chrXXI:59053_C/T 
                   # "1"                    "1"                    "1"

# Transpose genotypes and make a data frame for PCA

cat("\nTransposing genotypes array\n")
genotypes <- as.data.frame( t(genotypes), stringsAsFactors = FALSE)
gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  6353862 339.4   14644934  782.2  18683826  997.9
# Vcells 32324264 246.7  204580818 1560.9 362764926 2767.7

cat("\nConverting genotypes to 0 1 2 scores (for number of '1' alleles)\n")
# Count the number of alleles of the type "second allele" (stored in z)
genoScore <- mapply(genotypes, z, FUN = function(x,z){
	score <- strsplit(x, split = "/")
	score <- sapply(score, function(score){sum(score == z)})
	})
# genoScore[ , 1:3]
gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  6353894 339.4   25047118 1337.7  30641133 1636.5
# Vcells 60056889 458.2  219761608 1676.7 362764926 2767.7

rm(genotypes)

cat("\nPrincipal components analysis\n")

# library(pcaMethods)
# zscore <- prep(genoScore, scale="none", center=TRUE)
# resPCA <- pca(zscore, method="svd", center=FALSE, nPcs=5)
# gc()

# sDev(pcaRes)

# -----
# First, including all fish

pcaResList <- list() # save results but not genoScore for now

z <- prcomp(genoScore, center = TRUE)

cat("\nProportion of variance accounted for by each PC\n")
pcaVarProp <- 100*(z$sdev^2)/sum(z$sdev^2)
print( round( pcaVarProp, 2) )
 # [1] 25.81 17.40  7.49  4.60  3.35  2.81  1.98  1.70  1.44  1.31  1.23  1.07
# [13]  0.97  0.94  0.89  0.84  0.75  0.72  0.68  0.66  0.62  0.58  0.56  0.54
# [25]  0.53  0.50  0.48  0.48  0.47  0.47  0.47  0.45  0.45  0.45  0.44  0.43
# [37]  0.42  0.41  0.40  0.39  0.39  0.38  0.38  0.38  0.38  0.37  0.37  0.36
# [49]  0.36  0.36  0.35  0.35  0.35  0.34  0.34  0.33  0.33  0.32  0.32  0.32
# [61]  0.32  0.32  0.31  0.31  0.30  0.30  0.30  0.30  0.29  0.28  0.28  0.28
# [73]  0.28  0.27  0.27  0.27  0.26  0.26  0.25  0.25  0.23  0.22  0.22  0.21
# [85]  0.19  0.00

pcaRes <- cbind.data.frame(fish, groupnames[groups], stringsAsFactors = FALSE)
pcaRes <- cbind.data.frame(pcaRes, z$x)

pcaResList$nBiallelicSnpComplete <- nBiallelicSnpComplete

pcaResList$pcaRes <- pcaRes
pcaResList$pcaVarProp <- pcaVarProp

# Calculate Euclidean distances between all pairs of individuals
x <- z$x
dimnames(x) <- list(fish, fish)
pcaDist <- dist(x, method = "euclidean", diag = TRUE, upper = TRUE)

pcaResList$pcaDist <- pcaDist


# PCA axis plots and phylogram based on distances
pdf(file = paste(project, chrname, "pcaRes.pdf", sep = "."))

plot(PC2 ~ PC1, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
text(PC2 ~ PC1, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

plot(PC3 ~ PC2, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
text(PC3 ~ PC2, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

plot(PC4 ~ PC3, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
text(PC4 ~ PC3, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

plot(PC5 ~ PC4, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
text(PC5 ~ PC4, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)


# NJ tree using distances
library(ape)

pcaPhylo <- nj(pcaDist)

# plot(pcaPhylo, cex = 0.3)
# plot(pcaPhylo, cex = 0.3, type = "fan")
plot(pcaPhylo, cex = 0.3, type = "unrooted", lab4ut = "axial", main = chrname)
# plot(pcaPhylo, cex = 0.3, type = "radial")

# Improved NJ (?)
# pcaPhylo <- bionj(pcaDist/10) # can't have any values greater than 100
# plot(pcaPhylo, cex = 0.3)

dev.off()

save(pcaResList, file = pcaResFile) 
# load(file = pcaResFile) # object is pcaResList
# load("BenlimAllMarine.chrXXI.pcaRes.rdd") # object is pcaResList
