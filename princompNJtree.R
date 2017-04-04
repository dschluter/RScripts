#!/usr/bin/Rscript

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

# All files will be named "*.pcaRes.rdd" and "*.pcaRes.pdf" regardless of args
# So RENAME files after running to distinguish different runs using different subsets of fish.

# This code assumes that we are working only with ONE chromosome at a time
# Note: At most 3 ALT alleles permitted per snp in this version

# Muhua's VSQR filters are: 
#	QD < 2.00 || FS > 60.000 || MQ < 50.00 || MQRankSum < -12.500 || ReadPosRankSum < -8.000")
# GATK Best Practices recommends cutoffs: 
#	QD 2, FS 60, MQ 40, MQRankSum -12.5, ReadPosRankSum -8.000

# Based on chrM, most snp have QD > 2, FS < 60, MQ between 20 and 55, 
#	MQRankSum between -4 and 4, ReadPosRankSum between -4 and 4)

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# Filter based on FILTER
# NULL, LowQual, PASS, VQSRTrancheSNP99.00to99.50, VQSRTrancheSNP99.50to99.90, VQSRTrancheSNP99.50to99.90+

library(VariantAnnotation)

project 	<- NULL
chrInclude  <- NULL
removeIndels<- FALSE # Default
filter 		<- NULL 

args <- commandArgs(TRUE)
# args <- c("project=Benlim", "chrInclude=chrM,chrXXI", "removeIndels=FALSE", "filter=PASS")

# Parses the args into a data frame with two columns (V1=left and V2=right of each "=" sign)
# and then assigns V2 to variables whose names are in V1 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}

chrname <- unlist(strsplit(chrInclude, split = ","))

vcfparamFile <- paste(project, "vcfparam.rdd", sep = ".")
load(vcfparamFile) # object is vcfparam

Glazerize <- vcfparam$Glazerize
fishnames <- vcfparam$fishnames
groupcodes<- vcfparam$groupcodes
groupnames<- vcfparam$groupnames
print(groupnames)
 # [1] "paxl"       "paxb"       "pril"       "prib"       "qryl"      
 # [6] "qryb"       "ensl"       "ensb"       "marine-pac" "marine-atl"
# [11] "marine-jap" "solitary"   "sculpin"    "stream"    

if(is.null(chrInclude)) stop("Provide chrInclude= in arguments")
if(is.null(project)) stop("Provide project= in arguments")

for(i in chrname){
	# i <- "chrM"
	pcaResFile 	<- paste(project, chrname, "pcaRes.rdd", sep = ".")

	cat("\nProject, chrname, groupnames:\n")
	cat(project, chrname, "\n", paste(groupnames, collapse = ","), "\n")

	if(Glazerize){
		vcfresultsfile <- paste(project, ".", i, ".vcfNew.rdd", sep = "")
		} else {
		vcfresultsfile <- paste(project, ".", i, ".vcfresults.rdd", sep = "")}

	cat("\nLoading vcfresults file\n")
	load(file = vcfresultsfile)   # object is "vcf"
	
	gc()
	            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells  14201705  758.5   17756025  948.3  16864700  900.7
	# Vcells 146561218 1118.2  211583004 1614.3 159676022 1218.3

	# keep only cases corresponding to FILTER
	keep <- casefold(filt(vcf)) == casefold(filter)
	# table(keep)
	# keep
	# FALSE  TRUE 
	 # 4236   104
	
	vcf <- vcf[keep]
	
	cat("\nExtracting genotypes\n")
	genotypes <- geno(vcf)$GT
	pos <- info(vcf)$newPos
	
	# Don't need this object any more
	rm(vcf)


** CONTINUE HERE **

# For chrM test how many snp left if apply all hard filters:
# z <- cbind(QD,FS,MQ,MQRankSum,ReadPosRankSum)
# nrow(z)
# [1] 1154
# nrow(na.omit(z))
# [1] 443

if(filterFS){
	genotypes <- genotypes[FS < 60 | is.na(FS) , ]
	pos <- pos[FS < 60 | is.na(FS)]
	# FS <- FS[FS < 60 | is.na(FS)]
	}


gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  3916383 209.2   18306168  977.7  18683826  997.9
# Vcells 86723924 661.7  319657529 2438.8 362764926 2767.7

cat("\nNumber of snp\n")
print(nrow(genotypes))
# [1] 1148

# -----
# Keep only snp with all non-missing genotypes for PCA

z <- apply(genotypes, 1, function(x){sum(is.na(x))})
genotypes <- genotypes[ z == 0, ]
pos <- pos [z == 0]

cat("\nDropped all snp with any genotypes missing\n")
cat("\nNumber of snp remaining\n")
print(nrow(genotypes))
# [1] 336

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
  # 1   2   3 
# 100 225  11 

genotypes <- genotypes[ z == 2, ]
alleleFreq <- alleleFreq[ z == 2 ]
pos <- pos[ z == 2 ]

cat("\nNumber of snp remaining\n")
nBiallelicSnpComplete <- nrow(genotypes)
print(nBiallelicSnpComplete)
# [1] 225

# Grab the name of the second allele -- then count how many such alleles each genotype has as a score

z <- sapply(alleleFreq, function(x){dimnames(x)[[1]][2]})

# head(z)
# chrM:397_G/A chrM:480_C/A chrM:543_C/G chrM:607_T/A chrM:620_G/A chrM:648_C/T
         # "1"          "1"          "1"          "1"          "1"          "1" 


if(removeIndels){
	cat("\nremoveIndels not implemented yet\n")
	}
	
# Transpose genotypes and make a data frame for PCA
# rows are fish, columns are genes

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

# head(genoScore[ , 1:5])
     # chrM:397_G/A chrM:480_C/A chrM:543_C/G chrM:607_T/A chrM:620_G/A
# [1,]            2            2            0            0            0
# [2,]            2            2            0            0            0
# [3,]            2            2            0            0            1
# [4,]            2            2            0            0            0
# [5,]            2            2            0            0            0
# [6,]            2            2            0            0            0

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
  # [1] 46.19 32.29  2.51  2.19  1.36  1.20  1.12  0.95  0.85  0.66  0.61  0.52
 # [13]  0.49  0.49  0.40  0.39  0.38  0.35  0.33  0.33  0.32  0.30  0.28  0.28
 # [25]  0.28  0.27  0.26  0.25  0.24  0.22  0.21  0.20  0.20  0.19  0.19  0.17
 # [37]  0.16  0.16  0.15  0.15  0.14  0.13  0.10  0.10  0.10  0.10  0.10  0.10
 # [49]  0.10  0.09  0.09  0.07  0.07  0.07  0.07  0.06  0.06  0.06  0.05  0.05
 # [61]  0.05  0.04  0.02  0.02  0.02  0.02  0.02  0.02  0.01  0.01  0.00  0.00
 # [73]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
 # [85]  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
 # [97]  0.00  0.00  0.00  0.00

pcaRes <- cbind.data.frame(fishnames, groupnames[groupcodes], stringsAsFactors = FALSE)
pcaRes <- cbind.data.frame(pcaRes, z$x)

pcaResList$nBiallelicSnpComplete <- nBiallelicSnpComplete

pcaResList$pcaRes <- pcaRes
pcaResList$pcaVarProp <- pcaVarProp

# Calculate Euclidean distances between all pairs of individuals
x <- z$x
dimnames(x)[[1]] <- fishnames
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
