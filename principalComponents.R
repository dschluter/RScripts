#!/usr/bin/Rscript

# Makes vcfresults, containing all snp results

# Run in Unix as " Rscript readSaveDnpdByGroup.R ... "

# ALT alleles of "<*:DEL>" are set to variant type NA

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

# This code assumes that we are working only with ONE chromosome at a time
# Note: At most 3 ALT alleles permitted per snp in this version
# NOTE: windowmaskerSdust start position is 0-based, so added 1 (end position is 1-based!)

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn", "paxl", "paxb", "marine-pac")
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb", "marine-pac")
# args <- c("BenlimAllMarine", "chrXXI", "paxl","paxb","pril","prib","qryl","qryb","ensl","ensb","marine-pac","marine-atl","marine-jap","solitary")

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)]

cat("\nProject, chrname, groupnames:\n")
cat(project, chrname, groupnames,"\n")

dropRareAlleles	<- FALSE
saveBiAllelic 	<- FALSE # saves a second data set having exactly 2 snp per marker (not necessarily the REF), no indels
plotQualMetrics <- FALSE

# convert snp to Glazer assembly coordinates **** Not done yet here -- I'm using old chr ****
Glazerize 		<- TRUE # Requires file "glazerFileS4 NewScaffoldOrder.csv" in current directory

GTminFrac 		<- 2/3

# load "chrvec" for the current chromosome
chrno 				<- gsub("^chr", "", chrname)
chrmaskfile         <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd

fastaname		<- paste(chrname, "fa", sep = ".")
vcfname			<- paste(project, ".", chrname, ".var.vcf", sep="")
pcaResFile		<- paste(project, chrname, "pcaRes.rdd", sep = ".")

GTmissing		<- "."	# how GATK represents missing genotypes in the vcf file "./."
nMaxAlt			<- 3	# maximum number of ALT alleles

control <- list()
control$nMaxAlt <- nMaxAlt
control$snpOptions <- c(dropRareAlleles = dropRareAlleles, saveBiAllelic = saveBiAllelic)
control$GTminFrac <- GTminFrac

cat("\nControl settings on this run\n")
print(control)

library(VariantAnnotation, quietly = TRUE)

load(chrmaskfile) 	# object is named "chrvec"
which.chrvec.M <- which(chrvec == "M")

# object.size(chrvec)
#  93,740,176 bytes # chrXXI
# 500,401,968 bytes # chrUn

rm(chrvec)

# Garbage collection -- current memory consumption
# gcinfo(TRUE) # sends messages about memory garbage collection
gc()
          # used  (Mb) gc trigger  (Mb) max used  (Mb)
# Ncells 2998785 160.2    4418719 236.0  3493455 186.6
# Vcells 3675701  28.1   32139972 245.3 38033273 290.2

# ---------------------
# Alternatives to reading whole vcf file

# 1. Make tabix and use data ranges
# compressedVcfname <- paste(project, chrname, "var.vcf.bgz", sep = ".")
# zipped <- bgzip(vcfname, compressedVcfname)
# idx <- indexTabix(zipped, "vcf")
# vcfTabix <- TabixFile(zipped, idx)
# vcf <- readVcf(vcfTabix, chrname, IRanges(c(1, 100000))) # not working yet

# 2. Read a subset of fields
# vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("ALT", "QUAL"), geno="GT", info=NA))

# ---------------------
# Read variant VCF file

vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("ALT", "QUAL"), geno="GT", info=NA))

# object.size(vcf)
#  50,514,280 bytes # chrXXI
# 210,524,096 bytes # chrUn

cat("\nSuccessfully read vcf file\n")

# -----
# fish group codes 1, 2, ... Order is determined by sequence of groupnames, eg "paxl", "paxb", "marine-pac"
fishnames <- samples(header(vcf))
cat("\nfishnames:\n")
print(fishnames)

groupcodes <- vector(length = length(fishnames), mode = "integer")
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
cat("\ngroupcodes:\n")
print(groupcodes)  # 
 # [1]  8  8  8  8  8  8  7  7  7  7  7  7 10 10 10 11 11 11 11  9  9  9  9  9  9
# [26]  9  4  4  4  4  4  3  3  3  3  3  2  2  2  2  2  1  1  1  1  1  4  4  4  4
# [51]  3  3  3  3  6  6  6  6  6  6  5  5  5  5  5  5 12 12 12 12 12 12 12 12  2
# [76]  2  2  2  2  2  1  1  1  1  1  1
 

# ------
# Drop masked SNP (variants whose start position is masked in the chromosome)
# ** applied to CHRIV, this step removes base "12811481", the putative site of the Eda mutation **
# 	In windowsmaskerSdust.chrIV, the site is marked as masked.

keep <- !(start(ranges(vcf)) %in% which.chrvec.M) # i.e., includes only the good bases: only upper case, no "M"
vcf <- vcf[keep]
cat("\nCompleted removal of snp corresponding to masked bases\n\n")

rm(keep)
rm(which.chrvec.M) # remove later

# object.size(vcf)
# 36,360,752 bytes # chrXXI
gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  3695203 197.4    7700734  411.3   6910418  369.1
# Vcells 66303599 505.9  177407159 1353.6 177175701 1351.8

cat("\nNumber of snp remaining\n")
print(nrow(geno(vcf)$GT))
# [1] 648522

# ------
# set missing genotypes to NA instead of "."
geno(vcf)$GT[geno(vcf)$GT == GTmissing] <- NA

# -----
# Keep only snp with all non-missing genotypes
# "genotype" is all character data

genotypes <- na.omit( geno(vcf)$GT )
rm(vcf)
gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  3673555 196.2    7700734  411.3   6910418  369.1
# Vcells 33818067 258.1  230920621 1761.8 261518340 1995.3

cat("\nDropped all snp with any genotypes missing\n")
cat("\nNumber of snp remaining\n")
print(nrow(genotypes))
# [1] 310935

# -----
# Drop snp containing rare alleles if dropRareAlleles is TRUE
#  -- we are using ** RULE2 ** otherwise we get lots of rare alleles when many genotypes are missing
# RULE1: drop alleles less than 5% total, measured as percent of 2*sum(nInd), the maximum number of alleles possible
# RULE2: drop alleles less than 5% total, measured as percent of the total number of alleles in the sample.

if(dropRareAlleles){
	cat("\nSorry: drop snp with rare alleles not implemented yet\n")
	}


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
 # 80350 215258  14688    639

# Keep those with exactly 2 alleles

genotypes <- genotypes[ z == 2, ]
alleleFreq <- alleleFreq[ z == 2 ]

cat("\nNumber of snp remaining\n")
print(nrow(genotypes))
# [1] 215258

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
           # used  (Mb) gc trigger  (Mb)  max used   (Mb)
# Ncells  5608538 299.6    9040660 482.9   9040660  482.9
# Vcells 25609286 195.4   94585084 721.7 261518340 1995.3

cat("\nConverting genotypes to 0 1 2 scores (for number of '1' alleles)\n")
# Count the number of alleles of the type "second allele" (stored in z)
genoScore <- mapply(genotypes, z, FUN = function(x,z){
	score <- strsplit(x, split = "/")
	score <- sapply(score, function(score){sum(score == z)})
	})
# genoScore[ , 1:3]
gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  5608543 299.6   20248149 1081.4  24336408 1299.8
# Vcells 50594068 386.1  175086369 1335.9 261518340 1995.3

rm(genotypes)

cat("\nPrincipal components analysis\n")

# library(pcaMethods)
# zscore <- prep(genoScore, scale="none", center=TRUE)
# resPCA <- pca(zscore, method="svd", center=FALSE, nPcs=5)
# gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  5497641 293.7   16198519  865.1  24336408 1299.8
# Vcells 79760177 608.6  175086369 1335.9 261518340 1995.3

# sDev(pcaRes)

# -----
# First, including all fish

pcaResList <- list() # save results but not genoScore for now

z <- prcomp(genoScore, center = TRUE)

cat("\nProportion of variance accounted for by each PC\n")
pcaVarProp <- 100*(z$sdev^2)/sum(z$sdev^2)
print( round( pcaVarProp, 2) )
 # [1] 27.24 17.95  7.38  4.27  3.01  2.83  1.86  1.62  1.42  1.26  1.20  1.04
# [13]  0.94  0.90  0.84  0.81  0.74  0.72  0.67  0.65  0.59  0.56  0.54  0.53
# [25]  0.51  0.49  0.48  0.48  0.47  0.46  0.46  0.46  0.45  0.43  0.43  0.42
# [37]  0.41  0.40  0.39  0.39  0.38  0.38  0.37  0.36  0.36  0.35  0.35  0.34
# [49]  0.34  0.34  0.34  0.34  0.33  0.33  0.32  0.32  0.32  0.31  0.31  0.31
# [61]  0.31  0.30  0.30  0.29  0.29  0.29  0.28  0.28  0.28  0.28  0.27  0.27
# [73]  0.27  0.26  0.26  0.25  0.25  0.25  0.24  0.24  0.22  0.21  0.21  0.21
# [85]  0.19  0.00

pcaRes <- cbind.data.frame(fishnames, groupnames[groupcodes], stringsAsFactors = FALSE)
pcaRes <- cbind.data.frame(pcaRes, z$x)

pcaResList$pcaRes <- pcaRes
pcaResList$pcaVarProp <- pcaVarProp

# PCA axis plots and phylogram based on distances
pdf(file = paste(project, chrname, "pcaRes.pdf", sep = "."))

plot(PC2 ~ PC1, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
text(PC2 ~ PC1, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

plot(PC3 ~ PC2, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
text(PC3 ~ PC2, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

plot(PC4 ~ PC3, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
text(PC4 ~ PC3, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

plot(PC5 ~ PC4, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
text(PC5 ~ PC4, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)

# Calculate the distances between all pairs of individuals
# Euclidean
x <- z$x
dimnames(x) <- list(fishnames, fishnames)
pcaDist <- dist(x, method = "euclidean", diag = TRUE, upper = TRUE)

pcaResList$pcaDist <- pcaDist

# NJ tree using distances
library(ape)

pcaPhylo <- nj(pcaDist)

# plot(pcaPhylo, cex = 0.3)
# plot(pcaPhylo, cex = 0.3, type = "fan")
plot(pcaPhylo, cex = 0.3, type = "unrooted", lab4ut = "axial")
# plot(pcaPhylo, cex = 0.3, type = "radial")

# dev.off()

# Improved NJ (?)
# pcaPhylo <- bionj(pcaDist/10) # can't have any values greater than 100
# plot(pcaPhylo, cex = 0.3)
# dev.off()


# -----
# redo without Sea of Japan

if( "marine-jap" %in% groupnames ){
	groups <- groupnames[groupcodes]
	genoScore2 <- genoScore[ groups != "marine-jap", ]

	z <- prcomp(genoScore2, center = TRUE)
	cat("\nProportion of variance accounted for by each PC\n")
	pcaVarProp2 <- 100*(z$sdev^2)/sum(z$sdev^2)
	print( round( pcaVarProp2, 2) )
	 # [1] 29.54  9.82  5.63  3.96  3.73  2.45  2.13  1.87  1.67  1.58  1.37  1.24
	# [13]  1.18  1.11  1.06  0.98  0.95  0.88  0.85  0.78  0.74  0.72  0.70  0.67
	# [25]  0.65  0.63  0.61  0.60  0.60  0.57  0.56  0.55  0.54  0.52  0.52  0.51
	# [37]  0.50  0.50  0.49  0.48  0.47  0.47  0.46  0.45  0.45  0.45  0.45  0.44
	# [49]  0.44  0.43  0.43  0.42  0.42  0.41  0.41  0.41  0.40  0.40  0.39  0.38
	# [61]  0.38  0.38  0.37  0.37  0.37  0.36  0.36  0.36  0.35  0.35  0.34  0.34
	# [73]  0.33  0.32  0.32  0.31  0.29  0.28  0.27  0.27  0.24  0.00
	
	pcaRes2 <- cbind.data.frame(fishnames = fishnames[groups != "marine-jap"], 
				groupnames = groups[groups != "marine-jap"], stringsAsFactors = FALSE)
	pcaRes2 <- cbind.data.frame(pcaRes2, z$x)
	
	pcaResList$pcaRes2 <- pcaRes2
	pcaResList$pcaVarProp2 <- pcaVarProp2

	# Plots
	plot(PC2 ~ PC1, data = pcaRes2, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
	text(PC2 ~ PC1, data = pcaRes2, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)
	
	plot(PC3 ~ PC2, data = pcaRes2, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
	text(PC3 ~ PC2, data = pcaRes2, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)
	
	plot(PC4 ~ PC3, data = pcaRes2, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
	text(PC4 ~ PC3, data = pcaRes2, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)
	
	plot(PC5 ~ PC4, data = pcaRes2, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)) )
	text(PC5 ~ PC4, data = pcaRes2, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)


	# Calculate the distances between all pairs of individuals
	# Euclidean
	x <- z$x
	dimnames(x) <- list(fishnames[groups != "marine-jap"], fishnames[groups != "marine-jap"])
	pcaDist2 <- dist(x, method = "euclidean", diag = TRUE, upper = TRUE)
	
	pcaResList$pcaDist2 <- pcaDist2

	# NJ tree using distances
	library(ape)
	
	pcaPhylo <- nj(pcaDist2)
	
	# plot(pcaPhylo, cex = 0.3)
	# plot(pcaPhylo, cex = 0.3, type = "fan")
	plot(pcaPhylo, cex = 0.3, type = "unrooted", lab4ut = "axial")
	# plot(pcaPhylo, cex = 0.3, type = "radial"])
	
	# dev.off()
	
	# Improved NJ (?)
	# pcaPhylo <- bionj(pcaDist/10) # can't have any values greater than 100
	# plot(pcaPhylo, cex = 0.3)
	# dev.off()

	}

dev.off()

save(pcaResList, file = pcaResFile) 
# load(file = pcaResFile) # object is pcaResList
# load("BenlimAllMarine.chrXXI.pcaRes.rdd") # object is pcaResList
