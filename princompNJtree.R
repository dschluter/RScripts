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

project 		<- NULL
chrInclude  	<- NULL
removeIndels	<- FALSE
filter 			<- NULL 
includeJapSea 	<- FALSE # deletes from PCA and NJ trees

args <- commandArgs(TRUE)
# args <- c("project=Benlim", "chrInclude=chrM,chrXXI", "removeIndels=FALSE", "filter=PASS,.", "includeJapSea=FALSE")

# Parses the args into a data frame with two columns (V1=left and V2=right of each "=" sign)
# and then assigns V2 to variables whose names are in V1 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}

chrname <- unlist(strsplit(chrInclude, split = ","))
chrname <- g$chrOrderRoman(chrname)

filter	<- unlist(strsplit(filter, split=","))

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

if(includeJapSea == "FALSE"){
	noJapSea <- !grepl("JapSea", fishnames, ignore.case = TRUE)
	fishnames <- fishnames[noJapSea]
	groupcodes <- groupcodes[noJapSea]
	}

if(is.null(chrInclude)) stop("Provide chrInclude= in arguments")
if(is.null(project)) stop("Provide project= in arguments")

cat("\nProject, chrname, groupnames:\n")
cat(project, "\n", chrname, "\n", paste(groupnames, collapse = ","), "\n")

for(i in chrname){
	# i <- "chrM"
	pcaResFile 	<- paste(project, i, "pcaRes.rdd", sep = ".")

	if(Glazerize){
		vcfresultsfile <- paste(project, ".", i, ".vcfNew.rdd", sep = "")
		} else {
		vcfresultsfile <- paste(project, ".", i, ".vcfresults.rdd", sep = "")}

	cat("\nLoading vcfresults file\n")
	load(file = vcfresultsfile)   # object is "vcf"
	
	gc()

	# keep only cases corresponding to filter
	if(all(!is.na(filter) & filter != "NULL")){
		keep <- casefold(filt(vcf)) %in% casefold(filter)
		# table(keep)
		# FALSE  TRUE 
		 # 3847   493
		 
		vcf <- vcf[keep]
		rm(keep)
 		}
	
	cat("\nExtracting genotypes\n")
	genotypes <- geno(vcf)$GT
	pos <- info(vcf)$newPos
	
	# Don't need this object any more
	rm(vcf)

	cat("\nNumber of snp\n")
	print(nrow(genotypes))
	# [1] 493

	# -----
	# Keep only snp with all non-missing genotypes for PCA - no interpolation yet
	
	z <- apply(genotypes, 1, function(x){sum(is.na(x))})
	genotypes <- genotypes[ z == 0, ]
	pos <- pos [z == 0]
	
	cat("\nDropped all snp with any genotypes missing\n")
	cat("\nNumber of snp remaining\n")
	print(nrow(genotypes))
	# [1] 282


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
	# 274   6   2 

	genotypes <- genotypes[ z == 2, ]
	alleleFreq <- alleleFreq[ z == 2 ]
	pos <- pos[ z == 2 ]

	cat("\nNumber of snp remaining\n")
	nBiallelicSnpComplete <- nrow(genotypes)
	print(nBiallelicSnpComplete)
	# [1] 6

	# Grab the name of the second allele -- then count how many such alleles each genotype has as a score

	z <- sapply(alleleFreq, function(x){dimnames(x)[[1]][2]})
	
	# head(z)
	       # chrM:3106_ACCCTTG/A chrM:4449_T/TCAGGAAAGTACA           chrM:9312_G/GTA 
	                      # "1"                       "1"                       "1" 
	          # chrM:9314_ATG/A   chrM:11426_A/ATTGCAAGCC          chrM:14008_CTG/C 
	                      # "1"                       "1"                       "1" 


	if(removeIndels){
		cat("\nremoveIndels not implemented yet\n")
		}
	
	# Transpose genotypes and make a data frame for PCA
	# rows are fish, columns are genes
	
	cat("\nTransposing genotypes array\n")
	genotypes <- as.data.frame( t(genotypes), stringsAsFactors = FALSE)

	cat("\nConverting genotypes to 0 1 2 scores (for number of '1' alleles)\n")
	# Count the number of alleles of the type "second allele" (stored in z)
	genoScore <- mapply(genotypes, z, FUN = function(x,z){
		score <- strsplit(x, split = "/")
		score <- sapply(score, function(score){sum(score == z)})
		})
		

	# head(genoScore[ , 1:5])
	       # chrM:3106_ACCCTTG/A chrM:4449_T/TCAGGAAAGTACA chrM:9312_G/GTA
	# ENSB01                   0                         0               0
	# ENSB03                   0                         0               0
	# ENSB08                   0                         0               0
	# ENSB12                   0                         0               0
	# ENSB15                   0                         0               0
	# ENSB23                   0                         0               0
	       # chrM:9314_ATG/A chrM:11426_A/ATTGCAAGCC
	# ENSB01               0                       0
	# ENSB03               0                       0
	# ENSB08               0                       0
	# ENSB12               0                       0
	# ENSB15               0                       0
	# ENSB23               0                       0

	rm(genotypes)
	
	if(includeJapSea == "FALSE"){
		genoScore <- genoScore[noJapSea, ]
		}


	cat("\nPrincipal components analysis\n")

	# This is how you would do pcaMethods
	# library(pcaMethods)
	# zscore <- prep(genoScore, scale="none", center=TRUE)
	# resPCA <- pca(zscore, method="svd", center=FALSE, nPcs=5)
	# sDev(resPCA)
	
	pcaResList <- list() # save results but not genoScore for now
	
	z <- prcomp(genoScore, center = TRUE)
	
	cat("\nProportion of variance accounted for by each PC\n")
	pcaVarProp <- 100*(z$sdev^2)/sum(z$sdev^2)
	print( round( pcaVarProp, 2) )
	# [1] 61.41 32.43  2.07  2.07  2.01  0.00
	
	pcaRes <- cbind.data.frame(fishnames, groupnames=groupnames[groupcodes], stringsAsFactors = FALSE)
	pcaRes <- cbind.data.frame(pcaRes, z$x)
	
	pcaResList$nBiallelicSnpComplete <- nBiallelicSnpComplete
	
	pcaResList$pcaRes <- pcaRes
	pcaResList$pcaVarProp <- pcaVarProp
	
	# Calculate Euclidean distances between all pairs of individuals
	x <- z$x
	dimnames(x)[[1]] <- fishnames
	pcaDist <- dist(x, method = "euclidean", diag = TRUE, upper = TRUE)
	
	pcaResList$pcaDist <- pcaDist

	# # PCA axis plots and phylogram based on distances
	# pdf(file = paste(project, chrname, "pcaRes.pdf", sep = "."))
	
	# plot(PC2 ~ PC1, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
	# text(PC2 ~ PC1, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)
	
	# plot(PC3 ~ PC2, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
	# text(PC3 ~ PC2, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)
	
	# plot(PC4 ~ PC3, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
	# text(PC4 ~ PC3, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)
	
	# plot(PC5 ~ PC4, data = pcaRes, col = as.numeric(factor(groupnames)), pch = as.numeric(factor(groupnames)), main = chrname )
	# text(PC5 ~ PC4, data = pcaRes, labels = groupnames, pos = 4, offset = 0.5, cex = 0.7, xpd = NA)


	# NJ tree using distances
	library(ape)
	
	pcaPhylo <- nj(pcaDist)
	
	# # plot(pcaPhylo, cex = 0.3)
	# # plot(pcaPhylo, cex = 0.3, type = "fan")
	# plot(pcaPhylo, cex = 0.3, type = "unrooted", lab4ut = "axial", main = chrname)
	# # plot(pcaPhylo, cex = 0.3, type = "radial")
	
	# # Improved NJ (?)
	# # pcaPhylo <- bionj(pcaDist/10) # can't have any values greater than 100
	# # plot(pcaPhylo, cex = 0.3)

	# dev.off()
	
	pcaResList$pcaPhylo <- pcaPhylo

	save(pcaResList, file = pcaResFile) 
	# load(file = pcaResFile) # object is pcaResList
	# load("BenlimAllMarine.chrXXI.pcaRes.rdd") # object is pcaResList
	}
	
# Combine the results into a list of lists - read each file one by one and save in a single object
pcaList <- list()
for(i in chrname){
	# i <- "chrM"
	pcaResFile 	<- paste(project, i, "pcaRes.rdd", sep = ".")
	load(pcaResFile) # object is pcaResList
	pcaList[[i]] <- pcaResList
	file.remove(pcaResFile)
	}
	
# names(pcaList[[1]])
# [1] "nBiallelicSnpComplete" "pcaRes"                "pcaVarProp"           
# [4] "pcaDist"               "pcaPhylo"             

save(pcaList, file = paste(project, "pcaList.rdd", sep = "."))

# Create the combined distance matrix from all chromosomes
n <- nrow(as.matrix(pcaList[[1]]$pcaDist))
# [1] 105
pcaDistTot <- matrix(0, nrow=n, ncol=n)
for(i in 1:length(pcaList)){
	pcaDistTot <- pcaDistTot + as.matrix(pcaList[[i]]$pcaDist)^2
	}
pcaDistTot <- sqrt(pcaDistTot) # is a matrix
pcaNJtot <- nj(pcaDistTot) # neighbor joining again
save(pcaDistTot, file = paste(project, "pcaDist.rdd", sep = "."))
save(pcaNJtot, file = paste(project, "pcaNJtot.rdd", sep = "."))

pdf(file = paste(project, "pcaNJtot.pdf", sep = "."))
	plot(pcaNJtot, cex = 0.3, type = "unrooted", lab4ut = "axial", main = "all chr combined")
dev.off()

# Multidimensional scaling of the total distance matrix
# -------------------
x <- as.data.frame( cmdscale(pcaDistTot, k = 4) )
colnames(x) <- c("mds1","mds2", "mds3","mds4")
pdf(file = paste(project, "pcaCmdscale.pdf", sep = "."))
	plot(mds2 ~ mds1, data = x, main = "cmdscale of total euclidean distance matrix", pch = groupcodes)
	text(mds2 ~ mds1, data = x, labels = rownames(x), pos = 4, offset = 0.1, cex = 0.5, xpd = NA)
	plot(mds4 ~ mds3, data = x, main = "cmdscale of total euclidean distance matrix", pch = groupcodes)
	text(mds4 ~ mds3, data = x, labels = rownames(x), pos = 4, offset = 0.1, cex = 0.5, xpd = NA)
dev.off()

