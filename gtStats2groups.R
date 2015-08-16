#!/usr/bin/Rscript

# -------------------------------------------------------------------
# Analysis of differences between two groups (e.g., benthic-limnetic)

# This code assumes that you are working only with ONE chromosome
# Note: At most 3 ALT alleles permitted per snp in this version

# Code analyzes all snps, including indels. 
# Need to add step if want true snps only

# Q: How did Jones et al handle indels in the 21 genomes project?
# Answer: they did not call indels

# Obtains stats on one pair of populations
# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb")

includePfisher 	<- TRUE
includeFst     	<- TRUE
trueSnpOnly		<- FALSE
# dropRareAlleles <- FALSE # Set in "saveSnpsByGroup.R"

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)]
vcfname			<- paste(project, ".", chrname, ".var.vcf", sep="")
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

gtstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "-"), "rdd", sep = ".")
gtstats <- list()

if(length(groupnames) > 2 ) stop("Provide names of only two groups")

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcfresults"
# lapply(vcfresults, object.size)

# names(vcfresults)
# [1] "groupcodes"        "vcf"               "altUsedList"       "nAltUsed"          "snpTypeList"      
# [6] "alleleFreqByGroup"

****
*** Need to redo altUsedList, snpTypeList, for this pair of populations - make a function
*** If use only true snp, redo alleleFreqByGroup also
****


library(VariantAnnotation)
# library(GenomicFeatures)

# Group codes for thw two groups of interest here
# 1 and 2 will indicate the two groups of interest; all other fish are assigned a code of 0
fishnames <- samples(header(vcfresults$vcf)) # pulls fish names
groupcodes <- rep(0, length(fishnames))		 # initialize
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
cat("\ngroupcodes:\n")
print(groupcodes)
# [1] 0 0 0 0 0 0 0 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

# Grab genotypes of the two groups of interest only
groups <- groupcodes[groupcodes > 0]
# groups
 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

# ----------
# Pull out the genotypes and allele frequencies for just the two groups being analyzed

genotypes <- geno(vcfresults$vcf)$GT[ , groupcodes > 0] 
alleleFreqByGroup <- lapply(vcfresults$alleleFreqByGroup, function(x){x[groupnames,]})
# nrow(genotypes)
# [1] 281607

# This will save memory and not needed any more
vcfresults$alleleFreqByGroup <- NULL

# ----------
# Drop loci that have too few individuals
GTminFrac <- 2/3
nInd <- as.vector(table(groups)) 	# n individuals genotyped in each group eg 11 11  7
nMin <- floor(GTminFrac*nInd) 		# minimum number required in each group eg 7 7 4
nMin <- mapply(rep(5, length(groupnames)), nMin, FUN = min) # criterion is 5 or nMin, whichever is smaller, eg 5 5 4

# Identify and keep those loci with enough alleles in both groups
# The minimum number of alleles is nMin*2
test.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
	z <- rowSums(x, na.rm = TRUE)
	z[1] >= nMin[1] & z[2] >= nMin[2]
	})

genotypes <- genotypes[test.alleles.per.group, ]
alleleFreqByGroup <- alleleFreqByGroup[test.alleles.per.group]
altUsedList <- vcfresults$altUsedList[test.alleles.per.group] # not necessarily used by these two groups
snpTypeList <- vcfresults$snpTypeList[test.alleles.per.group] # not necessarily used by these two groups
# bring over coding annotations too

rm(test.alleles.per.group)
rm(vcfresults) # if we have extracted all the useful bits

# nrow(genotypes)
# [1] 274912

# -----------------
# Drop everything but true snp if TRUE
if(trueSnpOnly){
	tGT <- as.data.frame(t(genotypes), stringsAsFactors = FALSE)
	which.alt.not.snp <- lapply(snpTypeList, function(x){
		which(x != "snp")
		})
	# which.alt.not.snp["chrXXI:75476_ACAACT/AT"]
	# tGT[, "chrXXI:75476_ACAACT/AT"]
	z <- mapply(tGT, which.alt.not.snp, FUN = function(x, i){
		# x <- tGT[,"chrXXI:18599_AT/A"]; i <- 1
		if(sum(i) > 0) x[grep(i, x)] <- NA
		})
	}  

# ------------------
# Genotype frequencies
# if(getGenotypeFreq){
	# # Generate all possible genotypes for multi-allelic SNP 
	# # Allele 0 always refers to the REF allele, 1-3 to ALT alleles
	# z1 <- outer(c(0,1),c(0,1),paste, sep="/")
	# z2 <- outer(c(0,1,2),c(0,1,2),paste, sep="/")
	# z3 <- outer(c(0,1,2,3),c(0,1,2,3),paste, sep="/")
	# gtypes1 <- z1[upper.tri(z1,diag=TRUE)] # "0/0" "0/1" "1/1"
	# gtypes2 <- z2[upper.tri(z2,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2"
	# gtypes3 <- z3[upper.tri(z3,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2" "0/3" "1/3" "2/3" "3/3" 
	
	# # Genotype contingency table for every SNP: genotype by population
	# # Inside, change gtypes to a factor so that the same columns are present in every table
	
	# cat("\nCalculating genotype frequencies at every snp (maximum of 3 ALT alleles assumed)\n")
	# genotypeFreqTable <- lapply(geno,function(x){
					# z <- table(groups, factor(x, levels = gtypes3))
					# }) 
					
	# # head(genotypeFreqTable) # 0 is REF allele
	# # $`1`   
	# # pop 0/0 0/1 1/1 0/2 1/2 2/2 0/3 1/3 2/3 3/3
	  # # 1   0   1   1   0   0   0   0   0   0   0
	  # # 2   0   0   0   0   0   0   0   0   0   0
	# # $`2`
	# gtstats$genotypeFreqTable <- genotypeFreqTable
	# rm(genotypeFreqTable)
	# }
	

# Start saving key results
gtstats$groupnames <- groupnames 	# "paxl" "paxb"
gtstats$groups <- groups 			# 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
gtstats$genotypes <- genotypes	# columns are loci
gtstats$alleleFreqByGroup <- alleleFreqByGroup
gtstats$altUsedList <- altUsedList # needed?
gtstats$snpTypeList <- snpTypeList # needed?

# -------------------------------------------------------------
# Group differences in allele frequencies using Fisher exact test
# Note that one can use logistic regression to get equivalent of odds ratio when there are more than 2x2 categories.

if(includePfisher){ 
	# Fisher exact tests of allele and genotype frequencies using ALLELE frequencies
	
	cat("\nCalculating Fisher exact test p-values for group allele freq differences (maximum of 3 ALT alleles assumed)\n")
	# Loci that are not polymorphic will have a pfisher = 1
	pfisher <- lapply(alleleFreqByGroup, function(x){
		pfisher <- fisher.test(x)$p.value
		})
	
	# Write to a file in bedGraph format
	# g$write.bedGraph(paste(chrname,".pfisher.bedGraph", sep=""), chrinfo = chrname, vcfresults$POS, vcfresults$POS+1, round(-log10(pfisher),3),
		# graphname = "-log10Pfisher", col = "41,155,241", altcol = "100,100,100", ylinemark = 4,
		# description = "-log10Pfisher") 

	# Plot where the significant SNPs reside
	
	# pdf(file = paste(vcfdir, project, ".", chrname, ".pfisherPlot", ".pdf", sep=""))	
	# x <- start(vcfresults$rowdata)/(10^6)
	# y <- (-log10(vcfresults$pfisher) + rnorm(length(x), mean=0, sd=.005))
	# colour <- (round(-log10(vcfresults$pfisher)/1.2)+1)
	# plot( y ~ x, pch = ".", ylab="-log10 Pfisher", col = colour, xlab = chrname,
		# main = "Fisher exact test group allele freq differences (max 3 ALT alleles)") 
	# dev.off()
	
	# cat("\nTransition-transversion ratio for pfisher < 0.01 (uses first ALT allele if more than one)\n")
	# #snptype <- sapply(vcfresults$snptypeList, function(x){x[1]}) # first ALT allele: snp ins or del?

	# snp <- unlist(vcfresults$snpTypeList[nAltUsed > 0])
	# ref <- rep(as.vector(ref(vcfresults$vcf)), nAltUsed)
	# alt <- unlist(vcfresults$altUsedList)

	# snp <- snp[!is.na(snp) & snp == "snp"]
	# ref <- ref[!is.na(snp) & snp == "snp"]
	# alt <- alt[!is.na(snp) & snp == "snp"]
	# ref <- substr(ref, 1, 1) # keep the first letter of REF
	# alt <- substr(alt, 1, 1) # keep the first letter of ALT
	# print(g$tstv(ref[ref != alt], alt[ref != alt]))

	# REF <- as.character(vcfresults$rowdata$REF)
	# ALT <- sapply(vcfresults$ALTlist, function(x){x[1]})
	# print( g$tstv(REF[snptype == "snp" & vcfresults$pfisher <= 0.01], ALT[snptype == "snp" & vcfresults$pfisher <= 0.01]) )
	# $R
	# [1] 1.258801
	# $tstv
	     # alt
	# ref     pur   pyr
	  # pur 60301 47947
	  # pyr 47864 60306
	# $tot
	# [1] 216418
	}


# -----------------
# Prepare genotypes for Fst

# Subset genotypes to include only polymorphic sites
is.polymorphic <- sapply(alleleFreqByGroup, function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) >= 2
	})
	
genotypes.poly <- genotypes[is.polymorphic, ]

# Check that locus names are unique (for later matching)
# length(rownames(genotypes))
# [1] 274912
# length(unique(rownames(genotypes)))
# [1] 274912

# Transpose genotype array - lapply on SNPs easier this way
geno <- as.data.frame(t(genotypes.poly), stringsAsFactors = FALSE) 
names(geno) <- rownames(genotypes.poly)

#geno[1:22,1:5]

# table(unlist(geno), useNA = "always")
    # 0/0     0/1     0/2     0/3     1/1     1/2     1/3     2/2     2/3     3/3 
# 2118477  694045   16864     475  720186    4791     140   19681     112     673 
   # <NA> 
 # 100294 

gtstats$is.polymorphic <- is.polymorphic

# --------------

# library(Geneland)
# # Try Geneland method - this method below assumes that all genotypes are single-digit
# # Note: It was very fast but only gave a total Fis and Fst, not variance components per locus.
# #	so not as useful as hierfstat
# It also required the installation of a gigantic code development app

# genotypes <- as.data.frame(t(vcfresults$GT), stringsAsFactors = FALSE)
# GTsplit <- lapply(genotypes, function(x){
	# vcfresults <- substr(x, 1, 1)
	# x2 <- substr(x, 3, 3)
	# cbind(vcfresults, x2)
	# })
# GTsplit <- do.call("cbind", GTsplit)
# z <- Fstat(GTsplit, 2, groups)
# z
# $Fis
# [1] 0.02649774 0.05939996

# $Fst
          # [,1]      [,2]
# [1,] 0.0000000 0.4434629
# [2,] 0.4434629 0.0000000
	
# --------------

if(includeFst){

	library(hierfstat)
	
	cat("\nCalculating Fst (slowly) between groups\n")
	
	# Per-locus Fst using Weir and Cockrham estimates
	# Returns NA if there isn't at least 8 alleles for each pop
	# Fst is calculated as tsiga/sum(c(tsiga, tsigb, tsigw)), where
	# tsiga <- sum(lsiga)
	# tsigb <- sum(lsigb)
	# tsigw <- sum(lsigw)
	# lsiga, lsigb, and lsigw are the variance components: for each locus among populations,
	#   among individuals within populations, and within individuals, respectively
	# wc command does not like duplicate row names

	# Uses alleleFreqTable as calculated earlier
	# A snp needs at least 8 alleles for each pop to be analyzed
	# A snp needs >= 2 different alleles across all individuals (ie ALT snp not fixed in L and B) or wc will protest

	# groups
	# [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	
	# pop <- as.character(groups)

	# Convert genotypes to hierfstat format:
	# Columns are genotypes indicated as 11 (0/0), 12 (0/1), 22 (1/1), or NA (NA)

	# This next operation is redundant (did it already above) but can't use the duplicate names
	# Redoing it here so that names remain as "V1" "V2" etc

	# genotypes <- as.data.frame( t(vcfresults$GT), stringsAsFactors=FALSE) 

	#genotypes[,1:25]
	
	# dim(geno)
	# [1]     22 167079

	temp <- unlist(geno) # geno is the transposed genotypes.poly
	
	temp[temp == "0/0"] <- "11"
	temp[temp == "0/1"] <- "12"
	temp[temp == "1/1"] <- "22"
	temp[temp == "0/2"] <- "13"
	temp[temp == "1/2"] <- "23"
	temp[temp == "2/2"] <- "33"
	temp[temp == "0/3"] <- "14"
	temp[temp == "1/3"] <- "24"
	temp[temp == "2/3"] <- "34"
	temp[temp == "3/3"] <- "44"

	temp <- as.data.frame(matrix(as.integer(temp), nrow = length(pop)))
	#names(temp) <- names(geno)

	# dim(temp)
	# [1]     22 167079

	#head(names(temp)) 
	#[1] "V1" "V2" "V3" "V4" "V5" "V6"
	#tail(names(temp))
	#[1] "V167074" "V167075" "V167076" "V167077" "V167078" "V167079"

	#table(pop)
	# pop
	 # 1  2 
	# 11 11 
	# Must convert pop to a number -- using factor(pop) means that pops will be in alphabetical order
	# xpop <- as.numeric(factor(pop))

	# Calculate Fst - # took about 40 min
	# z <- wc(cbind(groups, temp))

	# Use this slightly faster version - took about 20 minutes
	z <- g$wc.revised(cbind(groups,temp))

	cat("\nFst\n")
	print(z$FST)
	#[1] 0.4020001
	# head(z$sigma.loc)

	# Check, per locus measures: ** this is the one we're keeping
	# tsiga <- sum(z$sigma.loc[,"lsiga"], na.rm=TRUE)
	# tsigb <- sum(z$sigma.loc[,"lsigb"], na.rm=TRUE)
	# tsigw <- sum(z$sigma.loc[,"lsigw"], na.rm=TRUE)
	# tsiga/sum(c(tsiga, tsigb, tsigw))

	# Another check, using per allele measures:
	# tsiga <- sum(z$sigma[,"siga"], na.rm=TRUE)
	# tsigb <- sum(z$sigma[,"sigb"], na.rm=TRUE)
	# tsigw <- sum(z$sigma[,"sigw"], na.rm=TRUE)
	# tsiga/sum(c(tsiga, tsigb, tsigw))

	z1 <- cbind.data.frame(z$sigma.loc, fst = z$per.loc$FST, fis = z$per.loc$FIS)
	# rownames(z1) <- rownames(vcfresults$GT) # error because of duplicate names
	
	# range of variance among groups at single loci - maximum is 1, same as psd
	# range(z1$fst$lsiga)
	# [1] -0.1666667  1.0000000
	
	
	# Need to MATCH result with larger genotype data set, set fst = NA for non-polymorphic loci
	
	gtstats$fst <- z1
	rm(z1)

	pdf(file = paste(vcfdir, project, ".", chrname, ".FstPlot", ".pdf", sep=""))

	# Plot Fst
	# add a tiny random number to reduce discrete bands
	y <- vcfresults$fst$fst + rnorm(length(start(vcfresults$rowdata)), mean = 0, sd = .005)
	x <- start(vcfresults$rowdata)/(10^6)
	colour <- abs(round(abs(vcfresults$fst$fst) * 4.8)) + 1
	plot( y ~ x, pch = ".", ylab="Fst", col = colour, xlab="MB",
		main = "Fst (max 3 ALT alleles)" ) 

	dev.off()
	
	# ----
	# Transition to transversion ratios for significant snp - snp only (not indels) first ALT allele only
	# Transition-transversion ratios - snp only (not indels) first ALT allele only

	cat("\nTransition-transversion ratio for true snps having Fst >= 0.95 snp (uses first ALT allele if more than one)\n")
	snptype <- sapply(vcfresults$snptypeList, function(x){x[1]}) # first ALT allele: snp ins or del?
	REF <- as.character(vcfresults$rowdata$REF)
	ALT <- sapply(vcfresults$ALTlist, function(x){x[1]})
	print(g$tstv(REF[snptype == "snp" & vcfresults$fst$fst >= .95], ALT[snptype == "snp" & vcfresults$fst$fst >= .95]))
	# $R
	# [1] 1.237113
	# $tstv
	     # alt
	# ref   pur pyr
	  # pur 469 396
	  # pyr 380 491
	# $tot
	# [1] 1736

	# g$tstv(vcfresults$REF[vcfresults$fst$fst >.95],vcfresults$ALT[fst$fst > 0.99])
	
	vcfresults <- vcfresults

	# cat("\nInterim save\n")
	# save(vcfresults, file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
	# load(file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
	# vcfresults <- vcfresults

	} # end if(includeFst}
	
# --------------

if(includePsd){
	
	vcfresults <- vcfresults

	cat("\nCalculating 'percent sequence divergence' between each pair of individuals at each variant
		(indels are treated simply as alleles)\n")
	
	psd <- g$psd(vcfresults$GT, groups)
	vcfresults$psd <- psd
	rm(psd)

	vcfresults <- vcfresults
	} # end if(includePsd)
	
# -------------


cat("\nSaving results\n")
save(vcfresults, file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
# load(file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = "")) # vcfresults
# vcfresults <- vcfresults
	
