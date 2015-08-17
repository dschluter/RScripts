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

GTminFrac <- 2/3
ancestor <- "marine-pac" # use to classify the most common ancestral allele for comparison

includePfisher 	<- TRUE
includeFst     	<- TRUE
includePsd		<- TRUE
trueSnpOnly		<- TRUE 
# note that an indel relative to reference genome, if shared between limnetic and benthic, 
# is not an indel between the latter. It is a fixed difference. 
# So maybe don't want to leave indels out, or at least count them differently (NOT DONE).
# However, a polymorphic REF indel between limnetic and benthic is an indel, so want to remove these. 

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
# [1] "groupcodes"        "vcf"               "altUsedList"       "snpTypeList"       "alleleFreqByGroup"

library(VariantAnnotation)
# library(GenomicFeatures)

# Group codes for the two groups of interest here
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
groups <- groupcodes[groupcodes > 0]
# groups
 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

nInd <- as.vector(table(groups)) 	# n individuals genotyped in each group eg 11 11  7
nMin <- floor(GTminFrac*nInd) 		# minimum number required in each group eg 7 7 4
nMin <- mapply(rep(5, length(groupnames)), nMin, FUN = min) # criterion is 5 or nMin, whichever is smaller, eg 5 5 4

# ----------
# Pull out the genotypes for just the two groups being analyzed
# Bring over coding annotations too, if present

genotypes <- geno(vcfresults$vcf)$GT[ , groupcodes > 0]
geno(vcfresults$vcf)$GT <- NULL
alleleFreqByGroup <- lapply(vcfresults$alleleFreqByGroup, function(x){x[groupnames,]})
vcfresults$alleleFreqByGroup <- NULL

# ----------
# Drop the rows with insufficient numbers of alleles

# Identify loci with enough alleles in both groups (the minimum number of alleles is nMin*2)
sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
	z <- rowSums(x, na.rm = TRUE)
	z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
	})
	
# length(sufficient.alleles.per.group)
# [1] 281607

# Drop the rows with insufficient numbers of alleles
genotypes <- genotypes[sufficient.alleles.per.group, ]
alleleFreqByGroup <- alleleFreqByGroup[sufficient.alleles.per.group]
snpTypeList <- vcfresults$snpTypeList[sufficient.alleles.per.group]
vcfresults$snpTypeList <- NULL
gtstats$vcf <- vcfresults$vcf[sufficient.alleles.per.group]

# nrow(genotypes)
# [1] 269349

# ---------------------------------
# Determine which loci are variable, ie have more than one allele across both groups.
# Those that are not polymorphic we want to keep as a record of an invariant.
# Note that even if an indel w.r.t. REF, an invariant between limnetic and benthic is not an indel.
is.polymorphic <- sapply(alleleFreqByGroup, function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) >= 2
	})

# length(is.polymorphic)
# [1] 269349

# -----------------
# Drop everything except *true snp and invariants* if trueSnpOnly = TRUE
if(trueSnpOnly){

	# tGT is transposed genotypers but does NOT include invariant sites between the two groups
	tGT <- as.data.frame(t(genotypes[is.polymorphic, ]), stringsAsFactors = FALSE)
	# ncol(tGT)
	# [1] 163572

	# Identify the ALT alleles that are indels (NAs are ignored)
	which.alt.not.snp <- lapply(snpTypeList[is.polymorphic], function(x){
		which(x != "snp")
		})
	# test <- c("chrXXI:75572_CACTG/C", "chrXXI:75793_CTTG/C", "chrXXI:76075_T/C")
	# genotypes[test, ]
	# snpTypeList[test]
	# which.alt.not.snp[test]
	
	# Sets tGT to NA if not a true SNP.
	z <- mapply(tGT, which.alt.not.snp, FUN = function(x, i){
		if(sum(i) > 0) x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		# x
		return(x)
		})
	z <- t(z) # is a matrix again, has the right rownames but colnames are absent
	colnames(z) <- rownames(tGT)
	z[test,]
	
	genotypes[is.polymorphic, ] <- z
	
	# groupcodes
	 # [1] 0 0 0 0 0 0 0 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	# groups
	 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

	# redo allele frequencies by group - use groups as argument not groupcodes. is.polymorphic saves a little time
	alleleFreqByGroup[is.polymorphic] <- g$tableAlleleFreqByGroup(genotypes[is.polymorphic, ], 
		groupnames, groups) 

	# Filter once more to drop rows having insufficient numbers of alleles after indels deleted
	sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
		z <- rowSums(x, na.rm = TRUE)
		z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
		})
	# table(sufficient.alleles.per.group)
	# sufficient.alleles.per.group
	 # FALSE   TRUE 
	 # 12702 256647

	genotypes <- genotypes[sufficient.alleles.per.group, ]
	alleleFreqByGroup <- alleleFreqByGroup[sufficient.alleles.per.group]
	snpTypeList <- snpTypeList[sufficient.alleles.per.group]
	gtstats$vcf <- gtstats$vcf[sufficient.alleles.per.group]
	
	# nrow(genotypes)
	# [1] 256647

	}  

# Do I still need this? Calculate it here after we have the right "genotype" matrix, however filtered above
altUsedList <- g$makeAltUsedList(genotypes, alt(gtstats$vcf))

rm(sufficient.alleles.per.group)
rm(vcfresults) # assuming we have extracted all the useful bits


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

gtstats$genotypes <- genotypes		# rows are loci
rm(genotypes)

gtstats$alleleFreqByGroup <- alleleFreqByGroup
rm(alleleFreqByGroup)

gtstats$altUsedList <- altUsedList # needed?
rm(altUsedList)

gtstats$snpTypeList <- snpTypeList # needed?
rm (snpTypeList)

is.polymorphic <- sapply(gtstats$alleleFreqByGroup, function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) >= 2
	})

# -------------------------------------------------------------
# Group differences in allele frequencies using Fisher exact test
# P = 1 if there are no differences or locus is monomorphic

if(includePfisher){ 
	# Fisher exact tests of allele and genotype frequencies using ALLELE (not genotype) frequencies

	# Set to NA non-polymorphic sites with sufficiently many genotypes
	gtstats$pfisher <- rep(NA, length(gtstats$alleleFreqByGroup))
	names(gtstats$pfisher) <- names(gtstats$alleleFreqByGroup)
	
	cat("\nCalculating Fisher exact test p-values for group allele freq differences (maximum of 3 ALT alleles assumed)\n")
	# Loci that are not divergent will have a pfisher = 1
	gtstats$pfisher[is.polymorphic] <- sapply(gtstats$alleleFreqByGroup[is.polymorphic], function(x){
		pfisher <- fisher.test(x)$p.value
		})
	
	# Write to a file in bedGraph format
	# g$write.bedGraph(paste(chrname,".pfisher.bedGraph", sep=""), chrinfo = chrname, vcfresults$POS, vcfresults$POS+1, round(-log10(pfisher),3),
		# graphname = "-log10Pfisher", col = "41,155,241", altcol = "100,100,100", ylinemark = 4,
		# description = "-log10Pfisher") 

	# Plot where the significant SNPs reside
	
	# pdf(file = paste(project, ".", chrname, ".pfisherPlot", ".pdf", sep=""))	
	# x <- start(rowData(gtstats$vcf))/(10^6)
	# y <- (-log10(gtstats$pfisher) + rnorm(length(x), mean=0, sd=.005))
	# colour <- (round(-log10(gtstats$pfisher)/1.2)+1)
	# plot( y ~ x, pch = ".", ylab="-log10 Pfisher", col = colour, xlab = chrname,
		# main = "Fisher exact test group allele freq differences (max 3 ALT alleles)") 
	# dev.off()
	
	# cat("\nTransition-transversion ratio for pfisher < 0.01\n")
	}
	
# --------------

if(includeFst){
	
	# Prepare genotypes for Fst

	# Set to NA non-polymorphic sites with sufficiently many genotypes
	gtstats$fst <- rep(NA, length(gtstats$alleleFreqByGroup))
	names(gtstats$fst) <- names(gtstats$alleleFreqByGroup)


	geno <- as.data.frame(t(gtstats$genotypes[is.polymorphic, ]), stringsAsFactors = FALSE) 

	library(hierfstat)
	
	cat("\nCalculating Fst (slowly) between groups\n")
	
	# Per-locus Fst using Weir and Cockrham estimates
	# Fst is calculated as tsiga/sum(c(tsiga, tsigb, tsigw)), where
	# tsiga <- sum(lsiga)
	# tsigb <- sum(lsigb)
	# tsigw <- sum(lsigw)
	# lsiga, lsigb, and lsigw are the variance components: for each locus among populations,
	#   among individuals within populations, and within individuals, respectively
	# wc command does not like duplicate row names

	# groups
	# [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

	pop <- as.character(groups)

	# Convert genotypes to hierfstat format:
	# Columns are genotypes indicated as 11 (0/0), 12 (0/1), 22 (1/1), or NA (NA)
	
	temp <- unlist(geno)
	
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

	temp <- as.data.frame(matrix(as.integer(temp), nrow = length(pop))) # original names lost
	#names(temp) <- names(geno)

	#head(names(temp)) 
	#[1] "V1" "V2" "V3" "V4" "V5" "V6"

	#table(pop)
	# pop
	 # 1  2 
	# 11 11 

	# Must use groups instead of pop in wc command because must be a number

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
	
