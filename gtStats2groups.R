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

# setwd("~/Desktop")
# git("genome.r")

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn", "paxl", "paxb")
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb")

GTminFrac <- 2/3
ancestor <- "marine-pac" # use to classify the most common ancestral allele for comparison

includePfisher 	<- FALSE
includeFst     	<- TRUE
includePsd		<- TRUE
trueSnpOnly		<- FALSE 
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

gcinfo(TRUE)
gc()

alleleFreqByGroup <- lapply(vcfresults$alleleFreqByGroup, function(x){x[groupnames,]})
vcfresults$alleleFreqByGroup <- NULL


# ----------
# Drop the rows with insufficient numbers of alleles - these will not be seen again, whether variant or invariant

# Identify loci with enough alleles in both groups (the minimum number of alleles is nMin*2)
sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
	z <- rowSums(x, na.rm = TRUE)
	z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
	})
	
# length(sufficient.alleles.per.group)
# [1] 281607

# Drop
genotypes <- genotypes[sufficient.alleles.per.group, ]

alleleFreqByGroup <- alleleFreqByGroup[sufficient.alleles.per.group]

# Still need this here if want to remove non-polymorphic indels
snpTypeList <- vcfresults$snpTypeList[sufficient.alleles.per.group]
vcfresults$snpTypeList <- NULL

gtstats$vcf <- vcfresults$vcf[sufficient.alleles.per.group] # contains rowData but not genotypes (GT empty)

rm(sufficient.alleles.per.group)
rm(vcfresults) # assuming we have extracted all the useful bits

gc()
# chrXXI
           # used  (Mb) gc trigger  (Mb) max used  (Mb)
# Ncells  6037587 322.5   10049327 536.7 10049327 536.7
# Vcells 17640420 134.6   39455907 301.1 39417064 300.8
# hermes chrUn
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells 12407797 662.7   25310187 1351.8  25310187 1351.8
# Vcells 50531600 385.6  128513897  980.5 128347380  979.3


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

# length(alleleFreqByGroup[is.polymorphic])
# [1] 153997

# Don't drop the non-polymorphic sites - these need to be counted with the invariants!

# -----------------
# Drop everything except *true snp and invariants* if trueSnpOnly = TRUE
# Remember: an indel w.r.t. REF might be invariant between limnetic and benthic, so is not an indel here.
# 	so we drop only the polymorphic loci that are classified as indels
if(trueSnpOnly){

	# tGT is transposed genotypes but does NOT include invariant sites between the two groups
	tGT <- as.data.frame(t(genotypes[is.polymorphic, ]), stringsAsFactors = FALSE)
	# ncol(tGT)
	# [1] 153997

	# Identify the ALT alleles that are indels (NAs are ignored)
	which.alt.not.snp <- lapply(snpTypeList[is.polymorphic], function(x){
		which(x != "snp")
		})
	# test <- c("chrXXI:75572_CACTG/C", "chrXXI:75793_CTTG/C", "chrXXI:76075_T/C")
	# genotypes[test, ]
	# snpTypeList[test]
	# which.alt.not.snp[test]
	
	# Sets tGT individual genotypes to NA if not a true SNP.
	z <- mapply(tGT, which.alt.not.snp, FUN = function(x, i){
		if(sum(i) > 0) x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		return(x)
		})
	z <- t(z) # is a matrix again, has the right rownames but colnames are absent
	colnames(z) <- rownames(tGT)
	# z[test,]
	
	genotypes[is.polymorphic, ] <- z
	
	# groupcodes
	 # [1] 0 0 0 0 0 0 0 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	# groups
	 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

	# redo allele frequencies by group - use groups as argument not groupcodes
	alleleFreqByGroup[is.polymorphic] <- g$tableAlleleFreqByGroup(genotypes[is.polymorphic, ], 
		groupnames, groups) 
		
	gc()
		       # used  (Mb) gc trigger  (Mb) max used  (Mb)
	# Ncells  6681616 356.9   11759451 628.1 11759451 628.1
	# Vcells 28811593 219.9   57298003 437.2 57298003 437.2


	# Filter once more to drop rows having insufficient numbers of alleles after indels deleted
	sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
		z <- rowSums(x, na.rm = TRUE)
		z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
		})
	# table(sufficient.alleles.per.group)
	 # FALSE   TRUE 
	 # 12629 256720

	genotypes <- genotypes[sufficient.alleles.per.group, ]
	
	alleleFreqByGroup <- alleleFreqByGroup[sufficient.alleles.per.group]
	
	gtstats$vcf <- gtstats$vcf[sufficient.alleles.per.group] # tossing the indels for good
	
	# snpTypeList <- snpTypeList[sufficient.alleles.per.group] # don't really need this anymore
	
	# nrow(genotypes)
	# [1] 256720

	}  

gc() # if trueSnpOnly = TRUE
           # used  (Mb) gc trigger  (Mb) max used  (Mb)
# Ncells  6796335 363.0   11759451 628.1 11759451 628.1
# Vcells 29382190 224.2   57298003 437.2 57298003 437.2


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

# We're running this again because... (might be needed only if removed indels above)
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
	
	gc()
	
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

gc()

# --------------

if(includeFst){
	
	# Prepare genotypes for Fst

	# Initialize so that we can set all stats for non-polymorphic sites to NA
	fst <- data.frame(row.names = rownames(gtstats$genotypes) )

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
	names(temp) <- names(geno)

	# head(names(temp)) # if not set names above
	#[1] "V1" "V2" "V3" "V4" "V5" "V6"

	#table(pop)
	# pop
	 # 1  2 
	# 11 11 

	gc() # when trueSnpOnly = FALSE
		       # used  (Mb) gc trigger  (Mb) max used  (Mb)
	# Ncells  6395311 341.6   14466089 772.6 12131744 648.0
	# Vcells 28639355 218.6   66338776 506.2 63103588 481.5

	rm(geno)
	
	# Must use groups instead of pop in wc command because must be a number

	# Calculate Fst - # took about 40 min
	# z <- wc(cbind(groups, temp))

	# Use this slightly faster version - took about 20 minutes
	z <- g$wc.revised(cbind(groups,temp))
	
	gc() # when trueSnpOnly = FALSE
		           # used  (Mb) gc trigger  (Mb) max used  (Mb)
	# Ncells  6860317 366.4   14466089 772.6 14466089 772.6
	# Vcells 32533361 248.3   85109004 649.4 73302178 559.3

	cat("\nFst\n")
	print(z$FST)
	#[1] 0.4076119
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
	
	fst$lsiga <- rep(NA, nrow(fst))
	fst[rownames(z$sigma.loc), 1] <- z$sigma.loc[,"lsiga"]
	fst$lsigb <- rep(NA, nrow(fst))
	fst[rownames(z$sigma.loc), 2] <- z$sigma.loc[,"lsigb"]
	fst$lsigw <- rep(NA, nrow(fst))
	fst[rownames(z$sigma.loc), 3] <- z$sigma.loc[,"lsigw"]
	fst$fst <- rep(NA, nrow(fst))
	fst[rownames(z$sigma.loc), 4] <- z$per.loc$FST
	fst$fis <- rep(NA, nrow(fst))
	fst[rownames(z$sigma.loc), 5] <- z$per.loc$FIS

	# head(fst)	
	                        # lsiga       lsigb      lsigw         fst         fis
	# chrXXI:16024_C/T           NA          NA         NA          NA          NA
	# chrXXI:16025_C/G           NA          NA         NA          NA          NA
	# chrXXI:16026_A/T -0.003954248  0.32956656 0.05263158 -0.01045423  0.86229243
	# chrXXI:18389_C/A -0.008046977 -0.01007063 0.19047619 -0.04668742 -0.05582218
	# chrXXI:18511_A/G  0.027777778  0.07222222 0.10000000  0.13888889  0.41935484
	# chrXXI:18549_C/T -0.009721618  0.09556847 0.04761905 -0.07283972  0.66743575	
	
	# Need to MATCH result with larger genotype data set, set fst = NA for non-polymorphic loci
	
	gtstats$fst <- fst
	
	gc()
		       # used  (Mb) gc trigger  (Mb) max used  (Mb)
	# Ncells  6865570 366.7   14466089 772.6 14466089 772.6
	# Vcells 33888695 258.6   85109004 649.4 73302178 559.3

	rm(fst)

	# pdf(file = paste(project, ".", chrname, ".FstPlot", ".pdf", sep=""))
	# Plot Fst
	# add a tiny random number to reduce discrete bands
	# y <- vcfresults$fst$fst + rnorm(length(start(vcfresults$rowdata)), mean = 0, sd = .005)
	# x <- start(vcfresults$rowdata)/(10^6)
	# colour <- abs(round(abs(vcfresults$fst$fst) * 4.8)) + 1
	# plot( y ~ x, pch = ".", ylab="Fst", col = colour, xlab="MB",
		# main = "Fst (max 3 ALT alleles)" ) 

	# dev.off()

	# Interim save
	# save(gtstats, file = gtstatsfile)
	
	} # end if(includeFst}

# load(file = gtstatsfile) # object is "gtstats"

gc()

# --------------

if(includePsd){
	
	cat("\nCalculating 'percent sequence divergence' between each pair of individuals at each variant
		(indels are treated simply as alleles)\n")
	
	psd <- g$psd(gtstats$genotypes, groups)
	
	gc()
	           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells  6866870 366.8   14466089  772.6  14466089  772.6
	# Vcells 96380275 735.4  312091343 2381.1 311994807 2380.4 ** whoa!

	gtstats$psd <- psd
	rm(psd)

	} # end if(includePsd)
	
# -------------

gc()

cat("\nSaving results\n")
save(gtstats, file = gtstatsfile)
# load(file = gtstatsfile) # object is "gtstats"
	
