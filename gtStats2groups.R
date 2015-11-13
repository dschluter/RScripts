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
# args <- c("BenlimAllMarine", "chrXXI", "marine-pac", "qryb")

GTminFrac <- 2/3

includePfisher 	<- FALSE
includeFst     	<- TRUE
includePsd		<- FALSE
trueSnpOnly		<- FALSE
# note: An indel relative to REF genome, if fixed and shared between limnetic and benthic, 
# is not an indel between limnetic and benthic, but a fixed difference. 
# However, REF indel that is polymorphic within or between limnetics and benthics is an indel.
# So handle indels carefully.
# Q: What if two alt alleles that are indels relative to REF are the same width. 
#	 Might they be true snps between benthic and limnetic? I have not investigated.

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)]
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

gtstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "."), "rdd", sep = ".")
gtstats <- list()

if(length(groupnames) > 2 ) stop("Provide names of only two groups")

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcfresults"
# lapply(vcfresults, object.size)

# names(vcfresults)
# [1] "groupnames"        "groupcodes"        "nInd" 	"control"     "vcf"       "altUsedList"      
# [7] "snpTypeList"       "alleleFreqByGroup"

control <- vcfresults$control
control$GTminFrac <- GTminFrac
control$trueSnpOnly <- trueSnpOnly
control$includePfisher <- includePfisher
control$includeFst <- includeFst
control$includePsd <- includePsd


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

cat("\nMinimum sample size criterion based on GTminFrac = ", GTminFrac, "*nInd or 5, whichever is smaller\n", sep = "")
print( data.frame(groupnames=groupnames, nInd, nMin) )

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
	
# length(sufficient.alleles.per.group[sufficient.alleles.per.group])
# [1] 269349

# Drop cases with insufficient numbers of individuals
# Keep snpTypeList too, need it if want to remove non-polymorphic indels later

genotypes <- genotypes[sufficient.alleles.per.group, ]

alleleFreqByGroup <- alleleFreqByGroup[sufficient.alleles.per.group]

gtstats$vcf <- vcfresults$vcf[sufficient.alleles.per.group] # contains rowData but not genotypes (GT empty)

snpTypeList <- vcfresults$snpTypeList[sufficient.alleles.per.group]

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

# table(genotypes[rownames(genotypes) == "chrXXI:54638_T/TG",])
# 0/0 1/1 
 # 20   1 

# ---------------------------------
# Determine which loci are variable, ie have more than one allele across both groups.
# Those that are not polymorphic we want to keep as a record of an invariant between these two groups

# "i" = invariant
# "v" = variant
# "n" = missing, indicates that the snp is dropped in a later step

status <- rep("v", length(alleleFreqByGroup)) # initialize
is.invariant <- sapply(alleleFreqByGroup, function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) < 2
	})
status[is.invariant] <- "i"
rm(is.invariant)

# table(status)
     # i      v
# 115352 153997

# status[rownames(genotypes) == "chrXXI:54638_T/TG"]
# [1] "v"

# Don't drop the invariant sites - these need to be counted with the invariants when doing block stats

# -----------------
# If trueSnpOnly = TRUE, drop everything except *true snp and invariants* 
# Sites have multiple alt alleles, and only some are true snps

# Remember: if an indel w.r.t. REF is fixed for the same allele in both limnetic and benthic, it is not an indel.
# 	so we drop only the polymorphic loci that are classified as indels
if(trueSnpOnly){

	# tGT is transposed genotypes but does NOT include invariant sites between the two groups
	tGT <- as.data.frame(t(genotypes[status == "v", ]), stringsAsFactors = FALSE)
	# ncol(tGT)
	# [1] 153997

	# Identify the ALT alleles that are indels (NAs are ignored)
	which.alt.not.snp <- lapply(snpTypeList[status == "v"], function(x){
		which(x != "snp")
		})
	# test <- c("chrXXI:75572_CACTG/C", "chrXXI:75793_CTTG/C", "chrXXI:76075_T/C")
	# genotypes[test, ]
	# snpTypeList[test]
	# which.alt.not.snp[test]
	
	# which.alt.not.snp["chrXXI:54638_T/TG"]
	# $`chrXXI:54638_T/TG`
	# [1] 1
	
	# Sets tGT individual genotypes to NA if not a true SNP.
	z <- mapply(tGT, which.alt.not.snp, FUN = function(x, i){
		if(sum(i) > 0) x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		return(x)
		})
	z <- t(z) # is a matrix again, has the right rownames but colnames are absent
	colnames(z) <- rownames(tGT)
	# z[test,]
	# table(z["chrXXI:54638_T/TG",])
	# 0/0 
	 # 20 
	
	genotypes[status == "v", ] <- z
	rm(which.alt.not.snp)
	
	# groupcodes
	 # [1] 0 0 0 0 0 0 0 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	# groups
	 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

	# redo allele frequencies by group - use groups as argument not groupcodes
	alleleFreqByGroup[status == "v"] <- g$tableAlleleFreqByGroup(genotypes[status == "v", ], 
		groupnames, groups) 

	# alleleFreqByGroup["chrXXI:54638_T/TG"]
	
	gc()
		       # used  (Mb) gc trigger  (Mb) max used  (Mb)
	# Ncells  6681616 356.9   11759451 628.1 11759451 628.1
	# Vcells 28811593 219.9   57298003 437.2 57298003 437.2


	# Filter once more to change the status of rows having insufficient numbers of alleles after indels deleted
	sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
		z <- rowSums(x, na.rm = TRUE)
		z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
		})
	# table(sufficient.alleles.per.group)
	 # FALSE   TRUE 
	 # 12629 256720
	
	status[!sufficient.alleles.per.group] <- "n"
	rm(sufficient.alleles.per.group)
	# table(status)
	     # i      n      v 
	# 115352  12629 141368
	
	# Who is still polymorphic (among those formerly with status "v")
	is.invariant <- sapply(alleleFreqByGroup[status == "v"], function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) < 2
	})
	status[status == "v"][is.invariant] <- "i"
	rm(is.invariant)

	# table(status)
	# status
	     # i      n      v 
	# 126841  12629 129879
		
	# From now on, analyze only the snp with status = "v"

	} # end if(trueSnpOnly)

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
gtstats$nInd <- nInd
gtstats$nMin <- nMin
gtstats$control <- control
gtstats$genotypes <- genotypes		# rows are loci
rm(genotypes)
gtstats$alleleFreqByGroup <- alleleFreqByGroup
rm(alleleFreqByGroup)
gtstats$status <- status
rm(status)

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
	gtstats$pfisher[gtstats$status == "v"] <- sapply(gtstats$alleleFreqByGroup[gtstats$status == "v"], function(x){
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

	geno <- as.data.frame(t(gtstats$genotypes[gtstats$status == "v", ]), stringsAsFactors = FALSE) 

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
	temp <- g$recode(temp, c("0/0","0/1","1/1","0/2","1/2","2/2","0/3","1/3","2/3","3/3"), 
							c("11","12", "22", "13", "23", "33", "14", "24", "34", "44"))
	temp <- as.data.frame(matrix(as.integer(temp), nrow = length(pop))) # original names lost
	colnames(temp) <- colnames(geno)

	# head(names(temp)) 

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
	#[1] 0.4078011
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
	
	# This saves back to the rows of "fst" where status == "v"
	fst$lsiga <- rep(NA, nrow(fst)) # initialize
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
	
	# head(cbind.data.frame(gtstats$status, fst))
		
	gtstats$fst <- fst
	rm(fst)
	
	gc()
		       # used  (Mb) gc trigger  (Mb) max used  (Mb)
	# Ncells  6865570 366.7   14466089 772.6 14466089 772.6
	# Vcells 33888695 258.6   85109004 649.4 73302178 559.3

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
	
	z <- g$psd(gtstats$genotypes[gtstats$status == "v", ], groups)
	
	# nrow(z)
	# head(z)
	
	# initialize
	psd <- as.data.frame( matrix(nrow = nrow(gtstats$genotypes), ncol = ncol(z)), stringsAsFactors = FALSE )
	rownames(psd) <- rownames(gtstats$genotypes)
	
	psd[ rownames(z) , ] <- z
	colnames(psd) <- colnames(z) # duplicate names were ok
	rm(z)
	
	gc()

	gtstats$psd <- psd
	rm(psd)

	} # end if(includePsd)
	
# -------------

gc()

cat("\nSaving results\n")
save(gtstats, file = gtstatsfile)
# load(file = gtstatsfile) # object is "gtstats"
	
