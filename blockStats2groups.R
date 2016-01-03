#!/usr/bin/Rscript

# Calculates preliminary results for sliding window analyses on a (repeat-masked) chromosome
# It breaks the chromosome into blocks or bins of size stepsize and calculates a value of interest in each bin
# Combines variant and invariant information in genome blocks of size 500 (default) nucleotides

# gtstats is the gtstats list containing the variables of interest for a single chromosome
# Method asks whether gtstats$fst and gtstats$psd exist and computes and saves blockstats accordingly.

# GATK 3.4 no longer gives multiple rows of snps at the same value of POS

# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimAllMarine", "chrXXI", "marine-pac", "paxb")
# args <- c("BenlimAllMarine", "chrXXI", "marine-pac", "paxl")
# args <- c("BenlimAllMarine", "chrXXI", 500, "paxl", "paxb")
# args <- c("BenlimAllMarine", "chrXX", 500, "pril", "prib") # one of the smaller ones that crashed because of memory

project <- args[1]
chrname <- args[2]
stepsize <- as.integer(args[3])
groupnames <- args[-c(1:3)]
chrno <- gsub("chr", "", chrname)

psdMissingAction <- "meanBW" 
	# psdMissingAction is for g$blockstats, how to average when psd values are missing. Must be one of the following:
	# 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
	# 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
	#   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
	#		identified by missingGroups.
	#		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
	#		"1,1" and "2,2" as the same (i.e., both are within-group).
if( !is.element(psdMissingAction, c("meanAll", "meanBW", "meanGroup")) )
	stop("psdMissingAction must be one of meanAll, meanBW, or meanGroup")

cat("\nProject, chrname, groupnames, stepsize", project, chrname, groupnames, stepsize, "\n")

gtstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "."), "gtstats.rdd", sep = ".")	# object is gtstats
blockstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "."), "blockstats", stepsize, "rdd", sep = ".")

load(gtstatsfile) # object is gtstats
# names(gtstats)
 # [1] "vcf"        "groupnames" "groups"     "nInd"       "nMin"      
 # [6] "control"    "genotypes"  "status"     "newPos"     "fst"       
# [11] "psd"      

# object.size(gtstats) 

control <- gtstats$control
Glazerize <- control$Glazerize
groups <- gtstats$groups
groupnames <- gtstats$groupnames
nInd <- gtstats$nInd 	# n individuals genotyped in each group eg 7 11
nMin <- gtstats$nMin 	# criterion is 5 or nMin, whichever is smaller, eg 4 5
# control$psdMissingAction <- psdMissingAction # not used right now
status <- gtstats$status
fst <- gtstats$fst # will be NULL if absent
psd <- gtstats$psd # will be NULL if absent

gc()

if(Glazerize){
	goodInvariantsFile <- paste(project, ".", chrname, ".goodInvNew.rdd", sep="")
	chrvecfile = paste("chrvec.", chrno, ".glazer.rdd", sep = "")
	pos <- gtstats$newPos # pos is for the variants, use POS for the invariants file
	
	} else {
		
	# *Need to confirm that this works if not glazerizing* (because gtstats$vcf is a vcf not a list)

	goodInvariantsFile <- paste(project, ".", chrname, ".goodInv.rdd", sep="") # object is goodInvariants
	chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")
	library(VariantAnnotation)
	pos <- start(rowData(gtstats$vcf)) # pos is for the variants, use POS for the invariants file
	}

rm(gtstats)

gc()

# head(pos)
# [1] 1794949 1794951 1794957 1794962 1794969 1794974

# All start positions are unique, snps and indels at the same start POS are just different alleles
# c( length(pos), length(unique(pos)) )
# [1] 871120 871120

load(goodInvariantsFile) # object is goodInvariants

# object.size(goodInvariants) 

# Fix the names in goodInvariants (change "marine.pac" to "marine-pac", etc)
names(goodInvariants) <- gsub("[.]", "-", names(goodInvariants))

if(Glazerize){
	goodInvariants <- goodInvariants[, c("newPos", groupnames)] # grab the columns of interest
	# rename "newPos" to "POS"
	names(goodInvariants)[names(goodInvariants) =="newPos"] <- "POS"
	} else{
	goodInvariants <- goodInvariants[, c("POS", groupnames)] # grab the columns of interest
	}

# Drop invariants that do not meet the minimum number of genotypes required
# This adds a lot of memory use, maybe try data.table instead
z <- apply(goodInvariants[, groupnames], 1, function(x){
	all(x >= nMin)
	})
goodInvariants <- goodInvariants[z, ]
rm(z)

# chrvec is a vector of the whole MASKED chromosome, all nucleotides as distinct elements, obtained by:
#       Masked bases are indicates with an "M"
load(chrvecfile) # object is chrvec

gc()

# -------------
# Get block window stats for group differences

	
# pop <- as.character(gtstats$groups)
pop <- as.character(groups)
k <- as.integer(stepsize) # this is the size of the block
		
# Establish the break points of the bins into which nucleotides will be grouped (e.g., k = 500 bases per bin)
# The last bin goes from the final bin of fully k nucleotides to the last nucleotide. 
# This last bin might be small.
	
nbases <- length(chrvec)
ibase <- c( seq(1, nbases, k), nbases + 1) # bases marking breaks between steps of size k
midbase <- ibase + (stepsize)/2 # If want something to indicate midpoint of blocks instead

# head(midbase)
# [1]  251  751 1251 1751 2251 2751
# head(ibase)
# [1]    1  501 1001 1501 2001 2501
# tail(ibase)
# [1] 17355501 17356001 17356501 17357001 17357501 17357773 # * last one is less than k

# ----
# 1) Count number of non-M bases in each nucleotide bin of stepsize "k". 
	
# Break nucleotide index of chr into bins of stepsize k
chrbins <- findInterval(1:nbases, ibase) # indicates in which "bin" of size k each base belongs
										 # if none missing, there should be k 1's, k 2's, etc.
# chrbins[499:505]
# [1] 1 1 2 2 2 2 2
	
# count up the number of unmasked nucleotides in each bin - These are the ones coded as ACGT
nUnmasked <- tapply(chrvec, chrbins, function(x){length(x[x %in% c("A","C","G","T")])})

# nUnmasked[1:20]
 # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
 # 0  0  0  9  0  0  0  0  0  0  0  0  0  0  4  0  0 11  0  0

rm(chrvec)

# ----
# 2) Count up the number of snp and number of invariants in each bin

# Need to split SNP counts into same bins of stepsize k and then sum up

# Count up the number of snp in each bin (cut works because ibase intervals are made as factors)
snpBins <- cut(pos, breaks = ibase, right = FALSE)
nSnp <- table( snpBins )

# head(nSnp)
	# [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001) 
	#       0           0           0           0           0           0 
    
# Data still includes some invariants, those sites not polymorphic in these two populations
# The following is now redundant -- I think this information is already in gtstats$status
		# is.polymorphic <- sapply(gtstats$alleleFreqByGroup, function(x){
		# z <- colSums(x, na.rm = TRUE)
		# sum( z > 0 ) >= 2
		# })
		# table(is.polymorphic, gtstats$status)
		# is.polymorphic      i      v
		         # FALSE 627449      0
		         # TRUE       0 243671
		# is.monomorphic <- split(!is.polymorphic, snpBins)

is.monomorphic <- split(status == "i", snpBins)
nMonomorphic <- sapply(is.monomorphic, sum)

# Count up the number of good invariants in each bin
# Remember: newPos has been renamed to POS if Glazerize is true
nInvariants <- table( cut(goodInvariants$POS, breaks = ibase, right = FALSE) ) 
	
# head(nInvariants)
    # [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001)  
    #       0           0           0           0           0           0 
    
rm(goodInvariants)

# all the following items have the same length

blockstats <- data.frame(
			ibase = ibase[-length(ibase)], 
			midbase = midbase[-length(midbase)], 
			nUnmasked = nUnmasked, 
			nSnp = as.vector(nSnp) - nMonomorphic, 
			nInvariants = as.vector(nInvariants) + nMonomorphic
			)

rm(ibase)
rm(midbase)
rm(nUnmasked)
rm(nSnp)
rm(nInvariants)
rm(nMonomorphic) 

gc()

# blockstats[110:130,]
	    # ibase midbase nUnmasked nSnp nInvariants
	# 110 54501   54751         0    0           0
	# 111 55001   55251         0    0           0
	# 112 55501   55751        28    0           0
	# 113 56001   56251         0    0           0
	# 114 56501   56751         9    0           0
	# 115 57001   57251         9    0           0
	# 116 57501   57751         5    0           0
	# 117 58001   58251        99    4          60
	# 118 58501   58751       193   33         160
	# 119 59001   59251        89    6          83
	# 120 59501   59751       411   19         389
	# 121 60001   60251       398    3         194
	# 122 60501   60751       343    5         222
	# 123 61001   61251       314    3         266
	# 124 61501   61751       456    5         180
	# 125 62001   62251       336    8          88
	# 126 62501   62751         4    0           0
	# 127 63001   63251         8    0           0
	# 128 63501   63751       106    0           0
	# 129 64001   64251       283    0           0
	# 130 64501   64751       363    0           6
	
# ----
# 3) Calculate Fst summary stats

if( !is.null(fst) ){
		
	cat("\nBeginning Fst calculations\n")

	# fst is a matrix, needs to be a data.frame of split will screw up
	fst <- as.data.frame(fst, stringsAsFactors = FALSE)
	# colnames(fst) 
	# [1] "lsiga" "lsigb" "lsigw" "fst"   "fis"

	# break the data frame into the bins
	fstBinned <- split(fst, snpBins) 
	
	# Sum Fst components for each bin
	z <- lapply(fstBinned, function(z){
			VARa <- sum(z$lsiga, na.rm = TRUE) # among populations
			VARb <- sum(z$lsigb, na.rm = TRUE) # among individuals within populations
			VARw <- sum(z$lsigw, na.rm = TRUE) # within individuals
			FST <- VARa/(VARa + VARb + VARw)
			return(c(VARa = VARa, VARb = VARb, VARw = VARw, fst = FST))
			})
	# z[1]
	# $`[1,501)`
	# VARa VARb VARw  fst 
	   # 0    0    0  NaN 
	   
	# z[20000]
	# $`[9.9995e+06,1e+07)`
	     # VARa      VARb      VARw       fst 
	# 2.6595041 0.3818182 2.6363636 0.4684134

	z <- as.data.frame(do.call("rbind", z))
	
	# Drop the NaN's
	z$fst[ is.nan(z$fst) ] <- NA
	
	# Put results into the blockstats data frame. 
	blockstats <- cbind.data.frame(blockstats, z)
	rm(z)
	rm(fstBinned)
	rm(fst)

	cat("\nDone Fst calculations\n")
	
	gc()
	
	} # end fst

  			
# ----
# 4) Calculate CSS quantities of interest within bins
  
if( !is.null(psd) ){
  		
	cat("\nBeginning CSS calculations\n")

	gtpairs <- names(psd)
	# gtpairs <- combn(pop, 2, FUN=function(x){paste(x, collapse=",")})
		
	# Eliminate missing pairwise differences by replacing with averages
			
	# 1. Figure out categories for "psdMissingAction" behavior
	
	wbPairs <- rep("w", length(gtpairs)) # initialize
	wbPairs[sapply(strsplit(gtpairs, split=","), function(x){length(unique(x))}) == 2] <- "b"
	# wbPairs
		  # [1] "w" "w" "w" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "b" "b" "b"
		 # [19] "b" "b" "b" "w" "w" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "b"
		 # [37] "b" "b" "b" "b" "b" "w" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w"
		 # [55] "b" "b" "b" "b" "b" "b" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w"
		 # [73] "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "b"
		 # [91] "b" "b" "b" "b" "b" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w"
		# [109] "w" "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w"
		# [127] "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "w" "b" "b" "b"
		# [145] "b" "b" "b" "w" "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w"
		# [163] "w" "w" "w" "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w"
		# [181] "b" "b" "b" "b" "b" "b" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "b"
		# [199] "b" "b" "b" "b" "b" "w" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"
		# [217] "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w"
	
	if(psdMissingAction == "meanBW") psdGroups <- wbPairs else
		if(psdMissingAction == "meanAll") psdGroups <- rep(1, length(wbPairs)) else
			if(psdMissingAction == "meanGroup") psdGroups <- gtpairs else
				stop("Unrecognizeable value for 'psdMissingAction'")


	# 2. Assign average pairwise distance to missing pairwise distance values at every marker separately
	
	cat("\nAssigning averages to missing pairwise distances\n")

	# This is slightly wasteful because it still includes monomorphic sites

	saveColNames <- colnames(psd) # maybe not needed - it was in case of problems with duplicate row names
	
	# transposing costs some memory
	psd <- as.data.frame(t(psd), rownames = NULL, stringsAsFactors = FALSE) # best if snps are columns not rows
	
	# levels(as.factor(psdGroups)) # this will be the order in which the means are given below
	# [1] "b" "w" 
	
	# z <- aggregate(psd, by = list(psdGroups), FUN = mean) # takes forever
	
	# the following seems to be faster than aggregate
	psdGroupMeans <- list()
	for(i in 1:length(levels( as.factor(psdGroups) ))){
		# i <- 1 # "b"
		psdGroupMeans[[i]] <- sapply( psd[psdGroups == levels( as.factor(psdGroups) )[i],], mean, na.rm = TRUE )
		psdGroupMeans[[i]][is.nan(psdGroupMeans[[i]])] <- NA
		}
	names(psdGroupMeans) <- levels(as.factor(psdGroups))
	psdGroupMeans <- as.data.frame(do.call("rbind", psdGroupMeans), stringsAsFactors = FALSE)
	# psdGroupMeans[,1:4]
	  # chrUn:5108405_C/A chrUn:5108407_C/G chrUn:5108413_C/T chrUn:5108418_C/A
	# b                NA        0.08333333                NA                NA
	# w                NA        0.04901961                NA                NA

	# psd[,2] # this column selected for testing because it isn't all NAs
	  # [1]  NA 0.0 0.0  NA 0.0 0.0  NA 0.0  NA 0.0 0.0 0.0 0.0 0.0 0.0  NA 0.5  NA
	 # [19]  NA 0.0 0.0  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA
	 # [37]  NA  NA  NA  NA  NA 0.0  NA 0.0 0.0  NA 0.0  NA 0.0 0.0 0.0 0.0 0.0 0.0
	 # [55]  NA 0.5  NA  NA 0.0 0.0  NA 0.0 0.0  NA 0.0  NA 0.0 0.0 0.0 0.0 0.0 0.0
	 # [73]  NA 0.5  NA  NA 0.0 0.0  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA
	 # [91]  NA  NA  NA  NA  NA 0.0  NA 0.0  NA 0.0 0.0 0.0 0.0 0.0 0.0  NA 0.5  NA
	# [109]  NA 0.0 0.0  NA 0.0  NA 0.0 0.0 0.0 0.0 0.0 0.0  NA 0.5  NA  NA 0.0 0.0
	# [127]  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA  NA 0.0 0.0 0.0
	# [145] 0.0 0.0 0.0  NA 0.5  NA  NA 0.0 0.0  NA  NA  NA  NA  NA  NA  NA  NA  NA
	# [163]  NA  NA  NA 0.0 0.0 0.0 0.0 0.0  NA 0.5  NA  NA 0.0 0.0 0.0 0.0 0.0 0.0
	# [181]  NA 0.5  NA  NA 0.0 0.0 0.0 0.0 0.0  NA 0.5  NA  NA 0.0 0.0 0.0 0.0  NA
	# [199] 0.5  NA  NA 0.0 0.0 0.0  NA 0.5  NA  NA 0.0 0.0  NA 0.5  NA  NA 0.0 0.0
	# [217]  NA  NA  NA  NA  NA  NA  NA 0.5 0.5  NA  NA  NA  NA  NA 0.0

	# position here indicates whether the rows of psd are matched to "b" (1) or "w" (2)
	position <- match(psdGroups, levels(as.factor(psdGroups))) 
	# position    
	  # [1] 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1
	 # [38] 1 1 1 1 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 2 1 1 1 1 1 2 2 2 2 2 2 1 1
	 # [75] 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2
	# [112] 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 2
	# [149] 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1
	# [186] 1 2 2 2 1 1 1 1 1 1 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2
	# [223] 2 2 2 2 2 2 2 2 2

	# Replace the NAs with the group means. This step is slow but hopefully doesn't require too much memory?
	# 	psd[2,][is.na(psd[2,])] <- z[position][is.na(psd[2,])] 
	z <- mapply(psd, psdGroupMeans, FUN = function(x, y){
		# x <- psd[,2]; y <- psdGroupMeans[,2]
		x[is.na(x)] <- y[position][is.na(x)]
		z <- x
		z
		})

	# z[,2] # The NAs above have been replaced with the means
	  # [1] 0.04901961 0.00000000 0.00000000 0.04901961 0.00000000 0.00000000
	  # [7] 0.08333333 0.00000000 0.08333333 0.00000000 0.00000000 0.00000000
	 # [13] 0.00000000 0.00000000 0.00000000 0.08333333 0.50000000 0.08333333
	 # [19] 0.08333333 0.00000000 0.00000000 0.04901961 0.04901961 0.04901961
	 # [25] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.04901961
	 # [31] 0.04901961 0.04901961 0.04901961 0.04901961 0.04901961 0.08333333
	 # [37] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.00000000
	 # [43] 0.04901961 0.00000000 0.00000000 0.08333333 0.00000000 0.08333333
	 # [49] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
	 # [55] 0.08333333 0.50000000 0.08333333 0.08333333 0.00000000 0.00000000
	 # [61] 0.04901961 0.00000000 0.00000000 0.08333333 0.00000000 0.08333333
	 # [67] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
	 # [73] 0.08333333 0.50000000 0.08333333 0.08333333 0.00000000 0.00000000
	 # [79] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.04901961
	 # [85] 0.04901961 0.04901961 0.04901961 0.04901961 0.04901961 0.08333333
	 # [91] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.00000000
	 # [97] 0.04901961 0.00000000 0.04901961 0.00000000 0.00000000 0.00000000
	# [103] 0.00000000 0.00000000 0.00000000 0.04901961 0.50000000 0.04901961
	# [109] 0.04901961 0.00000000 0.00000000 0.04901961 0.00000000 0.04901961
	# [115] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
	# [121] 0.04901961 0.50000000 0.04901961 0.04901961 0.00000000 0.00000000
	# [127] 0.04901961 0.04901961 0.08333333 0.08333333 0.08333333 0.08333333
	# [133] 0.08333333 0.08333333 0.04901961 0.04901961 0.04901961 0.04901961
	# [139] 0.04901961 0.04901961 0.04901961 0.00000000 0.00000000 0.00000000
	# [145] 0.00000000 0.00000000 0.00000000 0.04901961 0.50000000 0.04901961
	# [151] 0.04901961 0.00000000 0.00000000 0.08333333 0.08333333 0.08333333
	# [157] 0.08333333 0.08333333 0.08333333 0.04901961 0.04901961 0.04901961
	# [163] 0.04901961 0.04901961 0.04901961 0.00000000 0.00000000 0.00000000
	# [169] 0.00000000 0.00000000 0.08333333 0.50000000 0.08333333 0.08333333
	# [175] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
	# [181] 0.08333333 0.50000000 0.08333333 0.08333333 0.00000000 0.00000000
	# [187] 0.00000000 0.00000000 0.00000000 0.08333333 0.50000000 0.08333333
	# [193] 0.08333333 0.00000000 0.00000000 0.00000000 0.00000000 0.08333333
	# [199] 0.50000000 0.08333333 0.08333333 0.00000000 0.00000000 0.00000000
	# [205] 0.08333333 0.50000000 0.08333333 0.08333333 0.00000000 0.00000000
	# [211] 0.08333333 0.50000000 0.08333333 0.08333333 0.00000000 0.00000000
	# [217] 0.04901961 0.04901961 0.04901961 0.04901961 0.04901961 0.04901961
	# [223] 0.04901961 0.50000000 0.50000000 0.04901961 0.04901961 0.04901961
	# [229] 0.04901961 0.04901961 0.00000000
	
	psd <- as.data.frame(t(psd), stringsAsFactors = FALSE)

	cat("\nDone assigning averages to missing pairwise distances\n")
	
	gc()

	psdBins <- split(psd, snpBins)
	rm(psd)
  	cat("\nEnded assigning average pairwise distances\n")
		
  	cat("\nDropping monomorphic sites\n")
  	
	# 3. Drop the monomorphic sites to minimize confusion
	
	# Can't get this mapply to work
		# z <- mapply(is.monomorphic, psdBins, FUN = function(i, x){
			# # x <- psdBins[[10000]]; i <- is.monomorphic[[10000]]
			# z <- x#[!i,]
			# return(z)
			# }, SIMPLIFY = FALSE)
	# This worked instead
	z <- list()
	for(i in 1:length(psdBins)){
		z[[i]] <- psdBins[[i]][!is.monomorphic[[i]],]
		}
	psdBins <- z
	rm(z)
	
	gc()

 	cat("\nSumming the psd's within bins\n")
	
	# Sum the psd's within bins or blocks
	# colSum yields a value of 0 if there's no data! Fix? Or just pay attention to nSnp in blockstats
	psdSum <- lapply(psdBins, colSums)
	psdSum <- as.data.frame( do.call("rbind", psdSum), stringsAsFactors = FALSE )

	# psdSum[120:125,]
		          # 2,2       2,2       2,2       2,2      2,1       2,1       2,1
		# 120 4.1061111 3.8134125 5.6071429 4.6071429 6.000000 5.6071429 4.0000000
		# 121 0.7600000 0.6052842 0.6052842 0.6052842 0.750000 0.5238095 0.5863095
		# 122 0.3838589 0.2116101 0.2116101 1.3838589 0.225000 0.3678571 0.2250000
		# 123 0.4913219 2.0363636 0.4913219 0.3286713 1.100000 1.1000000 1.1000000
		# 124 0.3929746 0.9286889 0.9286889 0.3284585 1.963333 2.5633333 2.5633333
		# 125 3.5000000 2.8092589 2.8092589 2.5000000 2.000000 3.7000000 2.5000000
		         # 2,1      2,1       2,2       2,2       2,2        2,2       2,2
		# 120 6.500000 4.868056 6.5149063 4.2113717 5.6220492 6.00000000 4.6000000
		# 121 0.750000 0.750000 0.7600000 0.6052842 0.7600000 0.76000000 0.6052842
		# 122 1.225000 0.725000 0.4747680 0.1207010 0.4747680 0.12070099 1.1207010
		# 123 1.100000 1.100000 1.0363636 0.4913219 1.3286713 0.03636364 0.3286713
		# 124 1.563333 1.463333 0.3284585 0.3284585 0.3284585 0.32845850 0.9286889
		# 125 3.583692 0.500000 3.1202266 2.9264464 3.4200000 1.50000000 3.5490323
		# ....

	names(psdSum) <- paste("psd", names(psdSum), sep = ".")

	blockstats <- cbind.data.frame(blockstats, psdSum, stringsAsFactors = FALSE)
  	cat("\nEnding CSS calculations\n")
	}

# -------------

save(blockstats, file = blockstatsfile)
# load(blockstatsfile) # object is blockstats

