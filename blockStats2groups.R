#!/usr/bin/Rscript

# Calculates preliminary results for sliding window analyses on a (repeat-masked) chromosome
# It breaks the chromosome into blocks or bins of size stepsize and calculates a value of interest in each bin
# Combines variant and invariant information in genome blocks of size 500 (default) nucleotides

# gtstats is the gtstats list containing the variables of interest for a single chromosome
# Method asks whether gtstats$fst and gtstats$psd exist and computes and saves blockstats accordingly.

# GATK 3.4 no longer gives multiple rows of snps at the same value of POS

# Filter based on FILTER
# NULL means all, others can be chained, e.g., "PASS,." to include "PASS" and "."
# NULL, LowQual, PASS, VQSRTrancheSNP99.00to99.50, VQSRTrancheSNP99.50to99.90, VQSRTrancheSNP99.50to99.90+
# filter 	 	<- "PASS,." # PASSandDOT result files
# filter		<- "PASS"	# PASS result files - for now include everybody, and filter later

# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# fishPair must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

chrname 		<- NULL
project 		<- NULL
fishPair		<- NULL
stepsize 		<- NULL
psdMissingAction<- NULL 
# genomeDir		<- NULL 

args <- commandArgs(TRUE)
# args <- c("project=Benlim","chrname=chrXXI","psdMissingAction=meanBW","stepsize=500","filter=PASS,.","fishPair=paxl,paxb")

# Parses the args into a data frame with two columns and then assigns variables 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}
print(x)
                # V1        V2
# 1          project    Benlim
# 2          chrname    chrXXI
# 3 psdMissingAction    meanBW
# 4         stepsize       500
# 5           filter    PASS,.
# 6         fishPair paxl,paxb

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")
# if(is.null(genomeDir)) stop("Provide genomeDir= in arguments")
if(is.null(stepsize)) stop("Provide stepsize= in arguments")
if(is.null(psdMissingAction)) stop("Provide psdMissingAction= in arguments")
if(is.null(fishPair)) stop("Provide fishPair= in arguments")

# psdMissingAction is for g$blockstats, how to average when psd values are missing. Must be one of the following:
# 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
# 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
#   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
#		identified by missingGroups.
#		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
#		"1,1" and "2,2" as the same (i.e., both are within-group).

if( !is.element(psdMissingAction, c("meanAll", "meanBW", "meanGroup")) )
	stop("psdMissingAction must be one of meanAll, meanBW, or meanGroup")

chrno <- gsub("chr", "", chrname)

cat("\nProject, chrname, fishPair, stepsize, psdMissingAction\n", 
	project, chrname, fishPair, stepsize, psdMissingAction, "\n")

vcfparamFile <- paste(project, "vcfparam.rdd", sep = ".")
load(vcfparamFile) # object is vcfparam

Glazerize <- vcfparam$Glazerize
GTminFrac <- vcfparam$GTminFrac
fishnames <- vcfparam$fishnames
groupcodes<- vcfparam$groupcodes
groupnames<- vcfparam$groupnames
nMaxAlt <- vcfparam$nMaxAlt

fishPair <- unlist(strsplit(fishPair, split = ","))
if(length(fishPair) > 2 ) stop("Provide names of only two groups")

nInd <- vcfparam$nInd[fishPair]
nMin <- vcfparam$nMin[fishPair]

gtstatsfile 	<- paste(project, chrname, paste(fishPair, collapse = "."), "gtstats.rdd", sep = ".")
blockstatsfile 	<- paste(project, chrname, paste(fishPair, collapse = "."), "blockstats", stepsize, "rdd", sep = ".")

cat("\nLoading gtstats file\n")

load(gtstatsfile) # object is gtstats
# names(gtstats)
 # [1] "fishPair"    "groups"      "genotypes"   "status"      "newPos"     
 # [6] "altUsedList" "snpTypeList" "FILTER"      "fst"         "psd"        

# object.size(gtstats) 
# 3,127,131,704 bytes # 3Gb wow (chrXXI paxl paxb)

# library(data.table)

groups 	<- gtstats$groups
status 	<- gtstats$status
fst 		<- gtstats$fst # will be NULL if absent
psd 		<- gtstats$psd # will be NULL if absent
FILTER 	<- gtstats$FILTER


# library(data.table)
# psd <- as.data.table(gtstats$psd) # will be NULL if absent

gc()


if(Glazerize){
	goodInvariantsFile <- paste(project, chrname, "goodInvNew.rdd", sep=".")
	chrvecfile = paste("chrvec", chrno, "glazer.rdd", sep = ".")
	pos <- gtstats$newPos # pos is for the variants, use POS for the invariants file
	
	} else {
		
	# *Need to confirm that this works if not glazerizing
	goodInvariantsFile <- paste(project, ".", chrname, ".goodInv.rdd", sep="") # object is goodInvariants
	chrvecfile = paste("chrvec", chrno, "masked.rdd", sep = ".")
	pos <- gtstats$pos # pos is for the variants, use POS for the invariants file
	}

# Filter the SNPs using FILTER
table(FILTER)
                          # .                     LowQual 
                     # 204698                        9657 
                       # PASS  VQSRTrancheSNP99.00to99.50 
                     # 466788                       53132 
 # VQSRTrancheSNP99.50to99.90 VQSRTrancheSNP99.50to99.90+ 
                     # 120577                       31566 

filter	<- unlist(strsplit(filter, split=","))
filter
# [1] "PASS" "." 
if(all(!is.na(filter) & filter != "NULL")){
	keep <- casefold(FILTER) %in% casefold(filter)
	 # FALSE   TRUE 
	# 214932 671486 
	 
	status 	<- status[keep]
	fst 		<- fst[keep,]
	psd 		<- psd[keep,]
	FILTER 	<- FILTER[keep]
	pos 		<- pos[keep]
	rm(keep)
 	}

rm(gtstats)
cat("\nRemoving gtstats file from memory after extracting relevant elements\n")

# head(pos)
# [1] 1510759 1512807 1513124 1513272 1513334 1547865

# All start positions are unique, snps and indels at the same start POS are just different alleles
# c( length(pos), length(unique(pos)) )
# [1] 671486 671486


# -------------
# Get block window stats for group differences

pop <- as.character(groups)
k <- as.integer(stepsize) # this is the size of the block
		
cat("\nLoading masked genome\n")

# chrvec is a vector of the whole MASKED chromosome, all nucleotides as distinct elements, obtained by:
#       Masked bases are indicates with an "M"
load(chrvecfile) # object is chrvec
print(length(chrvec))
# [1] 17357772

# Establish the break points of the bins into which nucleotides will be grouped (e.g., k = 500 bases per bin)
# The last bin goes from the final bin of fully k nucleotides to the last nucleotide. 
# This last bin might be small.
	
nbases <- length(chrvec)
ibase <- c( seq(1, nbases, k), nbases + 1) # bases marking breaks between steps of size k
midbase <- ibase + (k)/2 # If want something to indicate midpoint of blocks instead

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
rm(chrbins)

blockstats <- data.frame(
			ibase = ibase[-length(ibase)], 
			midbase = midbase[-length(midbase)], 
			nUnmasked = nUnmasked
			)
rm(midbase)
rm(nUnmasked)

gc()


# ----
# 2) Count up the number of snp and number of invariants in each bin

cat("\nLoading good invariants file\n")

load(goodInvariantsFile) # object is goodInvariants

# object.size(goodInvariants) 
# 1652726128 bytes # 1.65 Gb

cat("\nFinished loading good invariants file\n")

gc()

# Fix the names in goodInvariants (eg change "marine.pac" to "marine-pac", etc)
names(goodInvariants) <- gsub("[.]", "-", names(goodInvariants))

if(Glazerize){
	goodInvariants <- goodInvariants[, c("newPos", fishPair)] # grab the columns of interest
	# rename "newPos" to "POS"
	names(goodInvariants)[names(goodInvariants) =="newPos"] <- "POS"
	} else{
	goodInvariants <- goodInvariants[, c("POS", fishPair)] # grab the columns of interest
	}
	
cat("\nDropping invariants not meeting the minimum number of genotypes\n")

# Drop invariants that do not meet the minimum number of genotypes required
# names(goodInvariants)
# [1] "POS"  "paxl" "paxb" # counts needed are in columns 2 and 3

i <- goodInvariants[,2] >= nMin[1] & goodInvariants[,3] >= nMin[2]
goodInvariants <- goodInvariants[i, ]
rm(i)

head(goodInvariants)
                  # POS paxl paxb
# chrUn.3918179 1793078    6   11
# chrUn.3918180 1793079    6   11
# chrUn.3918181 1793082    6   11
# chrUn.3918182 1793083    6   11
# chrUn.3918183 1793085    6   11
# chrUn.3918184 1793086    6   11

cat("\nBinning the good invariants\n")

# Count up the number of good invariants in each bin
# Remember: newPos has been renamed to POS if Glazerize is true

invariantBins <- cut(goodInvariants$POS, breaks = ibase, right = FALSE)
nInvariants <- table( invariantBins ) 
rm(goodInvariants)
rm(invariantBins)

# head(nInvariants)
    # [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001)  
    #       0           0           0           0           0           0 

# Count up the number of snp in each bin (cut works because ibase intervals are made as factors)
# "pos" refers to the positions of the snp

cat("\nBinning the good snps\n")

snpBins <- cut(pos, breaks = ibase, right = FALSE)
# is.factor(snpBins)
# [1] TRUE

nSnp <- table( snpBins )

# head(nSnp)
	# [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001) 
	#       0           0           0           0           0           0 

is.monomorphic <- split(status == "i", snpBins)
nMonomorphic <- sapply(is.monomorphic, sum)

# rm(status)

blockstats$nSnp = as.vector(nSnp) - nMonomorphic
blockstats$nInvariants = as.vector(nInvariants) + nMonomorphic

rm(ibase)
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
# 117 58001   58251        99    0          56
# 118 58501   58751       193   12         142
# 119 59001   59251        89    0          79
# 120 59501   59751       411    4         364
# 121 60001   60251       398    0         183
# 122 60501   60751       343    0         213
# 123 61001   61251       314    0         258
# 124 61501   61751       456    0         229
# 125 62001   62251       336    1         104
# 126 62501   62751         4    0           0
# 127 63001   63251         8    0           0
# 128 63501   63751       106    0           0
# 129 64001   64251       283    0           0
# 130 64501   64751       363    0           0
	
# ----
# 3) Calculate Fst summary stats

if( !is.null(fst) ){
		
	cat("\nBinning the Fst values\n")

	# fst is a matrix, needs to be a data.frame or split will screw up

	fst <- as.data.frame(fst, stringsAsFactors = FALSE)

	# colnames(fst) 
	# [1] "lsiga" "lsigb" "lsigw" "fst"   "fis"

	# break the data frame into the bins
	# Note: snpBins is a factor, so the resulting list is in the correct order
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
	
	# Change NaN's to NA
	z$fst[ is.nan(z$fst) ] <- NA
	
	# Put results into the blockstats data frame. 
	blockstats <- cbind.data.frame(blockstats, z)
	rm(z)
	rm(fstBinned)
	rm(fst)

	cat("\nDone Fst binning\n")
	
	gc()
	            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells   2388682  127.6    9320297  497.8  14562966  777.8
	# Vcells 218174273 1664.6  518283907 3954.2 518191814 3953.5

	} # end fst

  			
# ----
# 4) Calculate CSS quantities of interest within bins
  
if( !is.null(psd) ){
  		
	cat("\nCalculating quantities for CSS score\n")
	# psd is a data frame
	
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
	
	# Identify the columns of psd according to psdMissingAction
	# 
	if(psdMissingAction == "meanBW") psdGroups <- wbPairs else
		if(psdMissingAction == "meanAll") psdGroups <- rep("a", length(wbPairs)) else
			if(psdMissingAction == "meanGroup") psdGroups <- gtpairs else
				stop("Unrecognizeable value for 'psdMissingAction'")

	names(psd) <- psdGroups  
	# names(psd) # duplicate column names seem to be ok
		
	# head(psd)
	# Note that many loci have all NA's even though alleles with insufficient numbers of individuals genotyped
	# were dropped in "gtstats2groups.R". However, psd was calculated and saved only for loci with status == "v",
	# which means that psd for all the invariants are NA rather than 0.


	psdUnique <- unique(names(psd))
	# [1] "w" "b"
	
	# Try reducing psd just to variant rows
	# nrow(psd)

	cat("\nReducing psd just to variant cases\n")
	psd <- psd[status != "i", ]

	# nrow(psd)

	cat("\nAssigning averages to missing pairwise distances\n")

	for(p in psdUnique){
		# p <- psdUnique[1]
		
		doCols <- c(grep(p, names(psd))) # which columns are of interest on this iteration of loop?
		# doCols
		  # [1]   1   2   3   4  10  11  12  13  14  15  22  23  24  30  31  32  33  34  35  42  43  49  50  51  52  53
		 # [27]  54  61  67  68  69  70  71  72  84  85  86  87  88  89  96  97  98  99 106 107 108 109 110 111 112 113
		 # [53] 114 121 122 123 124 125 126 127 128 135 136 137 138 139 140 141 148 149 150 151 152 153 160 161 162 163
		 # [79] 164 165 166 167 168 169 170 177 178 179 180 187 188 189 196 197 204 217 218 219 220 221 222 223 224 225
		# [105] 226 227 228 229 230 231
		
		rowMean <- apply(psd[ , doCols], 1, mean, na.rm=TRUE)

		# THese steps carried out before psd reduced to variant rows:
		# Check that all NaN's correspond to status = "i"
		# all( status[is.nan(rowMean)] == "i" )
		# [1] TRUE
		# all( status[!is.nan(rowMean)] == "v" )
		# [1] TRUE
		
		cat("\nTest that no Nan's remain in rowMeans (TRUE is correct)\n")
		print( length(rowMean[is.nan(rowMean)]) == 0 )

		# rowMean[is.nan(rowMean)] <- NA
		
		for(k in doCols){
			# k <- doCols[1]
			# unname(rowMean[is.na(psd[ , k])])
			psd[is.na(psd[ , k]) , k] <- rowMean[is.na(psd[ , k])]
			}			
		}

	# Delete rowMean
	rm(rowMean)
	
	cat("\nDone assigning averages to missing pairwise distances\n")

	gc()
	           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells  1764518  94.3    5964989  318.6  14562966  777.8
	# Vcells 70026999 534.3  331701700 2530.7 518283862 3954.2
	
	
	cat("\nBinning psd by snpBins\n")
	# snpBins is a factor so the resulting list will be in the correct order:

	# Modify if we analyzed only rows with variants
	# psdBins <- split(psd, snpBins)

	psdBins <- split(psd, snpBins[status != "i"])
		
	rm(psd)
	
	print(gc())
	           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells  9957363 531.8   13263914  708.4  14562966  777.8
	# Vcells 90438778 690.0  265361360 2024.6 518283862 3954.2

	cat("\nDone binning psd values\n")

	psdSum <- lapply(psdBins, colSums, na.rm=TRUE)
	# length(psdSum)
	# [1] 34716
	# nrow(blockstats)
	# [1] 34716

	psdSum <- do.call("rbind", psdSum)
	colnames(psdSum) <- paste("psd", gtpairs, sep = ".")
	
	# psdSum[120:125,]
	                # psd.2,2   psd.2,2   psd.2,2   psd.2,2  psd.2,1   psd.2,1
	# [59501,60001) 4.1061111 3.8134125 5.6071429 4.6071429 6.000000 5.6071429
	# [60001,60501) 0.7600000 0.6052842 0.6052842 0.6052842 0.750000 0.5238095
	# [60501,61001) 0.3838589 0.2116101 0.2116101 1.3838589 0.225000 0.3678571
	# [61001,61501) 0.4913219 2.0363636 0.4913219 0.3286713 1.100000 1.1000000
	# [61501,62001) 0.3929746 0.9286889 0.9286889 0.3284585 1.963333 2.5633333
	# [62001,62501) 3.5000000 2.8092589 2.8092589 2.5000000 2.000000 3.7000000
	                # psd.2,1  psd.2,1  psd.2,1   psd.2,2   psd.2,2   psd.2,2
	# [59501,60001) 4.0000000 6.500000 4.868056 6.5149063 4.2113717 5.6220492
	# [60001,60501) 0.5863095 0.750000 0.750000 0.7600000 0.6052842 0.7600000
	# [60501,61001) 0.2250000 1.225000 0.725000 0.4747680 0.1207010 0.4747680
	# [61001,61501) 1.1000000 1.100000 1.100000 1.0363636 0.4913219 1.3286713
	# [61501,62001) 2.5633333 1.563333 1.463333 0.3284585 0.3284585 0.3284585
	# [62001,62501) 2.5000000 3.583692 0.500000 3.1202266 2.9264464 3.4200000

	blockstats <- cbind.data.frame(blockstats, psdSum, stringsAsFactors = FALSE)
	rm(psdSum)
	
  	cat("\nEnding CSS calculations\n")
  	
  	# gc()
	  	            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells  10588039  565.5   15375065  821.2  15375065  821.2
	# Vcells 244937212 1868.8  670525648 5115.8 578877838 4416.5

	}

# -------------

save(blockstats, file = blockstatsfile)
# load(blockstatsfile) # object is blockstats

