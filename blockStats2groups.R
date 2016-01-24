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
# args <- c("BenlimAllMarine", "chrXX", 500, "pril", "prib") # one of the smaller ones that crashed because of memory
# args <- c("BenlimAllMarine", "chrM", 500, "paxl", "paxb")
# args <- c("BenlimAllMarine", "chrXXI", 500, "paxl", "paxb")

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

cat("\nLoading gtstats file\n")

load(gtstatsfile) # object is gtstats
# names(gtstats)
 # [1] "vcf"        "groupnames" "groups"     "nInd"       "nMin"      
 # [6] "control"    "genotypes"  "status"     "newPos"     "fst"       
# [11] "psd"      

# object.size(gtstats) 
# 3111541328 bytes # 3Gb wow (chrXXI paxl paxb) - will be less when not including alleleFreqByGroups (newer version)

# library(data.table)

control <- gtstats$control
Glazerize <- control$Glazerize
groups <- gtstats$groups
groupnames <- gtstats$groupnames
nInd <- gtstats$nInd 	# n individuals genotyped in each group eg 7 11
nMin <- gtstats$nMin 	# criterion is 5 or nMin, whichever is smaller, eg 4 5
# control$psdMissingAction <- psdMissingAction # not used right now
status <- gtstats$status
fst <- gtstats$fst # will be NULL if absent
# library(data.table)
psd <- gtstats$psd # will be NULL if absent
# psd <- as.data.table(gtstats$psd) # will be NULL if absent


gc()
            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  10093132  539.1   15229393  813.4  10101556  539.5
# Vcells 255284629 1947.7  274105541 2091.3 255309030 1947.9

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
cat("\nRemoving gtstats file from memory after extracting relevant elements\n")

# head(pos)
# [1] 1794949 1794951 1794957 1794962 1794969 1794974

# All start positions are unique, snps and indels at the same start POS are just different alleles
# c( length(pos), length(unique(pos)) )
# [1] 871120 871120


# -------------
# Get block window stats for group differences

pop <- as.character(groups)
k <- as.integer(stepsize) # this is the size of the block
		
cat("\nLoading masked genome\n")

# chrvec is a vector of the whole MASKED chromosome, all nucleotides as distinct elements, obtained by:
#       Masked bases are indicates with an "M"
load(chrvecfile) # object is chrvec

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
rm(chrbins)

blockstats <- data.frame(
			ibase = ibase[-length(ibase)], 
			midbase = midbase[-length(midbase)], 
			nUnmasked = nUnmasked
			)
rm(midbase)
rm(nUnmasked)

gc()
            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells   2247472  120.1    6237958  333.2  10162262  542.8
# Vcells 214258536 1634.7  333521806 2544.6 317563341 2422.9

# ----
# 2) Count up the number of snp and number of invariants in each bin

cat("\nLoading good invariants file\n")

load(goodInvariantsFile) # object is goodInvariants

# object.size(goodInvariants) 
# 1560947856 bytes # 1.6 Gb

cat("\nFinished loading good invariants file\n")

gc()
            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  13724133  733.0   14562966  777.8  13724272  733.0
# Vcells 359852242 2745.5  518283907 3954.2 359852403 2745.5


# Fix the names in goodInvariants (change "marine.pac" to "marine-pac", etc)
names(goodInvariants) <- gsub("[.]", "-", names(goodInvariants))
	# setnames(goodInvariants, names(goodInvariants), gsub("[.]", "-", names(goodInvariants)))

if(Glazerize){
	goodInvariants <- goodInvariants[, c("newPos", groupnames)] # grab the columns of interest
	# rename "newPos" to "POS"
	names(goodInvariants)[names(goodInvariants) =="newPos"] <- "POS"
	} else{
	goodInvariants <- goodInvariants[, c("POS", groupnames)] # grab the columns of interest
	}

# gc()
            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells  13724116  733.0   15331114  818.8  13751976  734.5
# Vcells 273770656 2088.8  518283907 3954.2 359931884 2746.1
	
cat("\nDropping invariants not meeting the minimum number of genotypes\n")

# Drop invariants that do not meet the minimum number of genotypes required
# names(goodInvariants)
# [1] "POS"  "paxl" "paxb" # counts needed are in columns 2 and 3

i <- goodInvariants[,2] >= nMin[1] & goodInvariants[,3] >= nMin[2]
goodInvariants <- goodInvariants[i, ]
rm(i)

# head(goodInvariants)
                  # POS paxl paxb
# chrUn.3907275 1794692    5    7
# chrUn.3907276 1794693    5    7
# chrUn.3907277 1794694    5    7
# chrUn.3907278 1794695    5    7
# chrUn.3907279 1794697    5    7
# chrUn.3907280 1794698    5    7

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
            # used   (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells   2318481  123.9   12264891  655.1  14568721  778.1
# Vcells 223260813 1703.4  518283907 3954.2 515655976 3934.2


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


	# This is slightly wasteful because it still includes monomorphic sites
	
		# if(FALSE){ # Data table way - requires more memory than data frame way
		
		# library(data.table) # doesn't use as much memory if do here rather than use "as.data.table(gtstats$psd)" at start
		# psd <- as.data.table(psd)
		
		# cat("\nConverted psd to data table\n")
		# gc()
		            # # used   (Mb) gc trigger   (Mb)  max used   (Mb)
		# # Ncells   2411082  128.8    8262486  441.3  14189900  757.9
		# # Vcells 230305282 1757.1  551359436 4206.6 433308962 3305.9
		
		# psdUnique <- unique(names(psd)) # this assumes that we have many duplicate names -- haven't changed
		# # [1] "w" "b"
		
		# for(p in psdUnique){
			# # p <- psdUnique[1]
			
			# doCols <- c(grep(p, names(psd))) # which columns are of interest on this iteration of loop?
			# # doCols
			  # # [1]   1   2   3   4  10  11  12  13  14  15  22  23  24  30  31  32  33  34  35  42  43  49  50  51  52  53
			 # # [27]  54  61  67  68  69  70  71  72  84  85  86  87  88  89  96  97  98  99 106 107 108 109 110 111 112 113
			 # # [53] 114 121 122 123 124 125 126 127 128 135 136 137 138 139 140 141 148 149 150 151 152 153 160 161 162 163
			 # # [79] 164 165 166 167 168 169 170 177 178 179 180 187 188 189 196 197 204 217 218 219 220 221 222 223 224 225
			# # [105] 226 227 228 229 230 231
			
			# saveNames <- names(psd)[doCols]
			
			# # Give the columns a temporary standard name to help next step
	
			# setnames(psd, doCols, paste("column", doCols, sep = ""))
			# # names(psd)
			  # # [1] "column1"   "column2"   "column3"   "column4"   "b"         "b"         "b"         "b"        
			  # # [9] "b"         "column10"  "column11"  "column12"  "column13"  "column14"  "column15"  "b"        
			 # # [17] "b"         "b"         "b"         "b"         "b"         "column22"  "column23"  "column24" 
				# # ....		
			# useCols <- names(psd)[doCols] # [1] "column1"   "column2"   "column3"   "column4"   "column10" ....
			
			# psd[ , rowMean := apply(.SD, 1, mean, na.rm=TRUE),, .SDcols = doCols ]
			# psd[ is.nan(rowMean ), `:=`(rowMean = NA)]
	
			# # gc()
			            # # used   (Mb) gc trigger   (Mb)  max used   (Mb)
			# # Ncells   2426601  129.6    8262486  441.3  14189900  757.9
			# # Vcells 231207554 1764.0  670525648 5115.8 670253610 5113.7 # whoa! this is from p <- psdUnique[1]
	
			# for(k in useCols){
				# # k <- useCols[1]
				# i <- parse( text = paste("is.na(", k, ")") )
				# j <- parse( text = paste("`:=`(", k, "= rowMean)") )
				# psd[eval(i), eval(j)]
				# }
			# setnames(psd, doCols, saveNames)
			
			# }
	
		# # Delete rowMean
		# psd[ , c("rowMean") := NULL]
		
		# gc()
		            # # used   (Mb) gc trigger   (Mb)  max used   (Mb)
		# # Ncells   2429347  129.8    8262486  441.3  14189900  757.9
		# # Vcells 231211996 1764.1  704131930 5372.2 703966666 5370.9 # whoa!
		
		# psd$snpBins <- snpBins
	
		# psdSum <- psd[ , lapply(.SD, sum, na.rm = TRUE), keyby = snpBins] # keyby orders the snpBins but adds column
			# # head(psdSum)
		         # # snpBins         w         w         w         w          b          b
		# # 1: [11501,12001) 0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000
		# # 2: [12001,12501) 6.8913043 6.8913043 6.8913043 6.8913043 7.40000000 7.40000000
		# # 3: [12501,13001) 1.2434409 1.2434409 1.2434409 1.2434409 4.06031746 4.06031746
		# # 4: [15001,15501) 0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000
		# # 5: [16001,16501) 0.2585139 0.2585139 0.2585139 0.2585139 0.24735450 0.24735450
		# # 6: [42501,43001) 0.0000000 0.0000000 0.0000000 0.0000000 0.07142857 0.07142857
		# # .....
	
		# # nrow(psdSum) 
		# # [1] 30342
		# # nrow(blockstats)
		# # [1] 34716
		
		# ** Not finished yet **
		# ** Need to match column 1 of psdSum to blockstats **
	
		# blockstats <- cbind.data.frame(blockstats, psdSum, stringsAsFactors = FALSE)
		
		# rm(psd)
		# rm(psdSum)
		# gc()
	
		# } # End modified data.table way
	

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

