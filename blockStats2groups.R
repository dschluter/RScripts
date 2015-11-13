#!/usr/bin/Rscript

# Combines variant and invariant information in genome blocks of size 500 (default) nucleotides

# qsub -I -l walltime=02:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# stepsize <- 500     # Number of nucleotides in a genome block - now passed as an argument
psdMissingAction <- "meanBW" 
	# psdMissingAction is for g$blockstats, how to average when psd values are missing. Must be one of the following:
	# 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
	# 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
	#   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
	#		identified by missingGroups.
	#		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
	#		"1,1" and "2,2" as the same (i.e., both are within-group).

# ---

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimAllMarine", "chrXXI", "marine-pac", "paxb")
# args <- c("BenlimAllMarine", "chrXXI", "marine-pac", "paxl")
# args <- c("BenlimAllMarine", "chrXXI", 500, "marine-pac", "solitary")

project <- args[1]
chrname <- args[2]
stepsize <- as.integer(args[3])
groupnames <- args[-c(1:3)]

# open object containing genotype statistics
gtstatsfile 		<- paste(project, chrname, paste(groupnames, collapse = "."), 
							"rdd", sep = ".") 								# object gtstats
goodInvariantsFile 	<- paste(project, ".", chrname, ".goodInv.rdd", sep="") # object goodInvariants
blockstatsfile 		<- paste(project, chrname, paste(groupnames, collapse = "."), 
							"blockstats", stepsize, "rdd", sep = ".")

load(gtstatsfile) # object is gtstats
load(goodInvariantsFile) # object is goodInvariants

groupnames <- gtstats$groupnames

# process the invariants - keep rows meeting the nMin rule
nInd <- gtstats$nInd 	# n individuals genotyped in each group eg 7 11
nMin <- gtstats$nMin 	# criterion is 5 or nMin, whichever is smaller, eg 4 5

# Fix the names in goodInvariants (change "marine.pac" to "marine-pac", etc)
names(goodInvariants) <- gsub("[.]", "-", names(goodInvariants))

goodInvariants <- goodInvariants[, c("POS", "REF", groupnames)]
z<-apply(goodInvariants[, groupnames], 1, function(x){
	all(x >= nMin)
	})
goodInvariants <- goodInvariants[z, ]

# get the chrvec file name for blockstats
chrno <- gsub("chr", "", chrname)
chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")

gc()


# ----------
# Get block window stats for group differences

blockstats <- g$blockstats(gtstats, stepsize, goodInvariants = goodInvariants, 
					chrvecfile, psdMissingAction = psdMissingAction) 

save(blockstats, file = blockstatsfile)
# load(blockstatsfile)

