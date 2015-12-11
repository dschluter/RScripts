#!/usr/bin/Rscript

# Combines variant and invariant information in genome blocks of size 500 (default) nucleotides

# qsub -I -l walltime=02:00:00 -l mem=4gb 
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

project <- args[1]
chrname <- args[2]
stepsize <- as.integer(args[3])
groupnames <- args[-c(1:3)]
chrno <- gsub("chr", "", chrname)

cat("\nProject, chrname, groupnames, stepsize", project, chrname, groupnames, stepsize, "\n")

gtstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "."), "gtstats.rdd", sep = ".")	# object is gtstats
load(gtstatsfile) # object is gtstats
# names(gtstats)
 # [1] "vcf"               "groupnames"        "groups"           
 # [4] "nInd"              "nMin"              "control"          
 # [7] "genotypes"         "alleleFreqByGroup" "status"           
# [10] "newPos"            "fst"               "psd"

psdMissingAction <- "meanBW" 
	# psdMissingAction is for g$blockstats, how to average when psd values are missing. Must be one of the following:
	# 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
	# 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
	#   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
	#		identified by missingGroups.
	#		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
	#		"1,1" and "2,2" as the same (i.e., both are within-group).

control <- gtstats$control
Glazerize <- control$Glazerize
groupnames <- gtstats$groupnames
nInd <- gtstats$nInd 	# n individuals genotyped in each group eg 7 11
nMin <- gtstats$nMin 	# criterion is 5 or nMin, whichever is smaller, eg 4 5
control$psdMissingAction <- psdMissingAction

if(Glazerize){
	goodInvariantsFile <- paste(project, ".", chrname, ".goodInvNew.rdd", sep="")
	chrvecfile = paste("chrvec.", chrno, ".glazer.rdd", sep = "")
	} else {
	goodInvariantsFile <- paste(project, ".", chrname, ".goodInv.rdd", sep="") # object is goodInvariants
	chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")
	}
load(goodInvariantsFile) # object is goodInvariants
load(chrvecfile) # object is chrvec

blockstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "."), 
							"blockstats", stepsize, "rdd", sep = ".")

# Fix the names in goodInvariants (change "marine.pac" to "marine-pac", etc)
names(goodInvariants) <- gsub("[.]", "-", names(goodInvariants))

# Rename newPos to POS if using Glazer reassembly - no, keep them separate
# goodInvariants$POS <- goodInvariants$newPos

# Grab the columns of interest
goodInvariants <- goodInvariants[, c("newPos", "POS", "REF", groupnames)]

# Drop invariants that do not meet the minimum number of genotypes required
z <- apply(goodInvariants[, groupnames], 1, function(x){
	all(x >= nMin)
	})
goodInvariants <- goodInvariants[z, ]

gc()


# ----------
# Get block window stats for group differences

blockstats <- g$blockstats(gtstats, stepsize, goodInvariants = goodInvariants, chrvec = chrvec, control = control) 

save(blockstats, file = blockstatsfile)
# load(blockstatsfile)

