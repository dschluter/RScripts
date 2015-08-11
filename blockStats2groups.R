#!/usr/bin/Rscript
# Combines variant and invariant information in genome blocks of size 500 (default) nucleotides
# Files required in preparation for this script:
# 	genome.r									in scriptdir
#	chrvec.chrno.masked.rdd					in vcfdir
#	project.chrname.vcfresults.rdd			in vcfdir
#	project.chrname.processedInvariants.rdd	in vcfdir

# Run in Unix
# as " Rscript blockStats2groups.R chrno project blocksize psdMissingAction FSTtrueSnps CSStrueSnps & "
#    " Rscript ../scripts/blockStats2groups.R 21 "Paxton12fish" 500 meanAll FALSE TRUE >&Rscript.out & "

# Defaults - can be changed by including an alternative value as an argument IN THE RIGHT PLACE
blocksize <- 500     # Number of nucleotides in a genome block
FSTtrueSnps <- FALSE # if TRUE, Fst stats include only *true snp*, not indels; if FALSE, everything is included 
					 # (in which case a given POS might be counted twice, as a snp and an indel)
CSStrueSnps <- TRUE  # if TRUE, CSS stats only included *true snp*, not indels
psdMissingAction <- "meanAll"   
	# psdMissingAction is for g$blockstats, how to average when psd values are missing. Must be one of the following:
	# 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
	# 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
	#   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
	#		identified by missingGroups.
	#		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
	#		"1,1" and "2,2" as the same (i.e., both are within-group).


# Make sure that the latest version of genome.r is copied to scriptdir folder before sourcing this file

# ---
args <- commandArgs(TRUE)

# args <- c(21, "Paxton12fish", 500)

chrno <- as.character(as.roman(args[1]))
chrname <- paste("chr", chrno, sep="")
project <- args[2]
if(!is.na(args[3])) blocksize <- as.integer(args[3])
if(!is.na(args[4])) psdMissingAction <- args[4]
if(!is.na(args[5])) FSTtrueSnps <- as.logical(args[5])
if(!is.na(args[6])) CSStrueSnps <- as.logical(args[6])

cat("BlockStats settings:\n")
cat("\nChromosome name:", chrname, "\n")
cat("\npsdMissingAction:", psdMissingAction, "\n")
cat("\nFSTtrueSnps:", FSTtrueSnps, "\n")
cat("\nCSStrueSnps:", CSStrueSnps, "\n")


# westgrid is the default, otherwise this doesn't work in Rscript mode
#if( grepl("westgrid",Sys.info()["nodename"]) ){
	setwd("~")
	vcfdir <- ""
	genomedir <- ""
	invariantsummarydir <- vcfdir
	scriptdir <- ""
#	}

# To run on office computer (copy important files to the local genomics folder, don't access files on server - too slow)
if( grepl("greendrake",Sys.info()["nodename"]) ){
	setwd("/Users/schluter/Documents/genomics/")
	vcfdir <- ""
	genomedir <- ""
	invariantsummarydir <- vcfdir
	scriptdir <- "/Volumes/schluter/Dolph_data/R and S-Plus macros/"
	}

if( grepl("Wallacea",Sys.info()["nodename"]) ){
	setwd("/Users/schluter/Documents/Research/genomics/BenLim-q15/")
	vcfdir <- ""
	genomedir <- "/Users/schluter/Documents/Research/genomics/reference genomes/"
	invariantsummarydir <- vcfdir
	scriptdir <- "/Users/schluter/Documents/Research/R and SPlus macros/"
	}

if( grepl("SciBorg",getwd()) ){
	setwd("~/BenLim-q15/") 
	vcfdir <- ""
	genomedir <- "~/reference_genomes/"
	invariantsummarydir <- vcfdir
	scriptdir <- "~/scripts/"
	# folder <- recode(project, c("benlim", "paxbl", "priestbl", "qrybl", "enosbl"), 
							# c("BenLim-q15","Paxton-q15","Priest-q15","LQuarry-q15", "Enos-q15"))
	# setwd( paste("~", folder, "", sep="/") )  #
	# vcfdir <- paste("~", folder, "", sep="/") # 
	}
	
source(paste(scriptdir,"genome.r", sep = ""))

# Load masked "chrvec" object from .rdd file in vcfdir
load(file = paste(vcfdir, "chrvec.", chrno, ".masked.rdd", sep = "")) # chrvec
# head(chrvec)

# Load "vcfresults", the vcf variant results, from vcfdir
load(file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = "")) # vcfresults

# Load "processedInvariants", giving the chromosome positions meeting the DPmin criterion, from vcfdir
load(file = paste(vcfdir, project, ".", chrname, ".processedInvariants.rdd", sep="")) # processedInvariants
head(processedInvariants)
        # POS REF nGTgroup1 nGTgroup2
# 18308 18386   G         6         5
# 18309 18387   T         6         5
# 54659 54836   C         4         4
# 54660 54837   T         4         4
# 54661 54838   C         4         5
# 54662 54839   C         4         5

# names(vcfresults)
 # [1] "rowdata"      "info"         "ALTlist"      "GT"           "AD"           "DP"          
 # [7] "GQ"           "PL"           "snptypeList"  "is.trueSnp"   "loc"          "coding"      
# [13] "freq.table"   "allele.count" "pfisher"      "fst"          "psd"         

# Show the minimum number of genotypes required for each group
cat("\nSpecimens analyzed:\n")
colnames(vcfresults$GT)
cat("\nGroups:\n")
print(vcfresults$groups)
cat("\nVariant genotypes: Minimum number present in each group:\n")
print(vcfresults$nMin)
cat("\nInvariant genotypes: Minimum number present in each group:\n")
print(apply(processedInvariants[,-c(1:2)], 2, min))

# --------------------------------------------
# Get block window stats for group differences
# --------------------------------------------

blockstats <- g$blockstats(vcfresults,
					stepsize = blocksize,
					invariants = processedInvariants, 
					chrvec = chrvec,
					CSStrueSnps = CSStrueSnps, 
					FSTtrueSnps = FSTtrueSnps,
					psdMissingAction = psdMissingAction) 

# On output, "psd" refers to pairwise sequence divergence of a pair of specimens in each block
# and "np" refers to the number of non-missing psd values used in calculating the block average for each pair

# head(blockstats)					
  # ibase midbase nUnmasked nSnp nInvariants VARa VARb VARw fst nTrueSnp psd.1,1 psd.1,1 psd.1,1 psd.1,1
# 1     1     251        12    0           0    0    0    0 NaN        0       0       0       0       0
# 2   501     751         0    0           0    0    0    0 NaN        0       0       0       0       0
# 3  1001    1251         0    0           0    0    0    0 NaN        0       0       0       0       0
# 4  1501    1751         0    0           0    0    0    0 NaN        0       0       0       0       0
# 5  2001    2251         0    0           0    0    0    0 NaN        0       0       0       0       0
# 6  2501    2751         0    0           0    0    0    0 NaN        0       0       0       0       0
# ....
  # psd.2,2 psd.2,2 psd.2,2 psd.2,2 psd.2,2 psd.2,2 psd.2,2 psd.2,2 psd.2,2 psd.2,2 np.1,1 np.1,1 np.1,1
# 1       0       0       0       0       0       0       0       0       0       0      0      0      0
# 2       0       0       0       0       0       0       0       0       0       0      0      0      0
# 3       0       0       0       0       0       0       0       0       0       0      0      0      0
# 4       0       0       0       0       0       0       0       0       0       0      0      0      0
# 5       0       0       0       0       0       0       0       0       0       0      0      0      0
# 6       0       0       0       0       0       0       0       0       0       0      0      0      0


save(blockstats, file = paste(vcfdir, project, ".", chrname, ".blockstats.rdd", sep = ""))
# load(file = paste(vcfdir, project, ".", chrname, ".blockstats.rdd", sep = "")) # blockstats

