#!/usr/bin/Rscript

# Run in Unix as " Rscript readSaveDnpdByGroup.R ... "

# ALT alleles of "<*:DEL>" are set to variant type NA

# This reads the variant file, collects good bits, and saves the results. 
# Includes plots.
# No analyses are done by groups, but cases in which fewer than 2 groups have 1 genotype
# 	are dropped (this is a low threshold, use more stringent criteria later on)

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

# This code assumes that we are working only with ONE chromosome at a time
# Note: At most 3 ALT alleles permitted per snp in this version
# NOTE: windowmaskerSdust start position is 0-based, so added 1 (end position is 1-based!)

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn", "paxl", "paxb", "marine-pac")
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb", "marine-pac")

project <- args[1]
chrname <- args[2]
# DPmin <- as.integer(args[3])
# GTminFrac <- eval(parse(text=args[4]))
groupnames <- args[3:length(args)]

# load "chrvec" for the current chromosome
chrno 				<- gsub("^chr", "", chrname)
chrmaskfile         <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd
load(chrmaskfile) 	# object is named "chrvec"

fastaname		<- paste(chrname, "fa", sep = ".")
# chrgroupname		<- paste("group", chrno, sep="") # annotation database uses groupXXI not chrXXI
vcfname			<- paste(project, ".", chrname, ".var.vcf", sep="")
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

dropRareAlleles	<- FALSE
GTmissing		<- "."	# how GATK represents missing genotypes in the vcf file "./."
nMaxAlt			<- 3	# maximum number of ALT alleles

library(VariantAnnotation, quietly = TRUE)

# load(file = vcfresultsfile)
# vcf <- vcfresults$vcf
# GT <- vcfresults$GT
# altUsedList <- vcfresults$altUsedList
# snpTypeList <- vcfresults$snpTypeList
# nAltUsed <- sapply(altUsedList, length)

# ---------------------
# Read variant VCF file

# Use "vcfLong <- expand(vcf)" if want variants having multiple ALT alleles to appear on multiple rows, 
# one per ALT allele. This would be a straightforward way to get predictCoding for all alleles. 

vcf <- readVcf(file = vcfname, genome = fastaname) 
cat("Successfully read vcf file\n")
# You can import data subsets or fields if desired, see manual

# Useful accessor functions
# header(vcf) # Header information
# samples(header(vcf))		# names of fish in sample
# geno(header(vcf))			# info on GT, AD, DP, GQ, MIN_DP, PGT, PID PL RGQ SB
# # head(rowRanges(vcf), 3)	# Genomic positions. This command didn't work
# head(rowData(vcf), 3)		# Genomic positions. This command worked instead
# ref(vcf)					# extract the REF 
# alt(vcf)					# ALT alleles (DNAStringSetList)
# qual(vcf)	 				# SNP quality
# geno(vcf)					# Genotype data: List of length 10 names(10): GT AD DP GQ MIN_DP PGT PID PL RGQ SB
#							# 	RGQ is "Unconditional reference genotype confidence"
# geno(header(vcf))["DS",]	# didn't work - no DS in file
# geno(header(vcf))["SB",]	# Info on what a given genotype item is
# SB <-geno(vcf)$SB			# grab the genotype values of SB (Fisher's Exact Test to detect strand bias)
# info(vcf)					# Didn't give the same columns as in manual, gave
# names(info(vcf))    
 # [1] "AC"              "AF"              "AN"              "BaseQRankSum"   
 # [5] "ClippingRankSum" "DP"              "DS"              "END"            
 # [9] "FS"              "HaplotypeScore"  "InbreedingCoeff" "MLEAC"          
# [13] "MLEAF"           "MQ"              "MQRankSum"       "QD"             
# [17] "ReadPosRankSum"  "SOR"

# ---
# fish group codes 1, 2, ... Order is determined by sequence of groupnames, eg "paxl", "paxb", "marine-pac"
fishnames <- samples(header(vcf))
groupcodes <- vector(length = length(fishnames), mode = "integer")
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
cat("groupcodes:")
print(groupcodes)
 # [1] 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

# Not used here 
# nInd <- as.vector(table(groupcodes)) # n individuals genotyped in each group eg 11 11  7
# nMin <- floor(GTminFrac*nInd) # minimum number required in each group eg 7 7 4

# ---
# Drop masked SNP (variants whose start position is masked in the chromosome)

z1 <- chrvec[start(ranges(vcf))]  # hopefully this works even though start has duplicated entries
keep <- z1 %in% c("A","C","G","T") # i.e., includes only the good bases: only upper case, no "M"
vcf <- vcf[keep]
cat("\nCompleted removal of snp corresponding to masked bases\n")

# ---
# Decided to skip this DPmin step for vcf, it just causes problems with the allele set.
# Set SNP not meeting minimum DP criterion (DPmin = 1) to missing
# Set SNP with missing DP (i.e., DP is NA) to missing
# GT <- geno(vcf)$GT  # is a character matrix
# DP <- geno(vcf)$DP  # is a numeric matrix
# z <- mapply(GT, DP, FUN = function(gt, dp){
	# gt[dp < DPmin] <- GTmissing
	# gt[is.na(dp)] <- GTmissing
	# return(gt)
	# })
# z[z == GTmissing] <- NA # set "." to NA
# z <- matrix(z, nrow = nrow(GT))
# dimnames(z) <- dimnames(GT)
# GT <- z
# rm(z, DP)

GT <- geno(vcf)$GT  # is a character matrix
GT[ GT == GTmissing ] <- NA # set "." to NA
# head(GT)

# From now on use GT not geno(vcf)$GT

# ----------
# Grab ALT alleles and keep a copy that includes only alleles actually used in genotype calls
# Max of 3 ALT alleles per locus was used.

altallele <- alt(vcf) 			# is a compressed CharacterList, an IRanges object
altList <- as.list(altallele) 	# ordinary list
nAlt <- sapply(altList, length)

# table(nAlt, useNA = "always") # note there are no zeros
# nAlt
     # 1      2      3      4      5   <NA> 
# 268720  13263    763      9      1      0 

# Some ALT alleles are labeled as "<*:DEL>". 
# Broad says: "This means there is a deletion in that sample that spans this position, 
#	while in other samples there is another event at this site (typically a SNP)."
# table(unlist(altList))[1:4]
# <*:DEL>       A      AA     AAA  ...
   # 8352   74133      21       4  ...  

# Enumerate the alt alleles used in actual genotype calls
# whichAltAllelesUsed lists all the ALT alleles used in genotypes by their index (1, 2, ...)
cat("Determining which ALT alleles actually used in genotype calls; others ALT alleles set to NA\n")
whichAltAllelesUsed <- apply(GT, 1, function(x){
	# x <- GT[1,]
	x1 <- strsplit(x[!is.na(x)], split = "/")
	x1 <- sort( as.integer(unique(unlist(x1))) )
	whichAltAllelesUsed <- x1[x1 > 0] # drops the REF allele
	})
nAltUsed <- sapply(whichAltAllelesUsed, length) # can be 0 if only REF alleles are in the genotypes


# Set unused ALT alleles to NA, don't drop, so genotype integers remain indices of ALT alleles [1] [2] etc
# z <- which(nAltUsed < nAlt)
# altList[z[[1]]] 				#  "T" "C"
# whichAltAllelesUsed[z[[1]]] # integer(0) # ie, there are only REF genotypes, no ALT alleles are used
# altList[z[[2]]] 				#  "C"
# whichAltAllelesUsed[z[[2]]] # integer(0) # ie, there are only REF genotypes, no ALT alleles are used
# z <- which(nAltUsed < nAlt & nAltUsed > 0) integer(0) # ie here all the unused alleles involve only REF genotypes

altUsedList <- mapply(whichAltAllelesUsed, altList, FUN = function(x, y){
	# x <- whichAltAllelesUsed[[7]]; y <- altList[[7]]
	z <- y
	ialt <- 1:length(z)
	z[!(ialt %in% x)] <- NA
	return(z)
	})
altUsedList[nAltUsed == 0] <- NA 
# altList[101] 		# "T"
# altUsedList[101] 	# NA
# unique(GT[101,]) 	# "0/0" NA

# Note, there are still "<*:DEL>" alleles. 

# Note that NA has length 1, so a nonzero nAltUsed doesn't mean there are real ALT alleles there

# Compare the number of ALT alleles used vs number of ALT alleles called
# table(sapply(altList, length), useNA = "always")
     # 1      2      3      4      5 
# 268720  13263    763      9      1 

# table(nAltUsed, useNA = "always")
     # 0      1      2      3 
 # 26184 244465  11415    692 

# ---
# Drop variants in which fewer than 2 groups have at least one genotype

# Tally up the number of genotypes in each group
z <- split(colnames(GT), groupcodes) # Note that order of resulting groups is as in groupnames and groupcodes

# $`1`
 # [1] "PaxLim-PxCL09maleLM1-GS14" "PaxLim-PxLfemale6-GS18"   
 # [3] "PaxLim-PxLmale102-GS16"    "PaxLim-PxLmale106-GS15"   
 # [5] "PaxLim-PxLmale107-GS17"    "paxl01"                   
 # [7] "paxl05"                    "paxl09"                   
 # [9] "paxl10"                    "paxl13"                   
# [11] "paxl14"                   

# $`2`
 # [1] "PaxBen-PxBmale5-GS11"        "PaxBen-PxBmale6-GS12"       
 # [3] "PaxBen-PxBmale8-GS10"        "PaxBen-PxCL09femaleBF6-GS13"
 # [5] "PaxBen-RPxCL09maleBM2-GS9"   "paxb04"                     
 # [7] "paxb05"                      "paxb06"                     
 # [9] "paxb07"                      "paxb08"                     
# [11] "paxb09"                     

# $`3`
# [1] "Marine-Pac-BIGR-52_54_2008-02"  "Marine-Pac-Japan-01-Katie"     
# [3] "Marine-Pac-LittleCampbell-LC1D" "Marine-Pac-MANC_X_X05"         
# [5] "Marine-Pac-Oyster-06-Sara"      "Marine-Pac-Seyward-01-Sara"    
# [7] "Marine-Pac-WestCreek-01-Sara"  

nCalledGenotypes <- lapply(z, function(z){
	z1 <- apply(GT[ , z], 1, function(z){sum(!is.na(z))})
	})
names(nCalledGenotypes) <- groupnames # don't sort groupnames! This is the order that determined codes
z2 <- lapply(nCalledGenotypes, function(z1){z1 >= 1}) # indicates whether there's at least one genotype
z2 <- data.frame(z2)
z3 <- apply(z2, 1, sum)
keep <- z3 >= 2

cat("Keeping only snp that have at least one genotype in at least 2 groups\n")
vcf <- vcf[keep]
GT <- GT[keep, ]
altList <- altList[keep]
altUsedList <- altUsedList[keep]
nCalledGenotypes <- nCalledGenotypes[keep]
nAlt <- nAlt[keep]
nAltUsed <- nAltUsed[keep]


# Types of variants
# -----------------
# Variant type is defined for VCFtools as (http://vcftools.sourceforge.net/VCF-poster.pdf)
# SNPs
# Alignment  VCF representation#   ACGT     POS REF ALT
#   ATGT      2   C   T

# Deletions# Alignment  VCF representation
#   ACGT     POS REF ALT 
#   A--T      1  ACG  A

# Insertions
# Alignment  VCF representation
#  AC-GT     POS REF ALT
#  ACTGT      2   C   CT

# Complex events
# Alignment  VCF representation
#   ACGT     POS REF ALT 
#   A-TT      1  ACG  AT

cat("Determining variant types\n")
i <- nchar(ref(vcf)) # length of the reference sequence
i <- split(i, 1:length(i)) # split REF allele size into a list
j <- lapply(altUsedList, function(x){
	j <- nchar(x)            	# length of ALT alleles in altUsedList; 
								# is integer(0) if altUsedList[[i]] is character(0)
								# but is 2 if altUsedList[[i]] is NA
	j[x == "<*:DEL>"] <- NA  	# length of "<*:DEL>" set to NA
	j[is.na(x)] <- NA
	return(j)
	})

# Inspect some results from the first run
# Looking for cases where changes don't meet the simple criteria for snp, del, ins
# From the results just below, it looks like even when i and/or j isn't 1, 
#	i < j is always an "ins", i > j is always "del", and i = j is always "snp".
#y <- mapply(snpTypeList, altUsedList,
	#FUN = function(x, y){any(is.na(x)) & all(y[!is.na(y)] != "<*:DEL>") & any(nchar(y[!is.na(y)]) > 1) })
#ref <- as.vector(ref(vcf))
#z <- mapply(ref, snpTypeList, altUsedList, FUN = c)
#head(z[y])
# Here are some examples where snp type is NA
# "TCC"     -> "TC" "T";        NA  "del"     # is del, ref begins with alt
# "CGG"     -> "CG" "C";        NA  "del"     # is del, ref begins with alt
# "GAAGCACGTACT" -> "GACGTACT" "G"; NA "del"  # is del, 1st letter of alt same as ref, 2nd-5th bases dropped
# "ACAACT"  -> "AT"  "A";       NA  "del"     # is del? 1st letter of alt same as ref, 2nd-5th bases dropped
# "ACACTCACTCCGCATCCTATCG" -> "A" "ACACTCCGCATCCTATCG"; "del" NA; is del, 1st letter alt same as ref, 2nd-4th bases dropped
# "TAA"     -> "TTAAA" "T";     NA  "del"     # is ins, 1st letter of alt same as ref, inserted 2 bases, then rest is same
# "CTCCCA"  -> "C" "CGTGTCCCA"; "del" NA      # is ins, 1st letter of alt same as ref, inserted 3 bases, then rest is same
# "TCTCA"   -> "T" "TGCTCA";   "del"  NA      # is ins, 1st letter of alt same as ref, inserted 1 base, then rest is same
# "CGG"     -> "C"  "TGG";   "  del" NA       # is snp, first letter different
# "GGCCGGT" -> "CGCCGGT" "G";   NA  "del"     # is snp, first letter different
# "TC"      -> "CC"  "GC" "T";  NA   NA "del" # is snp, first letter different
# "GC"      -> "CC"  "G";       NA  "del"     # is snp, first letter different
# "GGT"     -> "AGT" "G";       NA  "del"     # is snp, first letter different
# "AG"      -> "A" "CG" "TG";  "del" NA NA    # is snp, first letter different
# "AAACTAC" -> "GAACTAC" "A";   NA  "del"     # is snp, first letter different
# "ATGAAAC" -> "A" "GTGAAAC";   "del" NA      # is snp, first letter different
# "CG"      -> "GG"  "C";       NA   "del"    # is snp, first letter different
# "ACT"     -> "TCT" "A";       NA   "del"    # is snp, first letter different

# 2nd run uses broader criteria to designate variant type
# Note: j is based on altUsedList, so snpTypeList is also
snpTypeList <- mapply(i, j, FUN = function(i, j){ # i and j refer to nchar of ref and alt
	snptype <- rep(NA, length(j)) # initialize with NAs
	# snptype[i == 1 & j == 1] <- "snp" # first run
	# snptype[i < j  & i == 1] <- "ins" # first run
	# snptype[i > j  & j == 1] <- "del" # first run
	snptype[i == j] <- "snp" # 2nd run
	snptype[i <  j] <- "ins" # 2nd run
	snptype[i >  j] <- "del" # 2nd run
	snptype[length(j) == 0] <- NA # snp type is NA if length of altUsedList is 0, ie no snp used
	return(snptype)
	})
names(snpTypeList) <- names(altUsedList)

# Check that all multi-based snp differ only by the first letter
# Remember that altUsedList still contains "<*:DEL>"
# snp check:
# snp <- unlist(snpTypeList[nAltUsed > 0])
# ref <- rep(as.vector(ref(vcf)), nAltUsed)
# alt <- unlist(altUsedList)
# x <- data.frame(ref, alt, snp, stringsAsFactors = FALSE)
# x <- x[x$alt != "<*:DEL>", ]
# x <- x[nchar(x$ref) == nchar(x$alt) & nchar(x$ref) > 1, ] 
# rownames(x) <- 1:nrow(x)
# head(x)
# table(x$snp) # confirm that they are all "snp"
# any( substr(x$ref, 1, 1) == substr(x$alt, 1, 1) ) # FALSE confirms that the first letter is always different
# all( substr(x$ref, 2, nchar(x$ref)) == substr(x$alt, 2, nchar(x$alt)) ) # TRUE confirms that rest is always the same

# del check:
# snp <- unlist(snpTypeList[nAltUsed > 0])
# ref <- rep(as.vector(ref(vcf)), nAltUsed)
# alt <- unlist(altUsedList)
# x <- data.frame(ref, alt, snp, stringsAsFactors = FALSE)
# x <- x[x$alt != "<*:DEL>", ]
# x <- x[x$snp == "del", ] # 
# rownames(x) <- 1:nrow(x)
# head(x)
# all( substr(x$ref, 1, 1) == substr(x$alt, 1, 1) ) # TRUE confirms that the first letter is always the same
# x <- x[nchar(x$alt) > 1, ] # 
# #      last n letters of alt             last n letters of ref
# all( substr(x$alt, 2, nchar(x$alt)) == substr(x$ref, nchar(x$ref) - nchar(x$alt) + 2, nchar(x$ref)) ) # TRUE!

# ins check:
# snp <- unlist(snpTypeList[nAltUsed > 0])
# ref <- rep(as.vector(ref(vcf)), nAltUsed)
# alt <- unlist(altUsedList)
# x <- data.frame(ref, alt, snp, stringsAsFactors = FALSE)
# x <- x[x$alt != "<*:DEL>", ]
# x <- x[x$snp == "ins", ] # 
# rownames(x) <- 1:nrow(x)
# head(x)
# all( substr(x$ref, 1, 1) == substr(x$alt, 1, 1) ) # TRUE confirms that the first letter is always the same
# x <- x[nchar(x$ref) > 1, ] # 
# #      last n letters of ref              last n letters of alt
# all( substr(x$ref, 2, nchar(x$ref)) == substr(x$alt, nchar(x$alt) - nchar(x$ref) + 2, nchar(x$alt)) ) # TRUE!


# --------------------------------------
# Make a table of the allele frequencies
cat("\nCalculating allele frequencies by group at every snp\n")
tGT <- as.data.frame(t(GT), stringsAsFactors = FALSE)
alleleFreqByGroup <- lapply(tGT, function(locus){  # columns of tGT are loci, so apply function locus by locus
	# locus <- tGT[,1]
	z1 <- split(locus, groupcodes) # the genotypes for each group at the locus
	z2 <- lapply(z1, function(x){ 
		unlist(strsplit(x, split = "/"))
		}) 		# the alleles for each group at the locus
	z3 <- z2	# initialize
	for(i in 1:length(z3)) z3[[i]] <- rep( groupnames[i], length(z3[[i]]) ) # replace z3 with corresponding group name
	table( unlist(z3), factor(unlist(z2), levels=0:nMaxAlt) ) # factor so that all alleles are counted in each table
	})

# alleleFreqByGroup[[1]]     
             # 0 1 2 3
  # marine-pac 2 0 0 0
  # paxb       0 0 0 0
  # paxl       1 3 0 0

# --------------------------------------
# Calculate the allele proportions and drop rare alleles if dropRareAlleles is TRUE

alleleProportions <- lapply(alleleFreqByGroup, function(x){
	z <- colSums(x)
	z1 <- z/sum(z)
	})

# alleleProportions[[1]]
  # 0   1   2   3 
# 0.5 0.5 0.0 0.0

if(dropRareAlleles){ # IN PROGRESS
	# Delete rare alleles (< 5%) based on vcfresults$alleleProportions
	# To do this, identify the allele and change the genotypes to missing that correspond to those rare alleles
	rareAlleles <- lapply(vcfresults$alleleProportions, function(x){
		names(x[x < .05])
		})
#	rareAlleles <- sapply(rareAlleles, function(x){paste(x, collapse = "")}) # very slow
	library(stringr)
	rareAlleles <- sapply(rareAlleles, function(x){str_c(c("[",x,"]"), collapse = "")}) # bit faster
	# head(rareAlleles)
	 # chrXXI:6316_C/T chrXXI:10364_T/A chrXXI:10365_G/T chrXXI:16024_C/T 
	          # "[23]"          "[023]"          "[023]"          "[123]" 
	# chrXXI:16025_C/G chrXXI:16026_A/T 
	         # "[123]"           "[23]"
	z <- lapply(genotypes, function(x){
		# x <- genotypes[,1]
		grepl( rareAlleles, x))
		})
		}


# --------------------------------------
# Save everything to a list for analyses

vcfresults <- list()
vcfresults$groupcodes <- groupcodes
vcfresults$vcf <- vcf
vcfresults$GT <- GT
vcfresults$nCalledGenotypes <- nCalledGenotypes
vcfresults$nAlt <- nAlt
vcfresults$altUsedList <- altUsedList
vcfresults$nAltUsed <- nAltUsed
vcfresults$snpTypeList <- snpTypeList # based on altUsedList
vcfresults$alleleFreqByGroup <- alleleFreqByGroup
vcfresults$alleleProportions <- alleleProportions

# cat("\nSaving results\n")
save(vcfresults, file = vcfresultsfile)
# load(file = vcfresultsfile) # saved object is "vcfresults"
	
# --------------------------------------
# Summary stats and quality plots

# nrow(vcfresults$GT)
# [1] 281607

# Expand the lists of ALT alleles, repeat REF as many times
snp <- unlist(snpTypeList[nAltUsed > 0])
ref <- rep(as.vector(ref(vcf)), nAltUsed)
alt <- unlist(altUsedList[nAltUsed > 0])

# c(length(snp), length(ref), length(alt))
# [1] 268287 268287 268287

# head(table(alt))                                # <*:DEL> are still included in alt
# table(snp[alt == "<*:DEL>"], useNA = "always")  # but they are not classified as snp, ins, or del, so ok

# Tabulate all snp types
cat("Table of variant types used in genotype calls (<*:DEL> is NA)\n")
print(table(snp, useNA="always"))
   # del    ins    snp   <NA> 
 # 24286  20909 216418   6674 

# Transition-transversion ratios - true snp only (not indels)
cat("\nTransition-transversion ratio - all true snp\n")
snp1 <- snp[!is.na(snp) & snp == "snp"]
ref1 <- ref[!is.na(snp) & snp == "snp"]
alt1 <- alt[!is.na(snp) & snp == "snp"]

# head(cbind(ref1,alt1))
# head(cbind(ref1,alt1)[nchar(ref1) > 1,])
# c(length(snp1), length(ref1), length(alt1))
# [1] 216418 216418 216418

ref1 <- substr(ref1, 1, 1) # keep the first letter of REF
alt1 <- substr(alt1, 1, 1) # keep the first letter of ALT

# c( length(ref1[ref1 != alt1]), length(alt1[ref1 != alt1]) )
# [1] 216418 216418

print(g$tstv(ref1[ref1 != alt1], alt1[ref1 != alt1]))
# $R
# [1] 1.258801
# $tstv
     # alt
# ref     pur   pyr
  # pur 60301 47947
  # pyr 47864 60306
# $tot
# [1] 216418

# Plots of quality metrics
cat("\nPlotting quality metrics\n")
pdf(file = paste(project, ".", chrname, ".vcfPlots", ".pdf", sep=""))

rowdata <- rowData(vcf)
# head(rowdata)
# GRanges object with 6 ranges and 5 metadata columns:
                   # seqnames         ranges strand | paramRangeID            REF
                      # <Rle>      <IRanges>  <Rle> |     <factor> <DNAStringSet>
   # chrXXI:6316_C/T   chrXXI [ 6316,  6316]      * |         <NA>              C
  # chrXXI:10364_T/A   chrXXI [10364, 10364]      * |         <NA>              T
  # chrXXI:10365_G/T   chrXXI [10365, 10365]      * |         <NA>              G
  # chrXXI:16024_C/T   chrXXI [16024, 16024]      * |         <NA>              C
  # chrXXI:16025_C/G   chrXXI [16025, 16025]      * |         <NA>              C
  # chrXXI:16026_A/T   chrXXI [16026, 16026]      * |         <NA>              A
                               # ALT      QUAL      FILTER
                   # <CharacterList> <numeric> <character>
   # chrXXI:6316_C/T               T    253.35           .
  # chrXXI:10364_T/A               A    517.44           .
  # chrXXI:10365_G/T               T    517.44           .
  # chrXXI:16024_C/T               T     77.92           .
  # chrXXI:16025_C/G               G     58.52           .
  # chrXXI:16026_A/T               T    748.99           .
  
# QUAL
# QUAL is probability of a polymorphism at the site, not a measure of quality of genotypes
# It is "Phred scaled Probability that REF/ALT polymorphism exists at this site given sequencing data"
# A value of 10 indicates p = 0.99, i.e., 1 - 0.01
# It depends on the quality of genotypes, but is not itself a measure of genotype quality.
# GATK2 drops snp with QUAL <= 30
hist(rowdata$QUAL[rowdata$QUAL <= 10000], breaks = 100, right = FALSE, col = "red", 
	main = "QUAL (prob of polymorphism) for called variants")

# Individual fish genotype quality scores GQ for called genotypes
# "if GT is 0/1, then GQ is really L(0/1) / (L(0/0) + L(0/1) + L(1/1))"
GQ <- matrix(geno(vcf)$GQ, ncol=1)
GTm <- matrix(GT, ncol=1)
hist(GQ[!is.na(GTm)], right = FALSE, col = "red", breaks = 50, 
	main = "Histogram of individual fish genotype quality score, GQ,\nfor called genotypes")

# Individual fish read depth DP for called genotypes
DP <- matrix(geno(vcf)$DP, ncol=1)
hist(DP[DP <= 50 & !is.na(GTm)], right = FALSE, col = "red", breaks = 50,
	main = "Histogram of individual fish read depth, DP,\nfor called genotypes")

# Relationship between GQ and DP for called genotypes
# GTm <- matrix(vcfresults$GT, ncol=1)
plot(GQ[DP <= 10 & !is.na(GTm)] ~ jitter(DP[DP <= 10 & !is.na(GTm)]), pch = ".", col = "red",
	main = "Relationship between GQ and DP for called genotypes,\n(note that GQ of called genotypes is often very low)")

# Strand bias
# From the GATK pages: "Higher SB values denote more bias (and therefore are more likely to indicate false positive calls)."
# SB <- geno(vcf)$SB # all NA

# Phred-scaled P-value for Fisher exact test of strand bias (higher is more biased)
FS <- info(vcf)$FS
hist(FS[FS <= 100], col="red", breaks = 100, right=FALSE, xlab = "Phred-scaled P",
	main = "Results of Fisher exact test of strand bias") # peaks at 0

# Strand bias vs Quality by Depth (see Fig 4a,b in de Pristo et al 2011)
QD <- info(vcf)$QD
plot(FS ~ QD, pch = ".", col = "red", ylab = "FS (Fisher test of strand bias)", 
	xlab = "QD (Quality by Depth)", main = "Strand Bias vs Quality by Depth (cf. Fig 4a,b in de Pristo et al 2011)")

# Total read depth
DPtotal <- info(vcf)$DP
hist(DPtotal[DPtotal <= length(groupcodes)*50], right = FALSE, col = "red", breaks = 200, 
	main = "Total read depth") 
dev.off()

