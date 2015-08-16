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
groupnames <- args[3:length(args)]

dropRareAlleles	<- FALSE
plotQualMetrics <- TRUE

# load "chrvec" for the current chromosome
chrno 				<- gsub("^chr", "", chrname)
chrmaskfile         <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd

fastaname		<- paste(chrname, "fa", sep = ".")
# chrgroupname		<- paste("group", chrno, sep="") # annotation database uses groupXXI not chrXXI
vcfname			<- paste(project, ".", chrname, ".var.vcf", sep="")
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

GTmissing		<- "."	# how GATK represents missing genotypes in the vcf file "./."
nMaxAlt			<- 3	# maximum number of ALT alleles

library(VariantAnnotation, quietly = TRUE)

load(chrmaskfile) 	# object is named "chrvec"
which.M <- which(chrvec == "M")
# object.size(chrvec)
#  93,740,176 bytes # chrXXI
# 500,401,968 bytes # chrUn
rm(chrvec)

# Garbage collection -- current memory consumption
gcinfo(TRUE) # sends messages about memory garbage collection
gc()

# ---------------------
# Alternatives to reading whole vcf file

# 1. Make tabix and use data ranges
# compressedVcfname <- paste(project, chrname, "var.vcf.bgz", sep = ".")
# zipped <- bgzip(vcfname, compressedVcfname)
# idx <- indexTabix(zipped, "vcf")
# vcfTabix <- TabixFile(zipped, idx)
# vcf <- readVcf(vcfTabix, chrname, IRanges(c(1, 100000))) # not working yet

# 2. Read a subset of fields
# vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("ALT", "QUAL"), geno="GT", info=NA))

# ---------------------
# Read variant VCF file

# Use "vcfLong <- expand(vcf)" if want variants having multiple ALT alleles to appear on multiple rows, 
# one per ALT allele. This would be a straightforward way to get predictCoding for all alleles. 

# Reads the whole vcf file
# vcf <- readVcf(file = vcfname, genome = fastaname) 
# object.size(vcf)
# 107,877,192 bytes # chrXXI

# Reads only a part of the vcf. ref(), alt(), qual(), geno()$GT still work

vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("ALT", "QUAL"), geno="GT", info=NA))

# object.size(vcf)
#  50,514,280 bytes # chrXXI
# 210,524,096 bytes # chrUn

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

# Q: What does QUAL = NA imply?

# ---
# Drop masked SNP (variants whose start position is masked in the chromosome)

keep <- !(start(ranges(vcf)) %in% which.M) # i.e., includes only the good bases: only upper case, no "M"
vcf <- vcf[keep]
cat("\nCompleted removal of snp corresponding to masked bases\n")
rm(keep)
# rm(which.M)
# object.size(vcf)
# 36,360,752 bytes # chrXXI
gc()

# ------
# set missing genotypes to NA
geno(vcf)$GT[geno(vcf)$GT == GTmissing] <- NA  # set "." to NA
gc()

# ---
# Drop variants in which fewer than 2 groups have at least one genotype

# Tally up the number of genotypes in each group
samplesByGroup <- split(samples(header(vcf)), groupcodes) # Order of resulting groups is as in groupnames and groupcodes

# $`1`
 # [1] "PaxLim-PxCL09maleLM1-GS14" "PaxLim-PxLfemale6-GS18"    "PaxLim-PxLmale102-GS16"   
 # [4] "PaxLim-PxLmale106-GS15"    "PaxLim-PxLmale107-GS17"    "paxl01"                   
 # [7] "paxl05"                    "paxl09"                    "paxl10"                   
# [10] "paxl13"                    "paxl14"                   

# $`2`
 # [1] "PaxBen-PxBmale5-GS11"        "PaxBen-PxBmale6-GS12"        "PaxBen-PxBmale8-GS10"       
 # [4] "PaxBen-PxCL09femaleBF6-GS13" "PaxBen-RPxCL09maleBM2-GS9"   "paxb04"                     
 # [7] "paxb05"                      "paxb06"                      "paxb07"                     
# [10] "paxb08"                      "paxb09"                     

# $`3`
# [1] "Marine-Pac-BIGR-52_54_2008-02"  "Marine-Pac-Japan-01-Katie"      "Marine-Pac-LittleCampbell-LC1D"
# [4] "Marine-Pac-MANC_X_X05"          "Marine-Pac-Oyster-06-Sara"      "Marine-Pac-Seyward-01-Sara"    
# [7] "Marine-Pac-WestCreek-01-Sara"  

nCalledGenotypes <- lapply(samplesByGroup, function(z){
	z1 <- apply(geno(vcf)$GT[ , z], 1, function(z){sum(!is.na(z))})
	})
names(nCalledGenotypes) <- groupnames # don't sort groupnames! This is the order that determined codes

z2 <- lapply(nCalledGenotypes, function(z1){z1 >= 1}) # indicates whether there's at least one genotype
z2 <- data.frame(z2)
z3 <- apply(z2, 1, sum)
keep <- z3 >= 2

cat("Keeping only snp that have at least one genotype in at least 2 groups\n")
vcf <- vcf[keep]
rm(nCalledGenotypes)
rm(z2)
rm(z3)
rm(keep)

gc()

# --------------------------------------
# Calculate the allele proportions and drop rare alleles if dropRareAlleles is TRUE

if(dropRareAlleles){ # IN PROGRESS
	# Delete rare alleles (< 5%) based on vcfresults$alleleProportions
	# This means 5% across populations? What if analysis includes wheatlandi eg, just a single individual?
	# Or should we just drop alleles with less than 5% within a group (but what if the allele is common in another group)
	# Maybe this rule: drop alleles less than 5% within groups only if also less than 5% across groups
	
	# To do this, identify the allele and change the genotypes to missing that correspond to those rare alleles

	alleleProportions <- lapply(alleleFreqByGroup, function(x){
		z <- colSums(x)
		z1 <- z/sum(z)
		})
	# alleleProportions[[1]]
	  # 0   1   2   3 
	# 0.5 0.5 0.0 0.0

	rareAlleles <- lapply(vcfresults$alleleProportions, function(x){
		names(x[x < .05])
		})
	# rareAlleles <- sapply(rareAlleles, function(x){paste(x, collapse = "")}) # very slow
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


# --------------------------
# Determine alleles actually used in genotype calls
# Max of 3 ALT alleles per locus was used.

# Some ALT alleles are labeled as "<*:DEL>". 
# Broad says: "This means there is a deletion in that sample that spans this position, 
#	while in other samples there is another event at this site (typically a SNP)."
# table(unlist(alt(vcf)))[1:4]
# <*:DEL>       A      AA     AAA  ...
   # 8338   73896      21       4  ...  

# Enumerate the alt alleles used in actual genotype calls
# whichAltAllelesUsed lists all the ALT alleles used in genotypes by their index (1, 2, ...)

cat("Determining which ALT alleles actually used in genotype calls; others ALT alleles set to NA\n")
# unused ALT alleles are set to NA -- they are NOT DELETED to preserve indices
altUsedList <- g$makeAltUsedList(geno(vcf)$GT, alt(vcf))
nAltUsed <- sapply(altUsedList, function(x){length(x[!is.na(x)])})

# Note, there are still "<*:DEL>" alleles. 

# Note that NA has length 1, so a nonzero nAltUsed doesn't mean there are real ALT alleles there

# Compare the number of ALT alleles used vs number of ALT alleles called
# table(sapply(alt(vcf), length), useNA = "always")
     # 1      2      3      4      5 
# 267591  13243    763      9      1 

# table(nAltUsed, useNA = "always")
     # 0      1      2      3   <NA>
 # 26102 243415  11398    692      0


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

# Examples from earlier analysis where snp type is NA
#  ref          altUsedList
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

snpTypeList <- g$makeSnpTypeList(REF = ref(vcf), ALTlist = altUsedList)

# Check that all multi-based snp differ only by the first letter
# Remember that altUsedList still contains "<*:DEL>"
# snp check:
# snp <- unlist(snpTypeList[nAltUsed > 0])
# ref <- rep(as.vector(ref(vcf)), nAltUsed)
# alt <- unlist(altUsedList[nAltUsed > 0])
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
# alt <- unlist(altUsedList[nAltUsed > 0])
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
# alt <- unlist(altUsedList[nAltUsed > 0])
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
# Not sure this is essential at this stage -- more useful when analyzing a subset of groups

alleleFreqByGroup <- g$tableAlleleFreqByGroup(geno(vcf)$GT, groupnames, groupcodes)

# alleleFreqByGroup[[1]]     
             # 0 1 2 3
  # marine-pac 2 0 0 0
  # paxb       0 0 0 0
  # paxl       1 3 0 0


# --------------------------------------
# Transition-transversion ratio

# nrow(geno(vcfresults$vcf)$GT)
# [1] 281607

g$vcfTsTv(ref(vcf), altUsedList, snpTypeList)

# Table of variant types used (<*:DEL> is NA)
# snp
   # del    ins    snp   <NA> 
 # 24286  20909 216418   6674 

# Transition-transversion ratio - all true snp
# $R
# [1] 1.258801

# $tstv
     # alt
# ref     pur   pyr
  # pur 60301 47947
  # pyr 47864 60306

# $tot
# [1] 216418


# --------------------------------------
# Save everything to a list for analyses

vcfresults <- list()
vcfresults$groupcodes <- groupcodes
vcfresults$vcf <- vcf
rm(vcf)
gc()
vcfresults$altUsedList <- altUsedList
rm(altUsedList)
gc()
vcfresults$nAltUsed <- nAltUsed
rm(nAltUsed)
vcfresults$snpTypeList <- snpTypeList # based on altUsedList
rm(snpTypeList)
vcfresults$alleleFreqByGroup <- alleleFreqByGroup
rm(alleleFreqByGroup)

# cat("\nSaving results\n")
save(vcfresults, file = vcfresultsfile)
# load(file = vcfresultsfile) # saved object is "vcfresults"
rm(vcfresults)

if(plotQualMetrics{
	# --------------------------------------
	# Plots of quality metrics
	# Current version drops only the rows corresponding to M in chrvec
	
	vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("QUAL"), geno=c("GT", "GQ", "DP"), 
		info=c("FS", "QD", "DP")))	
	gc()
	keep <- !(start(ranges(vcf)) %in% which.M) # i.e., includes only the good bases: only upper case, no "M"
	vcf <- vcf[keep]
	rm(keep)

	cat("\nPlotting quality metrics\n")
	pdf(file = paste(project, ".", chrname, ".vcfPlots", ".pdf", sep=""))
	
	# QUAL ---
	# QUAL is probability of a polymorphism at the site, not a measure of quality of genotypes
	# It is "Phred scaled Probability that REF/ALT polymorphism exists at this site given sequencing data"
	# A value of 10 indicates p = 0.99, i.e., 1 - 0.01
	# It depends on the quality of genotypes, but is not itself a measure of genotype quality.
	# GATK2 drops snp with QUAL <= 30
	
	hist(qual(vcf)[qual(vcf) <= 2000], breaks = 100, right = FALSE, col = "red", 
		main = "QUAL (prob of polymorphism) for called variants")
	
	# ---
	# Individual fish genotype quality scores GQ for called genotypes
	# "if GT is 0/1, then GQ is really L(0/1) / (L(0/0) + L(0/1) + L(1/1))"
	
	GQ <- matrix(geno(vcf)$GQ, ncol=1)
	GT <- matrix(geno(vcf)$GT, ncol=1)
	hist(GQ[!is.na(GT)], right = FALSE, col = "red", breaks = 50, 
		main = "Histogram of individual fish genotype quality score, GQ,\nfor called genotypes")

	# ---
	# Individual fish read depth DP for called genotypes
	
	DP <- matrix(geno(vcf)$DP, ncol=1)
	hist(DP[DP <= 50 & !is.na(GT)], right = FALSE, col = "red", breaks = 50,
		main = "Histogram of individual fish read depth, DP,\nfor called genotypes")
	
	# ---
	# Relationship between GQ and DP for called genotypes
	# Commented out - requires a lot of memory

	# plot(GQ[DP <= 10 & !is.na(GT)] ~ jitter(DP[DP <= 10 & !is.na(GT)]), pch = ".", col = "red",
		# main = "Relationship between GQ and DP for called genotypes,\n(note that GQ of called genotypes is often very low)")
	
	rm(GT)
	rm(DP)
	rm(GQ)

	# ---
	# Strand bias
	# From the GATK pages: "Higher SB values denote more bias (and therefore are more likely to indicate false positive calls)."
	# SB <- geno(vcf)$SB # all NA
	
	# Phred-scaled P-value for Fisher exact test of strand bias (higher is more biased)
	
	FS <- info(vcf)$FS
	hist(FS[FS <= 40], col="red", breaks = 100, right=FALSE, xlab = "Phred-scaled P",
		main = "Results of Fisher exact test of strand bias \n(higher P is more biased)") # peaks at 0
		
	# Strand bias vs Quality by Depth (see Fig 4a,b in de Pristo et al 2011)
	
	QD <- info(vcf)$QD
	plot(FS ~ QD, pch = ".", col = "red", ylab = "FS (Fisher test of strand bias)", 
		xlab = "QD (Quality by Depth)", main = "Strand Bias vs Quality by Depth (cf. Fig 4a,b in de Pristo et al 2011)")

	rm(FS)
	rm(QD)	

	# ---
	# Total read depth
	DPtot <- info(vcf)$DP
	hist(DPtot[DPtot <= length(groupcodes)*50], right = FALSE, col = "red", breaks = 200, 
		main = "Total read depth") 
	
	dev.off()
	}
