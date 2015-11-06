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

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn", "paxl", "paxb", "marine-pac")
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb", "marine-pac")
# args <- c("BenlimAllMarine", "chrXXI", "paxl","paxb","pril","prib","qryl","qryb","ensl","ensb","marine-pac","marine-atl","marine-jap","solitary")

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)]

dropRareAlleles	<- FALSE
plotQualMetrics <- FALSE
saveBiAllelic <- FALSE # saves a second data set having exactly 2 snp per marker (not necessarily the REF), no indels

# load "chrvec" for the current chromosome
chrno 				<- gsub("^chr", "", chrname)
chrmaskfile         <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd

fastaname		<- paste(chrname, "fa", sep = ".")
# chrgroupname		<- paste("group", chrno, sep="") # annotation database uses groupXXI not chrXXI
vcfname			<- paste(project, ".", chrname, ".var.vcf", sep="")
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")
vcfBiAllelicFile	<- paste(project, ".", chrname, ".vcfresultsBiAllelic.rdd", sep = "")

GTmissing		<- "."	# how GATK represents missing genotypes in the vcf file "./."
nMaxAlt			<- 3	# maximum number of ALT alleles

control <- list()
control$nMaxAlt <- nMaxAlt
control$snpOptions <- c(dropRareAlleles = dropRareAlleles, saveBiAllelic = saveBiAllelic)

cat("\nControl settings on this run\n")
print(control)

library(VariantAnnotation, quietly = TRUE)

load(chrmaskfile) 	# object is named "chrvec"
which.chrvec.M <- which(chrvec == "M")

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

# -----
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
 
nInd <- as.vector(table(groupcodes)) # number of individuals genotyped in each group
# [1] 11 11  7
names(nInd)<-groupnames

control$nInd <- nInd

# Q: What does QUAL = NA imply?

# ------
# Drop masked SNP (variants whose start position is masked in the chromosome)

keep <- !(start(ranges(vcf)) %in% which.chrvec.M) # i.e., includes only the good bases: only upper case, no "M"
vcf <- vcf[keep]
cat("\nCompleted removal of snp corresponding to masked bases\n")

rm(keep)
# rm(which.chrvec.M) # remove later

# object.size(vcf)
# 36,360,752 bytes # chrXXI
gc()

# ------
# set missing genotypes to NA instead of "."
geno(vcf)$GT[geno(vcf)$GT == GTmissing] <- NA

# -----
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

# Tabulate genotype frequencies by group
# genotypeFreqByGroup <- apply(geno(vcf)$GT, 1, function(x){
	# table(groupcodes, x)
	# # xtabs(~ groupcodes + x) # took twice as long!
	# })

# Indicate whether there's at least one genotype in at least two groups. 
# No need for a more stringent criterion here because one will need to be applied in gtStats2groups anyway

nCalledGenotypes <- lapply(samplesByGroup, function(z){
	z1 <- apply(geno(vcf)$GT[ , z], 1, function(z){sum(!is.na(z))})
	})
names(nCalledGenotypes) <- groupnames # don't sort groupnames! This is the order that determined codes

z2 <- lapply(nCalledGenotypes, function(z1){z1 >= 1})
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
# Make a table of the allele frequencies

alleleFreqByGroup <- g$tableAlleleFreqByGroup(geno(vcf)$GT, groupnames, groupcodes)
	# unname(geno(vcf)$GT[1, ])
	# groupcodes
	# alleleFreqByGroup[[1]]     
	             # 0 1 2 3
	  # paxl       1 3 0 0
	  # paxb       0 0 0 0
	  # marine-pac 2 0 0 0

gc()

# --------------------------------------
# Drop rare alleles if dropRareAlleles is TRUE

if(dropRareAlleles){
	# Delete rare alleles -- we are using ** RULE2 ** otherwise we get lots of rare alleles when many genotypes are missing
	# RULE1: drop alleles that are less than 5% total measured as a percent of 2*sum(nInd), the maximum number of alleles possible
	# or
	# RULE2: drop alleles that are less than 5% total measured as a percent of the total number of alleles in the sample.

	# Code below can also: identify alleles that are < 5% within every group
	#                      identify singleton alleles within groups
	# (currently commented out)
	
	# Note that afterward, some loci will become invariants

	cat("Deleting rare alleles and setting their genotypes to NA\n")

	# Test whether any alleles are rare globally
	whichRareTot <- lapply(alleleFreqByGroup, function(x){
		# x <- alleleFreqByGroup[["chrXXI:16024_C/T"]]
		              # 0  1  2  3
		  # paxl       20  0  0  0
		  # paxb       18  0  0  0
		  # marine-pac 11  1  0  0
		# p <- colSums(x)/(2*sum(nInd)) # RULE1
		p <- colSums(x)/sum(colSums(x)) # RULE2
		   # 0    1    2    3 
		# 0.98 0.02 0.00 0.00
		rareTot <- names(p[p > 0 & p < 0.05])
		# [1] "1"
		})
	# whichRareTot["chrXXI:16024_C/T"]
	
	# a slick way to comment something out
	if(FALSE){ 
		# Of those that are rare in total, figure out which are also rare (< 5%) within every group
		casesRareTot <- sapply(whichRareTot, length) > 0 # loci with at least one allele rare in total
		whichRareAll <- lapply(alleleFreqByGroup[casesRareTot], function(x){
			# x <- alleleFreqByGroup[["chrXXI:16024_C/T"]]
			# prop <- sweep(x, 1, rowSums(x), FUN = "/")
			                      # 0          1          2          3
			  # paxl       1.00000000 0.00000000 0.00000000 0.00000000
			  # paxb       1.00000000 0.00000000 0.00000000 0.00000000
			  # marine-pac 0.91666667 0.08333333 0.00000000 0.00000000
			# testAll <- apply(prop, 2, function(prop){all(prop < 0.05)})
			testAll <- apply(x, 2, function(x){all(x < 2)})
			rareAll <- names(testAll[testAll])
			# [1] "1" "2" "3"
			})
		# This next command finds the intersection of the alleles that are 
		#	both rare globally and within every group, and saves the result in "whichRareAll"
		whichRareAll <- mapply(whichRareTot[casesRareTot], whichRareAll, FUN = intersect)
		# head(whichRareAll)
		# $`chrXXI:16024_C/T`
		# [1] "1"
		# length(whichRareAll[[1]])
		# [1] 1
		whichRareTot[casesRareTot] <- whichRareAll

		rm(whichRareAll)
		rm(casesRareTot)

		} # END if FALSE

			  
	# View cases that have 3 rare alleles
	# whichRareTot[sapply(whichRareTot, length) > 2]
	# $`chrXXI:475792_G/GCGTCTGCTAAATGACC`
	# [1] "1" "2" "3"
	# $`chrXXI:488280_T/TAAAAATAAGTTCAACACTAATGCAAAATTGAGCAAAATGTTAC`
	# [1] "1" "2" "3"
	# $`chrXXI:1315110_G/GTCTGGCCGGCTCAATGGGTTGTAAGTCACAGCGCGCCGCACAT`
	# [1] "1" "2" "3"
	# head( alleleFreqByGroup[sapply(whichRareTot, length) > 2] )
	# $`chrXXI:475792_G/GCGTCTGCTAAATGACC`            
	              # 0  1  2  3
	  # paxl       22  0  0  0
	  # paxb        4  2  2  2 # all rare benthic alleles
	  # marine-pac 14  0  0  0
	
	# $`chrXXI:488280_T/TAAAAATAAGTTCAACACTAATGCAAAATTGAGCAAAATGTTAC`           
	              # 0  1  2  3
	  # paxl       22  0  0  0
	  # paxb       18  2  1  1 # all rare benthic alleles
	  # marine-pac 12  0  0  0
	
	# $`chrXXI:1315110_G/GTCTGGCCGGCTCAATGGGTTGTAAGTCACAGCGCGCCGCACAT`            
	              # 0  1  2  3
	  # paxl       21  0  0  1
	  # paxb       16  2  2  0 # rare benthic alleles
	  # marine-pac 13  0  0  1 # and one rare marine allele
	
	needFixing <- sapply(whichRareTot, length) > 0
	# length(whichRareTot[needFixing])
	# [1] 81429

	tGT <- as.data.frame(t(geno(vcf)$GT[needFixing, ]), stringsAsFactors = FALSE)
	# ncol(tGT)
	# [1] 81429
	
	z <- mapply(tGT, whichRareTot[needFixing], FUN = function(x, i){
		# x <- tGT[, "chrXXI:9878908_G/GTCGCCGGCCCT"]; i <- whichRareTot[["chrXXI:9878908_G/GTCGCCGGCCCT"]]
		# x
		 # [1] "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1"
		# [19] "1/1" "1/1" "1/2" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "0/3"
		x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		return(x)
		})
	
	# z is a matrix with cols corresponding to tGT
	# z[,"chrXXI:9878908_G/GTCGCCGGCCCT"]
		 # [1] "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1"
		# [19] "1/1" "1/1" NA    "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" "1/1" NA
	z <- t(z)
	geno(vcf)$GT[needFixing, ] <- z
	# geno(vcf)$GT["chrXXI:9878908_G/GTCGCCGGCCCT", ]
	
	rm(tGT)
	rm(z)
	
	# Need to recompute allele frequencies for cases needing fixing (rare alleles dropped)
	# 	(can't just edit the allele frequency tables because genotypes, not alleles, were dropped)
	# alleleFreqByGroup[["chrXXI:9878908_G/GTCGCCGGCCCT"]]
	              # 0  1  2  3
	  # paxl        1 20  0  1
	  # paxb        0 21  1  0
	  # marine-pac  0 14  0  0

	alleleFreqByGroup[needFixing] <- g$tableAlleleFreqByGroup(geno(vcf)$GT[needFixing, ], groupnames, groupcodes)
	# alleleFreqByGroup[["chrXXI:9878908_G/GTCGCCGGCCCT"]]
	              # 0  1  2  3
	  # paxl        0 20  0  0
	  # paxb        0 20  0  0
	  # marine-pac  0 14  0  0

	rm(whichRareTot)
	rm(needFixing)
	}

gc()


# --------------------------
# Determine alleles actually used in genotype calls
# Max of 3 ALT alleles per locus was used.

# Some ALT alleles are labeled as "<*:DEL>". 
# Broad says: "This means there is a deletion in that sample that spans this position, 
#	while in other samples there is another event at this site (typically a SNP)."
# table(unlist(alt(vcf)))[1:4]
# <*:DEL>       A      AA     AAA  ...
   # 8338   73896      21       4  ...  

cat("Determining which ALT alleles actually used in genotype calls; others ALT alleles set to NA\n")

# Enumerate the alt alleles used in actual genotype calls
# "whichAltAllelesUsed" lists all the ALT alleles used in genotypes by their index (1, 2, ...)
# unused ALT alleles are set to NA -- they are NOT DELETED, in order to preserve indices

altUsedList <- g$makeAltUsedList(geno(vcf)$GT, alt(vcf))

# Note, there are still "<*:DEL>" alleles. 

# Note that NA has length 1, so a nonzero nAltUsed doesn't mean there are real ALT alleles there

# Compare the number of ALT alleles used vs number of ALT alleles called
# Number called by GATK
# table(sapply(alt(vcf), length), useNA = "always")
     # 1      2      3      4      5   <NA>
# 267591  13243    763      9      1      0

# Number used in genotype calls (commands commented out)
if(FALSE){
	nAltUsed <- sapply(altUsedList, function(x){ length( x[!is.na(x)] ) })
	table(nAltUsed, useNA = "always")
	     # 0      1      2      3   <NA> # if dropRareAlleles = FALSE
	 # 26102 243415  11398    692      0
	rm(nAltUsed)
	} # End if FALSE

# Types of variants
# -----------------

snpTypeList <- g$makeSnpTypeList(REF = ref(vcf), ALTlist = altUsedList)

# Variant type is defined for VCFtools as (http://vcftools.sourceforge.net/VCF-poster.pdf)
# SNPs
# Alignment  VCF representation
#   ACGT     POS REF ALT
#   ATGT      2   C   T

# Deletions
# Alignment  VCF representation
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

# Examples from earlier analysis where snp type incomplete. 
# ** The way I'm typing them, snp and del at same POS are just alternative alleles, does this make sense????
# 	 It implies that both can't occur, yet they can in principle.
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

# Check that all multi-based snp differ only by the first letter
# Remember that altUsedList still contains "<*:DEL>"
# n1 <- sapply(altUsedList, length)
# n2 <- sapply(snpTypeList, length)
# all(n1 == n2) # should be TRUE
# snp check:
# snp <- unlist(snpTypeList)
# ref <- rep(as.vector(ref(vcf)), n1)
# alt <- unlist(altUsedList)
# snp1 <- snp[!is.na(snp) & snp == "snp"]
# ref1 <- ref[!is.na(snp) & snp == "snp"]
# alt1 <- alt[!is.na(snp) & snp == "snp"]
# table(snp1) # confirms they are all snp
# any( substr(ref1, 1, 1) == substr(alt1, 1, 1) ) # FALSE confirms that the first letter is always different
# all( substr(ref1, 2, nchar(ref1)) == substr(alt1, 2, nchar(alt1)) ) # TRUE confirms that rest is always the same

# del check:
# snp <- unlist(snpTypeList)
# ref <- rep(ref(vcf), sapply(snpTypeList, length))
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
# snp <- unlist(snpTypeList)
# ref <- rep(ref(vcf), sapply(snpTypeList, length))
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
# Transition-transversion ratio

# nrow(geno(vcfresults$vcf)$GT)
# [1] 281607

tstv <- g$vcfTsTv(ref(vcf), altUsedList, snpTypeList)
# tstv <- g$vcfTsTv(ref(vcfresults$vcf), vcfresults$altUsedList, vcfresults$snpTypeList)
# Results here are for dropRareAlleles = TRUE
# Table of variant types used (<*:DEL> and unused are NA)
# snp
   # del    ins    snp   <NA> 
 # 24286  20909 216418  34794 

print(tstv)
# $R
# [1] 1.258801
# $tstv
     # alt
# ref     pur   pyr
  # pur 60301 47947
  # pyr 47864 60306
# $tot
# [1] 216418

gc()

# --------------------------------------
# Save everything to a list for analyses

vcfresults <- list()
vcfresults$groupnames <- groupnames
vcfresults$groupcodes <- groupcodes
vcfresults$control <- control

vcfresults$vcf <- vcf
rm(vcf)

vcfresults$altUsedList <- altUsedList # remember: unused alt alleles are set to NA not deleted
rm(altUsedList)

vcfresults$snpTypeList <- snpTypeList # based on altUsedList
rm(snpTypeList)

vcfresults$alleleFreqByGroup <- alleleFreqByGroup
rm(alleleFreqByGroup)

gc()

# cat("\nSaving results\n")
save(vcfresults, file = vcfresultsfile)
# load(file = vcfresultsfile) # saved object is "vcfresults"

# ----
# Make a biallelic snp version of the data set * Warning: this might drop out whole populations at a marker if has unique alleles
if(saveBiAllelic){
	
	# Identify the ALT alleles that are true snp
	whichAltAreSnp <- lapply(vcfresults$snpTypeList, function(x){
		which(x == "snp")
		})

	# length(whichAltAreSnp)
	# [1] 281607

	# Keep only cases having at least one ALT snp
	keep <- sapply(whichAltAreSnp, length) >= 1
	
	vcf <- vcfresults$vcf[keep]
	tGT <- as.data.frame(t(geno(vcf)$GT), stringsAsFactors = FALSE)
	alleleFreqByGroup <- vcfresults$alleleFreqByGroup[keep]
	altUsedList  <- vcfresults$altUsedList[keep]
	whichAltAreSnp <- whichAltAreSnp[keep]
	
	rm(vcfresults)

	# length(whichAltAreSnp)
	# [1] 214283
	
	# ---
	# Erase the genotypes that are not pure snp
	whichAltAreNotSnp <- lapply(whichAltAreSnp, function(x){setdiff(c(1:nMaxAlt), x)})
	nNotSnp <- sapply(whichAltAreNotSnp, length)

	# table(nNotSnp) # the "0" values are cases in which all Alt alleles are snp
	     # 0      1      2 
	     # 9   2117 212157

	rm(whichAltAreSnp)

	z <- mapply(tGT[, nNotSnp > 0], whichAltAreNotSnp[nNotSnp > 0], FUN = function(x, i){
		# x <- tGT[, "chrXXI:63560_GCGC/G"]; i <- whichAltAreNotSnp["chrXXI:63560_GCGC/G"]
		# x
		 # [1] "0/0" NA    "2/2" NA    "2/2" NA    NA    NA    "1/1" NA    NA    "0/0" "2/2" "2/2" NA    NA    "0/0" NA    NA   
		# [20] NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
		# i
		# [1] 1 3
		x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		return(x)
		})

	# z is a matrix with cols corresponding to tGT
	# z[,"chrXXI:63560_GCGC/G"]
	 # [1] "0/0" NA    "2/2" NA    "2/2" NA    NA    NA    NA    NA    NA    "0/0" "2/2" "2/2" NA    NA    "0/0" NA    NA   
	# [20] NA    NA    NA    NA    NA    NA    NA    NA    NA    NA   

	geno(vcf)$GT[nNotSnp > 0, ] <- t(z)  # fast
	rm(z)
	
	# Redo alleleFreqByGroup
	alleleFreqByGroup[nNotSnp > 0] <- g$tableAlleleFreqByGroup(geno(vcf)$GT[nNotSnp > 0, ], groupnames, groupcodes)

	# alleleFreqByGroup["chrXXI:63560_GCGC/G"]
	             # 0 1 2 3
	  # paxl       2 0 4 0
	  # paxb       2 0 0 0
	  # marine-pac 2 0 4 0

	
	# -----
	# Keep only the two most common alleles (not necessarily the REF)
	
	allelesOrderedByFreq <- lapply(alleleFreqByGroup, function(x){
		# x <- alleleFreqByGroup[[1]]
		z <- colSums(x)
		z[order(z, decreasing = TRUE)]
		})
	# head(allelesOrderedByFreq)
	# $`chrXXI:6316_C/T`
	# 0 1 2 3 
	# 3 3 0 0 
	# $`chrXXI:10364_T/A`
	# 1 0 2 3 # 0 is second if it is equally the rarest
	# 8 0 0 0 
	
	# Identify all alleles other than the 2 most common alleles
	whichAllelesToDrop <- lapply(allelesOrderedByFreq, function(x){
		# x <- allelesOrderedByFreq[[1]]
		z <- names(x[3:(nMaxAlt + 1)])
		})
	
	# Count the number of alleles represented at each marker
	nUsedAlleles <- lapply(allelesOrderedByFreq, function(x){sum(x > 0)})
	# head(nUsedAlleles)
	
	# Drop those corresponding genotypes
	z <- mapply(tGT[, nUsedAlleles > 2], whichAllelesToDrop[nUsedAlleles > 2], FUN = function(x, i){
		# x <- tGT[, "chrXXI:18389_C/A"]; i <- whichAllelesToDrop[["chrXXI:18389_C/A"]]
		 # [1] "0/0" "0/2" "0/0" "0/2" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/1"
		# [17] "0/0" "0/0" "0/1" "0/0" "0/0" "0/1" "0/0" NA    "0/0" "0/0" "0/0" "0/0" "0/1"
		# i
		# [1] "2" "3"
		x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		return(x)
		})

	geno(vcf)$GT[nUsedAlleles > 2, ] <- t(z)
	# unname(geno(vcf)$GT["chrXXI:18389_C/A", ])
	 # [1] "0/0" NA    "0/0" NA    "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/0" "0/1"
	# [17] "0/0" "0/0" "0/1" "0/0" "0/0" "0/1" "0/0" NA    "0/0" "0/0" "0/0" "0/0" "0/1"
	
	# Redo alleleFreqByGroup again
	alleleFreqByGroup[nUsedAlleles > 2] <- g$tableAlleleFreqByGroup(geno(vcf)$GT[nUsedAlleles > 2, ], groupnames, groupcodes)
	
	# alleleFreqByGroup["chrXXI:18389_C/A"]
		
	vcfresultsBiAllelic <- list()
	vcfresultsBiAllelic$groupnames <- groupnames
	vcfresultsBiAllelic$groupcodes <- groupcodes
	vcfresultsBiAllelic$control <- control
	
	vcfresultsBiAllelic$vcf <- vcf
	rm(vcf)
	
	vcfresultsBiAllelic$altUsedList <- altUsedList
	rm(altUsedList)
		
	vcfresultsBiAllelic$alleleFreqByGroup <- alleleFreqByGroup
	rm(alleleFreqByGroup)
	
	gc()
	
	# cat("\nSaving results\n")
	save(vcfresultsBiAllelic, file = vcfBiAllelicFile)
	# load(file = vcfBiAllelicFile) # saved object is "vcfresultsBiAllelic"

	} else rm(vcfresults) 

if(plotQualMetrics){
	# --------------------------------------
	# Plots of quality metrics
	# Current version drops only the rows corresponding to M in chrvec
	# IE, rare alleles not excluded, etc.
	
	vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("QUAL"), geno=c("GT", "GQ", "DP"), 
		info=c("FS", "QD", "DP")))	
	gc()
	keep <- !(start(ranges(vcf)) %in% which.chrvec.M) # i.e., includes only the good bases: only upper case, no "M"
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