#!/usr/bin/Rscript

# Makes vcf, containing snp calling results
# Run in Unix as " Rscript readSaveByGroup.R ... "

# At most 3 ALT alleles permitted per snp in this version
# NOTE: windowmaskerSdust start position is 0-based, so added 1 (end position is 1-based!)

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# Expect to read these arguments from args
project 	<- NULL
chrname 	<- NULL

args <- commandArgs(TRUE)
# args <- c("chrname=chrXXI", "project=Benlim_99.0_SNP")

# Parses the args into a data frame with two columns (V1=left and V2=right of each "=" sign)
# and then assigns V2 to variables whose names are in V1 

x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){ assign(x[i,1], x[i,2]) }

# Check args
# print(chrname); print(project)
# [1] "chrM"                                                                                            
# [1] "Benlim_99.0_SNP" 

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")

cat("\nchrname is", chrname, "\n") # chrXXI

chrno 			<- gsub("^chr", "", chrname)
vcfparamFile <- paste(project, "vcfparam.rdd", sep = ".")
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

load(vcfparamFile) # object is vcfparam
Glazerize <- vcfparam$Glazerize

library(VariantAnnotation, quietly = TRUE)
# library(WhopGenome, quietly = TRUE)

# ------
# in progress
# Generate the new assembly vcf objects by reading selectively from old assembly vcf files,
# 	converting coordinates, and combining before saving.

# if(Glazerize){
	# if(chrno != "M" & !grepl("pitx1", chrno)){
		# glazer <- read.csv(Glazerfile, stringsAsFactors = FALSE)
		# glazer <- glazer[glazer$NewChr == chrNumeric, ] # all those rows that are converted to current chrno
		# chrFrom <- unique(glazer$OldChr) # [1] "Un" "21"
		# vcfList <- vector(mode = "list", length = 2) 
		# names(vcfList) <- chrFrom

# -------
# Or do this the old way:

# Read and process the vcf file for the specified chromosome

vcf <- g$readChrVcf(project, chrname, vcfparamFile)

	# Max of 3 ALT alleles per locus was used -- contained in alt(vcf)
	
	# Some ALT alleles are labeled as "<*:DEL>". 
	# Broad says: "This means there is a deletion in that sample that spans this position, 
	#	while in other samples there is another event at this site (typically a SNP)."
	# table(unlist(alt(vcf)))[1:4]
	# <*:DEL>       A      AA     AAA  ...
	   # 8338   73896      21       4  ...  

# Useful accessor functions

	# Use "vcfLong <- expand(vcf)" if want variants having multiple ALT alleles to appear on multiple rows, 
	# one per ALT allele. This would be a straightforward way to get predictCoding for all alleles. 
	
	# To print many rows, do this first:
	# options(showHeadLines=Inf)  # allows you to print many lines
	# options("showHeadLines"=NULL) # reverts back to default
	
	# header(vcf) 				# Header information
	# samples(header(vcf))		# names of fish in sample
	# geno(header(vcf))			# info on GT, AD, DP, GQ, MIN_DP, PGT, PID PL RGQ SB
	# rowData(vcf)				# paramRangeID, ALT, QUAL, FILTER
	# rowRanges(vcf)			# seqnames, ranges, strand, paramRangeID, REF, ALT QUAL FILTER
	# mcols(vcf)				# same
	# ref(vcf)					# extract the REF 
	# alt(vcf)					# ALT alleles (DNAStringSetList)
	# qual(vcf)	 				# SNP quality
	# filt(vcf)					# FILTER
	# geno(vcf)					# Genotype data: List of length 10 names(10): GT AD DP GQ MIN_DP PGT PID PL RGQ SB
	#							# 	RGQ is "Unconditional reference genotype confidence"
	# info(vcf)					# INFO field variables
	# info(vcf)$AC				# Allele count for each ALT allele in the same order as listed
	# info(vcf)$AF				# Allele frequency (fraction) for each ALT allele in the same order as listed
	# info(vcf)$AN				# Total number of alleles in called genotypes
	# info(vcf)$VQSLOD			# VQSLOD
	# info(vcf)$VariantType		# all NA, leave out of readVcf


# -----------------
# The alt alleles used in actual genotype calls
# info(vcf)$altUsedList

# -----------------
# Types of variants
# info(vcf)$snpTypeList

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

# -----------------
# Transition-transversion ratio
# tstv <- g$vcfTsTv(ref(vcf), info(vcf)$altUsedList, info(vcf)$snpTypeList)
# print(tstv)


# Save vcf file
save(vcf, file = vcfresultsfile)	# saved object is "vcf"
# load(file = vcfresultsfile) 

# Save parts of vcf needed to build new assembly vcf files

if(Glazerize){ # Requires conversion file "glazerFileS4 NewScaffoldOrder.csv" in current working directory

	# This doesn't work on the mac (returns "1 1" instead of POS)
	#   use rowRanges instead (different R version?)
	pos <- start(rowData(vcf)) 
	
	if(chrno != "M" & !grepl("pitx1", chrno) ){
		newCoords <- g$glazerConvertOld2New(chrname, pos)
		} else {
		newCoords <- data.frame(newChr = rep(chrno, length(pos)), newPos = pos)
		}

	# Put the new coordinates in the INFO field of vcf file
	newInfo <- DataFrame(Number=1, Type="Character",
	                      Description="Glazer assembly chromosome",
	                      row.names="newChr")
	info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
	info(vcf)$newChr <- newCoords$newChr

	newInfo <- DataFrame(Number=1, Type="Integer",
	                      Description="Glazer assembly chromosome positions",
	                      row.names="newPos")
	info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
	info(vcf)$newPos <- newCoords$newPos

	
	z <- unique(info(vcf)$newChr)
	# [1] "21" "Un"
	
	for(i in z){ # saved object is "vcfPart"
		vcfPart <- vcf[info(vcf)$newChr == i]
		save(vcfPart, file = paste(project, chrname, "vcfresultsPart", i, "rdd", sep = "."))
		}
	} # end if(Glazerize)

