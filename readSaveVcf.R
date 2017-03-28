#!/usr/bin/Rscript

# Makes vcf, containing snp calling results
# Run in Unix as " Rscript readSaveByGroup.R ... "

# In this revised function, only clean vcf is saved, not groupnames etc.
# Useful variables are added to the vcf object, not to a new list containing the vcf object.
# If Glazerize, then ranges are used to extract bits from chrUn etc before saving.

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

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
groupnames 	<- NULL
GTminFrac 	<- "2/3"
Glazerize 	<- "FALSE" 	# Requires file "glazerFileS4 NewScaffoldOrder.csv" in current directory
nMaxAlt		<- 3		# in GATK maximum number of ALT alleles set
dropRareAlleles	<- FALSE
plotQualMetrics <- FALSE
saveBiAllelic 	<- FALSE # saves second data set with 2 snp per marker (not necessarily the REF), removes indels
plotQualMetrics <- FALSE 
genome 			<- "gasAcu1pitx1new.fa"

args <- commandArgs(TRUE)
# args <- c("chrname=chrM", "project=Benlim_99.0_SNP", "groupnames=paxl,paxb,pril,prib,qryl,qryb,ensl,ensb,marine-pac,marine-atl,marine-jap,solitary,sculpin,stream", "Glazerize=TRUE", "GTminFrac=2/3", "postDrop=TRUE")

# Parses the args into a data frame with two columns (V1=left and V2=right of each "=" sign)
# and then assigns V2 to variables whose names are in V1 

x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){ assign(x[i,1], x[i,2]) }

# Check args
# c(chrname, project, groupnames, GTminFrac, Glazerize)
# [1] "chrM"                                                                                            
# [2] "Benlim_99.0_SNP"                                                                                
# [3] "paxl,paxb,pril,prib,qryl,qryb,ensl,ensb,marine-pac,marine-atl,marine-jap,solitary,sculpin,stream"
# [4] "2/3"                                                                                             
# [5] "TRUE"                                                                                         

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")
if(is.null(groupnames)) stop("Provide groupnames= in arguments (grounames separated by commas, no spaces)")

cat("\nchrname is", chrname, "\n")

groupnames <- unlist(strsplit(groupnames, split = ","))
cat("\ngroupnames\n")
cat(groupnames, sep = "\n")
# groupnames
 # [1] "paxl"       "paxb"       "pril"       "prib"       "qryl"      
 # [6] "qryb"       "ensl"       "ensb"       "marine-pac" "marine-atl"
# [11] "marine-jap" "solitary"   "sculpin"    "stream"    

GTminFrac <- eval(parse(text=GTminFrac)) # convert fraction to numeric

# load "chrvec" for the current chromosome
chrno 			<- gsub("^chr", "", chrname)
chrmaskfile     <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd

fastaname		<- paste(chrname, "fa", sep = ".")
# vcfname		<- paste(project, ".", chrname, ".var.vcf.gz", sep="")
vcfname			<- paste(project, ".", chrname, ".vcf.gz", sep="")
vcftbi 			<- paste(project, ".", chrname, ".vcf.gz.tbi", sep="")
vcffile			<- paste(project, ".", chrname, ".vcf.rdd", sep = "")

GTmissing		<- "."	# how GATK represents missing genotypes in the vcf file "./."

control <- list()
control$nMaxAlt <- nMaxAlt
control$snpOptions <- c(dropRareAlleles = dropRareAlleles)
control$GTminFrac <- GTminFrac

cat("\nControl settings on this run\n")
print(control)

library(VariantAnnotation, quietly = TRUE)
# library(WhopGenome, quietly = TRUE)

load(chrmaskfile) 	# object is named "chrvec"
which.chrvec.M <- which(chrvec == "M") # bases coded as "M" for masked

# object.size(chrvec)
#  93,740,176 bytes # chrXXI
# 500,401,968 bytes # chrUn

rm(chrvec)

# Garbage collection -- current memory consumption
gcinfo(TRUE) # sends messages about memory garbage collection
gc()

# -------------------------
# Build index, if not exist

if(!file.exists(vcftbi)){
	root <- sub("[.]gz$", "", vcfname) # remove the "gz" if present
	compressedVcfname <- paste(root, "bgz", sep = ".")
	# compress and save file name
	bgzname <- bgzip(vcfname, compressedVcfname)
	file.rename(bgzname, vcfname)
	# Create the index
	idx <- indexTabix(vcfname, "vcf")
	}

# ---------------------
# Read variant VCF file

# Use "vcfLong <- expand(vcf)" if want variants having multiple ALT alleles to appear on multiple rows, 
# one per ALT allele. This would be a straightforward way to get predictCoding for all alleles. 

# Useful accessor functions
# To print many rows, do this first:
# options(showHeadLines=Inf)  # allows you to print many lines
# options("showHeadLines"=NULL) # reverts back to default

# header(vcf) 				# Header information
# samples(header(vcf))		# names of fish in sample
# geno(header(vcf))			# info on GT, AD, DP, GQ, MIN_DP, PGT, PID PL RGQ SB
# rowData(vcf)				# paramRangeID, ALT, QUAL, FILTER
# mcols(vcf)				# same
# rowRanges(vcf)			# seqnames, ranges, strant, paramRangeID, REF, ALT QUAL FILTER
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

# Reads the whole vcf file
# vcf <- readVcf(file = vcfname, genome = fastaname) 

# object.size(vcf)

# test
# rngs <- GRanges("chrM", IRanges(c(1, 1000), c(2000, 3000))) # just these lines
# rng1 <- GRanges("chrM", IRanges(1, 1000)) 
# rng2 <- GRanges("chrM", IRanges(2000, 3000))

# vcf <- readVcf(file = vcfname, genome = fastaname, ScanVcfParam(fixed=c("ALT", "QUAL", "FILTER"), 
			# which = rngs, geno="GT", info=c("DP","FS","VQSLOD")))
vcf <- readVcf(file = vcfname, genome = genome, ScanVcfParam(fixed = c("ALT", "QUAL", "FILTER"), 
			# which = c(rng1, rng2), 
			geno="GT", info=c("AC", "DP","FS","VQSLOD")))

# object.size(vcf)

cat("\nSuccessfully read vcf file\n")

# -----
# fish group codes 1, 2, ... Order is determined by sequence of groupnames, eg "paxl", "paxb", "marine-pac"

fishnames <- samples(header(vcf))

# FUDGE TO FIX PRIEST in file names
# fishnames <- gsub("Priest", "Pri", fishnames, ignore.case = TRUE) 
# These two commands should change the names where they are stored in vcf (need both)
# header(vcf)@samples <- fishnames
# dimnames(vcf)[[2]] <- fishnames

cat("\nfishnames:\n")
print(fishnames)

# assign numbers corresponding to groups
groupcodes <- g$groupFishCodes(fishnames, groupnames)

cat("\ngroupcodes:\n")
print(groupcodes)  # 
  # [1]  8  8  8  8  8  8  7  7  7  7  7  7 10 10 10 11 11 11 11  9  9  9  9  9  9
 # [26]  9  9  4  4  4  4  4  4  3  3  3  3  3  3  3  3  2  2  2  2  2  1  1  1  1
 # [51]  1  1  1  4  4  4  4  3  3  6  6  6  6  6  6  5  5  5  5  5  5 13 13 13 13
 # [76] 13 13 13 13 12 12 12 12 12 12 12 12 12 14 14 14 14 14 14  2  2  2  2  2  2
# [101]  1  1  1  1  1  1  1  1  1

# table(groupcodes)
 # 1  2  3  4  5  6  7  8  9 10 11 12 13 
# 16 11 10 10  6  6  6  6  8  3  4  8  6 
 
nInd <- as.vector(table(groupcodes)) # number of individuals genotyped in each group
 # [1] 16 11 10 10  6  6  6  6  8  3  4  8  6

cat("\nNumber of individuals (nInd):\n")
names(nInd) <- groupnames

print(nInd) 
      # paxl       paxb       pril       prib       qryl       qryb       ensl 
        # 16         11         10         10          6          6          6 
      # ensb marine-pac marine-atl marine-jap   solitary    sculpin     stream 
         # 6          8          3          4          9          8          6 

nMin <- round(GTminFrac*nInd) # minimum number required in each group
# nMin
      # paxl       paxb       pril       prib       qryl       qryb       ensl 
        # 11          7          7          7          4          4          4 
      # ensb marine-pac marine-atl marine-jap   solitary    sculpin     stream 
         # 4          5          2          3          6          5          4 

nMin <- mapply(rep(5, length(groupnames)), nMin, FUN = min) # criterion is 5 or nMin, whichever is smaller, eg 5 5 4
names(nMin) <- groupnames

cat("\nMinimum no. genotyped individuals needed to use a snp when calculating Fst etc (nMin):\n")
print(nMin)  # 
      # paxl       paxb       pril       prib       qryl       qryb       ensl 
         # 5          5          5          5          4          4          4 
      # ensb marine-pac marine-atl marine-jap   solitary    sculpin     stream 
         # 4          5          2          3          5          5          4 

# Q: What does QUAL = NA imply?

# ------
# Drop masked SNP (variants whose start position is masked in the chromosome)
# ** applied to CHRIV, this step removes base "12811481", the putative site of the Eda mutation **
# 	In windowsmaskerSdust.chrIV, the site is marked as masked.

keep <- !(start(ranges(vcf)) %in% which.chrvec.M) # i.e., includes only the good bases: only upper case, no "M"

vcf <- vcf[keep]
cat("\nCompleted removal of snp corresponding to masked bases\n\n")

rm(keep)
# rm(which.chrvec.M) # remove later

# object.size(vcf)
# 36,360,752 bytes # chrXXI
gc()

# ------
# set missing genotypes to NA instead of "."
geno(vcf)$GT[geno(vcf)$GT == GTmissing] <- NA

# ------
print(table(filt(vcf)))
                          # .                     LowQual                        PASS 
                        # 389                           6                         104 
 # VQSRTrancheSNP99.00to99.50  VQSRTrancheSNP99.50to99.90 VQSRTrancheSNP99.50to99.90+ 
                         # 13                          89                        3739 

# -----
# Drop variants in which fewer than 2 groups have at least one genotype
*** and drop LowQual in FILTER field

# Tally up the number of genotypes in each group
samplesByGroup <- split(samples(header(vcf)), groupcodes) # Order of resulting groups is as in groupnames and groupcodes

samplesByGroup
	# $`1`
	 # [1] "PaxLim-PxCL09maleLM1-GS14"     "PaxLim-PxLfemale6-GS18"       
	 # [3] "PaxLim-PxLmale102-GS16"        "PaxLim-PxLmale106-GS15"       
	 # [5] "PaxLim-PxLmale107-GS17"        "PaxLim-formerlyPrLfemale1-GS8"
	 # [7] "PaxLim-formerlyPrLmale5-GS5"   "paxl01"                       
	 # [9] "paxl02-formerlyPril02"         "paxl05"                       
	# [11] "paxl09"                        "paxl10"                       
	# [13] "paxl13"                        "paxl14"                       
	# [15] "paxl20-formerlyPril20"         "paxl21-formerlyPril21"        
	
	# $`2`
	 # [1] "PaxBen-PxBmale5-GS11"        "PaxBen-PxBmale6-GS12"       
	 # [3] "PaxBen-PxBmale8-GS10"        "PaxBen-PxCL09femaleBF6-GS13"
	 # [5] "PaxBen-RPxCL09maleBM2-GS9"   "paxb04"                     
	 # [7] "paxb05"                      "paxb06"                     
	 # [9] "paxb07"                      "paxb08"                     
	# [11] "paxb09"                     
	
	# $`3`
	 # [1] "PRIL101"                     "PRIL102"                    
	 # [3] "PRIL104"                     "PRIL108"                    
	 # [5] "PRIL112"                     "PRIL16"                     
	 # [7] "PRIL17"                      "PRIL18"                     
	 # [9] "PriestLim-PrCL09maleLM3-GS7" "PriestLim-PrLmale2-GS6"     
	
	# $`4`
	 # [1] "PRIB02"                       "PRIB05"                      
	 # [3] "PRIB06"                       "PRIB07"                      
	 # [5] "PRIB11"                       "PRIB15"                      
	 # [7] "PriestBen-PrBfemale5-GS4"     "PriestBen-PrBmale3-GS2"      
	 # [9] "PriestBen-RPrCL09maleBM4-GS3" "PriestBen-RPrCL09maleBM6-GS1"
	
	# $`5`
	# [1] "QRYL04" "QRYL05" "QRYL07" "QRYL08" "QRYL09" "QRYL10"
	
	# $`6`
	# [1] "QRYB01" "QRYB06" "QRYB08" "QRYB11" "QRYB13" "QRYB25"
	
	# $`7`
	# [1] "ENSL172" "ENSL24"  "ENSL25"  "ENSL33"  "ENSL37"  "ENSL50" 
	
	# $`8`
	# [1] "ENSB01" "ENSB03" "ENSB08" "ENSB12" "ENSB15" "ENSB23"
	
	# $`9`
	# [1] "Marine-Pac-BIGR-52_54_2008-02"  "Marine-Pac-Japan-01-Katie"     
	# [3] "Marine-Pac-LITC_0_05_2008-FO"   "Marine-Pac-LittleCampbell-LC1D"
	# [5] "Marine-Pac-MANC_X_X05"          "Marine-Pac-Oyster-06-Sara"     
	# [7] "Marine-Pac-Seyward-01-Sara"     "Marine-Pac-WestCreek-01-Sara"  
	
	# $`10`
	# [1] "Marine-Atl-BITJ_X_X17"           "Marine-Atl-Denmark-BS27-Fuelner"
	# [3] "Marine-Atl-TYNE_1_2001-14"      
	
	# $`11`
	# [1] "Marine-JapSea-fem-Kitano-KIA31" "Marine-JapSea-fem-Kitano-KIA46"
	# [3] "Marine-JapSea-fem-Kitano-KIA48" "Marine-JapSea-male01-Katie"    
	
	# $`12`
	# [1] "Solitary-Black-BL4-Sara"  "Solitary-Bullock-19-Sara"
	# [3] "Solitary-Cranby-04-Sara"  "Solitary-Hoggan-13-Sara" 
	# [5] "Solitary-Kirk-12-Sara"    "Solitary-Stowell-04-Sara"
	# [7] "Solitary-Tom-01-Sara"     "Solitary-Trout-01-Sara"  
	
	# $`13`
	# [1] "SculpinLake-Ambrose-01-Sara" "SculpinLake-Brown05-Sara"   
	# [3] "SculpinLake-Cedar-01-Sara"   "SculpinLake-North-05-Sara"  
	# [5] "SculpinLake-Ormond08-Sara"   "SculpinLake-Pachena03-Sara" 
	# [7] "SculpinLake-Paq19-Sara"      "SculpinLake-Rosseau1F-Sara" 
	
	# $`14`
	# [1] "Stream-LittleCampbell-01-Sara"       
	# [2] "Stream-LittleCampbell_23_32_2008-306"
	# [3] "Stream-LittleCampbell_23_32_2008-324"
	# [4] "Stream-LittleCampbell_23_32_2008-347"
	# [5] "Stream-LittleCampbell_23_32_2008-356"
	# [6] "Stream-LittleCampbell_23_32_2008-744"


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

# names(nCalledGenotypes)
 # [1] "paxl"       "paxb"       "pril"       "prib"       "qryl"      
 # [6] "qryb"       "ensl"       "ensb"       "marine-pac" "marine-atl"
# [11] "marine-jap" "solitary"   "stream"    
# head(nCalledGenotypes$paxl)
 # chrM:62_G/A  chrM:90_G/A  chrM:91_A/G chrM:130_A/G chrM:140_G/T chrM:148_T/C 
          # 13           15           16           15           15           15

z2 <- lapply(nCalledGenotypes, function(z1){z1 >= 1})
z2 <- data.frame(z2)
z3 <- apply(z2, 1, sum)
keep <- z3 >= 2

cat("\nKeeping only snp that have at least one genotype in at least 2 groups\n\n")
vcf <- vcf[keep]

rm(nCalledGenotypes)
rm(z2)
rm(z3)
rm(keep)

gc()

# --------------------------------------
# Make a table of the allele frequencies by group
if(FALSE){
alleleFreqByGroup <- g$tableAlleleFreqByGroup(geno(vcf)$GT[, groupcodes > 0], groupnames, groupcodes[groupcodes > 0])
	}
	# unname(geno(vcf)$GT[1, ])
	# groupcodes
	# alleleFreqByGroup[[1]]     
	              # 0  1  2  3
	  # paxl       26  0  0  0
	  # paxb       18  0  0  0
	  # pril       16  2  0  0
	  # prib       16  0  0  0
	  # qryl        6  0  0  0
	  # qryb        8  0  0  0
	  # ensl        6  0  0  0
	  # ensb        0  0  0  0
	  # marine-pac 14  0  0  0
	  # marine-atl  0  0  0  0
	  # marine-jap  8  0  0  0
	  # solitary   12  0  0  0
	  # sculpin    12  0  0  0
	  # stream      4  0  0  0

gc()


# --------------------------------------
# Drop rare alleles if dropRareAlleles is TRUE
#  -- we are using ** RULE2 ** otherwise we get lots of rare alleles when many genotypes are missing
# RULE1: drop alleles less than 5% total, measured as percent of 2*sum(nInd), the maximum number of alleles possible
# RULE2: drop alleles less than 5% total, measured as percent of the total number of alleles in the sample.

# Note that we can get ALT allele counts from info(vcf)$AC (just need to get REF count)

# if(dropRareAlleles){
	# g$dropRareAlleles() # Not working yet as a function
	# }

# --------------------------
# Determine alleles actually used in genotype calls

# Max of 3 ALT alleles per locus was used -- contained in alt(vcf)

# Some ALT alleles are labeled as "<*:DEL>". 
# Broad says: "This means there is a deletion in that sample that spans this position, 
#	while in other samples there is another event at this site (typically a SNP)."
# table(unlist(alt(vcf)))[1:4]
# <*:DEL>       A      AA     AAA  ...
   # 8338   73896      21       4  ...  

# Enumerate the alt alleles used in actual genotype calls
# "whichAltAllelesUsed" lists all the ALT alleles used in genotypes by their index (1, 2, ...)
# Unused ALT alleles are set to NA -- they are NOT DELETED, in order to preserve indices
# Note, there are still "<*:DEL>" alleles. 
# Note that NA has length 1, so a nonzero nAltUsed doesn't mean there are real ALT alleles there

# Can also calculate this more quickly using info(vcf)$AC: a "0" implies no called genotypes use it
cat("\n\nDetermining which ALT alleles actually used in genotype calls; others ALT alleles set to NA\n")
altUsedList <- g$makeAltUsedList(geno(vcf)$GT[ , groupcodes > 0], alt(vcf))

newInfo <- DataFrame(Number=1, Type="CompressedList",
                      Description="ALT alleles actually used in genotype calls",
                      row.names="altUsedList")
info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
# rownames(info(header(vcf)))

info(vcf)$altUsedList <- as(unname(altUsedList), "CompressedList")

# This didn't work
# mcols(vcf)$altUsedList <- as(unname(altUsedList), "CompressedList")
	
# Compare the number of ALT alleles used vs number of ALT alleles called
# Number called by GATK
# table(sapply(alt(vcf), length), useNA = "always")
     # 1      2      3      4      5   <NA>
# 267591  13243    763      9      1      0

# Number of ALT alleles actually used in genotype calls ( if(FALSE) comments everything out )
	# nAltUsed <- sapply(altUsedList, function(x){ length( x[!is.na(x)] ) })
	# table(nAltUsed, useNA = "always")
	     # # 0      1      2      3   <NA> # if dropRareAlleles = FALSE
	 # # 26102 243415  11398    692      0
	# rm(nAltUsed)

# Types of variants
# -----------------

cat("\nMaking snp type list based on alt alleles actually used\n")
snpTypeList <- g$makeSnpTypeList(REF = ref(vcf), ALTlist = altUsedList)

newInfo <- DataFrame(Number=1, Type="CompressedList",
                      Description="Type of ALT allele",
                      row.names="snpTypeList")
info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
# rownames(info(header(vcf)))

info(vcf)$snpTypeList <- as(unname(snpTypeList), "CompressedList")
# info(vcf)

rm(snpTypeList)
rm(altUsedList)

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

tstv <- g$vcfTsTv(ref(vcf), altUsedList, snpTypeList)

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

# Save vcf file

save(vcf, file = vcffile)
# load(vcffile) # object is "vcf"