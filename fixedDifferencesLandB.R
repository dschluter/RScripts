#!/usr/bin/Rscript

# -------------------------------------------------------------------
# Pare down the vcf file to those cases where all limnetics and benthics are genotyped yet share no alleles

# All individuals from all populations are retained, only loci are dropped

# Obtains stats on one pair of populations
# qsub -I -l walltime=04:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# Filter based on FILTER
# NULL, LowQual, PASS, VQSRTrancheSNP99.00to99.50, VQSRTrancheSNP99.50to99.90, VQSRTrancheSNP99.50to99.90+

# fishPair must uniquely be substrings of the group names (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

project 		<- NULL
chrname  	<- NULL

args <- commandArgs(TRUE)
# args <- c("project=Benlim","chrname=chrXXI")

# Parses the args into a data frame with two columns and then assigns variables 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}
print(x)
       # V1     V2
# 1 project Benlim
# 2 chrname chrXXI

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")

vcfparamFile <- paste(project, "vcfparam.rdd", sep = ".")
load(vcfparamFile) # object is vcfparam

Glazerize <- vcfparam$Glazerize
# GTminFrac <- vcfparam$GTminFrac
fishnames <- vcfparam$fishnames
groupcodes<- vcfparam$groupcodes
groupnames<- vcfparam$groupnames
# nMaxAlt <- vcfparam$nMaxAlt

lims <- c("paxl","pril","qryl","ensl")
bens <- c("paxb","prib","qryb","ensb")

nLims <- sum(vcfparam$nInd[lims])
nBens <- sum(vcfparam$nInd[bens])

print( data.frame(nLims, nBens) )
	
# ----------------
# Manage group codes for this fish pair

library(VariantAnnotation)

# groupnames[groupcodes]
	  # [1] "ensb"       "ensb"       "ensb"       "ensb"       "ensb"      
	  # [6] "ensb"       "ensl"       "ensl"       "ensl"       "ensl"      
	 # [11] "ensl"       "ensl"       "marine-atl" "marine-atl" "marine-atl"
	 # [16] "marine-jap" "marine-jap" "marine-jap" "marine-jap" "marine-pac"
	 # [21] "marine-pac" "marine-pac" "marine-pac" "marine-pac" "marine-pac"
	 # [26] "marine-pac" "marine-pac" "prib"       "prib"       "prib"      
	 # [31] "prib"       "prib"       "prib"       "pril"       "pril"      
	 # [36] "pril"       "pril"       "pril"       "pril"       "pril"      
	 # [41] "pril"       "paxb"       "paxb"       "paxb"       "paxb"      
	 # [46] "paxb"       "paxl"       "paxl"       "paxl"       "paxl"      
	 # [51] "paxl"       "paxl"       "paxl"       "prib"       "prib"      
	 # [56] "prib"       "prib"       "pril"       "pril"       "qryb"      
	 # [61] "qryb"       "qryb"       "qryb"       "qryb"       "qryb"      
	 # [66] "qryl"       "qryl"       "qryl"       "qryl"       "qryl"      
	 # [71] "qryl"       "sculpin"    "sculpin"    "sculpin"    "sculpin"   
	 # [76] "sculpin"    "sculpin"    "sculpin"    "sculpin"    "solitary"  
	 # [81] "solitary"   "solitary"   "solitary"   "solitary"   "solitary"  
	 # [86] "solitary"   "solitary"   "solitary"   "stream"     "stream"    
	 # [91] "stream"     "stream"     "stream"     "stream"     "paxb"      
	 # [96] "paxb"       "paxb"       "paxb"       "paxb"       "paxb"      
	# [101] "paxl"       "paxl"       "paxl"       "paxl"       "paxl"      
	# [106] "paxl"       "paxl"       "paxl"       "paxl"      

# which(is.element(groupnames[groupcodes], lims))
 # [1]   7   8   9  10  11  12  34  35  36  37  38  39  40  41  47  48  49  50  51
# [20]  52  53  58  59  66  67  68  69  70  71 101 102 103 104 105 106 107 108 109

# Redo group codes for the two groups of interest here
# 1 and 2 will indicate the two groups of interest; all other fish are assigned a code of 0
z <- rep(0, length(fishnames))		 # initialize
z[is.element(groupnames[groupcodes], lims)] <- 1
z[is.element(groupnames[groupcodes], bens)] <- 2
groupcodes <- z

# cat("\ngroupcodes:\n")
# print(groupcodes)
  # [1] 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 1 1 1 1
 # [38] 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 1 1 2 2 2 2 2 2 1 1 1 1 1 1 0 0 0
 # [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1

# which(groupcodes > 0)

groups <- groupcodes[groupcodes > 0]
# groups
 # [1] 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1
# [39] 2 2 2 2 1 1 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1

if(Glazerize){
	vcfresultsfile <- paste(project, ".", chrname, ".vcfNew.rdd", sep = "")
	} else {
	vcfresultsfile <- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")}

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcf"

fixedLBfile <- paste(project, chrname, "fixedDiffLB.vcf.rdd", sep = ".")
gtstats <- list()

# ----------
# Count the number of non-missing genotypes and drop locus if not all genotypes
z <- apply( geno(vcf)$GT[ , groupcodes > 0], 1, function(x){sum(!is.na(x))} )

# Drop the loci without complete genotypes
	# nrow(vcf)
	# [1] 928484
vcf <- vcf[z == nLims + nBens]
	# nrow(vcf)
	# [1] 650463

# ----------
# Drop monomorphic loci
# one-way tables:
genoFreq <- apply( geno(vcf)$GT[, groupcodes > 0], 1, function(x){table(x)})
z <- sapply(genoFreq, dim)
# table(z)
# z
     # 1      2      3      4      5      6      7      8      9     10 
# 368097 100870 170942   5166   3795   1396    114     58      9     16 

vcf <- vcf[z > 1]
	# nrow(vcf)
	# [1] 282366

# ----------
# Drop cases where alleles are shared between L and B

# one-way genotype freq tables, separately for L and B:
genoL <- apply( geno(vcf)$GT[, groupcodes ==1], 1, function(x){table(x)})
genoB <- apply( geno(vcf)$GT[, groupcodes ==2], 1, function(x){table(x)})

allelesL <- lapply(genoL, function(x){unique(unlist(strsplit(names(x), split="/")))})
allelesB <- lapply(genoB, function(x){unique(unlist(strsplit(names(x), split="/")))})

z <- mapply(allelesL, allelesB, FUN = function(x,y){!any(x %in% y)})

vcf <- vcf[z]
	# nrow(vcf)
	# [1]


# ----------
cat("\nSaving results\n")
save(vcf, file = fixedLBfile)
# load(file = fixedLBfile) # object is "vcf"
	
