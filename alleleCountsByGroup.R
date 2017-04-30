#!/usr/bin/Rscript

# ----------------------------------------------------------------
# Tabulate the numbers of 0, 1, 2, 3 alleles for each SNP by group

# qsub -I -l walltime=04:00:00 -l mem=4gb # work interactively; use "exit" to exit
# cd ~/BenlimResults/
# module load R/3.1.2
# R

# NULL, LowQual, PASS, VQSRTrancheSNP99.00to99.50, VQSRTrancheSNP99.50to99.90, VQSRTrancheSNP99.50to99.90+

# fishPair must uniquely be substrings of the group names (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

project 	<- NULL
chrname  	<- NULL

args <- commandArgs(TRUE)
# args <- c("project=Benlim", "chrname=chrM", "groupsToDo=paxb,prib,qryb,ensb")
# args <- c("project=Benlim", "chrname=chrXXI", "groupsToDo=paxb,prib,qryb,ensb")

# Parses the args into a data frame with two columns and then assigns variables 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}
print(x)
          # V1                   V2
# 1    project               Benlim
# 2    chrname                 chrM
# 3 groupsToDo paxb,prib,qryb,enosb

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")
if(is.null(groupsToDo)) stop("Provide groupsToDo= in arguments")

groupsToDo <- unlist(strsplit(groupsToDo, split = "[ ,]"))

vcfparamFile <- paste(project, "vcfparam.rdd", sep = ".")
load(vcfparamFile) # object is vcfparam

Glazerize 	<- vcfparam$Glazerize
GTminFrac 	<- vcfparam$GTminFrac
fishnames 	<- vcfparam$fishnames
groupcodes	<- vcfparam$groupcodes
groupnames	<- vcfparam$groupnames
nMaxAlt 	<- vcfparam$nMaxAlt
nMin 		<- vcfparam$nMin
	
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

if(Glazerize){
	vcfresultsfile <- paste(project, ".", chrname, ".vcfNew.rdd", sep = "")
	} else {
	vcfresultsfile <- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")}

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcf"

for(i in groupsToDo){
	# i <- "paxb"
	
	z <- which(groupnames[groupcodes] == i)
	
	# alleles is a list
	alleles <- apply(geno(vcf)$GT[ , z], 1, function(x){unlist(strsplit(unname(x), split="[/|]"))})

	# tabulate the number of alleles from 0 to nMaxAlt (the preset maximum number of ALT alleles called)
	# tapply seems to be faster than table
	
	# alleleCount <- lapply(alleles, function(x){
		# table(factor(x, levels = 0:nMaxAlt))
		# })
	# system.time for chrXXI and "paxb"
	   # user  system elapsed 
	# 758.133   0.498 758.815 

	# This was fastest
	# prints "NA" instead of "0" when an allele has no 
	alleleCount <- lapply(alleles, function(x){
		tapply(x, factor(x, levels = 0:nMaxAlt), length)
		})
	
	# system.time for chrXXI and "paxb"
	   # user  system elapsed 
	# 264.930   0.211 265.291 


	# # Maybe not using mclapply properly (haven't specified resources)
	# alleleCount <- mclapply(alleles, FUN=function(x){
		# table(factor(x, levels = 0:nMaxAlt))
		# }, mc.cores=12)
	# system.time for chrXXI and "paxb"
	   # user  system elapsed 
	 # 23.776   1.385 812.254 
	
	# remove the snp names to save storage memory
	alleleCount <- unname(alleleCount)
	
	# To save memory, replace alleleCount list element with NA if total alleles less than nMin
	alleleCount <- lapply(alleleCount, function(x){
		if(sum(x, na.rm=TRUE) < nMin[i]*2) x <- NA
		x
		})

	# Put the new list in the INFO field of vcf file
	newInfo <- DataFrame(Number=1, Type="list",
	                      Description=paste("Allele count for", i),
	                      row.names = i)
	info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
	info(vcf)[[i]]<- alleleCount # will name new column the value of i (eg "paxb")

	}


# ----------
cat("\nSaving results back to vcfresultsfile\n")
save(vcf, file = vcfresultsfile)
# load(file = vcfresultsfile) # object is "vcf"
	
