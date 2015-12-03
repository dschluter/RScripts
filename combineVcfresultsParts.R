#!/usr/bin/Rscript

# Combines vcfresultsPart files having the same newChr (after Glazerized)
# chrname is the chromosome name whose vcfresultsPart objects are to be joined

# Run in Unix as " Rscript combineVcfResultsParts.R ... "

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname 
# args <- c("BenlimAllMarine", "chrXXI")

project <- args[1]
chrname <- args[2]

chrNumeric <- g$chrname2numeric(chrname)
# [1] 21

# Results file
vcfresultsNewName <- paste(project, chrname, "vcfresultsNew", "rdd", sep = ".")
# [1] "BenlimAllMarine.chrXXI.vcfresultsNew.rdd"

# glob to identify all vcfresultsPart files, they will have chrNumeric (e.g., "21") before the ".rdd"
vcfresultsPartNames <- paste(project, "chr*.vcfresultsPart", chrNumeric, "rdd", sep = ".")
# [1] "BenlimAllMarine.chr*.vcfresultsPart.21.rdd"

# The corresponding file names
z <- list.files( pattern=glob2rx( vcfresultsPartNames ), ignore.case=TRUE )
	# [1] "BenlimAllMarine.chrUn.vcfresultsPart.21.rdd" 
	# [2] "BenlimAllMarine.chrXXI.vcfresultsPart.21.rdd"

# Load the vcfresultsPart files, place in a temporary list: "parts"
parts <- vector("list", length(z)) # initiate
oldChrNames <- sapply( strsplit(z, split = "[.]"), function(x){x[2]})
# [1] "chrUn"  "chrXXI"
names(parts) <- oldChrNames

for(i in 1:length(z)){
	load(z[i])
	parts[[i]] <- vcfresultsPart
	}
	
# We can't simply combine the vcf VariantAnnotation parts because their sequence info is not the same 
# (chrUn.fa vs chrXXI.fa) so just keep as separate list objects of vcfresults,
#  even if there is just one part.
# Initiate
vcfresults <- parts[[1]]
vcfresults$vcf <- list(vcfresults$vcf) # converts vcf part into a list element
names(vcfresults$vcf)[1] <- oldChrNames[1]

if(length(z) > 1){
	
	# Check that names of all the elements are the same
	if(!all( unlist( lapply( parts, function(x){
		setequal( names(parts[[1]]), names(x) )
		}) ) )) stop("Names vary between vcfresultsPart objects")
		
	# Check that the group names are all the same
	if(!all( unlist( lapply( parts, function(x){
		setequal( parts[[1]]$groupnames, x$groupnames )
		}) ) )) stop("groupnames vary between vcfresultsPart objects")
		
	# Could also check that controls are the same, etc
	
	for(i in 2:length(parts)){
		# i <- 2
		vcfresults$vcf[[i]] <- parts[[i]]$vcf
		names(vcfresults$vcf)[i] <- oldChrNames[i]
		vcfresults$newChr <- c(vcfresults$newChr, parts[[i]]$newChr)
		vcfresults$newPos <- c(vcfresults$newPos, parts[[i]]$newPos)			
		vcfresults$altUsedList <- c(vcfresults$altUsedList, parts[[i]]$altUsedList)			
		vcfresults$snpTypeList <- c(vcfresults$snpTypeList, parts[[i]]$snpTypeList)			
		vcfresults$alleleFreqByGroup <- c(vcfresults$alleleFreqByGroup, parts[[i]]$alleleFreqByGroup)			
		}
	}
# names(vcfresults$vcf)
# [1] "chrUn"  "chrXXI"

save(vcfresults, file = vcfresultsNewName) # object is vcfresults
# load(vcfresultsNewName) # object is vcfresults
