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

# Find the corresponding files, chosen based on the "21" before the ".rdd"
z <- list.files( pattern=glob2rx( paste(project, "chr*.vcfresultsPart", chrNumeric, "rdd", sep = ".") ), ignore.case=TRUE )
	# [1] "BenlimAllMarine.chrUn.vcfresultsPart.21.rdd" 
	# [2] "BenlimAllMarine.chrXXI.vcfresultsPart.21.rdd"

oldChrNames <- sapply( strsplit(names(parts), split = "[.]"), function(x){x[2]})
# [1] "chrUn"  "chrXXI"

parts <- vector("list", length(z)) # initiate
names(parts) <- oldChrNames

# Load the vcfresultsPart files as list elements
for(i in 1:length(z)){
	load(z[i])
	parts[[i]] <- vcfresultsPart
	}
	
# We can't combine the vcf parts because sequence info is not the same (chrUn.fa vs chrXXI.fa)
# so just keep as a list of vcf objects, even if there is just one part.
x <- parts[[1]]
x$vcf <- list(x$vcf) # converts vcf part into a list element
names(x$vcf)[1] <- oldChrNames[1]

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
		x$vcf[[i]] <- parts[[i]]$vcf
		names(x$vcf)[i] <- oldChrNames[i]
		x$newChr <- c(x$newChr, parts[[i]]$newChr)
		x$newPos <- c(x$newPos, parts[[i]]$newPos)			
		x$altUsedList <- c(x$altUsedList, parts[[i]]$altUsedList)			
		x$snpTypeList <- c(x$snpTypeList, parts[[i]]$snpTypeList)			
		x$alleleFreqByGroup <- c(x$alleleFreqByGroup, parts[[i]]$alleleFreqByGroup)			
		}

vcfresultsList <- x
save(vcfresultsList, file = "") # object is vcfresultsList
# load("") # object is vcfresultsList
