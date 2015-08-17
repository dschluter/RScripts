g<-list()

 # [1] "makeAltUsedList"        "whichAltAllelesUsed"    "qsubRscriptPbs"         "makeChrvecMasked"      
 # [5] "qsubGenotypeGVCFsPbs"   "downloadSra"            "recode"                 "geno.diff"             
 # [9] "psd"                    "blockstats"             "slidewin"               "plot.reads"            
# [13] "peel"                   "flag.readable"          "gtf2thirdpositions"     "tstv"                  
# [17] "write.bedGraph"         "wc.revised"             "makeSnpTypeList"        "tableAlleleFreqByGroup"
# [21] "vcfTsTv"               

g$vcfTsTv <- function(REF, ALTlist, snpTypeList){
	# Calculates the raw transition-transversion ratio from the REF and list of ALT alleles, making
	# use of the snp types already provided in snpTypeList
	#
	nAltUsed <- sapply(ALTlist, function(x){length(x[!is.na(x)])})
	snp <- unlist(snpTypeList[nAltUsed > 0])
	ref <- rep(as.vector(REF), nAltUsed)
	alt <- unlist(ALTlist[nAltUsed > 0])
	cat("Table of variant types used (<*:DEL> is NA)\n")
	print(table(snp, useNA="always"))
	
	# Transition-transversion ratios - true snp only (not indels)
	cat("\nTransition-transversion ratio - all true snp\n")
	snp1 <- snp[!is.na(snp) & snp == "snp"]
	ref1 <- ref[!is.na(snp) & snp == "snp"]
	alt1 <- alt[!is.na(snp) & snp == "snp"]
	
	ref1 <- substr(ref1, 1, 1) # keep the first letter of REF
	alt1 <- substr(alt1, 1, 1) # keep the first letter of ALT
	
	# c( length(ref1[ref1 != alt1]), length(alt1[ref1 != alt1]) )
	# [1] 216418 216418
	
	print(g$tstv(ref1[ref1 != alt1], alt1[ref1 != alt1]))
	}

g$tableAlleleFreqByGroup <- function(GT = geno(vcf)$GT, groupnames, groupcodes, nMaxAlt = 3, split = "/"){
	# Generates a table of allele frequencies at each locus by group
	# Group is specified by the integer groupcodes, eg 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	# 
	groupnames <- factor(groupnames, levels = groupnames)
	tGT <- as.data.frame(t(GT), stringsAsFactors = FALSE) # annoying but apply won't return a list
	alleleFreqByGroup <- lapply(tGT, function(locus){
		# columns of tGT are loci, so apply function locus by locus
		# locus <- geno(vcf)$GT[6, ]
		z1 <- split(locus, groupcodes) # the genotypes for each group at the locus
		z2 <- lapply(z1, function(x){ 
			unlist(strsplit(x, split = split))
			}) 		# the alleles for each group at the locus
		z3 <- z2	# initialize
		for(i in 1:length(z3)) z3[[i]] <- rep( groupnames[i], length(z3[[i]]) ) # replace z3 with corresponding group name
		table( unlist(z3), factor(unlist(z2), levels=0:nMaxAlt) ) # factor so that all alleles are counted in each table
		})
	}

g$makeSnpTypeList <- function(REF = ref(vcf), ALTlist = altUsedList, AltDelSymbol = "<*:DEL>"){
	# Creates a list indicating what snp type is each ALT alleles: snp, ins, or del
	i <- nchar(REF) # length of the reference sequence
	i <- split(i, 1:length(i)) # split REF allele size into a list
	j <- lapply(ALTlist, function(x){
		j <- nchar(x)            	# length of ALT alleles in altList; 
									# is integer(0) if altUsedList[[i]] is character(0)
									# but is 2 if altUsedList[[i]] is NA
		j[x == AltDelSymbol] <- NA  # length of "<*:DEL>" set to NA
		j[is.na(x)] <- NA
		return(j)
		})
	snpTypeList <- mapply(i, j, FUN = function(i, j){ # i and j refer to nchar of ref and alt
		snptype <- rep(NA, length(j)) # initialize with NAs
		# snptype[i == 1 & j == 1] <- "snp" # first run
		# snptype[i < j  & i == 1] <- "ins" # first run
		# snptype[i > j  & j == 1] <- "del" # first run
		#
		# 2nd run below uses broader criteria to designate variant type
		# Note: j is based on ALTlist, so snpTypeList is also
		snptype[i == j] <- "snp" # 2nd run
		snptype[i <  j] <- "ins" # 2nd run
		snptype[i >  j] <- "del" # 2nd run
		snptype[length(j) == 0] <- NA # snp type is NA if length of altUsedList is 0, ie no snp used
		return(snptype)
		})
	names(snpTypeList) <- names(ALTlist)
	return(snpTypeList)
	}

g$makeAltUsedList <- function(GT, ALT = alt(vcf)){
	# Returns a list of the ALT alleles actually used in the genotypes in the sample
	# The list is based on alt(vcf) but unused ALT alleles are set to NA 
	# 	- NB they are NOT DELETED to preserve indices
	# GT is geno(vcf)$GT or a subset with rows as loci and columns as individuals
	#
	if(nrow(GT) != length(ALT)) stop("Number of rows of GT must equal length of ALT")
	z <- as.list(ALT) # mapply takes forever if use alt(vcf) instead of as.list(alt(vcf))
	whichAltAllelesUsed <- g$whichAltAllelesUsed(GT)
	altUsedList <- mapply(whichAltAllelesUsed, z, FUN = function(x, y){
		# x <- whichAltAllelesUsed[[7]]; y <- alt(vcf)[[7]]
		ialt <- 1:length(y)
		y[!(ialt %in% x)] <- NA
		return(y)
		})
	return(altUsedList)
	}

g$whichAltAllelesUsed <- function(GT, split = "/"){
	# returns a list indicating which alleles are actually used in the genotypes 
	# 	in the sample by their index (1, 2, ...
	# GT is geno(vcf)$GT or a subset with rows as loci and columns as individuals
	whichAltAllelesUsed <- apply(GT, 1, function(x){ 		# is a list
			# x <- GT[1, ]
			x1 <- strsplit(x[!is.na(x)], split = split)
			x1 <- sort( as.integer(unique(unlist(x1))) )
			whichAltAllelesUsed <- x1[x1 > 0] # drops the REF allele
			})
	}

g$qsubRscriptPbs <- function(Rscript = "", Rversion = "3.1.2", mem = 2, walltime = 24, run = TRUE){
	# Creates a *.pbs file to run the full Rscript command "Rscript"
	# The Rscript command executes a particular *.R file and provides any needed arguments
	# This version also downloads the .R file from github so that the latest version is available
	# This version does not include pbsfile as an argument, it will make from the .R root
	# This version automatically submits the pbs file to the queue
	# Example: (note that all the arguments to .R script are included in a single quotation
	# 	g$qsubRscriptPbs(Rscript = "countInvariantsByGroup.R BenlimPax22pacMar7 chrXXI 1 paxl paxb marine-pac")
	
	if(Rscript == "") stop("You need to provide an Rscript command")
		
	# Remove and repaste initial "Rscript" command if present in case it is missing from the Rscript text submitted
	# Also grab the name of the .R file argument, store in dotRfile
	Rscript <- sub("^[ ]*Rscript[ ]*|[ ]*", "", Rscript, ignore.case = TRUE)
	dotRfile <- strsplit(Rscript, split = " ")[[1]][1]
	dotRfile <- gsub("[.]R", "", dotRfile, ignore.case = TRUE)
	Rscript <- paste("Rscript", Rscript)
	
	# Attach date and time to name of file to make unique
	# pbsfile <- gsub(".pbs$", "", pbsfile)
	hour <- gsub("[ :]","-",Sys.time())
	pbsfile <- paste(dotRfile, "-", hour, ".pbs", sep = "")
	
	outfile <- file(pbsfile, "w")
	writeLines("#!/bin/bash", outfile)
	writeLines("#PBS -S /bin/bash", outfile)
	writeLines(paste("#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste("#PBS -l mem=", mem, "gb", sep = ""), outfile)
	#writeLines("#PBS -l epilogue=epilogue.script", outfile)

	#writeLines(paste("#PBS -N", jobname), outfile)
	writeLines("\n# Generic torque file to run Rscript", outfile)
	writeLines("# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	writeLines("# Run as \"qsub pbsfile.pbs\"\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	writeLines(paste("curl -o " , dotRfile, ".R", 
		" https://raw.githubusercontent.com/dschluter/RScripts/master/", dotRfile, ".R", sep = ""), outfile)

	writeLines(paste("module load R/", Rversion, sep = ""), outfile)
	writeLines(Rscript, outfile)

	writeLines("\necho \"Job finished with exit code $? at: \`date\`\"", outfile)
	
	close(outfile)
	if(run) system(paste("qsub", pbsfile))
	}

g$makeChrvecMasked <-  function(chrname, windowsmaskername){
	chr <- scan(paste(chrname, ".fa", sep = ""), what=character())  # get genome
	chr <- chr[-1] # drops the line with ">chr??"
	chr <- paste(chr, collapse="") # combined all the lines into a single word
	chrvec <- strsplit(chr, split="")[[1]] # break apart genome into individual bases
	
	#table(chrvec) # lower case refers to masked bases from the initial assembly (? maybe not in gasAcu1.fa?)
	#      A       C       G       N       T 
	#3133419 2580881 2578205  275208 3149774
	
	# Using repeatmasker output Felicity gave me on Oct 31, 2012
	chrmaskinfo <- read.table(windowsmaskername, stringsAsFactors = FALSE) 
	names(chrmaskinfo) <- c("chrno","from","to")
	
	#head(chrmaskinfo)
	   # chrno from   to
	# 1 chrXXI    0  423
	# 2 chrXXI  435 3277
	# 3 chrXXI 3285 6311
	# 4 chrXXI 6319 7888
	# 5 chrXXI 7976 8026
	# 6 chrXXI 8082 8874
	
	# ---
	# Get the positions of all the masked bases
	# NOTE: windowmasker start position is 0-based, so add 1 (end position is 1-based!)
	# Note: chrmaskinfo[ ,-1] used because otherwise each vector x is converted to all character and "seq" generates error
	
	z <- apply(chrmaskinfo[ ,-1], 1, function(x){seq(from=x[1] + 1, to=x[2])}) 
	m <- unlist(z) 
	
	chrvec[m] <- "M" # convert all these bases to "M" to indicate masked regions 
	
	cat("Summary of masked chromosome:")
	print(table(chrvec)) # Notice that "N" is not present anymore, so these were also masked in windowsmaskerSdust.chrXXI.txt
	# chrvec
	  # A       C       G       M       T 
	# 2464777 2130503 2128953 2521828 2471426 
	
	# Save to .rdd file in vcfdir
	chrno <- gsub("^chr", "", chrname)
	chrmaskfile <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd
	cat("\nSaving masked chrvec object to", chrmaskfile, "\n")
	save(chrvec, file = chrmaskfile)
	#return(chrvec)
	}        
 
g$qsubGenotypeGVCFsPbs <- function(gvcffiles, outvcfname, chromosome=NULL, GATK = "3.4.0", 
	mem = 4, walltime = 24, maxAltAlleles = 3, run = TRUE){
	# Generates qsub file "genotypeGVCFs.pbs" to carry out GATK genotypeGVCFs to call snps
	#	from multiple gvcf files inputted, one chromosome at a time
	# If chromosome is specified ("chrXXI" or "XXI") then only gvcf files having that chr are permitted.
	# If chromosome is not specified, it is extracted from gvcf file names.
	# Method assumes fasta file is named "chromosome.fa"
	# mem is in gb
	# walltime is in hours
	# Based on "makeGenotypeGVCFspbs.Rscript" but allows more flexibility in which files to analyze.
	
	cat("\ngvcffiles inputted:\n")
	for(i in 1:length(gvcffiles)){
		cat(gvcffiles[i],"\n")
		}

	# Get the chromosome number
	if(!is.null(chromosome)){
		chromosome <- gsub("chr","", chromosome)
		chromosome <- paste("chr", chromosome, sep = "")
		testfiles <- gvcffiles[grep(chromosome, gvcffiles)]
		if(length(testfiles) != length(gvcffiles)) stop("some filenames don't agree with chromosome name")
		}
	else{
		# Extract chromosome names
		z <- regexpr(text = gvcffiles, pattern = "[.]chr[A-z1]+[.]") # the "1" is for pitx1
		if( any(z < 1) ) stop("chromosome name missing from some filenames")
		chromosome <- substr(gvcffiles, z + 1, z + attr(z,"match.length") - 2)
		chromosome <- unique( chromosome )	
		# Check that there's only one chromosome included
		if( length(chromosome) != 1 ){
			stop("chr missing from some filenames, or more than one chromosome represented")
			}
		}
	pbsfilename <- paste("genotypeGVCFs", chromosome, "pbs", sep = ".")
	pbsfile <- file(pbsfilename, "w")
	cat("\ninstructions being written to", pbsfilename, "\n")
	
	writeLines("#!/bin/bash", pbsfile)
	writeLines("#PBS -S /bin/bash", pbsfile)
#	writeLines("#PBS -l walltime=24:00:00", pbsfile)
	writeLines(paste("#PBS -l walltime=", walltime, ":00:00", sep = ""), pbsfile)
#	writeLines("#PBS -l mem=4gb", pbsfile)
	writeLines(paste("#PBS -l mem=", mem, "gb", sep = ""), pbsfile)
	writeLines("#PBS -l epilogue=epilogue.script", pbsfile)
	
	jobname <- gsub("[.]vcf", "", outvcfname)
	writeLines(paste("#PBS -N", jobname), pbsfile)
	writeLines("\n# THIS FILE IS GENERATED BY AN RSCRIPT, DO NOT EDIT\n", pbsfile)
	writeLines("# This code uses multiple vgcf files from HaplotypeCaller to calls snps", pbsfile)
	writeLines("# Run as \"qsub genotypeGVCFs.pbs\"\n", pbsfile)

	if( !all(grepl("[.]vcf$", gvcffiles)) ) stop("Provide only .vcf files as arguments")

	# We assume that the fasta file for the reference chromosome is named chr[A-z1].fa	
	fastafile <- paste(chromosome, ".fa", sep = "")
	if( !file.exists(fastafile) ) stop("'chr.fa' fasta file missing")
	
#	writeLines("module load gatk/3.4.0", pbsfile)
	writeLines(paste("module load gatk/", GATK, sep = ""), pbsfile)
#	writeLines(paste("# Using GATK v ", GATK, sep = ""), pbsfile)
	cat("\nUsing GATK v", GATK, "\n")
	
	writeLines("# the -jar GenomeAnalysisTK.jar argument is included in the java call in gatk.sh so superfluous here", pbsfile)

	# Print the main java command to the pbs file (\\ is to get single backslash)
#	writeLines(paste("gatk.sh -Xmx4g -R", fastafile, "-T GenotypeGVCFs \\"), pbsfile)
	writeLines(paste("gatk.sh -Xmx", mem, "g -R ", fastafile, " -T GenotypeGVCFs \\", sep = ""), pbsfile)
	vcfarguments <- gvcffiles
	for(i in 1:length(gvcffiles)){
		vcfarguments[i] <- paste("     --variant", vcfarguments[i], "\\")
		writeLines(vcfarguments[i], pbsfile)
		}
#	writeLines(paste("     --includeNonVariantSites --max_alternate_alleles 3 -o", outvcfname), pbsfile)
#	cat("Setting max_alternate_alleles to 3\n")
	writeLines(paste("     --includeNonVariantSites --max_alternate_alleles", maxAltAlleles, "-o", outvcfname), 
		pbsfile)
	cat("Uses --max_alternate_alleles", maxAltAlleles, "\n")
	
	close(pbsfile)
	if(run) system(paste("qsub", pbsfilename))
	}
 
g$downloadSra <- function(SraSample, fishName){ # eg g$downloadSra("SRS639967", "KODK")
	# R commands for Westgrid to download SRA files from NCBI
	# based on "Using the SRAdb Package to Query the Sequence Read Archive", Zhu & Davis

	# SraSample <- "SRS639967"; fishName <- "KODK"
	# source("http://bioconductor.org/biocLite.R"); biocLite("SRAdb")
	library(SRAdb)
	sqlfile <- 'SRAmetadb.sqlite'
	if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
	
	z <- list.files(pattern = "^SRR.")
	if(length(z) > 0) stop("Method assumes there are no files beginning with 'SRR' in working directory")
	
	# Open a connection
	sra_con <- dbConnect(SQLite(), sqlfile) # opens a connection
	
	# Obtain all the experiment accessions, run accessions, etc, associated with the sample identifier
	sra_info <- sraConvert( SraSample, sra_con = sra_con )
	#sra_info
	      # sample submission     study experiment        run
	# 1  SRS639967  SRA164437 SRP042012  SRX621222 SRR1449383
	# 2  SRS639967  SRA164437 SRP042012  SRX608534 SRR1437855
	# 3  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462291
	# 4  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462120
	# 5  SRS639967  SRA164437 SRP042012  SRX621222 SRR1449714
	# 6  SRS639967  SRA164437 SRP042012  SRX608534 SRR1449659
	# 7  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462754
	# 8  SRS639967  SRA164437 SRP042012  SRX621222 SRR1449492
	# 9  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462185
	# 10 SRS639967  SRA164437 SRP042012  SRX608534 SRR1422274
	# 11 SRS639967  SRA164437 SRP042012  SRX608534 SRR1462825
	
	sra_list <- apply(sra_info, 2, unique)
	# $sample
	# [1] "SRS639967"
	# $submission
	# [1] "SRA164437"
	# $study
	# [1] "SRP042012"
	# $experiment
	# [1] "SRX621222" "SRX608534"
	# $run
	 # [1] "SRR1449383" "SRR1437855" "SRR1462291" "SRR1462120" "SRR1449714"
	 # [6] "SRR1449659" "SRR1462754" "SRR1449492" "SRR1462185" "SRR1422274"
	# [11] "SRR1462825"
	
	# Lists all the ftp addresses of the sra files
	# listSRAfile( SraSample, sra_con, fileType = 'sra' )
	
	# Gets the ftp addresses of the fastq files corresponding to the run
	# getFASTQinfo( sra_list$run, srcType = 'ftp' )
	# Downloads them!
	getSRAfile( sra_list$run, sra_con, fileType = 'fastq' )
	
	# Rename the files to match population for easy recognition using file.rename(from, to)
	z <- list.files(pattern = "^SRR.")
	for(i in 1:length(z)){
		file.rename(z[i], paste(fishName, z[i], sep="_"))
		}
	cat("Files downloaded and renamed:\n")
	list.files(pattern=paste("^", fishName,".", sep = ""))
	}

g$recode<-function (..., ret = c("numeric", "factor"), none = if (ret == "numeric") 0 else "none", na){
	# stolen from the Hmisc library
    ret <- match.arg(ret)
    w <- list(...)
    if (!is.logical(w[[1]]) && length(w) == 3) {
        z <- w[[3]][match(w[[1]], w[[2]])]
        if (!missing(none)) 
            z[if (is.numeric(none)) 
                is.na(z)
            else z == ""] <- none
        return(z)
    }
    nam <- names(w)
    if (missing(ret)) 
        ret <- if (all.is.numeric(nam)) 
            "numeric"
        else "factor"
    result <- rep(none, length(w[[1]]))
    for (i in 1:length(w)) result[w[[i]]] <- if (ret == "numeric") 
        numnam[i]
    else nam[i]
    if (ret == "factor") 
        result <- as.factor(result)
    if (!missing(na)) 
        result[is.na(na)] <- NA
    result
}

g$geno.diff <- function(GT1, GT2){
	# calculates the divergence between two genotype vectors of equal length (e.g., columns of GT matrix of vcf file)
	# Genotypes of individuals are in VCF "GT" format:
		# NA   "0/0"   "1/1"      NA   "0/1"   "0/1"   "0/0"      NA   "1/1"   "1/1"   "0/0"   "0/0"
	#
	# Tests that this is done correctly:
	# g$geno.diff("0/0", "0/1") # [1] 0.5
	# g$geno.diff("0/0", "1/1") # [1] 1
	# g$geno.diff("0/0", "0/0") # [1] 0
	# g$geno.diff("0/0", "0/1") # [1] 0.5
	# g$geno.diff("0/1", "0/0") # [1] 0.5
	# g$geno.diff("0/0", "1/1") # [1] 1
	# g$geno.diff("0/1", "2/3") # [1] 1
	# g$geno.diff("0/1", "1/2") # [1] 0.5
	# g$geno.diff(c("0/1","0/0"), c(NA,"0/0")) 	# [1] NA  0
	#
	# GT1 <- c(NA, "0/0", "0/0", "0/0", "0/1", "0/1", "0/1", "1/1", "1/1", "0/1", "0/1", "0/1")
	# GT2 <- c(NA, "0/0", "0/1", "1/1", "0/0", "0/1",    NA, "0/0", "0/1", "1/2", "2/3", "1/4")
	x <- strsplit(GT1, split = "/") 
	y <- strsplit(GT2, split = "/") 
	diff <- mapply(x, y, FUN = function(xi, yi){
		# xi <- x[[2]]; yi <- y[[2]]
		diff <- ceiling((4-sum(is.element(xi,yi), is.element(yi,xi)))/2)/2
		if(any(is.na(c(xi,yi)))) diff <- NA
		diff
		})
	diff
}

g$psd <- function(genotypes, groups){
	# Calculate the number of genotypic differences between all pairs of individuals, one pair at a time.
	# Pairwise sequence difference "psd" is 0, 0.5, and 1 according to alleles shared between individuals at the locus.
	# "genotype" is a matrix of genotypes of individuals in VCF "GT" format:
		# NA   "0/0"   "1/1"      NA   "0/1"   "0/1"   "0/0"      NA   "1/1"   "1/1"   "0/0"   "0/0"
	# "groups" indicates corresponding group membership of individuals.	
	# genotypes <- rbind( c(NA, "0/0", "0/0", "0/0", "0/1", "0/1", "0/1", "1/1", "1/1", "0/1", "0/1", "0/1"),
	#                     c(NA, "0/0", "0/1", "1/1", "0/0", "0/1",    NA, "0/0", "0/1", "1/2", "2/3", "1/4") )
	# groups <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)
	# genotypes <- x1$GT

	pop <- as.character(groups)  # population names

	gtpairs <- combn(pop, 2, FUN=function(x){paste(x, collapse=",")})
	# gtpairs
	 # [1] "1,1" "1,1" "1,1" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,1" "1,1" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2"
	# [20] "1,2" "1,2" "1,1" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2"
	# [39] "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2"
	# [58] "2,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2"
	
	colpairs <- combn(1:ncol(genotypes), 2)
	colpairs
	     # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
	# [1,]    1    1    1    1    1    1    1    1    1     1     1     2     2     2     2     2     2
	# [2,]    2    3    4    5    6    7    8    9   10    11    12     3     4     5     6     7     8
	     # [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
	# [1,]     2     2     2     2     3     3     3     3     3     3     3     3     3     4     4     4
	# [2,]     9    10    11    12     4     5     6     7     8     9    10    11    12     5     6     7
	     # [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49]
	# ...
	
	psd <- apply(colpairs, 2, function(z){
			psd <- g$geno.diff(genotypes[ ,z[1] ], genotypes[ ,z[2] ])
			})
	
	psd <- as.data.frame(psd, stringsAsCharacters = FALSE)
	names(psd) <- gtpairs
	psd
	}


g$blockstats <- function(vcfresults, stepsize = 500, invariants, 
		chrvec, CSStrueSnps = TRUE, FSTtrueSnps = FALSE, psdMissingAction = "none"){
	
	# This function calculates preliminary results for sliding window analyses on a (repeat-masked) chromosome
	# It breaks the chromosome into blocks or bins of size stepsize and calculates a value of interest in each bin

	# vcfresults is the VCF list containing all variables of interest for a single chromosome
	# Function assumes that vcfresults$fst and/or vcfresults$psd exist
	# Function assumes that the stats corresponding to masked nucleotides have already been dropped from vcfresults

	# yes it can handle multiple snps at the same value of POS

	if( !is.element(psdMissingAction, c("meanAll", "meanBW", "meanGroup")) )
		stop("psdMissingAction must be one of meanAll, meanBW, or meanGroup")

	# psdMissingAction is for g$blockstats, how to average when psd values are missing.
	# Must be one of the following if task include "css":
	# 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
	# 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
	#   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
	#		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
	#		"1,1" and "2,2" as the same (i.e., both are within-group).
	# On output, the "np" variables indicate the number of non-missing single-nucleotide psd values available
	# 	for that pair of fish to calculate the mean psd for that pair of fish across the 500 (eg) nucleotides 
	# 	in the block (in this case 500 minus "np" indicate the number of nucleotides in the block that were missing
	# 	for that pair and had been set to an average value (accurding to the rule specified by psdMissingAction)

	# chrvec is a vector of the whole MASKED chromosome, all nucleotides as distinct elements, obtained by:
	#       Masked bases are indicates with an "M"

	# invariants contains numbers indicating chromosome positions meeting the minimum-number-of genotypes criterion
	#		Probably in file "processedInvariants"

	# If FSTtrueSnps = TRUE, only true snps (not indels) are used here to calculate variance components.
	# If CSStrueSnps = TRUE, only true snps (not indels) are used here to calculate psd within blocks.

	# table(chrvec)
	      # A       C       G       M       T 
	# 2464777 2130503 2128953 2521828 2471426 
	
	x1 <- vcfresults
	stepsize <- as.integer(stepsize)	
	k <- stepsize # this is the size of the block
	nbases <- length(chrvec)

	# Establish the break points of the bins into which nucleotides will be grouped (e.g., k = 500 bases per bin)
	# The last bin goes from the final bin of fully k nucleotides to the last nucleotide. 
	# This last bin might be small.
	
	ibase <- c( seq(1, nbases, k), nbases + 1) # bases marking breaks between steps of size k
	
	# If want something to indicate midpoint of blocks instead:
	midbase <- ibase + (stepsize)/2
	#[1]  251  751 1251 1751 2251 2751

	# head(ibase)
	# [1]    1  501 1001 1501 2001 2501
	# tail(ibase)
	# [1] 11715001 11715501 11716001 11716501 11717001 11717488 * last one less than k

  # 1) Count number of non-M bases in each nucleotide bin of stepsize "k". 
	
	# Break nucleotide index of chr into bins of stepsize k

	chrbins <- findInterval(1:nbases, ibase) # indicates in which "bin" of size k each base belongs
											 # if none missing, there should be k 1's, k 2's, etc.
	# chrbins[499:505]
	# [1] 1 1 2 2 2 2 2
	
	# count up the number of unmasked nucleotides in each bin - These are the ones coded as ACGT
	nUnmasked <- tapply(chrvec, chrbins, function(x){length(x[x %in% c("A","C","G","T")])})

	# head(nUnmasked)
	#  1  2  3  4  5  6 
	# 12  0  0  0  0  0 

  # 2) Count up the number of snp and number of invariants in each bin
	# Need to split SNP counts into same bins of stepsize k and then sum up

	# Warning: some values for POS (or start(x1$rowdata)) are listed more than once
	# head(start(x1$rowdata)) # start positions of snp
	# [1] 18389 18389 54855 54866 54879 54979
		
	# Count up the number of good snp in each bin (cut works because ibase intervals are made as factors)
	library(GenomicRanges)
	library(Biostrings)
	snpBins <- cut(start(x1$rowdata), breaks = ibase, right = FALSE)
	nSnp <- table( snpBins )

	# head(nSnp)
    # [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001) 
    #       0           0           0           0           0           0 

	# Count up the number of good invariants in each bin
	# invariants is a large vector of numbers indicating bases that are ACGT not M
	#old nInvariants <- table( cut(invariants, breaks = ibase, right = FALSE) )
	nInvariants <- table( cut(invariants$POS, breaks = ibase, right = FALSE) )
	
	# head(nInvariants)
    # [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001)  
    #       0           0           0           0           0           0 

	# all the following items have the same length (23435)

	results <- data.frame(ibase = ibase[-length(ibase)], midbase = midbase[-length(midbase)], nUnmasked = nUnmasked, 
					nSnp = as.vector(nSnp), nInvariants = as.vector(nInvariants))

	#results[105:125,]
	    # ibase midbase nUnmasked nSnp nInvariants
	# 105 52001   52251         0    0           0
	# 106 52501   52751        25    0           0
	# 107 53001   53251       159    0           0
	# 108 53501   53751       372    0           0
	# 109 54001   54251       302    0           0
	# 110 54501   54751       480    4         128
	# 111 55001   55251       399    2          90
	# 112 55501   55751       500    0           0
	# 113 56001   56251       376    0           0
	# 114 56501   56751       500    0           0
	# 115 57001   57251       409    0           0
	# 116 57501   57751       257    0           0
	# 117 58001   58251       108    0           0
	# 118 58501   58751       440    1         157
	# 119 59001   59251       487    4         189
	# 120 59501   59751       443    0           0
	# 121 60001   60251       162    0           0
	# 122 60501   60751         0    0           0
	# 123 61001   61251       108    0           0
	# 124 61501   61751       203    0           0
	# 125 62001   62251       474    0           0
	
	
  # 3) Calculate Fst summary stats
	if( "fst" %in% names(vcfresults) ){
		
		fst <- as.data.frame(x1$fst, stringsAsFactors = FALSE)
		
		if(FSTtrueSnps){
			FstBins <- split(fst[x1$is.trueSnp, ], snpBins[x1$is.trueSnp])
			#
			# I haven't tested the following yet:
			results$nTrueSnp <- sapply(FstBins, nrow) 
			}
		else 
			FstBins <- split(fst, snpBins) # this works because snpBins is a factor created above
										   # by snpBins <- cut(start(x1$rowdata), breaks = ibase, right = FALSE)

		z <- lapply(FstBins, function(z){
				VARa <- sum(z$lsiga) # among populations
				VARb <- sum(z$lsigb) # among individuals within populations
				VARw <- sum(z$lsigw) # within individuals
				FST <- VARa/(VARa + VARb + VARw)
				return(c(VARa = VARa, VARb = VARb, VARw = VARw, fst = FST))
				})
		# z[1]
		# $`[1,501)`
		# VARa VARb VARw  fst 
		   # 0    0    0  NaN 
		   
		# z[110]
		# $`[54501,55001)`
		      # VARa       VARb       VARw        fst 
		# -0.3241358  1.3303872  0.7373737 -0.1858976

		z <- do.call("rbind", z)
		
		# Put results into the results data frame. 
		results <- cbind.data.frame(results, z)
		results$fst[is.nan(results$fst) & results$nSnp + results$nInvariants > 0] <- 0
		
		# head(results)
		  # ibase nUnmasked nSnp nInvariants VARa VARb VARw fst
		# 1     1        12    0           0    0    0    0 NaN
		# 2   501         0    0           0    0    0    0 NaN
		# 3  1001         0    0           0    0    0    0 NaN
		# 4  1501         0    0           0    0    0    0 NaN
		# 5  2001         0    0           0    0    0    0 NaN
		# 6  2501         0    0           0    0    0    0 NaN
		
		# results[105:115,]
		    # ibase midbase nUnmasked nSnp nInvariants       VARa     VARb      VARw        fst
		# 105 52001   52251         0    0           0  0.0000000 0.000000 0.0000000        NaN
		# 106 52501   52751        25    0           0  0.0000000 0.000000 0.0000000        NaN
		# 107 53001   53251       159    0           0  0.0000000 0.000000 0.0000000        NaN
		# 108 53501   53751       372    0           0  0.0000000 0.000000 0.0000000        NaN
		# 109 54001   54251       302    0           0  0.0000000 0.000000 0.0000000        NaN
		# 110 54501   54751       480    4         128 -0.3241358 1.330387 0.7373737 -0.1858976
		# 111 55001   55251       399    2          90 -0.0800000 0.300000 0.3250000 -0.1467890
		# 112 55501   55751       500    0           0  0.0000000 0.000000 0.0000000        NaN
		# 113 56001   56251       376    0           0  0.0000000 0.000000 0.0000000        NaN
		# 114 56501   56751       500    0           0  0.0000000 0.000000 0.0000000        NaN
		# 115 57001   57251       409    0           0  0.0000000 0.000000 0.0000000        NaN
		
		}

  			
  # 4) Calculate CSS quantities of interest within bins
  
  	# Note: include only true snps if CSStrueSnps = TRUE
  	
  	if( "psd" %in% names(vcfresults) ){
  		
  		# use "as.data.frame" not "data.frame", which puts an "X" at start of each colname
		psd <- as.data.frame(x1$psd, stringsAsFactors = FALSE)
		
		gtpairs <- names(psd)
		# gtpairs <- combn(pop, 2, FUN=function(x){paste(x, collapse=",")})

		# Eliminate missing pairwise differences by replacing with averages
			
		# 1. Save record of which pairwise psd values are NOT missing to array 'psdNonMissing'
		psdNonMissing <- as.data.frame(!is.na(psd))

		# 2. Figure out categories for "psdMissingAction" behavior
		wbPairs <- rep("w", length(gtpairs))
		wbPairs[sapply(strsplit(gtpairs, split=","), function(x){length(unique(x))}) == 2] <- "b"
		# wbPairs
		 # [1] "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w"
		# [25] "b" "b" "b" "b" "b" "b" "w" "w" "b" "b" "b" "b" "b" "b" "w" "b" "b" "b" "b" "b" "b" "b" "b" "b"
		# [49] "b" "b" "b" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w"
	
		if(psdMissingAction == "meanBW") psdGroups <- wbPairs else 
			if(psdMissingAction == "meanAll") psdGroups <- rep(1, length(wbPairs)) else 
				if(psdMissingAction == "meanGroup") psdGroups <- gtpairs else 
					stop("Unrecognizeable value for 'psdMissingAction'")
		
		# 3. Assign average pairwise distance to missing pairwise distance values at every marker separately
			# Function "ave" modified to allow missing values
		ave.rm <- function(x, ..., FUN = mean){
		    if(missing(...)) 
		        x[] <- FUN(x)
		    else{
		        g <- interaction(...)
		        split(x, g) <- lapply(split(x, g), FUN, na.rm = TRUE)
		    	}
		    x
			}

		z1 <- apply(psd, 1, function(x){
			z1 <- ave.rm(x, psdGroups) 
			x[is.na(x)] <- z1[is.na(x)]     # assign averages to missing pairs
			x
			})
		psd <- as.data.frame(t(z1))
		names(psd) <- gtpairs	
				
				
		if(CSStrueSnps){
			psdBins <- split(psd[x1$is.trueSnp, ], snpBins[x1$is.trueSnp])
			psdNonMissing <- split(psdNonMissing[x1$is.trueSnp, ], snpBins[x1$is.trueSnp])
			if(is.null(results$nTrueSnp)) results$nTrueSnp <- sapply(psdBins, nrow)
			}
		else{
			# genotypeBins <- split(GT, snpBins)
			psdBins <- split(psd, snpBins)
			psdNonMissing <- split(psdNonMissing, snpBins)
			}
		
		# psdNonMissing[1]
		# $`[1,501)`
		 # [1] 1,1 1,1 1,1 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,1
		# [25] 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2
		# [49] 1,2 1,2 1,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2
		# <0 rows> (or 0-length row.names)
		
		# psdNonMissing[110]
		# $`[54501,55001)`
		               # 1,1   1,1   1,1   1,1   1,1   1,2   1,2   1,2   1,2   1,2   1,2  1,1   1,1  1,1   1,1
		# chrXXI:54855 FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE TRUE FALSE TRUE  TRUE
		# chrXXI:54866 FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE TRUE  TRUE TRUE  TRUE
		# chrXXI:54879  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE TRUE  TRUE TRUE  TRUE
		# chrXXI:54979  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE TRUE  TRUE TRUE FALSE
		# chrXXI:54981  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE TRUE  TRUE TRUE FALSE
		# ...etc
		
		#colSums(psdNonMissing[[1]])
		# 1,1 1,1 1,1 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,1 1,2 
		  # 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
		# 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 
		  # 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
		# 1,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 
		  # 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
		
  		#colSums(psdNonMissing[[110]])
				# 1,1 1,1 1,1 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,1 1,2 
		  # 3   3   3   3   1   3   2   3   1   3   3   5   4   5   3   5   2   5   3   5   5   4   5   3   5 
		# 1,2 1,2 1,2 1,2 1,2 1,1 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,1 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 1,2 
		  # 2   5   3   5   5   4   2   4   2   4   2   4   4   3   5   2   5   3   5   5   3   0   3   3   3 
		# 1,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 2,2 
		  # 3   2   5   3   5   5   2   0   2   2   3   5   5   3   3   5 

		# Count up the number of non-missing psd values for each pair of fish in each bin
	
		psdCountByPair <- lapply(psdNonMissing, colSums)
		psdCountByPair <- as.data.frame( do.call("rbind", psdCountByPair) )
		names(psdCountByPair) <- paste("np", names(psdCountByPair), sep = ".")

		# Sum the psd's within bins or blocks

		psdSum <- lapply(psdBins, colSums)
		psdSum <- as.data.frame( do.call("rbind", psdSum), stringsAsFactors = FALSE )

		# psdSum[109:112,]
		                 # 1,1    1,1    1,1    1,1    1,1  1,2      1,2  1,2  1,2  1,2  1,2 1,1    1,1    1,1
		# [54001,54501) 0.0000 0.0000 0.0000 0.0000 0.0000 0.00 0.000000 0.00 0.00 0.00 0.00   0 0.0000 0.0000
		# [54501,55001) 1.6125 1.6125 1.6125 2.1125 1.7625 1.96 2.176667 1.96 1.56 1.46 1.46   3 2.5625 2.5000
		# [55001,55501) 1.0000 0.0000 0.0000 1.1875 0.5000 1.00 0.700000 0.25 0.00 0.50 1.00   1 1.0000 0.1875
		# [55501,56001) 0.0000 0.0000 0.0000 0.0000 0.0000 0.00 0.000000 0.00 0.00 0.00 0.00   0 0.0000 0.0000
		               # 1,1 1,2      1,2  1,2 1,2 1,2 1,2    1,1    1,1  1,1 1,2      1,2  1,2 1,2 1,2 1,2
		# [54001,54501) 0.00 0.0 0.000000 0.00 0.0 0.0 0.0 0.0000 0.0000 0.00 0.0 0.000000 0.00 0.0 0.0 0.0
		# [54501,55001) 1.65 1.5 2.676667 3.50 2.6 1.0 0.5 0.5625 1.5000 1.65 2.5 1.676667 0.50 0.6 3.0 2.5
		# [55001,55501) 0.50 2.0 0.700000 1.25 1.0 0.5 0.0 0.0000 1.1875 0.50 1.0 0.700000 0.25 0.0 0.5 1.0
		# [55501,56001) 0.00 0.0 0.000000 0.00 0.0 0.0 0.0 0.0000 0.0000 0.00 0.0 0.000000 0.00 0.0 0.0 0.0
		                 # 1,1    1,1 1,2      1,2  1,2 1,2 1,2 1,2    1,1  1,2      1,2  1,2  1,2  1,2  1,2 1,2
		# [54001,54501) 0.0000 0.0000   0 0.000000 0.00 0.0 0.0   0 0.0000 0.00 0.000000 0.00 0.00 0.00 0.00 0.0
		# [54501,55001) 1.5625 1.7125   2 1.676667 1.00 1.1 2.5   2 1.1500 1.00 1.676667 2.00 2.10 1.50 2.00 2.1
		# [55001,55501) 1.1875 0.5000   1 0.700000 0.25 0.0 0.5   1 0.6875 1.25 0.700000 1.25 1.25 0.75 0.25 1.5
		# [55501,56001) 0.0000 0.0000   0 0.000000 0.00 0.0 0.0   0 0.0000 0.00 0.000000 0.00 0.00 0.00 0.00 0.0
		                   # 1,2  1,2 1,2 1,2 1,2    2,2    2,2  2,2 2,2 2,2    2,2    2,2    2,2    2,2    2,2
		# [54001,54501) 0.000000 0.00 0.0 0.0 0.0 0.0000 0.0000 0.00 0.0   0 0.0000 0.0000 0.0000 0.0000 0.0000
		# [54501,55001) 1.776667 1.60 1.6 2.1 1.1 1.8325 3.0000 3.15 0.5   2 1.3325 1.9825 2.3325 2.8325 0.6500
		# [55001,55501) 0.700000 0.75 0.5 0.0 0.5 0.7275 0.1875 1.00 1.5   2 0.7275 0.7275 0.7275 0.7275 0.1875
		# [55501,56001) 0.000000 0.00 0.0 0.0 0.0 0.0000 0.0000 0.00 0.0   0 0.0000 0.0000 0.0000 0.0000 0.0000
		                 # 2,2    2,2  2,2  2,2 2,2
		# [54001,54501) 0.0000 0.0000 0.00 0.00 0.0
		# [54501,55001) 3.5000 3.0000 3.15 2.15 1.5
		# [55001,55501) 0.6875 1.1875 0.50 1.00 0.5
		# [55501,56001) 0.0000 0.0000 0.00 0.00 0.0

		names(psdSum) <- paste("psd", names(psdSum), sep = ".")

		# nTrueSnp has already been added to "results"
		results <- cbind.data.frame(results, psdSum, psdCountByPair, stringsAsFactors = FALSE)
		results
		}
}


g$slidewin <- function(blockstats, method = "FST", stepsize, nsteps.per.window, windowNmin){

	# Combines results of g$blockstats across nsteps.per.window of bins at a time in sliding window.
	# Averages and other stats are computed based only on the number of bases in the window regardless
	#		of whether CSStrueSnps and FSTtrueSnps were TRUE or FALSE.
	# "ibase" in blockstats is the centre of the sliding window of 5 blocks (e.g.).
	# rollapply adds padding to the left of the data range near the start of the data set and
	# 	to the right of the data range near the end; 
	#	padding is determined by "fill" option (NA is the default).
	# In all cases, windowNmin is the smallest number of good bases needed in each full sliding window to yield result, 
	# 	otherwise window's value is set to NA. 

	# FST: returns the Weir-Cockeram VARiances obtained by summing variance component across all bases in the window
	#		and then calculating Fst as the usual ratio of the variance components. Monomorphic
	#		bases in a window do not influence this Fst, because they have no variance. 
	#		However, the sliding windows for VARiance components include the 0's where there are data.
	#
	# CSS: returns Felicity's CSS score, following the formula in Jones et al 2012 Supplement Eq. 1.
	# 		I assume that she used metric multidimensional scaling on the pairwise sequences divergences between
	#			pairs of fish within each sliding window. Calculation of CSS uses Euclidean distances between 
	#			pairs of fish on just the first 2 axes from this MDS analysis.
	#		The formula for CSS given in the Jones et al Suppl is
	#				CSS <- ( sij/(m*n) ) - (m+n)*( sii/((m-1)*(m^2)) + sjj/((n-1)*(n^2)) )
	# 		but the correct formula is
	#				CSS <- ( sij/(m*n) ) - (1/(m+n))*( sii/((m-1)/2) + sjj/((n-1)/2) ) 
	# 		where sij/(m*n) is the mean of all pairwise distances between lim and ben
	# 			(1/(m+n))*( sii/((m-1)/2) ) here works out to be the mean pairwise distance within ben
	# 			(1/(m+n))*( sjj/((n-1)/2) ) here works out to be the mean pairwise distance within lim
	#					
	#		To handle missing values, Felicity et al calculated 21x21 pairwise distances
	# 		between the 21 populations of fish in each sliding window. Any pairwise distance based on fewer
	#		that 100 good bases (out of 2500 bases in the window, ie 4 %) was tossed and replaced with the
	#		matrix average of the remaining pairwise comparisons. If more than 50% of the pairwise distances
	#		had to be tossed, the whole result for the window was set to NA.		
	# 		My value for the minimum number of validated positions is blocksize * nblocks.per.window * windowFracMin
	# 		Note: I have already replaced NA's with the mean psd (separately for each locus.
	# 		So there should be no missing psd's within an informative locus by this point

	#		I handled missing values differently from Felicity et al. 
	#		1. For a given pair of fish I calculated psd (perc seq div) at every nucleotide (snp and non-snp).
	#		2. If a psd value was NA at a nucleotide, I used an average psd from the remaining pairs of fish
	#			(according to the rule specified by psdMissingAction)
	#		3. For the given pair of fish, I took the average psd across all nucleotides in the block of 500 (eg),
	#			(averaging over both the real and substituted "missing" psd values. 
	# 		4. Good bases already use a criterion for a minimum number of genotypes per group, for variant and invariant
	#			sites. This ensures that most pairwise divergence values are present for every good base.
	#		5. I also employ a minimum number of good bases per window to generate a result, otherwise the result
	#			from a window is set to NA.
	#		4. If needed, for every pair of fish, the corresponding "np" value in blockstats remembers the number of 
	#			real (non-missing) psd's used in the calculation of the average psd. This would allow me to set a 
	#			minimum number of non-missing psd's per pair of fish to yield a result. Meantime, previous criteria
	#			should suffice.

  	# Calculate sliding window values, weighting by n and N.
	library(zoo)
	
	k <- stepsize # this is the size of the block
	wink <- nsteps.per.window # number of blocks per sliding window
	
	if( casefold(method) == "fst" ){
		FST <- rollapply(blockstats, width = wink, by.column = FALSE, fill = NA, 
				FUN = function(x){    # remember: argument "x" is a window, so has multiple rows
					#print(x[,1:9])
					VARa <- x[, "VARa"] 
					VARb <- x[, "VARb"]
					VARw <- x[, "VARw"]
					
					nbases <- sum(x[, "nSnp"] + x[, "nInvariants"], na.rm = TRUE)
					if(nbases < windowNmin){
						totVARa <- NA
						totVARb <- NA
						totVARw <- NA
						FST <- NA
						}
					else{
						totVARa <- sum(VARa, na.rm=TRUE)
						totVARb <- sum(VARb, na.rm=TRUE)
						totVARw <- sum(VARw, na.rm=TRUE)
						FST <- totVARa/(totVARa + totVARb + totVARw)
						FST[is.nan(FST)] <- NA
						}
					#print(c(nbases, totVARa, totVARb, totVARw, FST))
					c(nbases = nbases, totVARa = totVARa, totVARb = totVARb, totVARw = totVARw, fst = FST)
					})
		#print(FST)
		results <- cbind.data.frame(ibase = blockstats$ibase, FST)		
		}


	if( casefold(method) == "css" ){
				
		# indicate identity of pairs, ie "1,1", "1,2", "2,2"
		# gtpairs <- combn(pop, 2, FUN = function(x){paste(x, collapse=",")}) # 1,1 1,2 2,2 indicates groups

		is.psd <- grep("psd.", names(blockstats))
		gtpairs <- names(unlist(blockstats[1,]))[is.psd] # only way I found to stop R adding suffix to duplicate names
		gtpairs <- gsub( "psd.", "", gtpairs )
		npop <- ceiling(sqrt(2*length(gtpairs)))
		
		# This next bit assumes just two populations types
		#poptype <- unique(pop)
		poptype <- unique(unlist(strsplit(gtpairs, split = ",")))
		if(length(poptype) != 2) stop("CSS only works for two population types")
		
		ii <- gtpairs == paste(poptype[1], poptype[1], sep = ",") # cases of "1,1"
		jj <- gtpairs == paste(poptype[2], poptype[2], sep = ",") # cases of "2,2"
		ij <- gtpairs == paste(poptype[1], poptype[2], sep = ",") # cases of "1,2" (there are no "2,1")
		m <- ceiling(sqrt(2 * length(gtpairs[ii])))
		n <- ceiling(sqrt(2 * length(gtpairs[jj])))
		
		# check
		if( length(gtpairs[ij]) != m*n ) stop("Number of pairs of individuals doesn't match up")
			
		# A <- grep(poptype[1],colnames(euclid)) # e.g., benthics
		# B <- grep(poptype[2],colnames(euclid)) # e.g., limnetics
		# m <- length(A)
		# n <- length(B)

		library(gdata) # needed for lowerTriangle()
		CSS <- rollapply(blockstats, width = wink, by.column = FALSE, fill = NA, 
				FUN = function(x){
				# FUN = function(x,  # remember: argument "x" is a window, so has multiple rows
						# is.psd = is.psd,
						# gtpairs = gtpairs,
						# npop = npop,
						# poptype = poptype,
						# ii = ii,
						# jj = jj,
						# m = m,
						# n = n
						# ){   
				# x <- blockstats[106:110, ]
					
					# xpairs <- x[, gtpairs]
					xpairs <- x[ ,is.psd ]
					
					sumdist <- apply(xpairs, 2, sum, na.rm = TRUE) # sum of summed pairwise dist over blocks in window
					names(sumdist) <- gtpairs
					
					# total number of bases in window, assume that CSStrueSnps = TRUE
					nbases <- sum(x[, "nTrueSnp"] + x[, "nInvariants"], na.rm = TRUE) 
										
					if(nbases < windowNmin) CSS <- NA
					else{
						# Here I'm just calculating the sum of psd for all pairs, for now ignoring
						# snpCountByFish, ie whether some individuals may have too few genotypes

						# Put sumdist into matrix
						# first get [i,j] of the pairs for the matrix d
												
						d <- matrix(NA, nrow = npop, ncol = npop)
						# lowerTriangle(d) <- gtpairs # a test
						lowerTriangle(d) <- sumdist
						d <- t(d)
						lowerTriangle(d) <- sumdist
						diag(d) <- rep(0, nrow(d)) # the diagonal should all be zeros
						
						# Following the formula in Jones et al 2012 Supplement Eq. 1
						# Assume that she used metric multidimensional scaling
						# The following command generates a lot of warnings of the sort:
						#  "1: In cmdscale(d, k = 2, add = TRUE) : only 0 of the first 2 eigenvalues are > 0"
						# 	so I have suppresed them in case they halt execution.
						suppressWarnings(
						z <- cmdscale(d, k = 2, add = TRUE)$points  # it crashed if add=FALSE used
						)
						
						# rownames(z) <- pop
						
						# calculate pairwise distance between observations on the cmdscale axes
						euclid <- as.matrix(dist(z, method = "euclidean")) 	# matrix of distances; pop is row and col names
						euclid.lower <- euclid[lower.tri(euclid)] 			# pairwise distances, order matches gtpairs
						
						# sij <- sum(euclid[A,B])    # is the same as when use euclid[B,A] instead
						# sii <- sum(euclid[A,A][lower.tri(euclid[A,A])])
						# sjj <- sum(euclid[B,B][lower.tri(euclid[B,B])])
						
						sii <- sum(euclid.lower[ii])
						sjj <- sum(euclid.lower[jj])
						sij <- sum(euclid.lower[ij])
						
						CSS <- ( sij/(m*n) ) - (1/(m+n))*( sii/((m-1)/2) + sjj/((n-1)/2) ) 
						}
						
					return(c(nbases = nbases, CSS = CSS))
					})
		results <- cbind.data.frame(ibase = blockstats$ibase, CSS)			
		detach(2) # to avoid conflicts among libraries
		}
		
	# if("OR" %in% method){
		# # Note that the value of meanOR is already averaged over all sites within a block
		# # with nbases = N - (nsnp - n)
		# meanOR <- g$blockstats(x1, pop = substr(names(x1$GT),1,4), method = "OR", 
						# zeros = TRUE, stepsize = 500, chrvec = chrvec)
		
		# OR <- rollapply(meanOR, width=wink, by.column=FALSE, fill = NA, 
				# FUN=function(x){    # remember: argument "x" is a window, so has multiple rows
					# mean.OR = x[, 3]
					# # n = x[, 2]
					# # nsnp = x[,4] # n < nsnp if some snp did not meet criteria to calculate OR
					# # N = x[, 5] # good bases
					# # nbases =   N - (  nsnp   -   n)
					# nbases <- x[, 5] - (x[,4] - x[, 2])
					# if(sum(nbases) < windowNmin) OR <- NA
					# else OR <- ( sum(nbases * mean.OR, na.rm=TRUE) )/sum(nbases, na.rm=TRUE)
					# #return(OR)
					# })
		# #meanOR$ORwin <- OR
		# results <- cbind.data.frame(ibase = meanOR$ibase, meanOR = OR)
		# }

		
	results
}

g$plot.reads <- function(x1, from = 1, to = nrow(x1), scalefactor = 1, minqual = 30, keepunmappedmate=TRUE){
	# x1 is a data frame made from scanBam object
	# If they don't already exist, the following LOGICAL variables are created:
	# "unmapped" "unmappedmate", "bothforward","bothreverse"
	# Positive insert size: 
	# 	draws a line connecting the start of a read (top) to the end of a read (bottom)
	# 	draws only reads whose position (pos) at the START of read is between "from" and "to"
	# Negative insert size: 
	#	draws a line connecting the end of a read (bottom) backwards to the start of a read (top)
	# 	draws only reads whose position (pos) at the END of read is between "from" and "to"
	# scalefactor is the amount that the x-axis in plot is extended beyond "for" and "to", in multiples of to-from

	if(!is.data.frame(x1)) stop("Error: object must be a data frame\n")
	interval <- to - from
	if(interval < 0) stop("Error: to must be greater than from\n")
	library(bitops)
	if(!exists("unmapped")) unmapped <- bitAnd(x1$flag, 0x0004) == 0x0004
	if(!exists("unmappedmate")) unmappedmate <- bitAnd(x1$flag, 0x0008) == 0x0008
	if(!exists("bothforward")) bothforward <- 
			bitAnd(x1$flag, 0x0010) == 0x0000 & bitAnd(x1$flag, 0x0020) == 0x0000 & !unmapped & !unmappedmate
	if(!exists("bothreverse")) bothreverse <- 
			bitAnd(x1$flag, 0x0010) == 0x0010 & bitAnd(x1$flag, 0x0020) == 0x0020 & !unmapped & !unmappedmate

	# regular isizes, positive insert sizes
	if(keepunmappedmate){
		i <- !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0} 
	else{
		i <- !unmappedmate & !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0}
	dat <- x1[i & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="positive inserts",ylab="")
	segments(dat$pos, rep(1,nrow(dat)), dat$pos+71+dat$isize, rep(0,nrow(dat)), lwd=0.2)

	# regular isizes, negative insert sizes
	if(keepunmappedmate){
		j <- !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0} 
	else{
		j <- !unmappedmate & !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0}
	datj <- x1[j & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="negative inserts",ylab="")
	segments(datj$pos, rep(0,nrow(datj)), datj$pos+datj$isize, rep(1,nrow(datj)), lwd=0.2)

	# both reads point same way, positive insert sizes
	i <- bothforward & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0
	j <- bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0
	dati <- x1[i & x1$pos>from & x1$pos<to,]
	datj <- x1[j & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="positive inserts, both reads point same direction",ylab="")
	segments(dati$pos, rep(1,nrow(datj)), dati$pos+71+dati$isize, rep(0,nrow(dati)), lwd=0.5, col="red")
	segments(datj$pos, rep(1,nrow(datj)), datj$pos+71+datj$isize, rep(0,nrow(datj)), lwd=0.5, col="blue")

	# both reads point same way, negative insert sizes
	i <- bothforward & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0
	j <- bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0
	dati <- x1[i & x1$pos>from & x1$pos<to,]
	datj <- x1[j & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="negative inserts, both reads point same direction",ylab="")
	segments(dati$pos, rep(0,nrow(datj)), dati$pos+dati$isize, rep(1,nrow(dati)), lwd=0.5, col="red")
	segments(datj$pos, rep(0,nrow(datj)), datj$pos+datj$isize, rep(1,nrow(datj)), lwd=0.5, col="blue")

	invisible()
	}


g$peel <- function(x,i){
	# Peels off the i'th elements from a comma-separated character-based set of numbers
	# 	e.g., "peel(PV4,1)" peels the first element "0.34" of PV4 "0.34,1.5e-20,1,1"
	# Repeats for every entry if x is a vector
	# Returns a data frame if i is a vector
	# 	e.g., "peel(PV4,c(1,3))" returns the first and third elements of PV4
	z1 <- strsplit(x, split=",")
	z2 <- sapply(z1,function(x,i){x[i]},i[1])
	if(length(i) > 1){
		z3 <- data.frame(z2)
		for(j in 2:length(i)){
				z3[j] <- sapply(z1,function(x,i){x[i]},i[j])
			}
		names(z3) <- paste("peel",1:length(i), sep="")
		return(z3)
		}
	else return(z2)
	}


g$flag.readable <- function(i){
	# i must be a vector of integers, which are samtools flags
	# i is converted to a human-readable version of the flag
	# For example, i-69 is converted to:
	# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    #   01   00   01   00   00   00   01   00
	#    p    P    u    U    r    R    1    2  samtools flag codes
	# interpretation of the 8 binary numbers (1=yes, 0=no)
	# 1	read paired?
	# 2	read mapped in a proper pair?
	# 3   query sequence unmapped?
	# 4	mate is unmapped?
	# 5 	query is reverse strand?
	# 6   mate is reverse strand?
	# 7	read is first in pair?
	# 8	read is second in pair?
	# Note also that it is possible to query integer i using the bitops library
	# E.g., to query if mate is unmapped, use "bitAnd(x1$flag, 0x0004) == 0x0004"
	# (not implemented)
	# library(bitops)
	# 0x0001	p	the read is paired in sequencing
	# 0x0002	P	the read is mapped in a proper pair
	# 0x0004	u	the query sequence itself is unmapped
	# 0x0008	U	the mate is unmapped
	# 0x0010	r	strand of the query (1 for reverse)
	# 0x0020	R	strand of the mate
	# 0x0040	1	the read is the first read in a pair
	# 0x0080	2	the read is the second read in a pair
	z <- sapply(i, function(i){
		z <- as.logical(intToBits(i)[1:8])
		z1 <- c("p","P","u","U","r","R","1","2")[z]
		z2 <- paste(z1,collapse="")
		})
	return(z)
	}


g$gtf2thirdpositions <- function(gtffilename, refseqname, directory = "", check=TRUE, dropNonTriplets=TRUE,
	collapse=TRUE){
	# Analyzes a single chromosome
	# gtffilename is name of "group.gtf" file (in quotations)
	# refseqname is name of "chr.fa" file (in quotations)
	# if collapse=FALSE, returns a list by transcript, otherwise the list is collapsed into a single
	#	vector with unique bases
	# If collapse=FALSE, result names are a paste of transcript name : gene name
	# use as "z <- g$gtf2thirdpositions("groupXXI.gtf", "chrXXI.fa", check=TRUE, dropNonTriplets=TRUE, collapse=TRUE)"
		
	library(GenomicFeatures)
	library(Biostrings)
	
	# setwd("/Users/schluter/Documents/Research/genomics/reference genomes")
	# gtffilename <- "groupXXI.gtf"
	# refseqname <- "chrXXI.fa"
	
	# make the transcript data base
	gtf <- makeTranscriptDbFromGFF(file= gtffilename, format="gtf") # has class "TranscriptDb"
	
	# Pull out gene names corresponding to each transcript
	z <- select(gtf, keys = keys(gtf), keytype="GENEID",
		cols = c("CDSSTART", "CDSEND", "EXONID", "EXONSTRAND", "EXONSTART","EXONEND", "TXNAME"))
	# Warning message:
	# In .generateExtraRows(tab, keys, jointype) :
  	# 'select' resulted in 1:many mapping between keys and return rows
	transcriptGeneName <- tapply(z$GENEID, z$TXNAME, unique)


	# Extract GRanges from transcript data base
	grange <- transcripts(gtf)
	
	#head(grange)
      # seqnames           ranges strand |     tx_id            tx_name
   	     # <Rle>        <IRanges>  <Rle> | <integer>        <character>
  	# [1] groupXXI [240732, 252486]      + |         1 ENSGACT00000002162
  	# [2] groupXXI [280690, 284889]      + |         2 ENSGACT00000002166
  	# [3] groupXXI [308111, 318180]      + |         3 ENSGACT00000002174
  	# [4] groupXXI [308262, 317945]      + |         4 ENSGACT00000002176
  	# [5] groupXXI [614128, 678875]      + |         5 ENSGACT00000002201
  	# [6] groupXXI [698931, 733565]      + |         6 ENSGACT00000002210
	#length(grange)
	# [1] 599
	if(check) cat("Number of transcripts = ", length(grange), "\n")

	# Pull out all the coding regions by transcript
	# NOTE: expect different transcripts of the same gene to overlap?
	# The method uses an internal ordering to ensure that transcripts with the same name are unique, just in case,
	# but this makes it hard to identify the gene and transcript in the end. 
	# Can do "use.names=TRUE" if no duplicate names - this worked for groupXXI
	cdslist <- cdsBy(gtf, by = "tx", use.names=TRUE)

	#head(cdslist)
	# $ENSGACT00000002162   # these are transcript names
    	  # seqnames           ranges strand |    cds_id    cds_name exon_rank
    	     # <Rle>        <IRanges>  <Rle> | <integer> <character> <integer>
  	# [1] groupXXI [240732, 240840]      + |         1        <NA>         1
  	# [2] groupXXI [241074, 241090]      + |         2        <NA>         2
  	# [3] groupXXI [242286, 242437]      + |         3        <NA>         3
  	# [4] groupXXI [242593, 242612]      + |         4        <NA>         4
  	# [5] groupXXI [243731, 243789]      + |         5        <NA>         5
  	# [6] groupXXI [243941, 244036]      + |         6        <NA>         6
  	# [7] groupXXI [244147, 244257]      + |         7        <NA>         7
  	# [8] groupXXI [244463, 244599]      + |         8        <NA>         8
  	# [9] groupXXI [251813, 251912]      + |         9        <NA>         9
	#
	# $ENSGACT00000002166 
    	  # seqnames           ranges strand | cds_id cds_name exon_rank
  	# [1] groupXXI [281806, 281887]      + |     10     <NA>         2
  	# [2] groupXXI [281980, 282131]      + |     11     <NA>         3
  	# [3] groupXXI [282801, 283077]      + |     12     <NA>         4
  	# [4] groupXXI [284468, 284502]      + |     13     <NA>         5
	# ...

	# Number of transcripts represented
	# This is not 599 beause some transcripts are missing data?
	#length(cdslist) 
	# [1] 586
	#setdiff( grange$tx_name,  names(cdslist)) # these are the missing transcripts
 	# [1] "ENSGACT00000029571" "ENSGACT00000029368" "ENSGACT00000002915" "ENSGACT00000029474" "ENSGACT00000029188"
 	# [6] "ENSGACT00000029179" "ENSGACT00000029436" "ENSGACT00000028331" "ENSGACT00000029353" "ENSGACT00000029344"
	# [11] "ENSGACT00000028203" "ENSGACT00000029513" "ENSGACT00000028077"
	if(check) cat("Number of transcripts with coding region groups = ", length(cdslist), "\n")

	# Check that all cda's belonging to each transcript have the same strand ("+" or "-")
	z1 <- sapply(cdslist, function(x){
		length(unique(strand(x)))
		})
	if(check){
		cat("Frequency table indicating the number of unique strand directions per CDS group","\n",
			"(should be all 1)","\n")
		print(table(z1))
		}
		# z1
  		# 1
		# 586 

	# Check that the CDS are divisible by 3
	z1 <- sapply(cdslist, function(x){
		sum(width(x))/3
		})
		#head(z1)
		# ENSGACT00000002162 ENSGACT00000002166 ENSGACT00000002174 ENSGACT00000002176 ENSGACT00000002201 ENSGACT00000002210 
	               # 267                182                552                510                573               1039 
	z2 <- z1[z1!=round(z1,0)]
	if(check){
		cat("These transcripts do not have an integer number of triplets","\n")
		if(dropNonTriplets) cat("they have been removed","\n")
		print(z2)
		}	
		# ENSGACT00000003716 ENSGACT00000006614 
        		  # 441.3333           174.6667 
	
	#cdslist['ENSGACT00000003716'] # nothing unusual about this transcript
	
	if(dropNonTriplets & length(z2) > 0) cdslist <- cdslist[-which(is.element(names(cdslist),names(z2)))]

	# Check whether exons are always arranged along chr in increasing order (or decreasing in the case of strand = "-")
	if(check){
		 z1 <- sapply(cdslist, function(z){
			if(as.logical(strand(z)[1]=="+")) z1 <- all(start(z) == sort(start(z))) else 
			z1 <- all(start(z) == sort(start(z), decreasing=TRUE)) 
			})
	cat("Test of whether exons occur in sequential order along chromosome","\n")
	print(table(z1))
		# z1
		# TRUE 
 		# 584
		}

	# Identify bases in first, second, or third base positions of transcripts
	chromseq <- readDNAStringSet(refseqname, "fasta")[[1]]

	# Find the base pairs corresponding to third positions
	thirdpositions <- lapply(cdslist, function(z){
	# z <- cdslist[['ENSGACT00000002162']] # strand is "+"
	# z <- cdslist[['ENSGACT00000003737']] # strand is "-"
		z1 <- start(z) # start position for each range
		z2 <- end(z)   # end position of each range
		if(as.logical(strand(z)[1]=="+")){ 
			z3 <- mapply(FUN = seq, z1, z2) # seq(from=z1, to=z2) to get all nucleotides for each range
			z3 <- unlist(z3) 				# collapse nucleotide positions from all ranges to a single vector
			# all(z3==sort(z3))
			z4 <- z3[seq(3, length(z3), by = 3)] # start at the 3rd base and increment to end by 3 bases at a time
			} else{
			z3 <- mapply(FUN = seq, z2, z1) # seq(from=z2, to=z1) as above but reversing the order
			z3 <- unlist(z3) 				# collapse nucleotide positions from all ranges to a single vector
			# all(z3==sort(z3, decreasing=TRUE))
			# remember: z3 starts high and counts down, 
			# 	i.e., if head(z3) is (7102135 7102134 7102133 7102132 7102131 7102130), 
			# 	then the first instance of a third base position is 7102133
			# so this is the correct sequence of third-position bases:
			z4 <- z3[seq(3, length(z3), by = 3)] # starts at the 3rd from the end and increments down by 3
			}
		return(z4)
		})
	# thirdpositions['ENSGACT00000002162']
	# thirdpositions['ENSGACT00000003737']
	if(collapse) thirdpositions <- unique(unlist(thirdpositions))
	else{
		# Change names of the list to include the gene name, so that this information is retained
		z1 <- transcriptGeneName[match(names(thirdpositions), names(transcriptGeneName))]
		names(thirdpositions) <- paste(names(thirdpositions), z1, sep=":")
		}
	return(thirdpositions)
	}
	


g$tstv <- function(REF,ALT){
	# calculates the raw transition to transversion ratio R
	if(!(length(REF)==length(ALT))) stop("Vectors not of same length")
	ref <- rep(NA,length(REF))
	alt <- rep(NA,length(ALT))
	ref[is.element(REF,c("A","G"))] <- "pur"
	ref[is.element(REF,c("C","T"))] <- "pyr"
	alt[is.element(ALT,c("A","G"))] <- "pur"
	alt[is.element(ALT,c("C","T"))] <- "pyr"
	z <- table(ref,alt)
	R <- list(R = (z[1,1]+z[2,2]) / (z[1,2]+z[2,1]), tstv = z, tot = sum(z))
	return(R)
}
	
g$write.bedGraph <- function(filename, chrinfo = "chrXXI:1-11,717,487", from, to, stat,
	graphname = "DolphsGraph", col = "41,155,241", altcol = "100,100,100", ylinemark = 5,
	description = "chrXXI sliding window")
	{
	# bedGraph allows scientific notation in the stat variable, but the base lmits must be integers.
	z <- cbind.data.frame(chr = rep(gsub("[:][0-9,-]*", "",chrinfo), length(stat)), from = from, to = to, stat = stat)
	z <- z[!is.na(z$stat),]
	lim <- round(range(stat, na.rm=TRUE),5)
	lim <- paste(lim[1],lim[2],sep=":")
	write(file=filename, paste("browser position",chrinfo, sep=" ")) 
	write(file=filename, append=TRUE, 
		paste("track type=bedGraph name=", graphname, ' description="', description, '"',
	" color=",col, " altColor=",altcol, " viewLimits=",lim, " yLineOnOff=on yLineMark=", ylinemark, 
	" windowingFunction=mean autoScale=off visibility=full priority=20", sep=""))
	# "format(z)" suppresses scientific notation in the integer bases
	write.table(format(z), file=filename, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
	invisible()	
}

g$wc.revised <- function(ndat){
	# modified "wc" from hierfstat version 0.04-14 to speed it up
	library(hierfstat)
	# library(parallel)
	# ndat <- cbind(xpop,geno[,1:10000]) # test

    diploid <- TRUE # no longer an argument
    pol <- 0        # no longer an argument
    pop <- ndat[, 1]
    ni <- length(pop)
    dat <- ndat
    loc.names <- names(dat)[-1]

	# ---
	# The following code to count number of genotyped individuals per population is too slow, 
	# replace with subsequent block of code

    # system.time({
    # n <- t(ind.count(dat)) 
	# }) # system time 4-7 sec with 10000 markers analyzed on 2 groups size 6 each
    
    # system.time({
    popf <- factor(pop)
    z <- split(pop, popf)
    z1 <- lapply(z, function(z){
		z1 <- apply(dat[ which(pop %in% z), -1], 2, function(x){sum(!is.na(x))})
		})
	n <- do.call("cbind", z1)
	# }) # system time is about 1 second
	# (n[n[,1]==0 | n[,2]==0,])[1:100,] # many lines with one of the columns having a 0
	        # 1 2
	# V1      2 0
	# V22     1 0
	# V25     1 0
	# V42     0 4
	# V43     0 6
# ---
    
    nt <- apply(n, 1, sum, na.rm = TRUE)
    z <- table(nt)
    if(!is.na(z["0"])) warning("some markers have 0 individuals, which will affect length of output")
    untyped.loc <- which(nt == 0) # integer(0)
    typed.loc <- which(nt != 0)
    if (length(untyped.loc) > 0) {
        dat <- dat[, -(untyped.loc + 1)]
        n <- t(ind.count(dat))
        nt <- apply(n, 1, sum, na.rm = TRUE)
    	}
    
	# -------
	# The following code to count total number of alleles is too slow, 
	# replace with subsequent block of code

    # system.time({
    # alploc <- nb.alleles(cbind(rep(1, ni), dat[, -1]))
    # }) # takes about 10 seconds
    	
    # system.time({
    # Break up genotypes into alleles and keep the result -- slows things here but speeds up a later step
    # Missing alleles are set to "9" and then back to NA so I can split them
    allAlleles <- lapply(dat[ ,-1], function(x){
    		x[is.na(x)] <- 99
    		z <- unlist(strsplit(as.character(x), split = ""))
    	z[z == "9"] <- NA
    	return(z)
    		})
    	nalleles <- lapply(allAlleles, function(x){length(unique(na.omit(x)))})
    	alploc <- matrix(nalleles, ncol = 1, dimnames = list(loc.names, NULL))
    # }) # takes about 2 seconds
	# -------
    
    np <- dim(n)[2]
    npl <- apply(n, 1, tempfun <- function(x) sum(!is.na(x)))
    nl <- dim(n)[1]
    
    # -----
    # system.time({
    # p <- pop.freq(dat, diploid=diploid)
    # pb <- pop.freq(cbind(rep(1, length(pop)), dat[, -1]), diploid)
    # }) # took 22 seconds
    
    # system.time({
    popf <- factor(pop) 							# [1] 1 1 1 1 1 1 2 2 2 2 2 2
    popalleles <- rep(popf, rep(2, length(popf)))	# [1] 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2
    alleleFreqTable <- lapply(allAlleles, function(x){ 
#    alleleFreqTable <- mclapply(allAlleles, FUN=function(x){ # parallel lapply - does it work on hermes?
    		# x <- allAlleles[[175]]
		z1 <- table(x, popalleles)
		})	
    p <- lapply(alleleFreqTable, function(x){
		z1 <- sweep(x, 2, colSums(x), FUN = "/")
		})
	pb <- lapply(alleleFreqTable, function(x){
		z <- rowSums(x)
		z2 <- matrix(z/sum(z), ncol = 1, dimnames = list(names(z), "1"))
		})
    # }) # took 8 seconds
    
    # -----
    
    n <- matrix(unlist(n), ncol = np)
    nal <- n[rep(1:nl, alploc), ]
    nc <- (nt - apply(n^2, 1, sum, na.rm = TRUE)/nt)/(npl - 1)
    ntal <- rep(nt, alploc)
    ncal <- rep(nc, alploc)
    p <- matrix(unlist(lapply(p, t)), ncol = np, byrow = TRUE)
    pb <- matrix(unlist(pb), ncol = 1)
	dum <- getal.b(dat[, -1]) # slow
	all.loc <- apply(dum, 2, tempfun1 <- function(y) as.numeric(dimnames(table(y))[[1]])) # slow
	hetpl <- apply(dum, 2, fun <- function(z) { # slow
		lapply(as.numeric(dimnames(table(z))[[1]]), who.is.het <- function(y) apply(z == 
			y, 1, ind.is.het <- function(x) xor(x[1], x[2])))
			})
	mho <- lapply(hetpl, tempfun2 <- function(x) matrix(unlist(lapply(x, 
			tempfun3 <- function(y) tapply(y, pop, sum, na.rm = TRUE))), 
			ncol = np))
	mho <- matrix(unlist(mho), ncol = np, byrow = TRUE)
	mhom <- (2 * nal * p - mho)/2
    SSG <- apply(nal * p - mhom, 1, sum, na.rm = TRUE)
    dum <- nal * (p - 2 * p^2) + mhom
    SSi <- apply(dum, 1, sum, na.rm = TRUE)
    dum1 <- nal * (sweep(p, 1, pb))^2
    SSP <- 2 * apply(dum1, 1, sum, na.rm = TRUE)
    ntalb <- rep(npl, alploc)
    MSG <- SSG/ntal
    MSP <- SSP/(ntalb - 1)
    MSI <- SSi/(ntal - ntalb)
    sigw <- MSG
    sigb <- 0.5 * (MSI - MSG)
    siga <- 1/2/ncal * (MSP - MSI)
    Fxy <- function(x) x[1]/sum(x, na.rm = TRUE)
    FST.pal <- apply(cbind(siga, sigb, sigw), 1, Fxy)
    FIS.pal <- apply(cbind(sigb, sigw), 1, Fxy)
    loc <- rep(1:nl, alploc)
    lsiga <- tapply(siga, loc, sum, na.rm = TRUE)
    lsigb <- tapply(sigb, loc, sum, na.rm = TRUE)
    lsigw <- tapply(sigw, loc, sum, na.rm = TRUE)
    sigloc <- cbind(lsiga, lsigb, lsigw)
    if (length(untyped.loc) > 0) {
        x <- order(c(typed.loc, untyped.loc))
        dum1 <- matrix(numeric(3 * length(untyped.loc)), ncol = 3)
        dum <- rbind(sigloc, dum1, deparse.level = 0)[x, ]
        sigloc <- dum
    	}
    sigloc <- data.frame(sigloc)
    names(sigloc) = c("lsiga", "lsigb", "lsigw")
    rownames(sigloc) <- loc.names
    lFST <- apply(cbind(lsiga, lsigb, lsigw), 1, Fxy)
    lFIS <- apply(cbind(lsigb, lsigw), 1, Fxy)
    tsiga <- sum(siga, na.rm = TRUE)
    tsigb <- sum(sigb, na.rm = TRUE)
    tsigw <- sum(sigw, na.rm = TRUE)
    tFST <- Fxy(c(tsiga, tsigb, tsigw))
    tFIS <- Fxy(c(tsigb, tsigw))
    # res <- list(call = cl, sigma = cbind(loc, siga, sigb, sigw), 
        # sigma.loc = sigloc, per.al = list(FST = FST.pal, FIS = FIS.pal), 
        # per.loc = list(FST = lFST, FIS = lFIS), FST = tFST, FIS = tFIS)
    res <- list(sigma = cbind(loc, siga, sigb, sigw), 
        sigma.loc = sigloc, per.al = list(FST = FST.pal, FIS = FIS.pal), 
        per.loc = list(FST = lFST, FIS = lFIS), FST = tFST, FIS = tFIS)
    #}
    class(res) <- "wc"
    res
}

