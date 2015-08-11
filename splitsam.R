#!/usr/bin/Rscript
args <- commandArgs(TRUE) # first argument will be samfile name

# Run in Unix as " /Linux/R-3.0.1/bin/Rscript splitsam.R subset.sam & "
#             or " /global/software/R-3.1.2/bin/Rscript splitsam.R subset.sam & "
#             or "qsub -v samfile=ensb_01.sam splitsam.R.pbs"

# To practice, extract a few lines from the samfile and save in temp.sam
# head -50000 paxl_01.noclip.sam > temp.sam

# This version keeps singletons (reads in which one of the pair doesn't map or maps to another chromosome)
# in the same chromosome file as the paired-end reads.

# setwd("/Users/schluter/Documents/Research/genomics/r-scripts")

# Define function that does the work (command comes afterward)

splitsam <- function(samfile, nFirstLines = 10^2, nLinesAtaTime = 10^3){
	# Open file and read nFirstLines lines from vcf file
	# --------------------------------------------------
	INFILE <- file(samfile, open="r")
	#close(INFILE)

	# according to http://www.r-bloggers.com/faster-files-in-r/ the following reads much faster
	# I'd need to fiddle with the last line if it is not all the way to the "\n" at the end
	# buf = readChar( fname, s, useBytes=T)
	# strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]

	x <- readLines(con = INFILE, n = nFirstLines) 

	# get number of header lines ('@')
	nhead <- length(grep('^@', x))
	nSQ <- nhead - 1
	nPG <- 1
	headerlines <- x[1:nhead]
	
	# Keep the non-header part of x
	x <- x[-c(1:nhead)]
	nlines <- length(x)
		
	samname <- gsub("(^[A-Za-z0-9_.\\/-]+)[.]sam", "\\1", samfile)
	#outdir <- paste("splitsam",samname, sep=".")
	
	chrnames <- gsub("^@[SQPG]+\tSN:([A-Za-z0-9.-_]+)\t[LN:0-9]+", "\\1", headerlines[-nhead]) # awkward but faster
	#chrnames <- sapply( strsplit(headerlines[-nhead], split = "[\t:]+"), function(z){z[3]}) # takes twice as long
	
	# Open multiple write connections
	outfiles <- vector("list", length(chrnames))
	for(i in chrnames) outfiles[[i]] <- file(paste(samname, i, "sam", sep="."), "w")
	badoutfile <- file(paste(samname, "chrBAD", "sam", sep="."), "w")
	
	# Write header lines to files
	for(i in seq_along(headerlines[-nhead])){
		writeLines(headerlines[i], outfiles[[chrnames[i]]])
		writeLines(headerlines[nhead], outfiles[[chrnames[i]]])
		}
	
	#for(i in chrnames) close(outfiles[[i]])
	#close(badoutfile)

	# Loop
	while(nlines > 0){

		# Need to extract third column from x (identifies chromosome)
		# One way is to strsplit the strings for this to work easily (it is much faster with fixed = TRUE)
		chrom <- sapply( strsplit(x, split="\t", fixed = TRUE), function(x){x[3]})
		
		# Another way is to grab the third column using gsub - this takes twice the time as the above
		#chrom <- gsub("^[A-z0-9_:-]+\t[0-9]+\t([^\t]+)+.+", "\\1", x) # "[^\t]" is not \t, "." stands for anything
	
		# Split x by chromosome
		xsplit <- split(x, chrom)
	
		# Put the bad reads into badoutfile and then delete them from xsplit
		if(length(xsplit[["*"]]) > 0){
			writeLines(xsplit[["*"]], badoutfile)
			xsplit["*"] <- NULL
			}
	
		# print to files
		for(i in names(xsplit)) writeLines(xsplit[[i]], outfiles[[i]])

		x <- readLines(con = INFILE, n = nLinesAtaTime)
		nlines <- length(x)
	
		} # end while loop

	close(INFILE)
	for(i in chrnames) close(outfiles[[i]])
	close(badoutfile)

	} # end of splitsam function

splitsam(samfile = args[1])
