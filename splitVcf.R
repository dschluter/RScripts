#!/usr/bin/Rscript
# Code to extract variant sites from a .vcf file 
#	and to save the invariant DP and RGQ values to a second file named etc.DP.inv

# Assumes that VCF file is named "project.chr.vcf"
# Run on hermes as:
#	/global/software/R-3.1.2/bin/Rscript splitVcf.R Paxton12fish.chrXXI.vcf

args <- commandArgs(TRUE) # get arguments

vcffile <- args
#vcffile <- "Paxton12fish.chrXXI.vcf"
#vcffile <- "temp.vcf"

# Extract project and chr name
# project <- gsub("(^[A-z0-9_-]*[.]chr[A-z1]+)[.]vcf", "\\1", vcffile)
project <- gsub("[.]vcf$", "", vcffile)

# OPTIONS
nFirstLines <- 1000    # first number of lines to read to extract header
nLinesAtaTime = 100000

# Open file and read nFirstLines lines from vcf file
INFILE <- file(vcffile, open = "r")

x <- readLines(con = INFILE, n = nFirstLines)

# Header line operations (do just once)
nhead <- length(grep('^#', x)) # get header lines ('#')
if(nhead == nFirstLines) stop("Increase nFirstLines to get whole header")
headerlines <- x[1:nhead]
x <- x[-c(1:nhead)] # keep the non-header part of x
nlines <- length(x)

# Open 2 output file connections and write headers
outfile1 <- file(paste(project, "var", "vcf", sep="."), "w")
outfile2 <- file(paste(project, "DP", "inv", sep="."), "w")

writeLines(headerlines, outfile1)

# only the last line of the header goes into the invariant file, but also keep POS REF and fish names
invheader <- gsub("(^#)CHROM\t(POS)\tID(\tREF)\tALT\tQUAL\tFILTER\tINFO\tFORMAT(\t.+$)", "\\1\\2\\3\\4",headerlines[nhead])
writeLines(invheader, outfile2)

# Loop
while(nlines >0){

	# Split according to ALT (5th column) and write to separate files
	snp <- gsub("^[a-zA-Z0-9._]+[\t][0-9]+[\t][.][\t][ACGTNacgtn]+[\t]([ACGTacgtn.]+)[\t].+$", "\\1", x) # extract 5th column
	x1 <- split(x, snp == ".")
	# names(x1) # [1] "FALSE" "TRUE" ( TRUE = invariants )
	
	# Write the variants
	if(!is.null(x1[["FALSE"]])) writeLines(x1[["FALSE"]], outfile1) # the variants
	
	# Simplify the invariants, keep only POS, REF and simplified genotype information (DP and RGQ)
	if(!is.null(x1[["TRUE"]])){
		x2 <- strsplit(x1[["TRUE"]], split = "\t")  # splits up rows of invariants. 
		startCol <- 10             # Genotypes start in column 10
		endCol <- length(x2[[1]])  # number of columns in invariant file
		# Harvest the DP and RGQ values from the genotype values in x2
		x3 <- sapply(x2, function(z){ z[9] }) # 9th column has format "GT:AD:DP:PGT:PID:RGQ" "GT:AD:DP:RGQ" "GT:DP:RGQ" etc
		format <- strsplit(x3, split = ":")   # Splits up the format "GT:DP:RGQ" to "GT" "DP" "RGQ"
		is.DP <- lapply(format, function(z){ which(z == "DP") })   # indicates which element of format is "DP"
		is.RGQ <- lapply(format, function(z){ which(z == "RGQ") }) # indicates which element of format is "RGQ"
		x4 <- mapply(is.DP, is.RGQ, x2, FUN = function(i, j, z){
			# i <- is.DP[[1]]; j <- is.RGQ[[1]]; z <- x2[[1]]
			z1 <- strsplit( z[startCol:endCol], split = ":" ) # grab genotype columns, split by colon
			z2 <- sapply(z1, function(z1){ paste(z1[i], z1[j], sep = ":") }) # grab only DP and RGQ elements
			z2 <- gsub("0:0", "0", z2) # simplify double zero to single
			return( paste( c(z[c(2,4)], z2), collapse = "\t") ) # keep only POS and REF from row
			}, SIMPLIFY = TRUE)
	
		writeLines(x4, outfile2)
		}
		
	x <- readLines(con = INFILE, n = nLinesAtaTime)
	nlines <- length(x)

	} # end while loop

close(INFILE)
close(outfile1)
close(outfile2)

