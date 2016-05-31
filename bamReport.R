#!/usr/bin/Rscript
# Creates a pdf file with a visual summary of the contents of one or more bam file.
# Requires bai files of the same name in the same directory as the bam file
# Each bam file is read into memory, so for large bam files pay attention to RAM available

# Requires a chromosome name as last argument or embedded in name of first argument
# 	because otherwise memory requirements too great for all chromosomes at once.

# Run in Unix as " Rscript bamReport.R bamfile1 bamfile2, ... & "
# 				 " Rscript bamReport.R paxb*.chrXXI.realigned.bam & "

# setwd("/Users/schluter/Documents/Research/genomics/scripts")
# bamfile <- "Oyster12.chrXXI.realigned.bam"; baifile <- "Oyster12.chrXXI.realigned.bai"

# qsub -I -l walltime=04:00:00 -l mem=10gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# module load R/3.1.2
args <- commandArgs(TRUE) # only arguments will be stored
# args <- "paxb*.chrXXI.realigned.bam"
# args <- c("PRIL17.recal.bam", "chrXXI")

# test for regular expression as first argument
fnames <- list.files(pattern = glob2rx(args[1])) 
args <- c(fnames, args[-1])

# grab chromosome name: provided as last argument or is embedded in file name of first argument
chrname <- args[length(args)]
if(grepl("chr", chrname)){ 		# if last arg is a chr
	args <- args[-length(args)] # drop chr from vector of args
	} else{
	if(!grepl("^.*(chr[XVIUnMpitx1new]+).*$",args[1])) stop("Missing chr argument") else
		chrname <- grepl("^.*(chr[XVIUnMpitx1new]+).*$",args[1])
	}

if(any(grep("bai$", args) > 0)) stop("Provide only bam files as arguments (bai files assumed to have same name)")

cat("Loading RSamtools ...\n")
library(Rsamtools, quietly = TRUE)
library(bitops)

nargs <- length(args)
maxLen <- 10^8 # a number larger than the largest chromosome

for(i in 1:nargs){
	# i <- 1
	bamfile <- args[i]
	root <- gsub("(.)[.]bam$", "\\1", bamfile)
	bamfile <- paste(root, "bam", sep = ".")
	# baifile <- paste(root, "bai", sep = ".")
	baifile <- root


	cat("Reading bamfile. See ?scanBam for explanation of columns.\n")
	# x <- scanBam(bamfile, index = baifile)
	x <- scanBam(bamfile, index = baifile, param = ScanBamParam(
			what=c("pos","qwidth","flag","rname","strand","qual","isize","mapq","isize","seq"),
			which=GRanges(seqnames = chrname, IRanges(1, maxLen)) ))
	# x1 <- x[[1]]
	# rm(x)
	
	# names(x[[1]])
	 # [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"   
	# [13] "qual"
	
	# Calculate coverage (code by Michael Dondrup at https://www.biostars.org/p/5165/)
  	# filter reads without matching positions
  	# Uses: pos, qwidth, qname
  	ind <- !is.na(x[[1]]$pos)
	bam <- lapply(x[[1]], function(x) x[ind]) # strips the bam object of unaligned reads

	# This next command assembles all the reads, with start, end, and width into an IRanges object. 
  	ranges1 <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))

  	# This next command assembles the reads into a GRanges object
  	# This turns out not to be necessary if you don't want mean coverage calculated separately for each chromosome
	#ranges2 <- GRanges(seqnames=Rle(bam$rname), ranges=ranges1, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
	#mean.coverage <- mean(coverage(ranges2))
	cover <- as.vector(coverage(ranges1))
	mean.coverage <- mean(cover)
	median.coverage <- median(cover)

	medianwidth <- median(x[[1]]$qwidth, na.rm = TRUE)

	# Query the flag for specific items, eg is mate unmapped?
	# These might now be redundant since x[[1]]$pos is NA if read didn't align
	#	and x[[1]]$mpos is NA if mate didn't align
	# codes are below
	# 0x0001	p	the read is paired in sequencing
	# 0x0002	P	the read is mapped in a proper pair*
	# 0x0004	u	the query sequence itself is unmapped
	# 0x0008	U	the mate is unmapped
	# 0x0010	r	strand of the query (1 for reverse)
	# 0x0020	R	strand of the mate
	# 0x0040	1	the read is the first read in a pair
	# 0x0080	2	the read is the second read in a pair

	x[[1]]$mapped <- bitAnd(x[[1]]$flag, 0x0004) != 0x0004  
	x[[1]]$mappedmate <- bitAnd(x[[1]]$flag, 0x0008) != 0x0008
	x[[1]]$mapped.in.properpair <- bitAnd(x[[1]]$flag, 0x0002) == 0x0002

	cat("Mean coverage, mapped reads:", round(mean.coverage,2), "\n")
	cat("Median coverage, mapped reads:", median.coverage, "\n")
	cat("Total number of reads in bamfile:", length(x[[1]]$pos), "\n")
	cat("Chromosomes included:"); print(table(x[[1]]$rname, useNA="ifany"))
	cat("Number of mapped reads:", sum(x[[1]]$mapped), "\n") # same as length(na.omit(x[[1]]$pos))
	cat("Number of reads whose mate also mapped:", sum(x[[1]]$mapped & x[[1]]$mappedmate), "\n")
	cat("Number of reads mapped in a proper pair:", sum(x[[1]]$mapped.in.properpair), "\n")
	cat("Strands to which reads are aligned:"); print(table(x[[1]]$strand))
	cat("Median query width:", medianwidth, "\n")

	cat("\nGenerating plots\n")
	cat("pos is the genomic coordinate at the start of the alignment\n")
	cat("Coordinates are left-most, at the 5' end of read on '+' strand\n") 
	cat("     and at the 3' end on '-' strand; and are 1-based.\n")
	cat("Position excludes clipped nucleotides, though soft-clipped nucleotides\n")
	cat("     are included in seq.\n")
	cat("...........\n")

	pdf(paste(root, chrname, "Report", "pdf", sep = "."))
	# Plot base quality by cycle
	# These plots give guidance on where to trim (-q option in bwa aln)
	# Mean quality scores >=30 are very good
	# Convert ASCII quality scores to integer values
	# Need to separate those on - and + strands because reads oriented differently
	# Analyze just the reads that are "meadianwidth" bases 
	#     (conversion to matrix won't work unless all are equal width)
	plus <- x[[1]]$qual[ which(x[[1]]$strand=="+" & width(x[[1]]$seq) == medianwidth) ]
	minus <- x[[1]]$qual[ which(x[[1]]$strand=="-" & width(x[[1]]$seq) == medianwidth) ]
	mplus <- as(plus, "matrix")
	mminus <- as(plus, "matrix")
	#plot(colMeans(mplus), type="l", ylim=c(0,max(mplus)), xlim=c(0, medianwidth))
	#lines(colMeans(mminus), col = "red", lty = 2)
	
	plot(colMeans(rbind(mplus, mminus)), type="l", ylim=c(0,max(30, mplus)), xlim=c(0, medianwidth),	
		xlab = "Base position", ylab = "Base quality", main = "Base quality by cycle", col = "red", lwd = 2)
	mtext(bamfile, side = 3)
	lines(c(0, medianwidth), c(30,30), lty = 2, lwd = 1.5)

	# reads <- as(x[[1]]$qual[width(x[[1]]$seq) == medianwidth], "matrix") # this doesn't work because of read strands
	# plot(colMeans(reads), type="l", ylim=c(0,max(mplus)), xlim=c(0, medianwidth),	
	# xlab = "Base position", ylab = "Base quality", main = "Base quality by cycle")

	text(0,24,paste("Mean coverage of mapped reads:", round(mean.coverage, 2)), pos=4)
	text(0,22,paste("Median coverage of mapped reads:", round(median.coverage, 2)), pos=4)
	text(0,20,paste("Total number of reads in bamfile:", length(x[[1]]$pos)), pos=4)
	text(0,18,paste("Number of mapped reads:", sum(x[[1]]$mapped)), pos=4)
	text(0,16,paste("Number of reads whose mate also mapped:", sum(x[[1]]$mapped & x[[1]]$mappedmate)), pos=4)
	text(0,14,paste("Number of reads mapped in a proper pair:", sum(x[[1]]$mapped.in.properpair)), pos=4)
	text(0,12,paste("Median query width:", medianwidth), pos=4)

	# Mapping quality histogram
	hist(x[[1]]$mapq, right=FALSE, main="Mapping quality, all reads", col="red")
	mtext(bamfile, side = 3)
	
	# Coverage frequency histogram
	hist(cover[cover <= 50], right=FALSE, breaks = 50, main="Coverage", col="red")
	mtext(bamfile, side = 3)

	# Median coverage across the chromosome
	# This code involves fooling R into plotting a barplot using "hist"
	k <- 50000
	nbases <- length(cover)
	ibase <- c( seq(1, nbases, k), nbases + 1) # bases marking breaks between steps of size k
	z <- hist(1:nbases, breaks = ibase, plot=FALSE, right = FALSE)
	coverbins <- cut(1:nbases, breaks = ibase, right = FALSE)
	y <- tapply(cover, coverbins, median)
	z$counts <- y
	z$breaks <- z$breaks/10^6
	plot(z, freq = FALSE, col = "red", ylab = "Median coverage", xlab="Position (million bases)", 
		main="Median coverage") # ignore warnings!
	mtext(bamfile, side = 3)

	# Insert sizes. Negatives refer to outward facing? Or to the direction of the placement? They seem symmetric.
	# range(x[[1]]$isize, na.rm=TRUE)
	hist(x[[1]]$isize[abs(x[[1]]$isize) <= 700], right = FALSE, col = "red", breaks=100, 
		main = "Insert sizes (those 700 bases or less)")
	mtext(bamfile, side = 3)

	# Read position
	hist(x[[1]]$pos/10^6, right=FALSE, main="Positions of all mapped reads", col="red",
		breaks = 200, xlab = "Position (million bases)")
	mtext(bamfile, side = 3)

	# Placement of reads whose mate did not map
	hist(x[[1]]$pos[!x[[1]]$mappedmate]/10^6, right=FALSE, main="Reads whose mate did not map", col="red",
		breaks = 200, xlab = "Position (million bases)")
	mtext(bamfile, side = 3)

	z<-quantile(x[[1]]$mapq, prob=c(0.25,0.5,0.75), na.rm = TRUE)
	# Placement of reads with mapping quality in lowest quartile
	hist(x[[1]]$pos[x[[1]]$mapq <= z["25%"]]/10^6, right=FALSE, 
		main="Reads with mapping quality in LOWEST quartile", col="red",
		breaks = 200, xlab = "Position (million bases)")
	mtext(bamfile, side = 3)

	# Placement of reads with mapping quality in highest quartile
	hist(x[[1]]$pos[x[[1]]$mapq >= z["75%"]]/10^6, right=FALSE, 
		main="Reads with mapping quality in HIGHEST quartile", col="red",
		breaks = 200, xlab = "Position (million bases)")
	mtext(bamfile, side = 3)

	# Placement of reads with large isize
	hist(x[[1]]$pos[abs(x[[1]]$isize) > 500]/10^6, right=FALSE, 
		main="Reads with excessive abs insert sizes (> 500)", col="red",
		breaks = 200, xlab = "Position (million bases)")
	mtext(bamfile, side = 3)

	dev.off()
	cat("Finished generating plots\n")
	}

