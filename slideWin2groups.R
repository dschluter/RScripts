#!/usr/bin/Rscript
# Grabs fst or pairwise seguence divergence data from blockstats.rdd and
#	 calculates sliding window statistics in windows of nblocks per window.
# Current methods are "FST" or "CSS" (case insensitive)
# Files required in preparation for this script:
# 	genome.r									in scriptdir
#	project.chrname.blockstats.blocksize.rdd 	in vcfdir

# Run in Unix
# as " Rscript slideWin2groups.R infile method nblocks.per.window windowFracMin & "
#    " Rscript slideWin2groups.R "Paxton12fish.chrXXI.blockstats.rdd" "fst" 5 0.4 >&Rscript.out & "

# Defaults 
# These can be changed by including an alternative value as an argument IN THE RIGHT PLACE
nblocks.per.window <- 5  # the number of blocks included in each sliding window
windowFracMin <- 0.4     # the min fraction of good nucleotides in a window required to yield a result, 	else return NA

# These can't yet be changed
bedGraph <- FALSE							 


# Make sure that the latest version of genome.r is copied to scriptdir folder before sourcing this file

# ---
args <- commandArgs(TRUE)

# args <- c("Paxton12fish.chrXXI.blockstats.rdd", "css", 5, 0.4)

z <- unlist(strsplit(args[1], split = "[.]"))
project <- z[1]
chrname <- z[2]
method <- args[2]
if( !(casefold(method) %in% c("fst","css")) ) stop("Only methods FST and CSS available")
if(!is.na(args[3])) nblocks.per.window <- as.integer(args[3])
if(!is.na(args[4])) windowFracMin <- as.numeric(args[4])

cat("SlideWin settings:\n")
cat("\nChromosome name:", chrname, "\n")
cat("\nnblocks.per.window:", nblocks.per.window, "\n")
cat("\nwindowFracMin:", windowFracMin, "\n")

# westgrid is the default, otherwise this doesn't work in Rscript mode
#if( grepl("westgrid",Sys.info()["nodename"]) ){
	setwd("~")
	vcfdir <- ""
	genomedir <- ""
	invariantsummarydir <- vcfdir
	scriptdir <- ""
#	}

# To run on office computer (copy important files to the local genomics folder, don't access files on server - too slow)
if( grepl("greendrake",Sys.info()["nodename"]) ){
	setwd("/Users/schluter/Documents/genomics/")
	vcfdir <- ""
	genomedir <- ""
	invariantsummarydir <- vcfdir
	scriptdir <- "/Volumes/schluter/Dolph_data/R and S-Plus macros/"
	}

if( grepl("Wallacea",Sys.info()["nodename"]) ){
	setwd("/Users/schluter/Documents/Research/genomics/BenLim-q15/")
	vcfdir <- ""
	genomedir <- "/Users/schluter/Documents/Research/genomics/reference genomes/"
	invariantsummarydir <- vcfdir
	scriptdir <- "/Users/schluter/Documents/Research/R and SPlus macros/"
	}

if( grepl("SciBorg",getwd()) ){
	setwd("~/BenLim-q15/") 
	vcfdir <- ""
	genomedir <- "~/reference_genomes/"
	invariantsummarydir <- vcfdir
	scriptdir <- "~/scripts/"
	# folder <- recode(project, c("benlim", "paxbl", "priestbl", "qrybl", "enosbl"), 
							# c("BenLim-q15","Paxton-q15","Priest-q15","LQuarry-q15", "Enos-q15"))
	# setwd( paste("~", folder, "", sep="/") )  #
	# vcfdir <- paste("~", folder, "", sep="/") # 
	}
	
source(paste(scriptdir,"genome.r", sep = ""))

# Load blockstats.rdd file in vcfdir
load(file = args[1]) # blockstats
# head(blockstats)

# names(blockstats)
  # [1] "ibase"       "midbase"     "nUnmasked"   "nSnp"        "nInvariants" "VARa"        "VARb"       
  # [8] "VARw"        "fst"         "nTrueSnp"    "psd.1,1"     "psd.1,1"     "psd.1,1"     "psd.1,1"    
 # [15] "psd.1,1"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"    
 # [22] "psd.1,1"     "psd.1,1"     "psd.1,1"     "psd.1,1"     "psd.1,2"     "psd.1,2"     "psd.1,2"    
 # [29] "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,1"     "psd.1,1"     "psd.1,1"     "psd.1,2"    
 # [36] "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,1"     "psd.1,1"    
 # [43] "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,1"    
 # [50] "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"    
 # [57] "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.1,2"     "psd.2,2"     "psd.2,2"    
 # [64] "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"    
 # [71] "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"     "psd.2,2"     "np.1,1"     
 # [78] "np.1,1"      "np.1,1"      "np.1,1"      "np.1,1"      "np.1,2"      "np.1,2"      "np.1,2"     
 # [85] "np.1,2"      "np.1,2"      "np.1,2"      "np.1,1"      "np.1,1"      "np.1,1"      "np.1,1"     
 # [92] "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,1"     
 # [99] "np.1,1"      "np.1,1"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"     
# [106] "np.1,2"      "np.1,1"      "np.1,1"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"     
# [113] "np.1,2"      "np.1,2"      "np.1,1"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"     
# [120] "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"      "np.1,2"     
# [127] "np.1,2"      "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"     
# [134] "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"      "np.2,2"     
# [141] "np.2,2"      "np.2,2"     

blocksize <- blockstats$ibase[2] - blockstats$ibase[1]
cat("\nBlock size:", blocksize, "\n")
windowNmin <- blocksize * nblocks.per.window * windowFracMin
cat("\nMinimum number of good bases in a window:", windowNmin, "\n")

# --------------------------------------------
# Sliding window analyses of group differences
# --------------------------------------------
resultWin <- g$slidewin(blockstats, 
						stepsize = blocksize, 
						nsteps.per.window = nblocks.per.window,
						method = method, 
						windowNmin = windowNmin
						)

# Fst
if(casefold(method) == "fst"){
	pdf( file = paste(project, chrname, "FSTslidewin.pdf", sep = ".") )
	
	plot(fst ~ I(ibase/10^6), data = resultWin, type = "l", col = "red", 
		xlab = "million bases", las = 1, ylab = "Fst", main = paste("block size = ", blocksize,
		", no. blocks per window = ", nblocks.per.window, " (", blocksize * nblocks.per.window, " bases)",
		"\nmin fraction good bases per window = ", windowFracMin, " (", windowNmin, " bases)", sep = "") )

	plot(I(100*totVARa/nbases) ~ I(ibase/10^6), data = resultWin, type = "l", col = "red", 
		xlab = "million bases", las = 1, ylab = "totVARa per 100 bases")
	plot(I(100*(totVARb+ totVARw)/nbases) ~ I(ibase/10^6), data = resultWin, pch = ".", col = "red", 
		xlab = "million bases", las = 1, ylab = "totVARb + totVARw per 100 bases")
	plot(I(100*totVARb/nbases) ~ I(ibase/10^6), data = resultWin, type = "l", col = "red", 
		xlab = "million bases", las = 1, ylab = "totVARb per 100 bases")
	plot(I(100*totVARw/nbases) ~ I(ibase/10^6), data = resultWin, type = "l", col = "red", 
		xlab = "million bases", las = 1, ylab = "totVARw per 100 bases")
	
	dev.off()
	
	if(bedGraph){
		# write VARa to bedgraph - make sure chromosome positions are "0-based half open"
		g$write.bedGraph( "VARa.bedGraph", chrinfo = chrname, from = resultWin$ibase,
			to = (resultWin$ibase + blocksize-1), stat = 100*resultWin$totVARa/resultWin$nbases, 
			graphname = "totVARa", col = "41,155,241", altcol = "100,100,100", ylinemark = 2, 
			description = paste(chrno, "totVARa per 100 bp sliding window") ) 
		}
	}


# CSS		
if(casefold(method) == "css"){

	pdf( file = paste(project, chrname, "CSSslidewin.pdf", sep = ".") )

	plot(I(100*CSS/nbases) ~ I(ibase/10^6), data = resultWin, type = "l", col = "red",
		xlab = "million bases", las = 1, ylab = "CSS per 100 bases", main = paste("block size = ", blocksize,
		", no. blocks per window = ", nblocks.per.window, " (", blocksize * nblocks.per.window, " bases)",
		"\nmin fraction good bases per window = ", windowFracMin, " (", windowNmin, " bases)", sep = "") )

	dev.off()
	# range(100*resultWin$CSS/resultWin$nbases, na.rm=TRUE)
	# [1] -0.1011879  2.6064735
	
	if(bedGraph){
		g$write.bedGraph( "CSSwin.bedGraph", chrinfo = "chrXXI", from = resultWin$ibase,
			to = (resultWin$ibase + blocksize-1), stat = 100*resultWin$CSS/resultWin$nbases, 
			graphname = "CSSwin", col = "41,155,241", altcol = "100,100,100", ylinemark = 2, 
			description = paste(chrno, "CSS per 100 pb sliding window") ) 
			}
	}