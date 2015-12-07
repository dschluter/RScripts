#!/usr/bin/Rscript

# Converts chrvec for each chromosome (which contains windowsmaskerSdust info) 
# from old (Jones et al) to new assembly (Glazer et al).

# qsub -I -l walltime=02:00:00 -l mem=4gb 
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # chrname job ("split" or "join")
# args <- c("chrXXI")

chrname <- args[1]
job <- args[2]
chrno <- gsub("^chr", "", chrname)
chrNumeric <- g$chrname2numeric(chrname)
# [1] 21

if(job == "split"){
	cat("\nChromosome being split is", chrname, "\n")
	chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")
	load(chrvecfile) # object name is "chrvec"
	POS <- seq(1:length(chrvec))
	
	if(chrno != "M" & chrno != "VIIpitx1" ){
		newCoords <- g$glazerConvertOld2New(chrname, POS)
		} else {
		newCoords <- data.frame(newChr = rep(chrno, length(pos)), newPos = pos)
		}
	
	# head(newCoords)
	  # newChr  newPos
	# 1     21 1494736
	# 2     21 1494737
	# 3     21 1494738
	# 4     21 1494739
	# 5     21 1494740
	# 6     21 1494741
	
	z <- unique(newCoords$newChr)
	# [1] 21 Un
	
	for(i in z){ # saved object is "chrvecPart"
		chrvecPart <- chrvec[newCoords$newChr == i]
		chrvecPart <- data.frame(newPos = newCoords$newPos[newCoords$newChr == i], base = chrvecPart)
			   # newPos base
		# 1 1494736    M
		# 2 1494737    M
		# 3 1494738    M
		# 4 1494739    M
		# 5 1494740    M
		# 6 1494741    M
		save(chrvecPart, file = paste("chrvec.", chrno, ".Part.", i, ".rdd", sep = ""))
		}
	} else if(job == "join"){
	} else stop("Job argument must be 'join' or 'split'")
	
	

