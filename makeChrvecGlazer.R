#!/usr/bin/Rscript

# Converts chrvec for each chromosome (which contains windowsmaskerSdust info) 
# from old (Jones et al) to new assembly (Glazer et al).

# qsub -I -l walltime=02:00:00 -l mem=4gb 
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # chrname job ("splitup" or "rejoin")
# args <- c("chrXXI", "splitup")

chrname <- args[1]
job <- args[2]
chrno <- gsub("^chr", "", chrname)
chrNumeric <- g$chrname2numeric(chrname)
# [1] 21

if(job == "splitup"){
	# SPLITUP
	
	cat("\nChromosome being split is", chrname, "\n")
	chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")
	load(chrvecfile) # object name is "chrvec"
	POS <- seq(1:length(chrvec))
	
	if(chrno != "M" & chrno != "VIIpitx1" ){
		newCoords <- g$glazerConvertOld2New(chrname, POS)
		} else {
		newCoords <- data.frame(newChr = rep(chrno, length(POS)), newPos = POS)
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
		save(chrvecPart, file = paste("chrvec.", chrno, ".Part.", i, ".rdd", sep = "")) # object is chrvecPart
		}
	} else if(job == "rejoin"){

	# REJOIN
	
	cat("\nChromosome being rejoined is", chrname, "\n")
	# Results file
	chrvecfileNewName <- paste("chrvec.", chrno, ".glazer.rdd", sep = "") 
	# [1] "chrvec.XXI.glazer.rdd"
	
	# glob to identify all chrvecPart files, they will have chrNumeric (e.g., "21") before the ".rdd"
	chrvecfilePartNames <- paste("chrvec.*.Part", chrNumeric, "rdd", sep = ".")
	# [1] "chrvec.*.Part.21.rdd"

	# The corresponding file names
	z <- list.files( pattern=glob2rx( chrvecfilePartNames ), ignore.case=TRUE )
	# [1] "chrvec.Un.Part.21.rdd"  "chrvec.XXI.Part.21.rdd" 

	# Load the chrvecPart files, place in a temporary list: "parts"
	parts <- vector("list", length(z)) # initiate
	oldChrNames <- sapply( strsplit(z, split = "[.]"), function(x){x[2]})
	# [1] "chrUn"  "chrXXI"
	names(parts) <- oldChrNames

	for(i in 1:length(z)){
		load(z[i]) # name of object is chrvecPart
		parts[[i]] <- chrvecPart
		}
	
	chrvecFrame <- do.call("rbind", parts)
	# head(chrvecFrame)
		      # newPos base
	# Un.1 1792806    T
	# Un.2 1792807    A
	# Un.3 1792808    T
	# Un.4 1792809    C
	# Un.5 1792810    A
	# Un.6 1792811    G
	
	chrvec <- as.character( chrvecFrame$base[order(chrvecFrame$newPos)] )

	save(goodInvariants, file = chrvecfileNewName) # object is chrvec
	# load(chrvecfileNewName) # object is chrvec
	
	# length(chrvec)
	# [1] 17349772

	} else stop("Job argument must be 'rejoin' or 'splitup'")
	
	

