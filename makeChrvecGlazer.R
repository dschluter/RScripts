#!/usr/bin/Rscript

# Converts chrvec for each chromosome (which contains windowsmaskerSdust info) 
# from old (Jones et al) to new assembly (Glazer et al).
# Requires the file "" in the working directory

# qsub -I -l walltime=02:00:00 -l mem=4gb 
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # chrname job ("splitup" or "rejoin")
# args <- c("chrXXI", "splitup")
# args <- c("chrXXI", "rejoin")

chrname <- args[1]
job <- args[2]
chrno <- gsub("^chr", "", chrname)
chrNumeric <- g$chrname2numeric(chrname)
# [1] 21

scafFile = "glazerFileS4 NewScaffoldOrder.csv"

if(job == "splitup"){
	# SPLITUP
	
	cat("\nChromosome being split is", chrname, "\n")
	chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")
	load(chrvecfile) # object name is "chrvec"
	POS <- seq(1:length(chrvec))
	
	# *keep only the non-M to avoid overlap later*
	POS <- POS[chrvec != "M"]
	
	if(chrno != "M" & chrno != "VIIpitx1" ){
		newCoords <- cbind.data.frame( g$glazerConvertOld2New(chrname, POS), oldPos = POS, stringsAsFactors = FALSE )
		} else {
		newCoords <- data.frame(newChr = rep(chrno, length(POS)), newPos = POS, oldPos = POS)
		}
	
	# NewStart	NewEnd	OldChr	OldStart		OldEnd	# 1			  517387		Un	22323409	 	22840795	# 518388		 1018815		Un	22841796	 	23342223	# 1019816	 1276953		Un	29838944	 	30096081	# 1277954	 1493735		Un	31561396	 	31777177	# 1494736	 1791805		21	       1		  297070	# 1792806	 4441218		Un	 5106262		 7754674	# 4442219	13638880		21	  298071		 9494732	# 13639881	15542716		21	 9495733		11398568	# 15543717	16437086		Un	15742364	 	16635733	# 16438087 	16704853		21	11399569	 	11666335	# 16705854 	16996730		Un	28000487	 	28291363	# 16997731	17222261		Un	31335865 	31560395	# 17223262	17357772		Un	37062399 	37196909
	
	# Grab the bases corresponding to oldPos, so that they are associated with the newPos and newChr
	newCoords$base <- chrvec[newCoords$oldPos]

	# head(newCoords)
	  # newChr  newPos oldPos base
	# 1     21 1495159    424    T
	# 2     21 1495160    425    C
	# 3     21 1495161    426    T
	# 4     21 1495162    427    G
	# 5     21 1495163    428    A
	# 6     21 1495164    429    C
	
	# length(chrvec)
	# [1] 11717487

	# range(newCoords$oldPos)
	# [1]      424 11717487
	
	z <- unique(newCoords$newChr)
	# [1] 21 Un
	
	for(i in z){ # saved object is "chrvecPart"
		# i <- z[1]
		chrvecPart <- newCoords[newCoords$newChr == i, ]
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
	           # newChr  newPos  oldPos base
	# Un.4134155     21 1792806 5106262    T
	# Un.4134156     21 1792807 5106263    A
	# Un.4134157     21 1792808 5106264    T
	# Un.4134158     21 1792809 5106265    C
	# Un.4134159     21 1792810 5106266    A
	# Un.4134160     21 1792811 5106267    G
	
	# load("chrvec.Un.masked.rdd")
	# head( chrvec[chrvecFrame$oldPos ] )
	# [1] "T" "A" "T" "C" "A" "G"
	
	# table(chrvecFrame$base)
	      # A       C       G       T 
	# 3360296 2902636 2896769 3371368

	# Check that there are no redundant positions	
	if( length(unique(chrvecFrame$newPos)) != length(chrvecFrame$newPos) ) stop("POS not unique")

	# Determine the size of the new chromosome
	# max(chrvecFrame$newPos) # not infallible method
	# [1] 17357772

	x <- read.csv(scafFile, stringsAsFactors = FALSE) # better method
	chrvecMax <- max( x$NewEnd[as.character(x$NewChr) == as.character(chrNumeric)] ) 

	# chrvecMax
	# [1] 17357772 # size of the new chromosome
	
	chrvec <- rep("M", chrvecMax) # initiate
	chrvec[chrvecFrame$newPos] <- chrvecFrame$base
	
	# table(chrvec)
	      # A       C       G       M       T 
	# 3360296 2902636 2896769 4826703 3371368 


	save(chrvec, file = chrvecfileNewName) # object is chrvec
	# load(chrvecfileNewName) # object is chrvec
	
	} else stop("Job argument must be 'rejoin' or 'splitup'")
	
	

