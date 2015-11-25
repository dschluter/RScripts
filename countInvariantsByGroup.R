#!/usr/bin/Rscript

# Run in Unix as " Rscript countGoodInvariantsByGroup.R ... "

# One strategy is to use a very low threshold when deciding whether to drop a base,
#	eg keep any base that has at least one genotype in at least 2 of the groups

# qsub -I -l walltime=03:00:00 -l mem=2gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname DPmin groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrXXI", "1", "paxl", "paxb", "marine-pac")
# args <- c(project, chrname, DPmin, groupnames)
# args <- c("BenlimAllMarine", "chrXXI", "paxl","paxb","pril","prib","qryl","qryb","ensl","ensb","marine-pac","marine-atl","marine-jap","solitary")

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)]

# convert snp to Glazer assembly coordinates
Glazerize 		<- TRUE # Requires file "glazerFileS4 NewScaffoldOrder.csv" in current directory

GTminFrac	<- 2/3
DPmin 		<- 1

# DPmin <- as.integer(args[3])
# GTminFrac <- eval(parse(text=args[4]))

GTmissing <- "."                     # how GATK represents missing genotypes in the vcf file "./."

# load "chrvec" for the current chromosome
chrno 					<- gsub("^chr", "", chrname)
chrmaskfile             <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd
load(chrmaskfile) 		# object is named "chrvec"

invariantsummaryname    <- paste(project, ".", chrname, ".DP.inv", sep="")
tempfile				<- paste(project, ".", chrname, ".goodInv.txt", sep="")
goodInvariantsFile 		<- paste(project, ".", chrname, ".goodInv.rdd", sep="")

nLinesAtaTime <- 1000000

INFILE <- file(invariantsummaryname, open = "r")
OUTFILE <- file(tempfile, "w")

x <- readLines(con = INFILE, n = nLinesAtaTime)
x1 <- strsplit(x, split = "\t")
headerline <- x1[[1]]
headerline[1] <- "POS"
x1 <- x1[-1]
nlines <- length(x1)

writeLines(paste(c(headerline[1:2], groupnames), collapse = "\t"), OUTFILE)

# assign numbers corresponding to groups
fishnames <- headerline[-c(1,2)]
groupcodes <- vector(length = length(fishnames), mode = "integer")
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
# groupcodes
 # [1]  8  8  8  8  8  8  7  7  7  7  7  7 10 10 10 11 11 11 11  9  9  9  9  9  9
# [26]  9  4  4  4  4  4  3  3  3  3  3  2  2  2  2  2  1  1  1  1  1  4  4  4  4
# [51]  3  3  3  3  6  6  6  6  6  6  5  5  5  5  5  5 12 12 12 12 12 12 12 12  2
# [76]  2  2  2  2  2  1  1  1  1  1  1
 
# Masked rows in the reference genome
maskedPOS <- which(chrvec == "M")

# Loop
while(nlines >0){
	
	x2 <- as.data.frame(do.call("rbind", x1), stringsAsFactors = FALSE)
	names(x2) <- headerline

	z <- x2$POS %in% maskedPOS
	x <- x2[!z,]
	
	# Elements of invariantsummaryname should be "0" or like "3:6" but some may be like ".:3"
	# Note: the columns haven't been divided into groups yet
	y <- lapply(x[, -c(1:2)], function(x){
	#		y <- as.integer(sub("([0-9]*)[:/][0-9]+", "\\1", x))
		y <- sub("([0-9.]*)[:][0-9.]+", "\\1", x)
		y[y == "."] <- "0"
		y <- ( as.integer(y) >= DPmin ) # TRUE if DP of a given genotype is bigger than DPmin
		})
	# Reunite columns
	y <- data.frame(y)
	
	# Split the data frame by group and take the sum of each row to get the total number of individuals that 
	# meet the DP minimum depth of coverage criterion DPmin within each group.
	z <- split(names(y), groupcodes)

	# cat("\nIdentified groups:\n")
	# print(z)
	
	# $`1`
	 # [1] "PaxLim.PxCL09maleLM1.GS14" "PaxLim.PxLfemale6.GS18"   
	 # [3] "PaxLim.PxLmale102.GS16"    "PaxLim.PxLmale106.GS15"   
	 # [5] "PaxLim.PxLmale107.GS17"    "paxl01"                   
	 # [7] "paxl05"                    "paxl09"                   
	 # [9] "paxl10"                    "paxl13"                   
	# [11] "paxl14"                   
	
	# $`2`
	 # [1] "PaxBen.PxBmale5.GS11"        "PaxBen.PxBmale6.GS12"       
	 # [3] "PaxBen.PxBmale8.GS10"        "PaxBen.PxCL09femaleBF6.GS13"
	 # [5] "PaxBen.RPxCL09maleBM2.GS9"   "paxb04"                     
	 # [7] "paxb05"                      "paxb06"                     
	 # [9] "paxb07"                      "paxb08"                     
	# [11] "paxb09"                     
	
	# $`3`
	# [1] "Marine.Pac.BIGR.52_54_2008.02"  "Marine.Pac.Japan.01.Katie"     
	# [3] "Marine.Pac.LittleCampbell.LC1D" "Marine.Pac.MANC_X_X05"         
	# [5] "Marine.Pac.Oyster.06.Sara"      "Marine.Pac.Seyward.01.Sara"    
	# [7] "Marine.Pac.WestCreek.01.Sara"  
	
	# Count how many genotypes per group
	z1 <- lapply(z, function(z){
		z1 <- apply(y[, z], 1, sum) # summing logicals treats TRUE as 1's
		})
		 
	goodInvariants <- data.frame(x[, c(1:2)], z1)
	
	# Drop rows not having at least 1 good genotype in at least 2 groups
	# Count number of groups having at least 1 good genotype
	z2 <- lapply(z1, function(z){z >= 1})
	z2 <- data.frame(z2)
	z3 <- apply(z2, 1, sum)
	
	# table(z3)
	
	keep <- z3 >= 2
	goodInvariants <- goodInvariants[keep, ]
	goodInvariants <- apply(goodInvariants, 1, function(x){paste(x, collapse = "\t")})
	
	if(length(goodInvariants) > 0) writeLines(goodInvariants, OUTFILE)
	
	x <- readLines(con = INFILE, n = nLinesAtaTime)
	x1 <- strsplit(x, split = "\t")
	nlines <- length(x)

	} # end while loop

close(INFILE)
close(OUTFILE)

# Read the whole thing into memory and save
goodInvariants <- read.table(file = tempfile, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

if(Glazerize){ # Requires conversion file "glazerFileS4 NewScaffoldOrder.csv" in current working directory
	pos <- goodInvariants$POS
	if(chrno != "M" & chrno != "VIIpitx1" ){
		chrNumeric <- chrno
		chrNumeric[chrno != "Un"] <- as.numeric( as.roman( chrNumeric[chrno != "Un"] ) )
		newCoords <- g$glazerConvertCoordinate(rep(chrNumeric, length(pos)), pos)
		
		} else {
		
		newCoords <- data.frame(newChr = rep(chrno, length(pos)), newPos = pos)
		
		}

	newChr <- unname(unlist(newCoords[ ,"newChr"]))
	newPos <- unname(unlist(newCoords[ ,"newPos"]))
	
	goodInvariants <- cbind.data.frame( data.frame(newChr = newChr, newPos = newPos), goodInvariants)
	
	z<- unique(newChr)
	# [1] "21" "Un"
	
	goodList <- split(goodInvariants, goodInvariants$newChr) # make sure that the list elements are named "21" and "Un"
	print(names(goodList)) # checking names of split data set
	
	for(i in z){ # saved object is "goodInvariantsPart"
		goodInvariantsPart <- goodList[[i]]
		save(goodInvariantsPart, file = paste(project, chrname, "goodInvPart", i, "rdd", sep = "."))
		# load(goodInvariantsPart)  # object named "goodInvariantsPart"
		}
		
	} # end if(Glazerize)

save(goodInvariants, file = goodInvariantsFile)
# load(goodInvariantsFile)  # object named "goodInvariants"

