# Methods to read and write DNA sequences in R

setwd("/Users/schluter/Documents/Research/genomics/reference genomes")

# ---------------------------------------------------------
# Read BAC pitx1 sequence from GenBank using accession number
# According to Chan et al, the deletion in PAXB is 1868bp

library(ape)

# accession number for the Pitx1 BAC from Salmon River
z <- read.GenBank("GU130435.1") # "DNAbin" object 
z
# 1 DNA sequence in binary format stored in a list.
# Sequence length: 377852 
# Label: GU130435.1 
# Base composition:
    # a     c     g     t 
0.265 0.236 0.237 0.263 

length(z[[1]])
#[1] 377852
names(z)
[1] "GU130435.1"
# names(z) <- "chrVIIpitx1"
# write.dna(z, "chrVIIpitx1.fa", format="fasta", colsep = "") # beware of blanks between columns

# is all lower case
write.dna(z, "GU130435.1.fa", format="fasta", colsep = "") # beware of blanks between columns

# Read in as a list - is all lower case
z <- read.dna("GU130435.1.fa", format="fasta", as.character=TRUE, as.matrix=FALSE) # reads as a vector

head(z[[1]], 20)
# [1] "g" "a" "a" "t" "t" "c" "t" "t" "t" "c" "c" "t" "t" "c" "c" "t" "c" "g" "t" "t"
length(z[[1]])
# [1] 377852
table(z[[1]], useNA = "always")
     # a      c      g      k      m      r      s      t      w      y   <NA> 
# 100189  88982  89440      4      4      5      2  99216      7      3      0 

# Convert to uppercase
z[[1]] <- toupper(z[[1]])

# eliminate the k's etc
z[[1]] <- recode(z[[1]], c("G","A","T","C","S","W","K","M","Y","R"), c("G","A","T","C","N","N","N","N","N","N") )
table(z[[1]], useNA = "always")
     # A      C      G      N      T   <NA> 
# 100189  88982  89440     25  99216      0 

# rename
names(z)[1] <- "chrVIIpitx1"

# write to new fasta file
write.dna(z, "test.fa", format="fasta", colsep = "") # beware of blanks between columns


# z <- read.FASTA("chrVIIpitx1.fa") # always returns a list of class "DNAbin".

# ---------------------------------------------------------
# Read Felicity's new pitx1 contig

# When running these analyses in January [2016], I realised that Jane Grimwood had updated the BAC sequences 
# for this region.  As a result I stitched together a new contig and used that in my analyses 
# (hence the unfortunate flip in Pitx1 orientation). FYI, I have uploaded this new contig onto the ftp site too: 
# You will see it as the file named:  SALR.Pitx1.118G22-164F21.Combined.fa

# Feb 2, 2016, email: "Pel enhancer is located at 248kb. Pitx1 216->206kb"
# 	(this is presumably the reverse direction position - this BAC is reversed compared to GU130435.1
#	WARNING: converts to lower case!
z1 <- read.dna("SALR.Pitx1.118G22-164F21.Combined.fa", format="fasta", as.character=TRUE, 
		as.matrix=FALSE)
names(z1)
# [1] "CH213-118G22.fulllength+CH213-164F21.pitx1.5'truncated"
head(z1[[1]], 20)
 # [1] "g" "a" "a" "t" "t" "c" "a" "g" "t" "a" "c" "g" "a" "t" "t" "c" "t" "a" "a" "a"
length(z1[[1]])
# [1] 377804 # is shorter than GU130435.1 by 48 bases, no "n"

# Convert to uppercase
z1[[1]] <- toupper(z1[[1]])

table(z1[[1]], useNA = "always")
     # A      C      G      T   <NA> 
 # 99205  89429  88981 100189      0 
 
# Try sequinr and ape commands
# library(seqinr)
# # ignore warning
# z1 <- read.fasta("SALR.Pitx1.118G22-164F21.Combined.fa", forceDNAtolower = FALSE) # returns list
# table(z1[[1]])
     # # A      C      G      T 
 # # 99205  89429  88981 100189
# # rename and write to new fasta file
# names(z1)[1] <- "chrVIIpitx1new"
# write.dna(z1, "chrVIIpitx1new.fa", format="fasta", colsep = "") # beware of blanks between columns

# # Combine with whole genome
# y1 <- read.dna("gasAcu1.fa", format="fasta", as.character=TRUE, as.matrix=FALSE)
# object.size(y1)
# # 3,706,843,648 bytes # wow! Much larger than Biostrings
# # To upper case
# y1 <- lapply(y1, toupper)
# names(y1)
 # # [1] "chrI"     "chrII"    "chrIII"   "chrIV"    "chrIX"    "chrUn"    "chrV"     "chrVI"   
 # # [9] "chrVII"   "chrVIII"  "chrX"     "chrXI"    "chrXII"   "chrXIII"  "chrXIV"   "chrXIX"  
# # [17] "chrXV"    "chrXVI"   "chrXVII"  "chrXVIII" "chrXX"    "chrXXI"   "chrM"    
# z1 <- c(y1, z1)
# # reshape to put chrVIIpitx1new after chrVII
# names(z1)[c(1:9, 24, 10:23)]
# z1 <- z1[c(1:9, 24, 10:23)]
# # Write it to a new file - contains gasAcu1 with pitx1new BAC. 
# # Didn't work on my Mac! Use Biostrings
# # write.dna(z1, "gasAcu1pitx1new.fa", format="fasta", colsep = "") # beware of blanks between columns

library(Biostrings)
# This is much faster than read.fasta from seqinr

z1 <- readDNAStringSet("chrVIIpitx1new.fa", "fasta") # has class "DNAStringSet"
names(z1)
# [1] "chrVIIpitx1new"
alphabetFrequency(z1)
         # A     C     G      T M R W S Y K V H D B N - + .
# [1,] 99205 89429 88981 100189 0 0 0 0 0 0 0 0 0 0 0 0 0 0

y1 <- readDNAStringSet("gasAcu1.fa", "fasta")
object.size(y1)
# 463,359,952 bytes

# combine them ( c() works!)
z <- c(y1,z1)

names(z)
 # [1] "chrI"           "chrII"          "chrIII"         "chrIV"          "chrIX"         
 # [6] "chrUn"          "chrV"           "chrVI"          "chrVII"         "chrVIII"       
# [11] "chrX"           "chrXI"          "chrXII"         "chrXIII"        "chrXIV"        
# [16] "chrXIX"         "chrXV"          "chrXVI"         "chrXVII"        "chrXVIII"      
# [21] "chrXX"          "chrXXI"         "chrM"           "chrVIIpitx1new"

# reshape to put chrVIIpitx1new after chrVII
z <- z[c(1:9,24,10:23)]

# Write it to a new file - contains gasAcu1 plus pitx1 BAC
writeXStringSet(z, file="gasAcu1pitx1new.fa") # very fast


# ---------------------------------------------------------
library(Biostrings)
# Read fasta file as a set of DNAString's
# This is much faster than read.fasta from seqinr
# This also converts all letters to upper case

z <- readDNAStringSet("GU130435.1.fa", "fasta") # has class "DNAStringSet"

z
  # A DNAStringSet instance of length 1
     # width seq                                                                      names               
# [1] 377852 GAATTCTTTCCTTCCTCGTTCAGGGAATCAACCTG...ATTTTCAAAAGGGTTTTAGAATCGTACTGAATTC chrVIIpitx1

# rename
names(z)[1] <- "chrVIIpitx1"

alphabetFrequency(z)
          # A     C     G     T M R W S Y K V H D B N - +
# [1,] 100189 88982 89440 99216 4 5 7 2 3 4 0 0 0 0 0 0 0

# Result is in an XStringViews format and hard to access contents
matchPattern("M", as.character(z))
     # start    end width
# [1] 166492 166492     1 [M]
# [2] 166495 166495     1 [M]
# [3] 167072 167072     1 [M]
# [4] 168244 168244     1 [M]

z1 <- as.character(z) # this causes it to lose the name attribute

# Convert "M" to "N"
z2 <- gsub("M","N",z1)
z2 <- gsub("R","N",z2)
z2 <- gsub("W","N",z2)
z2 <- gsub("S","N",z2)
z2 <- gsub("Y","N",z2)
z2 <- gsub("K","N",z2)

z2 <- DNAStringSet(z2) # This is better than DNAString() because behaves like a list and can name elements

z <- z2

alphabetFrequency(z)
          # A     C     G     T M R W S Y K V H D B  N - +
# [1,] 100189 88982 89440 99216 0 0 0 0 0 0 0 0 0 0 25 0 0
    
# Combine with whole genome
y <- readDNAStringSet("gasAcu1.fa", "fasta")

z <- c(y,z)

names(z)
 [1] "chrI"        "chrII"       "chrIII"      "chrIV"       "chrIX"      
 [6] "chrUn"       "chrV"        "chrVI"       "chrVII"      "chrVIII"    
[11] "chrX"        "chrXI"       "chrXII"      "chrXIII"     "chrXIV"     
[16] "chrXIX"      "chrXV"       "chrXVI"      "chrXVII"     "chrXVIII"   
[21] "chrXX"       "chrXXI"      "chrM"        "chrVIIpitx1"

# reshape to put chrVIIpitx1 after chrVII
z <- z[c(1:9,24,10:23)]

# Write it to a new file - contains gasAcu1 plus pitx1 BAC
writeXStringSet(z, file="gasAcu1pitx1.fa")

# Write after removing chrM
z1<- z[-c(which(names(z)=="chrM"))]
names(z1)
writeXStringSet(z, file="gasAcuWpitx1noM.fa")

# Write individual chromosomes to separate files
for(i in names(z)){
	# i <- "chrM"
	writeXStringSet(z[i], file=paste(i,".fa", sep=""))
	}
