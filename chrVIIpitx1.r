# Methods to read and write DNA sequences in R

setwd("/Users/schluter/Documents/Research/genomics/reference genomes")

# ---------------------------------------------------------
# Read BAC pitx1 sequence from GenBank using accession number

library(ape)

z <- read.GenBank("GU130435.1") # accession number for the Pitx1 BAC from Salmon River
length(z[[1]])
#[1] 377852
names(z)
[1] "GU130435.1"
names(z) <- "chrVIIpitx1"
write.dna(z, "chrVIIpitx1.fa", format="fasta", colsep = "") # beware of blanks between columns

z <- read.dna("chrVIIpitx1.fa", format="fasta", as.character=TRUE, as.matrix=FALSE) # reads as a vector
z[[1]][1:10]
 [1] "g" "a" "a" "t" "t" "c" "t" "t" "t" "c"

# or
z <- read.FASTA("chrVIIpitx1.fa") # always returns a list of class "DNAbin".
 
# ---------------------------------------------------------
library(Biostrings)
# Read fasta file as a set of DNAString's
# This is much faster than read.fasta from seqinr
# This also converts all letters to upper case

z <- read.DNAStringSet("chrVIIpitx1.fa", "fasta") # has class "DNAStringSet"
# z <- readFASTA("chrVIIpitx1.fa") deprecated

z
  # A DNAStringSet instance of length 1
     # width seq                                                                      names               
# [1] 377852 GAATTCTTTCCTTCCTCGTTCAGGGAATCAACCTG...ATTTTCAAAAGGGTTTTAGAATCGTACTGAATTC chrVIIpitx1

names(z)
# [1] "chrVIIpitx1"

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
names(z2)<-"chrVIIpitx1"

z <- z2

alphabetFrequency(z)
          # A     C     G     T M R W S Y K V H D B  N - +
# [1,] 100189 88982 89440 99216 0 0 0 0 0 0 0 0 0 0 25 0 0
    
# Combine with whole genome
y <- read.DNAStringSet("gasAcu1.fa", "fasta")

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
write.XStringSet(z, file="gasAcu1pitx1.fa")

# Write after removing chrM
z1<- z[-c(which(names(z)=="chrM"))]
names(z1)
write.XStringSet(z, file="gasAcuWpitx1noM.fa")

# Write individual chromosomes to separate files - NOTE: this overwrites the original chrVIIpitx1.fa
for(i in names(z)){
	# i <- "chrM"
	write.XStringSet(z[i], file=paste(i,".fa", sep=""))
	}
