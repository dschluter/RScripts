#!/usr/bin/Rscript

# Rscript ~/Documents/RScripts/testArgs.R chrname=chrXXI project=Benlim groups=paxl,paxb,pril,prilb drop=Sarafish

args <- commandArgs(TRUE)
# print(args)

# args<-unlist(strsplit("chrname=chrXXI project=Benlim groups=paxl,paxb,pril,prilb drop=Sarafish", split=" "))

chrname <- NULL
project <- NULL
groups <- NULL
drop <- NULL
x <- read.table(text = args, sep = "=", colClasses = "character")
# x
       # V1                   V2
# 1 chrname               chrXXI
# 2 project               Benlim
# 3  groups paxl,paxb,pril,prilb
# 4    drop             Sarafish

for(i in 1:nrow(x)){
	assign(x[i,1], x[i,2])
	}

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")
if(is.null(groups)) stop("Provide groups= in arguments (grounames separated by commas, no spaces)")

groups <- unlist(strsplit(groups, split = ","))

# ----
test <- function(...){
	x <- c(...)
	for(i in 1:length(x)) print(x[i])
	}
test("chrname","project","groups","drop")