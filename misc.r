curl<-function(Rfile){
	system(paste("curl -o " , Rfile, " https://raw.githubusercontent.com/dschluter/RScripts/master/", Rfile, sep = ""))
	}

showq<-function(){system("showq -u schluter")}

month2integer <- function(x){
	z1 <- casefold(substr(x, 1, 3)) # round to 3 digits
	months <- c("january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december")
	months.abb <- casefold(substr(months, 1, 3))
	z2 <- recode(z1, months.abb, 1:12)
	}

nice <- function(x){
	if(is.data.frame(x)) write.table(x, sep = " ", row.names = FALSE, quote = FALSE)
	else if(is.list(x)) for(i in 1:length(x)){
		cat("$",names(x)[i],"\n", sep="")
		print(x[[i]])
		}
	else if(is.vector(x)) sapply(1:length(x), function(i){ cat(x[i], "\n") })
	else stop("Not a data frame, list, or vector")
	return(invisible())
	}

nicelist <- function(z){
	# Function to print out list output (e.g., bio300$anova.oneway) for easy copy and paste to .r file
	for(i in 1:length(z)){
		cat("$",names(z)[i],"\n", sep="")
		print(z[[i]])
		}
	}
	
nicevector <- function(x){
	if(!is.vector(x)) stop("not a vector\n")
	z <- sapply(1:length(x), function(i){
		cat(x[i], "\n")
		})
	return(invisible())
	}
	
niceframe <- function(x){
	if(!is.data.frame(x)) stop("not a data frame\n")
	write.table(x, sep = " ", row.names = FALSE, quote = FALSE)
	return(invisible())
	}

names.quotes <- function(...){
	# use as names.quotes(c(a,b,c))
	#names.quotes <- function(x = c(a,b,c)){  # alternative
	x <- deparse(substitute(...), width.cutoff=500)
	x <- gsub("c(", "", x, fixed=TRUE)
	x <- gsub(")", "", x, fixed=TRUE)
	x <- unlist(strsplit(x, split="[ ,]+"))
	return(x)
	}
addquotes <- names.quotes # synonym

integer2binary <- function(x){
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base2 <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base2[, ndigits-i+1] <- (x %% 2)
                x <- (x %/% 2)
        }
        if(N ==1) Base2[1, ] else Base2
} 
