is_number <- function(A){ !grepl("[^0-9]+", A) }

curl <- function(Rfile){
	system(paste("curl -o " , Rfile, " https://raw.githubusercontent.com/dschluter/RScripts/master/", Rfile, sep = ""))
	}

showq <- function(){system("showq -u schluter")}

qdel <- function(x){system( paste("qdel", x) )}

qstat <- function(x){system( paste("qstat -f", x) )}

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

logit<-function(x){
	log(x/(1. - x))
}

antilogit<-function(x){
	exp(x)/(exp(x) + 1)
}

tally<-function(x, weights=NULL){
	# Counts up cases in different categories as a prelude to log-linear analysis
	# "x" is a data frame of categorical variables ONLY (or discrete variables)
	# The function counts the number of cases in each combination of categories
	# "weights" refers to a frequency variable if the data are already partly tallied by category,
	# otherwise the weight is assumed to be 1 for each row.
	x<-na.omit(as.data.frame(x))
	if(is.null(weights)){
		weights<-rep(1,nrow(x))
		}

	y<-as.list(x)
	# adjust the levels of y in case some are missing
	for(i in 1:length(y)){
		y[[i]]<-as.factor(as.character(y[[i]]))
		}

	z<-tapply(weights,y,sum)
	z<-as.vector(z)
	z[is.na(z)]<-0
	
	vars<-list(length(y))
	n<-length(z)
	levs<-unlist(lapply(y,function(x){length(levels(x))}))
	levs<-c(1,levs[-length(levs)])
	for(i in 1:length(y)){
		k<-prod(levs[1:i])
		vars[[i]]<-rep(levels(y[[i]]),each=k,length.out=n)
		}
	do.call("cbind.data.frame",vars)
	z<-cbind.data.frame(vars,z)
	names(z)<-c(names(y),"frequency")
	return(z)
}

matchcol<-function(matchto, matchfrom, fromvar){
	# Gets the position in "matchfrom" of the elements in the single vector variable "matchto"
	# Then "fromvar" is the variable/array in "matchfrom" corresponding to those positions
	# At this point matchto and fromvar will have the same length
	# normally after running this macro you will want to one of the following:
	#  1) bind "matchto" to result
	#  2) replace the variable "matchfrom" in result with the elements of "matchto" to eliminate NA's
	position <- match(matchto, matchfrom)
	if(is.null(ncol(fromvar))){
		result <- fromvar[position]
		}
	else{
		result <- fromvar[position,]
		}
	return(result)
}

recode <- function(..., ret = c("numeric", "factor"), none = if (ret == "numeric") 0 else "none", na){
	# stolen from the Hmisc library
    ret <- match.arg(ret)
    w <- list(...)
    if (!is.logical(w[[1]]) && length(w) == 3) {
        z <- w[[3]][match(w[[1]], w[[2]])]
        if (!missing(none)) 
            z[if (is.numeric(none)) 
                is.na(z)
            else z == ""] <- none
        return(z)
    }
    nam <- names(w)
    if (missing(ret)) 
        ret <- if (all.is.numeric(nam)) 
            "numeric"
        else "factor"
    result <- rep(none, length(w[[1]]))
    for (i in 1:length(w)) result[w[[i]]] <- if (ret == "numeric") 
        numnam[i]
    else nam[i]
    if (ret == "factor") 
        result <- as.factor(result)
    if (!missing(na)) 
        result[is.na(na)] <- NA
    result
}

