# .Rprofile on hermes
git <- function(githubfile){
	library(RCurl, quietly = TRUE)
	script <- getURL( paste("https://raw.githubusercontent.com/dschluter/genomeScripts/master/", githubfile, sep=""), 
				.opts = list(ssl.verifypeer = FALSE) )
	eval(parse(text = script), envir = .GlobalEnv)
	}
.First <- function(){
	git("genome.r")
	
	git("misc.r")
	
	}
