# get new Journal of Ecology paper metadata
library(rplos); library(RMendeley)

######
getdois <- function (input) {
	dois <- sapply(input, function(x) x$doi) # parse out DOIs
	dois <- if(any(sapply(dois, is.null)) == TRUE){
		dois[-(which(sapply(dois,is.null),arr.ind=TRUE))] } else
		{dois}
	dois[sapply(dois,str_length) > 0] # limit to DOIs > zero character length
}
je <- msearch('published_in:\'Journal of Ecology\'', numItems=20)
mc <- mendeley_auth()
je <- msearch('published_in:\'Journal of Ecology\'', numItems=20)
je
je <- msearch('published_in:\'Journal of Ecology\' year:2012', numItems=20)
je
je <- msearch('year:2012', numItems=20)
je
je <- msearch('published_in:\'Journal of Ecology\' year:\'2012\'', numItems=20)
'published_in:\'Journal of Ecology\' year:\'2012\''
msearch
query='published_in:\'Journal of Ecology\' year:2012'
query <- gsub(":", "%3A", query)  # use html symbol for colons
query
query <- gsub(" ", "%20", query)  # use html symbol for spaces
query
query='published_in:\'Journal of Ecology\'year:2012'
query
query <- gsub(":", "%3A", query)  # use html symbol for colons
query <- gsub(" ", "%20", query)  # use html symbol for spaces
query <- gsub("\"", "%22", query)  # use html symbol for quotes
query

# 

library(XML)
setwd("/Users/ScottMac/github/SChamberlain/work/jecol")
tt <- htmlParse("currentissue.html")
ttt <- xpathSApply(tt, "//a")
tttt <- sapply( ttt, function(x) xmlAttrs(x)["href"] )
urls <- tttt[grep("http://www.journalofecology.org/details/journalArticle/", tttt)]

urls[[1]]
content(GET(urls[[1]]))