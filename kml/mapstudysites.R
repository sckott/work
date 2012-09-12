##  map of study sites in pollinator-seed predator study
require(XML); require(plyr)
temp <- xmlToList("helianthussites_gmaps.kml")[25:46]
getit <- function(x) {
  a <- x[[1]]
  dd <- as.numeric(str_split(x[[length(x)]], ",")[[1]][1:2])
  c(a, dd)  
}
ldply(temp, getit)