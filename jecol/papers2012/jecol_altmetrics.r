# Journal of Ecology top 2012 papers - Altmetrics
dois <- as.character(read.csv("~/github/sac/work/jecol/papers2012/jecol_dois.csv")[,1])
dois

# trim stuff from DOIs
library(stringr)
dois <- str_trim(str_replace_all(dois, "DOI: ", ""), "both")

# Search for altmetrics
library(rImpactStory); library(rAltmetric); library(plyr)

out <- llply(dois, function(x) altmetrics(doi=x), .progress="text")
out <- out[!sapply(out, is.null)] # remove NULLs
outdf <- ldply(out, altmetric_data)
head(outdf)

# Output results to csv file
write.csv(outdf, file = "/Users/scottmac2/Mac/JEcol_blog/DOIs/altmetrics_data.csv")


# Plot results