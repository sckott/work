# Journal of Ecology top 2012 papers - Altmetrics
# dois <- read.csv("~/github/sac/work/jecol/papers2012/jecol_dois.csv")
# commented this ^ line because my path is different.
library(stringr)
library(rImpactStory)
library(rAltmetric)
library(plyr)
# read only the first column
dois <- read.csv('jecol_dois.csv')[,1]
dois_alt <- str_trim(str_replace_all(dois, "DOI: ", "doi/"), "both")

# Metrics from rAltmetric
out <- llply(dois_alt, altmetrics, .progress="text")
out <- compact(out) # remove NULLs
outdf <- ldply(out, altmetric_data)
head(outdf)
# Output results to csv file
write.csv(outdf, file = "/Users/scottmac2/Mac/JEcol_blog/DOIs/altmetrics_data.csv")


# Plot results
