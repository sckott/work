# Journal of Ecology top 2012 papers - Altmetrics
source("/Users/scottmac2/Dropbox/CANPOLIN_treeshape_ms/data/simulations/balance_cutoff_values/theme_myblank.r")
dois <- as.character(read.csv("~/github/sac/work/jecol/papers2012/jecol_dois.csv")[,1])
dois

# trim stuff from DOIs
library(stringr)
dois <- str_trim(str_replace_all(dois, "DOI: ", "doi/"), "both")

# Search for altmetrics
# install("/Users/scottmac2/github/ropensci/rAltmetric")
library(rAltmetric); library(plyr)

## Altmetric.com
out <- llply(dois, altmetrics, .progress="text")
out <- compact(out) # remove NULLs
outdf <- ldply(out, altmetric_data)
str(outdf)

# Output results to csv file
write.csv(outdf, file = "~/github/sac/work/jecol/papers2012/altmetrics_data.csv")

# Plot results
library(reshape2); library(ggplot2); library(stringr); library(reshape)
outdf <- read.csv("~/github/sac/work/jecol/papers2012/altmetrics_data.csv")
outdf_limited <- outdf[,c(2,11:13,15,22,23,25,47)]
names(outdf_limited)[2:9] <- c("g+","FB","Tw","Feeds","Alt","M","CiteUL","WoS")
outdf_limited$record <- 1:nrow(outdf_limited)
outdf_limited <- sort_df(outdf_limited, "Alt")
outdf_limited$Alt2 <- outdf_limited$Alt
str(outdf_limited)

outdf_limited_m <- melt(outdf_limited, id.vars=c(1,10,11))
outdf_limited_m$value[is.na(outdf_limited_m$value)] <- 0
outdf_limited_m$record <- as.factor(outdf_limited_m$record)
outdf_limited_m$variable <- factor(outdf_limited_m$variable, levels=c("Alt","M","Tw","WoS","FB","CiteUL","g+","Feeds"))
str(outdf_limited_m)
theme_myblank <- function(...){
	ggplot2::theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
				panel.background = element_blank(),
		# 		plot.background = element_blank(),
		axis.title.x = element_text(size=18),
		axis.title.y = element_blank(),
		axis.text.x = element_text(size=18),
		axis.text.y = element_blank(),
		# 		axis.line = element_blank(),
				axis.ticks.y = element_blank(),
		...
	)
}

# Heat Map
# png("~/github/sac/work/jecol/papers2012/toppapers2012_2.png", width=7, height=5, units="in", res=150)
# jpeg("~/github/sac/work/jecol/papers2012/toppapers2012_2.jpeg", width=7, height=5, units="in", res=150)
ggplot(outdf_limited_m, aes(variable, reorder(record, Alt2))) + 
	#   stat_bin(aes(fill=value), geom="tile", binwidth=5, position="identity") +
	geom_tile(aes(fill = log10(value+1)), colour = "white", binwidth=3) +
	scale_fill_gradient("Score", low = "white", high = "steelblue", breaks=seq(0,2,0.5)) +
	theme_bw(base_size = 16) + 
	labs(x = "", y = "") + 
	scale_x_discrete(expand = c(0, 0)) +
	scale_y_discrete(expand = c(0, 0)) +
	theme(
		axis.ticks = element_blank(), 
		axis.text.x = element_text(size = 12, hjust = 0.6),
# 		axis.text.y = element_text(size = 8, hjust = 0, colour = "grey50"),
		axis.text.y = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank()
)
ggsave("~/github/sac/work/jecol/papers2012/toppapers2012_3.png")
# dev.off()