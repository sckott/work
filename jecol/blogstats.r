library(XML); library(stringr); library(ggplot2); library(RCurl)
temp <- htmlParse(getURL("https://raw.github.com/SChamberlain/work/master/jecol/sitestats.html"))
temp2 <- xpathSApply(temp, "//td", xmlValue)
temp3 <- temp2[sapply(temp2, function(x) nchar(x), USE.NAMES=F) > 0]

titles <- temp3[seq(1, length(temp3), 2)]
views <- temp3[seq(2, length(temp3), 2)]
dat <- data.frame(titles=titles, views=views)
dat <- dat[!dat$titles %in% c("Home page / Archives","About","Contributors"),] # remove home page views
dat$views <- as.numeric(str_replace_all(as.character(dat$views), ",", ""))
dat <- dat[order(dat$views, decreasing=T),]

ggplot(dat, aes(reorder(titles, views), views)) +
	geom_point(size = 3) +
	coord_flip() +
	theme_bw() +
	labs(x = "", y = "Views") + 
	theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), 
				axis.text.x = element_text(size=18), axis.title.x = element_text(size=18))