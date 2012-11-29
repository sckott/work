##### 
### Vamosi-Heard project
## SCOTT !!! Always edit the version of this file in "~/github/SChamberlain/work/vamosi_heard/"
## Just for Scott: sync code files on git repo to dropbox whenever necessary
# system("cp -r ~/github/SChamberlain/work/vamosi_heard/ ~/Dropbox/Vamosi_Heard/code/")

#####
# load libraries
# install_github('taxize_', 'ropensci') # if not installed already
library(gdata); library(plyr); library(ggplot2); library(stringr); library(taxize)
require(phangorn); require(ape)

#####
# Read in data, clean
# setwd("~/Dropbox/Vamosi_Heard/data")
# data_nb <- read.xls("new brunswick fish list.xlsx")[-c(50:51),]
# data_s <- read.xls("sudbury fish list.xlsx", sheet="taxa")
# 
# # clean data 
# data_nb <- colwise(str_trim)(data_nb)
# data_nb$spname <- paste(data_nb$genus, data_nb$species)
# 
# data_s <- colwise(str_trim)(data_s)
# data_s$spname <- paste(data_s$genus, data_s$species)
# 
# # make unique species names vector, combining from two data sets
# splist <- sort(unique(c(data_nb$spname, data_s$spname)))

# New species list, after name corrections
setwd("~/Dropbox/Vamosi_Heard/data")
splist <- as.character(read.csv("splist.csv")[,1])

#####
# Search for and download sequences
# 	genes used by NearEtal2012PNAS: Glyt, myh6, plagl2, rag1, zic1, Ptr, sreb2, SH3PX3, tbr1
library(multicore); library(doMC) # which parallel pkg you use depends on your oper. system
df_towrite <- data.frame(taxon=NA, gene_desc=NA, gi_no=NA, acc_no=NA, length=NA, sequence=NA, spused=NA)

# CO1
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_coi.txt", row.names=F)
do <- llply(splist, get_seqs, gene = c("coi", "co1"), 
	seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_coi.txt")
fish_coi_seqs_coi_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_coi.txt", header=T)[-1,]

# RAG1
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_rag1.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "rag1", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_rag1.txt")
fish_coi_seqs_rag1_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_rag1.txt", header=T)[-1,]

# Glyt
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_Glyt.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "Glyt", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_Glyt.txt")
fish_coi_seqs_Glyt_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_Glyt.txt", header=T)[-1,]
fish_coi_seqs_Glyt_df

# myh6
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_myh6.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "myh6", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_myh6.txt")
fish_coi_seqs_myh6_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_myh6.txt", header=T)[-1,]

# plagl2
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_plagl2.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "plagl2", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_plagl2.txt")
fish_coi_seqs_plagl2_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_plagl2.txt", header=T)[-1,]

# Ptr
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_Ptr.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "Ptr", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_Ptr.txt")
fish_coi_seqs_Ptr_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_Ptr.txt", header=T)[-1,]

# SH3PX3
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_SH3PX3.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "SH3PX3", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_SH3PX3.txt")
fish_coi_seqs_SH3PX3_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_SH3PX3.txt", header=T)[-1,]

# sreb2
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_sreb2.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "sreb2", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_sreb2.txt")
fish_coi_seqs_sreb2_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_sreb2.txt", header=T)[-1,]

# tbr1
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_tbr1.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "tbr1", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_tbr1.txt")
fish_coi_seqs_tbr1_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_tbr1.txt", header=T)[-1,]

# zic1
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_zic1.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "zic1", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_zic1.txt")
fish_coi_seqs_zic1_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_zic1.txt", header=T)[-1,]

# cytb
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_cytb.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "cytb", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_cytb.txt")
fish_coi_seqs_cytb_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_cytb.txt", header=T)[-1,]

# NADH5
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_NADH5_new.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "NADH5|ND5", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_NADH5_new.txt")
fish_coi_seqs_NADH5_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_NADH5.txt", header=T)[-1,]

# twelveS
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_twelveS.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "12S ribosomal|12S large subunit", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_twelveS.txt")
fish_coi_seqs_twelveS_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_twelveS.txt", header=T)[-1,]

# sixteenS
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_sixteenS.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "16S ribosomal|16S large subunit", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_sixteenS.txt")
fish_coi_seqs_sixteenS_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_sixteenS.txt", header=T)[-1,]

# tweightS
registerDoMC(cores=4)
write.table(df_towrite, "~/Dropbox/Vamosi_Heard/data/fishseqsout_tweightS.txt", row.names=F)
dp <- llply(splist, get_seqs, gene = "28S ribosomal|28s large subunit", 
						seqrange = "1:3000", getrelated=T, writetodf=T, .parallel=T, filetowriteto="fishseqsout_tweightS.txt")
fish_coi_seqs_tweightS_df <- read.table("~/Dropbox/Vamosi_Heard/data/fishseqsout_tweightS.txt", header=T)[-1,]

#####
# Make summary table of genes available for each species
library(reshape); library(stringr); library(plyr)
setwd("~/Dropbox/Vamosi_Heard/data/")
prepdf <- function(x, column) {
	x_ <- read.table(x, header=T)[-1,c("taxon",column)]
	x__ <- na.omit(x_)
	if(identical(str_detect(x__$taxon, "Chrosomus eos$"), logical(0)))
		{NULL} else { x__ <- x__[!str_detect(x__$taxon, "Chrosomus eos$"),] }
	if(identical(str_detect(x__$taxon, "Acipenser oxyrhynchus$"), logical(0)))
		{NULL} else { x__ <- x__[!str_detect(x__$taxon, "Acipenser oxyrhynchus$"),] }
	if(identical(str_detect(x__$taxon, "Catostomus commersoni$"), logical(0)))
		{NULL} else { x__ <- x__[!str_detect(x__$taxon, "Catostomus commersoni$"),] }
	if(identical(str_detect(x__$taxon, "Micropterus dolomieui$"), logical(0)))
		{NULL} else { x__ <- x__[!str_detect(x__$taxon, "Micropterus dolomieui$"),] }
	dols <- x__[str_detect(x__$taxon, "Micropterus dolomieu$"),]
	if(length(dols)>1){
			temp <- dols[sample(1:nrow(dols), 1),]
			x__ <- x__[!str_detect(x__$taxon, "Micropterus dolomieu$"),]
			x__ <- rbind(x__, temp)
	} else {NULL}
	names(x__)[2] <- strsplit(strsplit(x, ".txt")[[1]], "_")[[1]][[2]]
	x__
}
filenames <- c("fishseqsout_zic1.txt","fishseqsout_tbr1.txt","fishseqsout_sreb2.txt",
							"fishseqsout_SH3PX3.txt","fishseqsout_Ptr.txt","fishseqsout_coi.txt",
							"fishseqsout_rag1.txt","fishseqsout_plagl2.txt","fishseqsout_myh6.txt",
							"fishseqsout_Glyt.txt","fishseqsout_tweightS.txt","fishseqsout_sixteenS.txt",
							 "fishseqsout_twelveS.txt","fishseqsout_NADH5.txt","fishseqsout_cytb.txt")
dfs <- lapply(filenames, prepdf, column="length")
summary_length <- merge_recurse(dfs, by="taxon")
summary_length <- sort_df(summary_length, "taxon")
write.csv(summary_length, "summarydf.csv", row.names=F)

dfs <- lapply(filenames, prepdf, column="spused")
summary_names <- merge_recurse(dfs, by="taxon")
summary_names <- sort_df(summary_names, "taxon")
write.csv(summary_names, "summary_names.csv", row.names=F)

#####
# Read in sequence files
setwd("/Users/ScottMac/Dropbox/Vamosi_Heard/data")

# fxn to read in file, removing first row of NA's, other NA rows, and duplicates
prepfasta <- function(x){
	dir <- "/Users/ScottMac/Dropbox/Vamosi_Heard/data/"
	temp <- read.table(paste0(dir, "fishseqsout_", x, ".txt"), header=T)[-1,] # header TRUE, and remove first NA row
	temp2 <- temp[!sapply(split(temp, as.numeric(row.names(temp))), function(x) all(is.na(x[,-1])) ), ]
	temp3 <- temp2[,c("spused","sequence")]
	doit <- function(x){ if(nrow(x)>1) { x[1,] } else { x } }
	out <- ldply(split(temp3, temp3$spused), doit)[,-1]
	seqs <- as.character(out$sequence)
	names(seqs) <- out$spused
	seqs
}
# coi <- prepfasta("coi")

#####
# Multiple sequence alignment using CLUSTAL
# Download Clustalw from http://www.clustal.org/download/current/
#   and I just put the clustal folder in my /Applications folder on my Mac
# Uses prepfasta function above
setwd("/Applications/clustalw-2.1-macosx/") # set to your directory that has clustalw
readalign <- function(genename) {
	file_ <- prepfasta(genename)	
	write_fasta(file_, paste0(genename, ".fas")) # from package sacbox at https://github.com/schamberlain/sacbox
	system(paste0('"./clustalw2" ', paste0(genename, ".fas"))) # run clustal multiple alignment
	read.dna(paste0(genename, ".aln"), format="clustal") # read aligned sequences
}
# (out <- readalign("zic1"))
genenames <- c("zic1","tbr1","sreb2","SH3PX3","Ptr","coi","rag1","plagl2","myh6",
							 "Glyt","tweightS","sixteenS","twelveS","NADH5","cytb")
alignments <- llply(genenames, readalign, .progress="text")

# Or read in if alignments already done
library(phangorn)
# setwd("/Applications/clustalw-2.1-macosx/")
# tt <- c("cytb.aln", "NADH5.aln", "twelveS.aln", "sixteenS.aln", "tweightS.aln", "Glyt.aln",
# 	"myh6.aln", "plagl2.aln", "rag1.aln", "coi.aln", "Ptr.aln", "SH3PX3.aln", 
# 	"sreb2.aln", "tbr1.aln", "zic1.aln")
# alignments <- llply(tt, function(x) read.dna(x, format="clustal"))
setwd("~/Dropbox/Vamosi_Heard/data/alignments/")
# tt <- c("myh6.aln.phy", "Glyt.aln.phy", "plagl2.aln.phy", "NADH5.aln.phy", "rag1.aln.phy", "coi.aln.phy")
tt <- c("sreb2.aln.phy", "SH3PX3.aln.phy", "Ptr.aln.phy", "tbr1.aln.phy", "zic1.aln.phy")
alignments <- llply(tt, function(x) read.dna(x))

#####
# Fit various evol models to sequence alignment, pick best model for each gene
fitevolmods <- function(x) {
	# convert ape format to phangorn format
	phang_ <- phyDat(x, return.index = T) 
	# fit different evolutionary models to the dataset using fxn modelTest
	njtree <- NJ(dist.logDet(phang_)) # create neighbor joining tree for use in modelTest
	evolmodels <- modelTest(phang_, tree = njtree) # comparison of evolutionary models
	evolmodels[order(evolmodels$AIC), ] # sorted by minimum AIC
}
fitevolmods_safe <- plyr::failwith(NULL, fitevolmods, quiet=T)
evolmods_results <- llply(alignments, fitevolmods_safe, .progress="text")
# names(evolmods_results) <- 
# 	c("cytb", "NADH5", "twelveS", "sixteenS", "tweightS", "Glyt","myh6", "plagl2",
# 	"rag1", "coi", "Ptr", "SH3PX3", "sreb2", "tbr1", "zic1")
# names(evolmods_results) <- c("myh6", "Glyt", "plagl2", "NADH5", "rag1", "coi")
names(evolmods_results) <- c("sreb2", "SH3PX3", "Ptr", "tbr1", "zic1")
evolmods_results_2 <- evolmods_results[!sapply(evolmods_results, is.null)]
evolmods_results_df <- ldply(evolmods_results_2, function(x) as.data.frame(x))
write.csv(evolmods_results_df, "~/Dropbox/Vamosi_Heard/data/evolmods_results_df_new_new.csv")

#####
# Concatenate sequence alignments into one alignment, noting where each partition starts/stops
# genenames <- c("zic1","tbr1","sreb2","SH3PX3","Ptr","coi","rag1","plagl2","myh6",
# 							 "Glyt","tweightS","sixteenS","twelveS","NADH5","cytb")
genenames <- c("coi","rag1","plagl2","Glyt","NADH5","myh6","zic1", "Ptr", "sreb2", "SH3PX3", "tbr1")

readaligns <- function(genename) {
	ttt <- read.table(paste0("~/Dropbox/Vamosi_Heard/data/alignments/", genename, ".aln.phy"), skip=1, fill=T)
	names(ttt) <- c("species", "seq")
	ttt
}
all <- llply(genenames, readaligns, .progress="text")

out <- merge(all[[1]], all[[2]], by="species", all.x=T)
names(out)[2:3] <- c("coi","rag1")
out <- merge(out, all[[3]], by="species", all.x=T)
names(out)[4] <- "plagl2"
out <- merge(out, all[[4]], by="species", all.x=T)
names(out)[5] <- "Glyt"
out <- merge(out, all[[5]], by="species", all.x=T)
names(out)[6] <- "NADH5"
out <- merge(out, all[[6]], by="species", all.x=T)
names(out)[7] <- "myh6"
out <- merge(out, all[[7]], by="species", all.x=T)
names(out)[8] <- "zic1"
out <- merge(out, all[[8]], by="species", all.x=T)
names(out)[9] <- "Ptr"
out <- merge(out, all[[9]], by="species", all.x=T)
names(out)[10] <- "sreb2"
out <- merge(out, all[[10]], by="species", all.x=T)
names(out)[11] <- "SH3PX3"
out <- merge(out, all[[11]], by="species", all.x=T)
names(out)[12] <- "tbr1"
str(out)
head(out)

rep_ <- apply(out[3:12], 2, function(x) max(nchar(as.character(x))))
# alllengths <- apply(out[2:7], 2, function(x) max(nchar(as.character(x))))
# write.csv(data.frame(gene = names(alllengths), length=alllengths, row.names=NULL), "~/dat.csv")
coilength <- max(nchar(as.character(out[,2])))
doit <- function(x) {
	if(is.na(x[,3])) { temp <- paste(rep("-", rep_[[names(x)[3]]]),collapse="") } else { temp <- as.character(x[,3]) }
	if(is.na(x[,4])) { temp1 <- paste(rep("-", rep_[[names(x)[4]]]),collapse="") } else { temp1 <- as.character(x[,4]) }
	if(is.na(x[,5])) { temp2 <- paste(rep("-", rep_[[names(x)[5]]]),collapse="") } else { temp2 <- as.character(x[,5]) }
	if(is.na(x[,6])) { temp3 <- paste(rep("-", rep_[[names(x)[6]]]),collapse="") } else { temp3 <- as.character(x[,6]) }
	if(is.na(x[,7])) { temp4 <- paste(rep("-", rep_[[names(x)[7]]]),collapse="") } else { temp4 <- as.character(x[,7]) }
	if(is.na(x[,8])) { temp5 <- paste(rep("-", rep_[[names(x)[8]]]),collapse="") } else { temp5 <- as.character(x[,8]) }
	if(is.na(x[,9])) { temp6 <- paste(rep("-", rep_[[names(x)[9]]]),collapse="") } else { temp6 <- as.character(x[,9]) }
	if(is.na(x[,10])) { temp7 <- paste(rep("-", rep_[[names(x)[10]]]),collapse="") } else { temp7 <- as.character(x[,10]) }
	if(is.na(x[,11])) { temp8 <- paste(rep("-", rep_[[names(x)[11]]]),collapse="") } else { temp8 <- as.character(x[,11]) }
	if(is.na(x[,12])) { temp9 <- paste(rep("-", rep_[[names(x)[12]]]),collapse="") } else { temp9 <- as.character(x[,12]) }
	itt <- paste0(x[,2], temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9)
	names(itt) <- x$species
	itt
}
outout <- lapply(split(out, out$species), doit)
outoutdf <- ldply(outout)
# sapply(outoutdf[,2], function(x) nchar(as.character(x)), USE.NAMES=F)
# write.csv(outoutdf, "final_aligned_concat.csv")
justseqs <- outoutdf[,2]
names(justseqs) <- outoutdf[,1]
write_fasta(justseqs, "~/final_aligned_concat_new.fas")


#####
# Write to file to use in MrBayes (Bayesian) and RaxML (ML) methods


#####
# Run ML tree recontsruction protocol


#####
# Read in specieslist and prune trees to each set of species
library(ape)
splist <- read.csv("~/Dropbox/Vamosi_Heard/data/splist.csv")
tree <- read.nexus(file="final_aligned_concat_tomrbayes2.nex.con.tre")
plot(tree[[1]])
nodelabels(tree[[1]]$node.label)