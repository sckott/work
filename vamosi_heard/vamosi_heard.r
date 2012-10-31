##### 
### Vamosi-Heard project
## Just for Scott: sync code files on git repo to dropbox whenever necessary
# system("cp -r ~/github/SChamberlain/work/vamosi_heard/ ~/Dropbox/Vamosi_Heard/code/")

#####
# load libraries
# install_github('taxize_', 'ropensci') # if not installed already
library(gdata); library(plyr); library(ggplot2); library(stringr); library(taxize)

#####
# Read in data, clean
setwd("~/Dropbox/Vamosi_Heard/data")
data_nb <- read.xls("new brunswick fish list.xlsx")[-c(50:51),]
data_s <- read.xls("sudbury fish list.xlsx", sheet="taxa")

# clean data 
data_nb <- colwise(str_trim)(data_nb)
data_nb$spname <- paste(data_nb$genus, data_nb$species)

data_s <- colwise(str_trim)(data_s)
data_s$spname <- paste(data_s$genus, data_s$species)

# make unique species names vector, combining from two data sets
splist <- sort(unique(c(data_nb$spname, data_s$spname)))

#####
# Search for and download sequences
# 	genes used by NearEtal2012PNAS: Glyt, myh6, plagl2, Ptr, rag1, SH3PX3, sreb2, tbr1, and zic1
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

#####
# Make summary table of genes available for each species
setwd("~/Dropbox/Vamosi_Heard/data/")
prepdf <- function(x, column) {
	x_ <- read.table(x, header=T)[-1,c("taxon",column)]
	x__ <- na.omit(x_)
	names(x__)[2] <- strsplit(strsplit(x, ".txt")[[1]], "_")[[1]][[2]]
	x__
}
filenames <- c("fishseqsout_zic1.txt","fishseqsout_tbr1.txt","fishseqsout_sreb2.txt",
							"fishseqsout_SH3PX3.txt","fishseqsout_Ptr.txt","fishseqsout_coi.txt",
							"fishseqsout_rag1.txt","fishseqsout_plagl2.txt","fishseqsout_myh6.txt",
							"fishseqsout_Glyt.txt")
dfs <- lapply(filenames, prepdf, column="length")
summary_length <- merge_recurse(dfs, by="taxon")
summary_length <- sort_df(summary_length, "taxon")
write.csv(summary_length, "summarydf.csv", row.names=F)

dfs <- lapply(filenames, prepdf, column="spused")
summary_names <- merge_recurse(dfs, by="taxon")
summary_names <- sort_df(summary_names, "taxon")
write.csv(summary_names, "summary_names.csv", row.names=F)

#####
# Read in sequence files and get data.frame for each regional set of species


#####
# Alignment using CLUSTAL


#####
# Fit various evol models to sequence alignment


#####
# 