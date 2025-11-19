
#instrall package 
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("Biostrings")
#install.packages("reshape2") 

library(Biostrings)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

inputfile = read.table(args[1])

#HW aggrument|inputfile <- read.table("Stat_raw/barcode49/outfilemaping/mapwithraf.tsv")
colnames(inputfile) <- c("qseqid","sseqid", "pident" ,"length","qlen")


#boxplot
#ggplot(inputfile, aes(x=reorder(qseqid,pident), y=pident,fill = qlen)) + 
 #geom_boxplot() +
  #xlab("contig") + coord_flip()


#Fasta Read 
fa <- readDNAStringSet("Reference/sequences.fasta")

##get name file input 
fa_n <- names(fa)

out <- data.frame()

for (i in 1:length(fa_n)) {
  faa <-fa_n[i]
  faa2 <- str_split(faa,"\\|")[[1]]
  df<-as.data.frame(t(faa2))
  
  out <- rbind(out,df)

}
#colnames() showname
#chang colname 
colnames(out) <- c("sseqid", "fullname" ,"contry","seqment")

out$sseqid <- str_remove(out$sseqid," ")


###merge zone
#inputfile
#out 
mdf <- dplyr::inner_join(inputfile,out,by = "sseqid")
sum_mary <- mdf |> group_by(seqment,qseqid) |> summarise(me=mean(pident),Sd=sd(pident),cnt=n(),lc=max(qlen))  
write.csv(x = sum_mary,file = args[2])


#######MapNameBlast
Reffile <- readDNAStringSet("Reference/sequences.fasta")
Reffile_name <- names(Reffile)
outfile <- data.frame()

for (i in 1:length(Reffile_name)) {
  Reffile_alinname <-Reffile_name[i]
  Reffile_alinname2 <- str_split(Reffile_alinname,"\\|")[[1]]
  dff<-as.data.frame(t(Reffile_alinname2))
  
  outfile <- rbind(outfile,dff)
}
colnames(outfile) <- c("sseqid", "fullname" ,"contry","seqment")
outfile$sseqid <- str_remove(outfile$sseqid," ")

inputfileblast <- read.table("Stat_raw/barcode49/outfileblastref/blastwithraf.tsv")
colnames(inputfileblast) <- c("qseqid","sseqid", "pident" ,"length","qlen")

mdf <- dplyr::inner_join(inputfileblast,outfile,by = "sseqid")
sum_mary2 <- mdf |> group_by(seqment,sseqid) |> summarise(me=mean(pident),Sd=sd(pident),cnt=n(),lc=max(qlen))  
write.csv(x = sum_mary,file = "Stat_raw/barcode49/outfileblastref/Sblastwithraf_final")



