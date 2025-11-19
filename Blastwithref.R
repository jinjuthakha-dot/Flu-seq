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
