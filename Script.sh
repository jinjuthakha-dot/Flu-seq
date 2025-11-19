###Day-1 20oct2025 Consept
###Day-2 31oct2025 Install package,step 1.Summary statisc 2.Trim Primer and cut lenght

Assembly() { 
  #DownlodeReference from Database NCBI viralnlm influenA+B Complete genmoe & human host
  #build database 1 time > makeblastdb -in Reference/sequences.fasta -dbtype nucl 
  #input is QC fastq 
  #output is Scaffolds
  
  type=$1
  if [ $type == "Denovo"];
  then 
    echo "3 Denovo Assembly"
    spades.py -s $outfile/trimprimer.fq.gz --only-assembler --careful -o $outfile/Denoass
    
    #stat after assembly
    seqkit stat $outfile/Denoass/scaffolds.fasta |tail -n +2 >> $outfile/statba49.txt
    
    sed -i '' 's/,//g' $outfile/statba49.txt
  
    map="$outfile/outfilemaping"
    mkdir -p $map
    blastn -query $outfile/Denoass/scaffolds.fasta -db Reference/sequences.fasta -outfmt '6 qseqid sseqid pident length qlen' > $map/mapwithraf.tsv  
    
    #Annnotation 
    Rscript ./Annotate_contig.R $map/mapwithraf.tsv $map/mapwithraf_final.tsv
  else
    echo "1243"
  fi
} 


#code for run > ./Script.sh RawSample/b49.fq.gz 300
echo "1 Summary statisc"

inputfile=$1
min_lenght=$2
outfile="Stat_raw/barcode49"
mkdir -p $outfile

#stat raw data
rawb49="RawSample/b49.fq.gz"
seqkit stat $rawb49 > $outfile/statba49.txt

echo "2 Trim Primer (Real Read length) & minimun lenght &quality read form nanopore "
f="AGCAAAAGCAGG"
r="AGTAGAAACAAGG"
cpr="CCTTGTTTCTACT"
cutadapt -g $f...$cpr -m $min_lenght --revcomp -e 0.2 -q 10 -o $outfile/trimprimer.fq.gz $inputfile 

#stat after QC
seqkit stat $outfile/trimprimer.fq.gz|tail -n +2 >> $outfile/statba49.txt
sed -i '' 's/,//g' $outfile/statba49.txt 

#Blast
seqkit fq2fa $outfile/trimprimer.fq.gz | gzip -9 > $outfile/trimprimer.fasta.gz 



blastref="$outfile/outfileblastref"
mkdir -p $blastref
gunzip -c $outfile/trimprimer.fasta.gz |\
  blastn -query - -db Reference/sequences.fasta -outfmt '6 qseqid sseqid pident length qlen' > $blastref/blastwithraf.tsv 
  



#Select the est ref. 
#Mapping 
#samtools consensus sample.bam -f thebestref.fa -m simple -d 30 > $blastref/consensus.fa







#echo "3 QC (Filter Read Length >Qulity Read)"
#echo "4 Reference Assembly"
#echo "5 Qulity Consensus (Depth >30x ,Coverage)"
#echo "6 Annotation (Genotype FluA or B,Stain H1N1 or H3N2,Clade,Sub clade,Antigenicity Prediction,Evolotion)"


