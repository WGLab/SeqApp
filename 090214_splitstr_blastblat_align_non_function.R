######################Starting with BLAST##########################
###################BLAST in WINDOWS#######################

#Before function
#Need to have: BLAST location,seq input (format check),hg19 ref location,output format, output name
#Need to let user know: web access needed for remote. Otherwise download database.
#*need to check if user input is at least 36 bps, so that upstream and downstram of the breakpoints will have ~18bps to be aligned
###################remote BLAST#######################
###################Database download BLAST##################
setwd("C:/Users/lieberlab5/Dropbox/ma_sharefolder_Lieberlab/breakpt_project/v1.0/function1_blastblat/")

#blastn
#note the default task is megablast in this case
start.time1 = Sys.time()
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db nt -query C:/Users/lieberlab5/Dropbox/ma_sharefolder_Lieberlab/microRNA_project/bcl1mtc_fragile.fasta -out testt.out -outfmt 6 -remote")
end.time1 = Sys.time()
time1=end.time1-start.time1
#trial1 Time difference of 2.574291 mins
#trial2 Time difference of 1.748892 mins
#trial3 Time difference of 1.520219 mins

#GRCh38:
 # GPIPE/9606/106/GCF_000001405.26_top_level
start.time2 = Sys.time()
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db GPIPE/9606/106/GCF_000001405.26_top_level -query C:/Users/lieberlab5/Dropbox/ma_sharefolder_Lieberlab/microRNA_project/bcl1mtc_fragile.fasta -out testt.out -outfmt 6 -remote")
end.time2 = Sys.time()
time2=end.time2-start.time2
#trial1 Time difference of 41.39214 secs
#trial2 Time difference of 1.53862 mins
#trial3 Time difference of 2.073691 mins

#GRCh37p.13:
 # GPIPE/9606/105/GCF_000001405.25_top_level
start.time3 = Sys.time()
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db GPIPE/9606/105/GCF_000001405.25_top_level -query C:/Users/lieberlab5/Dropbox/ma_sharefolder_Lieberlab/microRNA_project/bcl1mtc_fragile.fasta -out testt.out -outfmt 6 -remote")
end.time3 = Sys.time()
time3=end.time3-start.time3
#trial1 Time difference of 39.76398 secs
#trial2 Time difference of 26.05661 secs
#trial3 Time difference of 4.26681 mins


###################USE grch37########################
start.time3 = Sys.time()
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db GPIPE/9606/105/GCF_000001405.25_top_level -query test.fasta -out test.out -outfmt 6 -remote")
end.time3 = Sys.time()
time3=end.time3-start.time3
#trial 1 Time difference of 33.51035 secs

###########set parameters
start.time3 = Sys.time()
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db GPIPE/9606/105/GCF_000001405.25_top_level -query test.fasta -out test_word_size.out -outfmt 6 -remote -word_size=15")
end.time3 = Sys.time()
time3=end.time3-start.time3
#trial 1 Time difference of 1.721122 mins, set win_size - fail
#culling_limit: good parameter, but not feeding our purpose because we don't know if upstream is the best hit and then downstream is the second best or third or what.
#word_size 15: Time difference of 2.2875 mins


###########try blastn
start.time3 = Sys.time()
#have to use \ before " to make it still a ". It cannot be switched by ' in blast format sub options
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db GPIPE/9606/105/GCF_000001405.25_top_level -query test.fasta -out blastn_final_format.out -task blastn -outfmt \"6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -remote -word_size=15")
end.time3 = Sys.time()
time3=end.time3-start.time3
#Time difference of 59.77098 secs
#Time difference of 2.793861 mins

###########try default megablast
start.time3 = Sys.time()
#have to use \ before " to make it still a ". It cannot be switched by ' in blast format sub options
shell("C:/Users/lieberlab5/Desktop/Users/shuangcm/software/BLAST/blast-2.2.29+/bin/blastn -db GPIPE/9606/105/GCF_000001405.25_top_level -query test.fasta -out megablast_final_format.out -outfmt \"6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -remote -word_size=15")
end.time3 = Sys.time()
time3=end.time3-start.time3
#Time difference of 1.779506 mins

#following is with early tries
#trial 1 Time difference of 53.23175 secs
#trial 2 with format 7: Time difference of 5.878225 mins
#trial 3 with format 3: Time difference of 1.756226 mins
#trial 4 with format 1: Time difference of 1.721155 mins
###########function related
InputFile="test.fasta"
fasta=read.fasta(InputFile,as.string=T)

seq=list()
seqname=list()
for (i in 1:length(fasta)){
  seq=c(seq,fasta[i][[1]][1])
  seqname=c(seqname,attributes(fasta[i][[1]])$name[1])
}

minimum.length=30

for (i in 1:length(seq)){
  seqlength=nchar(seq[[i]])
  if (seqlength<minimum.length){
    warning(paste(seqname[[i]],"has only",seqlength,"bps, but we would recommand at least", minimum.length,"bps for an accurant alignment. Please refet to the user manual for details.",sep=" "))
  }
}

###################Database download BLAST##################


#read.blast in RFLPtools will help import BLAST result