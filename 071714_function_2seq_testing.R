##
##Version 1.0, Spring 2014. By Shuangchao Ma. Mentors: Dr. Michael Lieber，Dr. Kai Wang and Dr. Zhengfei Lu
##
##The goal can be found in the most recent report in Lieberlab dropbox
##

#v1.0 working directory
setwd("/Users/lyuan0388/Library/Dropbox/ma_share_folder/breakpt_project/tools translocations/breakpoints_txt/upsteam")
setwd("C:/Users/lieberlab5/Dropbox/ma_sharefolder_Lieberlab/breakpt_project/tools translocations/breakpoints_txt/upstream")


library("Biostrings","IRanges")

ref = "cacatcggtgtgcaggtttttgcgtggacgtctcaactcctttggataaaggcgaggagcataattgctgcactgcatattcggttagactgtgattagctttctaaaaagtggttttgttagatgtaaaaaatgaatatgacattctgaaacagaaaaaaataacttactctttatctgcgtgggatgagattaaactgcgtcttcttcgtggtttgaacgcaagagctccctgaacacctggcgctgccattggcgtgaacgaggggaagcccctcctgacagctggatggtaggacaaagccctctaagccccctctccccgtcacatccccccgaccctgcccacaagggaacctggggcactgggtgttcacctgcctcccactaggtgagatctttcctttttggctcctctatcctaatcctcacccacagcgttcccacggtggcctttgcagctaatatccagcgtccagagagtccctgggcttggttagcgttcctcccaaggctgaccctgagctaaagatgtagggacaggggaggtggtcccaggaagctcaggtgagacagagagaaaagcccaaagtgtg"
#refer to re-format file, or the original txt
subref = "ctttatctgcgtgggatgagattaaactgcgtcttcttcgtggtttgaacgcaagagctccctgaacacctggcgctgccattggcgtgaacgaggggaagcccctcctgacagctggatggtaggacaaagccctctaagccccctctc"

motif.distance = function(refseq,subref,upstreamseq,motif){
  upstreamseq = read.table(upstreamseq)
  
  align=apply(upstreamseq["upstreamseq"],1,FUN=pairwiseAlignment,subject=refseq,type = "local")#*--------apply function with other names. with all requirement of such function followed in the same argument
  #need to know how to make Views(align) work for more info w/o loop
  #"local" is essential here
  
  view.align = lapply(align,FUN=Views)
  
  align.subref = pairwiseAlignment(subref,refseq,type = "local")
  view.align.subref = Views(align.subref)
  
  subref.start = view.align.subref@ranges@start
  subref.end = subref.start + nchar(subref) -1
  print (subref.start)
  print (subref.end)
  
  start.upstreamseq = list()

  for (i in 1:length(view.align)){
    start.upstreamseq = c(start.upstreamseq,view.align[[i]]@ranges@start)
  }

  print (start.upstreamseq)
  
  brpt.pos = list()
  
  for (i in 1:length(upstreamseq[,1])){
    brpt.pos = c(brpt.pos,unlist(start.upstreamseq[i]) + nchar(as.character(upstreamseq[i,])) - 1)
  }

  #print(brpt.pos)
  #need to check
  #unlist(brpt.pos)>unlist(start.upstreamseq)
  
  
  brpt.tb = cbind(upstreamseq,t(data.frame(brpt.pos)))

  refseq = DNAString(refseq)
  subref = DNAString(subref)
  
  index.motif.info = matchPattern(motif,subref,fixed = F) 

  index.motif.start = c(index.motif.info@ranges@start) + subref.start -1
  index.motif.end = c(index.motif.start + nchar(motif) -1) + subref.start -1
  
  print(index.motif.start)
  
  index.motif.interval = IRanges(index.motif.start,index.motif.end)

  index.brpt = list()

  for (i in 1:length(brpt.pos)){
    index.brpt = c(index.brpt,IRanges(brpt.pos[[i]],brpt.pos[[i]]))
  }

  #print(index.brpt)
  
  

  distance = data.frame()
  for (i in 1:length(index.brpt)){
    distance= rbind(distance,distance(index.motif.interval,index.brpt[[i]]))
  }

  
  distance.min = data.frame()
  
  for (i in 1:length(index.brpt)){
      distance.min = rbind(distance.min,min(distance(index.motif.interval,index.brpt[[i]])))
  }

  distance.summary = data.frame(distance.min,distance)
  
  print(distance.summary)
  
  return(mean(unlist(distance.min)))

}

motif.distance(refseq = ref,subref = subref,upstreamseq = "bcl-1 mtc _exclude18.txt.seq.txt",motif = "cg")



#++++++++++#
