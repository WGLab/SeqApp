#v1.0 working directory
setwd("/Users/lyuan0388/Library/Dropbox/ma_share_folder/breakpt_project/tools translocations/breakpoints_txt/upstream")
setwd("C:/Users/lieberlab5/Dropbox/ma_sharefolder_Lieberlab/breakpt_project/tools translocations/breakpoints_txt/upstream")

#library
library("Biostrings","IRanges","seqinr")

#sequences the function will request
#for scl deletions (sil 25) (file name scl-sil.js for seq, file name scl 19 for reference)
ref = "tgtgaaaccttctcccagacaccgtggcgatggtactatactgtctgtactcacacgcatagcattctcgcactttgctgtttacactttacttccagtagtccatgaacttcaaacattaactgatcccttactacgtagcagagtgtgagccagaaattgtgaggaatgcaagagcgagcattcaaggtttcagtcctcaagaaattcccaccgtaaaaacaactttaaatggcacgttatcgactgcaaacctctttcacataccttagctcagatgatacccaaggatttcatggactgcttaatttctcccactagctccccttctcgtcatctgctccaggctcctggagggcggaacccactacagggctttgaatgttacccaccaaccttcccagggccccgctgccgcctctcaaggggaaaccaggagcacaaagctccacttaccgcaacttccgcggagctgaggtctgtttgcagggtaggagcgggagccggctccaaggagcgccaccgccgactccggccccgcctctgggacgttggggtcgcgggaactgaggcggcaaacacaagctcgcgaaactgaaggccgcgccacgcctgctgactggtcagaaacctcgccagtcagggagcgggggtggtgctagagagccgcgcacggtcgccgttacgtattggtgtg"
#refer to re-format file, or the original txt
subref = "tgagccagaaattgtgaggaatgca"

#ref = c2s(rev(comp(s2c(ref))))
#subref = c2s(rev(comp(s2c(subref))))
#**#?note this file has ref in negative strand, but for the brpt seq, > is used for positive strand, so have to find rev.comp for ref

#for scl deletions (sil 19) (file name scl-sil.js for seq, file name scl 19 for reference)
#ref = "cacaggcgggcctccgggcagagcagccgccgaccgggcgctgtccgcccacccaagccaacaactggctcccgaatacatcataatttggaataaaatgtgaaatcccgctccccgcccccatgccgccgcccccaccagcgcctcgatctctcgctcgccctccccccactcccgcccccagcgatttgcaaacgcacctctaaaggacacaggcacacaggcatacaactcagtgcggacaggaccacacagggtccagccccacagaagggcagcaaacaaacaccacctagcactgccccagaagccgacttggtcagccccgcacactgcccaacaggacacaacgaaatcagtcaaacgcagcggctcacggacacacaatccagcacagtcgggatcacacacgcccgacataacgacgacaccacccaaacacagtcgcaggggccacacacccccacacacagaatcagatccctgctgagaacaccaacggagcacaatcgtacgcaacgaagaaaacgcagaagggcctcgaagggtccacatctacacaccccaaccgcgcaggcacaccacactcggacacagagcctgtcgccaagaagaccacacttagaagcagccaacgccgcccacagttctcatgaacgcactctcacaatcccaccgcatgcacacaaccacgaagaagaaatgaaaaccaaccacagcctcgcgcatttctgtatattgcgtaaggaaaagggggaaggaaggaagagagtctccgaggcgggaggggcggcggcagccggggcgggggcgtccgtggaaaatgcccccccaccgccccccccggccatcgaaaggaaccgagggaaagggagcgaagacctcttttgcggacgtggggaaaaagaggaggaaattgatgaagaattcggtgggaggccgcgggttgcttttcccctagaaaaaggccagaatgctcacgtttttccgctacgggggtacgcgtgtg"
#refer to re-format file, or the original txt
#subref = "gaagaaatgaaaaccaacc"

#for bcl-1 mtc
#ref = "cacatcggtgtgcaggtttttgcgtggacgtctcaactcctttggataaaggcgaggagcataattgctgcactgcatattcggttagactgtgattagctttctaaaaagtggttttgttagatgtaaaaaatgaatatgacattctgaaacagaaaaaaataacttactctttatctgcgtgggatgagattaaactgcgtcttcttcgtggtttgaacgcaagagctccctgaacacctggcgctgccattggcgtgaacgaggggaagcccctcctgacagctggatggtaggacaaagccctctaagccccctctccccgtcacatccccccgaccctgcccacaagggaacctggggcactgggtgttcacctgcctcccactaggtgagatctttcctttttggctcctctatcctaatcctcacccacagcgttcccacggtggcctttgcagctaatatccagcgtccagagagtccctgggcttggttagcgttcctcccaaggctgaccctgagctaaagatgtagggacaggggaggtggtcccaggaagctcaggtgagacagagagaaaagcccaaagtgtg"
#refer to re-format file, or the original txt
#subref = "ctttatctgcgtgggatgagattaaactgcgtcttcttcgtggtttgaacgcaagagctccctgaacacctggcgctgccattggcgtgaacgaggggaagcccctcctgacagctggatggtaggacaaagccctctaagccccctctc"

#for bcl-2 mbr (175)
#ref = "cacaaatcctaaaagaagcattgaagtgaggtgtcatggattaattgacccctgtctatggaattacatgtaaaacattatcttgtcactgtagtttggttttatttgaaaacctgacaaaaaaaaagttccaggtgtggaatatgggggttatctgtacatcctggggcattaaaaaaaaaatcaatggtggggaactataaagaagtaacaaaagaagtgacatcttcagcaaataaactaggaaatttttttttcttccagtttagaatcagccttgaaacattgatggaataactctgtggcattattgcattatataccatttatctgtattaactttggaatgtactctgttcaatgtttaatgctgtggttgatatttcgaaagctgctttaaaaaaatacatgcatctcagcgtttttttgtttttaattgtatttagttatggcctatacactatttgtgagcaaaggtgatcgttttctgtttgagatttttatctcttgattcttcaaaagcattctgagaaggtgagataagccctgagtctcagctacctaagaaaaacctggatgtcactggccactgaggagctttgtttcaaccaagtcatgtgcatttccacgtcaacagaattgtttattgtgacagttatatctgttgtccctttgaccttgtttcttgaaggtttcctcgtccctgggcaattccgcatttaattcatggtattcaggattacatgcatgtttggttaaacccatgagattcattcagttaaaaatccagatggcaaatgaccagcagattcaaatctatggtggtttgacctttagagagttgctttacgtggcctgtttcaacacagacccacccagagccctcctgccctccttccgcgggggctttctcatggctgtccttcagggtcttcctgaaatgcagtggtgcttacgctccaccaagaaagcaggaaacctgtggtatgaagccagacctccccggcgggcctcagggaacagaatgatcagacctttgaatgattctaatttttaagcaaaatattattttatgaaaggtttacattgtcaaagtgatgaatatggaatatccaatcctgtgctgctatcctgccaaaatcattttaatggagtcagtttgcagtatgctccacgtggtaagatcctccaagctgctttagaagtaacaatgaagaacgtggacgtttttaatataaagcctgttttgtcttttgttgttgttcaaacgggattcaca"
#subref = "ggcctgtttcaacacagacccacccagagccctcctgccctccttccgcgggggctttctcatggctgtccttcagggtcttcctgaaatgcagtggtgcttacgctccaccaagaaagcaggaaacctgtggtatgaagccagacctccccggcgggcctcagggaacagaatg"

#for scl deletions (file name scl-sil.js for seq, file name scl 19 for reference)

#*NOTE: motif can be in ref, but brpt is only in subref(fragile zone)



#within the function: read breakpoints sequences
upstreamseq = apply(read.table("scl-sil.txt.seq.txt"),2,factor)


#upstreamseq = read.table("bcl-1 mtc _exclude18.txt.seq.txt")
#upstreamseq = read.table("bcl-2.txt.seq.txt")

#*---------motif.represent.start and moif.represent.end consider which positions in motif will be used to represent the motif------#
#e.g. CG will be |C|G|, but CAC will be only |CAC or GTG|
#situation like |G|YW will be considered with start and end given - so first find all GYW, but for distance, it will 
#only consider|G| in GYW; this is a same case like |CAC, except the end is 1.

motif.represent.start = 0
motif.represent.end = 0

#can be 0N1N2N3N4N5...

motif = "CAC"
motif.partial = substr(motif,motif.represent.start,motif.represent.end-1)
print (motif.partial)
#index.motif.partial = IRanges(motif.represent.start,motif.represent.end)

#*--------upordown consider if only look to the 5'or 3' of the motif--------#
#e.g. CAC is for 5' only
# it is already the up or down for the partial motif selected, so |G|YW will be all because it is all of G as the partial motif.

upordown = "up"
#others are down and all - either up/down, or all, all is default
#*so not considering super complicated cases like |CACA|C|G| #*it might be helpful to look into "precede" and "follow" functions in Iranges
#*nor case like |CA|C, when it is 1,3 but all except mid of partial motif. 
#basically, has to be up or down or all of partial motif. - this is doable though#*
#we will NOT consider any | within a partial motif, if not an all case


#*---------comp consider complementary situation, upon user's request--------#
#comp can only be t or f

comp = T

if(comp == F){
  motif = motif
}else {
  comp.motif = c2s(rev(comp(s2c(motif))))
  #comp.motif.represent.start = nchar(motif) +1 - motif.represent.end
  #comp.motif.represent.end = nchar(motif) +1 - motif.represent.start
}
#?no need of comp.motif.partial
#*note that can't use motif = c(motif, c2s(rev(comp(s2c(motif))))) as there are too many motif defined in one line and confused the system
#*note ese should be with the previous { so that R knows more codes are coming
#comp.motif

#*---------for palindromic motif only: e.g. CG. ---------#
#CAC not considered

upstreamseq.revcomp = t(data.frame(lapply(upstreamseq[,"upstreamseq"],function(x){c2s(rev(comp(s2c(x))))})))
upstreamseq = cbind(upstreamseq,upstreamseq.revcomp)

colnames(upstreamseq) = c("names","dir","unstreamseq.ori","upstreamseq","upstreamseq.revcomp")

#Step 0: preprocess to ensure:
#1. alignment of brpt seq and all breakpoints seq not excluded are 100% match
align=lapply(upstreamseq[,"upstreamseq.revcomp"],FUN=pairwiseAlignment,subject=ref,type = "local")
#*------apply function with other names. with all requirement of such function followed in the same argument
#?------need to know how to make Views(align) work for more info w/o loop
#*------"local" is essential here
#?#*----it only include exact match here


align.start = list()
align.width = list()
for (i in 1:length(align)){   
  align.start = c(align.start,align[[i]]@subject@range@start)
  align.width = c(align.width,align[[i]]@subject@range@width) 
} 

align.tb = cbind(upstreamseq,t(data.frame(align.start)),t(data.frame(align.width)))
colnames(align.tb) = c("names","dir","unstreamseq.ori","upstreamseq","upstreamseq.revcomp","start","width")

align.tb = apply(align.tb,2,factor)

print(align.tb)

#**subset dataframe based on a comparison of two columns
align.exrefout.tb=align.tb[apply(align.tb,1,function(x)x["width"]==nchar(x["upstreamseq"])),]


length(align.exrefout.tb[,"upstreamseq.revcomp"])
#?so far only those upstreamseq can completely match with ref 

#Step I: Alignment for subref
align.subref = pairwiseAlignment(subref,ref,type = "local")

#------change ref and subref into DNAString for matchpattern
refseq = DNAString(ref)
subref = DNAString(subref)
#Step II: Index to relative postions we need for start/end of breakpoint seq, subref, and motifs
#1. for breakpoints seq, only the bkpt.pos. It considers end for >, start for <

#####**NOTE dir function got switched for sil/scl samples
dir.f=function(x){
  if(x["dir"]=="<"){
    brpt.pos=as.numeric(x["start"])+as.numeric(nchar(x["upstreamseq"])-1)
  }
  else{
    brpt.pos=as.numeric(x["start"])-1
  }
}

brpt.pos=data.frame(apply(align.exrefout.tb,1,dir.f))

brpt.exrefout.tb=apply(cbind(align.exrefout.tb,brpt.pos),2,factor)
colnames(brpt.exrefout.tb) = c("names","dir","unstreamseq.ori","upstreamseq","upstreamseq.revcomp","start","width","brpt.pos")
print(brpt.exrefout.tb)

length(brpt.exrefout.tb[,"brpt.pos"])

#?*should include a check ensure no brpt.pos is duplicated

#2. start and end point of subref
view.align.subref = Views(align.subref)

subref.start = view.align.subref@ranges@start -1
subref.end = subref.start + nchar(subref) -1 +1
#*------the reasone to include 1bp up and downstream is in case brpt cut right before subref

print (subref.start)
print (subref.end)

#2.5 exclude any brpt that is in ref, but not in subref.
#*------we only interested in brpt falls in frigle zone.
filter.brpt=subset(brpt.exrefout.tb, brpt.exrefout.tb[,"brpt.pos"] >= subref.start & brpt.exrefout.tb[ ,"brpt.pos"] <= subref.end)  
filter.brpt = apply(filter.brpt,2,factor)

length(filter.brpt[,"brpt.pos"])

#3. start and end point of motif
#------consider rev.comp situation------#
if (comp == F){
  index.motif.info = matchPattern(motif,refseq,fixed = F)
  index.motif.partial.info.start = index.motif.info@ranges@start + motif.represent.start
  index.motif.partial.info.end = index.motif.info@ranges@start +motif.represent.end
}else{
  index.motif.info = matchPattern(motif,refseq,fixed = F)
  index.motif.partial.info.start = index.motif.info@ranges@start + motif.represent.start
  index.motif.partial.info.end = index.motif.info@ranges@start +  motif.represent.end
  
  index.comp.motif.info = matchPattern(comp.motif,refseq,fixed = F)
  index.comp.motif.partial.info.start = index.comp.motif.info@ranges@start + nchar(motif) - motif.represent.end -1
  index.comp.motif.partial.info.end = index.comp.motif.info@ranges@start + nchar(motif) - motif.represent.start -1
  }
#*------again, the subref got extend for 1 bps and has to consider for the above calculation
print(index.motif.partial.info.start)
print(index.motif.partial.info.end)
print(index.comp.motif.partial.info.start)
print(index.comp.motif.partial.info.end)

#Step III: Range index for brpt and motifs introduced
index.filter.brpt = mapply(IRanges,as.numeric(filter.brpt[,"brpt.pos"]),as.numeric(filter.brpt[,"brpt.pos"]))
index.motif.partial.interval = IRanges(index.motif.partial.info.start,index.motif.partial.info.end)
index.comp.motif.partial.interval = IRanges(index.comp.motif.partial.info.start,index.comp.motif.partial.info.end)

print(index.filter.brpt)
print(index.motif.partial.interval)
print(index.comp.motif.partial.interval)
#Step IV: calculate distances and min in all distances

#*if it is upstream, but not adjacent (3 situation: NN(bk)N)
#*probably need to simplify the loop in loop structure

#dir.f=function(x){
#  ifelse(index.motif.interval[[i]][2] < as.numeric(x["brpt.pos"]),
#         distance(index.motif.interval[i],IRanges(as.numeric(x["brpt.pos"]),as.numeric(x["brpt.pos"])))+1,
#         distance(index.motif.interval[i],IRanges(as.numeric(x["brpt.pos"]),as.numeric(x["brpt.pos"]))))
#}

#dis=list()
#for (i in 1:length(index.motif.interval)){
#  dis = c(dis,apply(filter.brpt,1,dir.f))
#  }
#20
#do.call(rbind.data.frame,dis)

#how to simply, the above method would take longer time

#min.distance = "unknown"

#080914checked everything above

#distance.tb=cbind(filter.brpt,min.distance)

if(upordown=="up"){
  
  min.distance = list()
  
  for (i in 1:length(filter.brpt[,"brpt.pos"])){
    now.index.filter.brpt.num = as.numeric(index.filter.brpt[[i]]) 
    now.index.filter.brpt=index.filter.brpt[[i]]
    #consider |CAC
    index.motif.partial.interval.up = index.motif.partial.interval[index.motif.partial.interval@start>now.index.filter.brpt.num]
    min.distance.motif = distanceToNearest(now.index.filter.brpt,index.motif.partial.interval.up)@elementMetadata$distance
    
    #consider GTG|,but not when breakpoint is at |
    if (any(index.comp.motif.partial.interval@start==now.index.filter.brpt.num)){
    index.comp.motif.partial.interval.down = index.comp.motif.partial.interval[index.comp.motif.partial.interval@start==now.index.filter.brpt.num]

    min.distance.comp.motif.down = distanceToNearest(now.index.filter.brpt,index.comp.motif.partial.interval.down)@elementMetadata$distance
    }else if(any(index.comp.motif.partial.interval@start<now.index.filter.brpt.num)){
    #consider GTG|, when breakpoint is at |
    
    index.comp.motif.partial.interval.down = index.comp.motif.partial.interval[index.comp.motif.partial.interval@start<now.index.filter.brpt.num]
    min.distance.comp.motif.down = distanceToNearest(now.index.filter.brpt,index.comp.motif.partial.interval.down)@elementMetadata$distance+1
    }
    
    min.distance = c(min.distance,min(min.distance.motif,min.distance.comp.motif.down))
    
  }

}else if (upordown=="down"){
 #*#?not tested yet 
  min.distance = list()
  
  for (i in 1:length(filter.brpt[,"brpt.pos"])){
    now.index.filter.brpt.num = as.numeric(index.filter.brpt[[i]]) 
    now.index.filter.brpt=index.filter.brpt[[i]]
    #consider CAC|
    index.comp.motif.partial.interval.up = index.comp.motif.partial.interval[index.comp.motif.partial.interval@start>now.index.filter.brpt.num]
    min.distance.comp.motif = distanceToNearest(now.index.filter.brpt,index.comp.motif.partial.interval.up)@elementMetadata$distance
    
    #consider |GTG,but not when breakpoint is at |
    if (any(index.motif.partial.interval@start==now.index.filter.brpt.num)){
      index.motif.partial.interval.down = index.motif.partial.interval[index.motif.partial.interval@start==now.index.filter.brpt.num]
      
      min.distance.motif.down = distanceToNearest(now.index.filter.brpt,index.motif.partial.interval.down)@elementMetadata$distance
    }else if(any(index.motif.partial.interval@start<now.index.filter.brpt.num)){
      #consider GTG|, when breakpoint is at |
      
      index.motif.partial.interval.down = index.motif.partial.interval[index.motif.partial.interval@start<now.index.filter.brpt.num]
      min.distance.motif.down = distanceToNearest(now.index.filter.brpt,index.motif.partial.interval.down)@elementMetadata$distance+1
    }
    
    min.distance = c(min.distance,min(min.distance.motif,min.distance.comp.motif.down))
  }

}else {
  index.combine.motif.partial.interval = c(index.comp.motif.partial.interval,index.motif.partial.interval)
  
  min.distance = list()
  
  for (i in 1:length(filter.brpt[,"brpt.pos"])){
    #consider |GTG,but not when breakpoint is at |
    if (any(index.combine.motif.partial.interval@start>=now.index.filter.brpt.num)){
      index.combine.motif.partial.interval.down = index.combine.motif.partial.interval[index.motif.partial.interval@start>=now.index.filter.brpt.num]
      
      min.distance.combine.motif.down = distanceToNearest(now.index.filter.brpt,index.combine.motif.partial.interval.down)@elementMetadata$distance
    }else{
      #consider GTG|, when breakpoint is at |
      
      index.combine.motif.partial.interval.up = index.combine.motif.partial.interval[index.combine.motif.partial.interval@start<now.index.filter.brpt.num]
      min.distance.combine.motif.up = distanceToNearest(now.index.filter.brpt,index.combine.motif.partial.interval.up)@elementMetadata$distance+1
      }
    
    min.distance = c(min.distance,min(min.distance.combine.motif.up,min.distance.combine.motif.down))
    }
  }
}

min.distance.num = list()
for (i in 1:length(min.distance)){
  min.distance.num = c(min.distance.num,min.distance[[i]]@elementMetadata$distance)
}

print (min.distance.num)

#distance.tb = data.frame(matrix(NA,ncol = length(index.motif.interval),nrow= length(index.filter.brpt)))
#for (i in 1:nrow(distance.tb)){
#  for (u in 1:length(index.motif.interval)){
#    if (index.motif.interval[[u]][length(index.motif.interval[[u]])] < index.filter.brpt[[i]]@start ){
#      #length(index.motif.interval[[u]]) because in this s4 obj, after [[]], all position
#      #for each motif letters will show
      
#      min.distance.tb[i,u]= distance(index.motif.interval[u],index.filter.brpt[[i]])+1
#    }
#    else{
#      min.distance.tb[i,u]= distance(index.motif.interval[u],index.filter.brpt[[i]])
#    }
#  }
#}

#*A warning is printed if the desired length is not an even multiple of the original vector's length.

print(distance.tb)

distance.min = data.frame()
for (i in 1:length(index.filter.brpt)){
  distance.min = rbind(distance.min,min(distance.tb[i,]))
}
#?not sure how to simply the above
#Step V: present result in a df

distance.summary = cbind(distance.min,distance.tb)

print(distance.summary)

meanvalue = mean(unlist(distance.min))
stdv = sd(unlist(distance.min))

print(meanvalue)
print(stdv)

df.distance.final = cbind (filter.brpt,distance.min)

print(df.distance.final)

