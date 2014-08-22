
  upstreamseq = read.table("bcl-1 mtc _exclude18.txt.seq.txt")
  
  align=apply(upstreamseq["upstreamseq"],1,FUN=pairwiseAlignment,subject=ref,type = "local")#*--------apply function with other names. with all requirement of such function followed in the same argument
  #need to know how to make Views(align) work for more info w/o loop
  #"local" is essential here
  
  
  align.exclude = list()
  
  for (i in 1:length(align)){   
    if(align[[i]]@subject@range@width >= nchar(as.character(upstreamseq["upstreamseq"][i,]))) {
    align.exclude = c(align.exclude,align[i])
    }
  }
  
  upstreamseq.exclude = list()
  upstreamseq.names = list()
  upstreamseq.dir = list()
  
  for (i in 1:length(upstreamseq[,"upstreamseq"])){   
    if(align[[i]]@subject@range@width >= nchar(as.character(upstreamseq["upstreamseq"][i,]))) {
      upstreamseq.exclude = c(upstreamseq.exclude,as.character(upstreamseq[i,"upstreamseq"]))
      upstreamseq.names = c(upstreamseq.names,as.character(upstreamseq[i,"names"]))
      upstreamseq.dir = c(upstreamseq.dir,as.character(upstreamseq[i,"dir"]))
    }
  }
  
  view.align = lapply(align.exclude,FUN=Views)
  
  
  align.subref = pairwiseAlignment(subref,ref,type = "local")
  view.align.subref = Views(align.subref)
  
  #the reasone to include 1bp up and downstream is in case brpt cut right before subref
  
  subref.start = view.align.subref@ranges@start -1
  subref.end = subref.start + nchar(subref) -1 +1
  print (subref.start)
  print (subref.end)
  
  start.upstreamseq.exclude = list()
  
  for (i in 1:length(view.align)){
    start.upstreamseq.exclude = c(start.upstreamseq.exclude,view.align[[i]]@ranges@start)
  }
  
  print (start.upstreamseq.exclude)
  
  brpt.pos = list()
  
  for (i in 1:length(upstreamseq.exclude)){
    if (upstreamseq.dir[[i]] == "<"){
      brpt.pos = c(brpt.pos,unlist(start.upstreamseq.exclude[i]))
    }
    else{
      brpt.pos = c(brpt.pos,unlist(start.upstreamseq.exclude[i]) + nchar(as.character(upstreamseq.exclude[[i]])) - 1)
    }
    
  }
  
  print(brpt.pos)
  
  brpt.tb = cbind(t(data.frame(upstreamseq.names)),upstreamseq.exclude,t(data.frame(brpt.pos)),t(data.frame(upstreamseq.dir)))
  colnames(brpt.tb) = c("names","seq","brpt.pos","dir")
  
  print(brpt.tb)
  
  refseq = DNAString(ref)
  subref = DNAString(subref)
  
  index.motif.info = matchPattern("CG",refseq,fixed = F) 
  
  index.motif.start = c(index.motif.info@ranges@start)
  index.motif.end = index.motif.start + nchar("CG") -1
  
  #again, the subref got extend for 1 bps and has to consider for the above calculation
  
  print(index.motif.start)
  print(index.motif.end)
  
  index.motif.interval = IRanges(index.motif.start,index.motif.end)
  
  
  
  filter.brpt=subset(brpt.tb, brpt.tb[,"brpt.pos"] >= subref.start & brpt.tb[ ,"brpt.pos"] <= subref.end)  
  
  index.filter.brpt = list()
  
  for (i in 1:length(filter.brpt[,"brpt.pos"])){
    index.filter.brpt = c(index.filter.brpt,IRanges(filter.brpt[[i,"brpt.pos"]],filter.brpt[[i,"brpt.pos"]]))
  }
  
  print(index.filter.brpt)
  
  
  distance = data.frame()
  for (i in 1:length(index.filter.brpt)){
    distance= rbind(distance,distance(index.motif.interval,index.filter.brpt[[i]]))
  }
  
  print(distance)
  
  
  distance.min = data.frame()
  
  for (i in 1:length(index.filter.brpt)){
    distance.min = rbind(distance.min,min(distance(index.motif.interval,index.filter.brpt[[i]])))
  }
  
  distance.summary = data.frame(distance.min,distance)
  
  print(distance.summary)
  
  meanvalue = mean(unlist(distance.min))
  sd = sd(unlist(distance.min))
  
  df.distance.final = cbind (brpt.tb,distance.min)
  
  print(df.distance.final)
