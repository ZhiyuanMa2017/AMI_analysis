setwd("C:/Users/mzy/Documents/学习/CHD病理/CHD病理/新建文件夹")
library(omicade4)
a = list.files()
mrna = read.table(a[1],header = T,row.names=1)
protein = read.csv(a[2],header = T,row.names=1)
protein1=protein[,c(16,17,18,13,14)]
list1=list(mrna=mrna,protein=protein1)
mcoin <- mcia(list1)
class(mcoin)
cancer_type <- colnames(list1$protein)
cancer_type <- sapply(strsplit(cancer_type, split="_"), function(x) x[1])
plot(mcoin, axes=1:2, phenovec=cancer_type, sample.lab=FALSE, df.color=1:2)
chd_gene1=selectVar(mcoin, a1.lim=c(2,Inf))
chd_gene2=selectVar(mcoin, a1.lim=c(-Inf,-1))
gene=rbind(chd_gene1,chd_gene2)
setwd("C:/Users/mzy/Documents")
#mcia结果
proteingene=gene[which(gene[,3]==TRUE),]
proteingene=proteingene[,1]
mrnagene=gene[which(gene[,2]==TRUE),]
mrnagene=mrnagene[,1]
for(i in 1:63){proteingene[i]=sub('^.','',proteingene[i])}
for(i in 1:2264){mrnagene[i]=sub('^.','',mrnagene[i])}
geneallname=gene[,1]
for(i in 1:2316){geneallname[i]=sub('^.','',geneallname[i])}
unique(geneallname) 

mrnag=data.frame()
for (i in 1:10034){
  if(rownames(mrna[i,])%in%mrnagene==T){
    mrnag=rbind(mrnag,mrna[i,])
  }
}
proteing=data.frame()
for (i in 1:567){
  if(rownames(protein1[i,])%in%proteingene==T){
    proteing=rbind(proteing,protein1[i,])
  }
}
nametemp1=row.names(proteing)
nametemp2=row.names(mrnag)

library(WGCNA)
mrnaallpvalue=corPvalueStudent(cor(t(mrnag),method=("kendall")),n=5)
proteinallpvalue=corPvalueStudent(cor(t(proteing),method=("kendall")),n=5)

mrnaallpvaluefdr=p.adjust(mrnaallpvalue,method="fdr")
mrnaallpvaluefdr=matrix(mrnaallpvaluefdr,2264,2264)
row.names(mrnaallpvaluefdr)=nametemp2
colnames(mrnaallpvaluefdr)=nametemp2

proteinallpvaluefdr=p.adjust(proteinallpvalue,method="fdr")
proteinallpvaluefdr=matrix(proteinallpvaluefdr,63,63)
row.names(proteinallpvaluefdr)=nametemp1
colnames(proteinallpvaluefdr)=nametemp1

mrnaallpvaluefdr[as.numeric(mrnaallpvaluefdr)>=0.05]=2
mrnaallpvaluefdr[as.numeric(mrnaallpvaluefdr)<0.05]=1
mrnaallpvaluefdr[mrnaallpvaluefdr==2]=0
for(i in 1:2264){mrnaallpvaluefdr[i,i]=0}


proteinallpvaluefdr[as.numeric(proteinallpvaluefdr)>=0.05]=2
proteinallpvaluefdr[as.numeric(proteinallpvaluefdr)<0.05]=1
proteinallpvaluefdr[proteinallpvaluefdr==2]=0
for(i in 1:63){proteinallpvaluefdr[i,i]=0}

library(igraph)
graphmrna=graph_from_adjacency_matrix(mrnaallpvaluefdr,mode =  "undirected")
graphprotein=graph_from_adjacency_matrix(proteinallpvaluefdr,mode =  "undirected")
uniongraph=union(graphmrna,graphprotein,byname=T)
allpagerank=page.rank(uniongraph)$vector

tempdf<- data.frame(Object = 1:2313, PageRank =allpagerank)
tempdf[,2]=as.numeric(tempdf[,2])
pagerankmax=max(tempdf[,2])
pagerankmin=min(tempdf[,2])
for(i in 1:2313){
  tempdf[i,3]=(tempdf[i,2]-pagerankmin)/(pagerankmax-pagerankmin)
}
tempdf=  tempdf[order(  tempdf[,3],decreasing=T),]
#所有的pagerank值

library(KEGGgraph)
pathwayid=read.csv("C:/Users/mzy/Documents/pathway.csv")
keggid=c()
for(i in 1:64){
  x=pathwayid[i,1]
  x=sub('^.','',x)
  x=sub('^.','',x)
  x=sub('^.','',x)
  keggid=c(keggid,x)
}
# sapply(keggid,function(x)  retrieveKGML(x, organism="rno", destfile=paste("C:/Users/mzy/Documents/xml/",x, ".xml"),method = "internal"))
filelist=list.files(path="C:/Users/mzy/Documents/xml")
for(i in 1:64){
  filelist[i]=paste("C:/Users/mzy/Documents/xml/",filelist[i],sep="")
}
kegggraphg=list()
for( i in 1:64){
  kegggraphg[[i]]=parseKGML2Graph(filelist[i],expandGenes=TRUE)
}

mergedKEGG=mergeKEGGgraphs(kegggraphg)
geneid=translateKEGGID2GeneID(nodes(mergedKEGG)) 
selectmrnaname=row.names(mrna)
selectproteinname=row.names(protein1)
selectid=c()
for(i in 1:10034){
  if((selectmrnaname[i]%in%geneid)&&(all(mrna[i,1:5]==0)==FALSE)){
    selectid=c(selectid,selectmrnaname[i])
  }
}
for(i in 1:567){
  if((selectproteinname[i]%in%geneid)&&(all(protein1[i,1:5]==0)==FALSE)){
    selectid=c(selectid,selectproteinname[i])
  }
}
selectid=unique(selectid)

chongfugene=intersect(row.names(tempdf),selectid)

selectidpr=data.frame()
tempid=c()
for (i in 1:2313){
  if(row.names(tempdf)[i]%in%chongfugene){
    tempid=c(row.names(tempdf)[i],tempdf[i,2],tempdf[i,3])
    selectidpr=rbind(selectidpr,tempid, stringsAsFactors = F)
    tempid=c()
  }
}
selectidpr[,2]=as.numeric(selectidpr[,2])
selectidpr[,3]=as.numeric(selectidpr[,3])
selectidpr=selectidpr[order(selectidpr[,2],decreasing=T),]
colnames(selectidpr)=c('genename','pagerank','biaozhunhuapagerank')


selectkeggid=translateGeneID2KEGGID(selectid,organism="rno")
finalmergedgraph=subGraph(selectkeggid,mergedKEGG)

fgraph=igraph.from.graphNEL(finalmergedgraph)
finalgraph=fgraph-names(degree(fgraph)[degree(fgraph)==0])
aa=edge.betweenness.community(finalgraph)
bb=walktrap.community(finalgraph)
cc=leading.eigenvector.community(finalgraph)

aascore=data.frame()

for(i in 1:82){
  print(i)
  ccc=c()
  
  for(j in 1:length(bb[[i]])){
    x=bb[[i]][j]
    x=sub('^.','',x)
    x=sub('^.','',x)
    x=sub('^.','',x)
    x=sub('^.','',x)
    ccc=c(ccc,x)}
  
  ndiff=length(ccc[ccc%in%intersect(ccc,chongfugene)==F])
  missvalue=0
  missvalue=ndiff/(length(ccc)*507)
  tempp=data.frame()
  phit=0
  pmiss=0
  pscore=0
  score=c()
  
  for(k in 1:507){
    if(selectidpr[k,1]%in%ccc){
      phit=phit+selectidpr[k,3]
      
    }
    if(selectidpr[k,1]%in%ccc==F){
      pmiss=pmiss+missvalue
    }
    pscore=phit-pmiss
    p=c(phit,pmiss,pscore)
    tempp=rbind(tempp,p)
  }
  
  score=c(paste("bb", i, sep = ""),max(tempp[,3]))
  aascore=rbind(aascore,score,stringsAsFactors = F )
  
}
for(i in 1:82){aascore[i,3]=length(bb[[i]])}

randomscore=data.frame()


for(i in 42:82){
  scorelist=c()
  for(j in 1:10000){
    changdu=length(bb[[i]])
    genelist=c()
    genelist=sample(geneid,changdu)
    ccc=c()
    
    for(f in 1:changdu){
      x=genelist[f]
      ccc=c(ccc,x)}
    ndiff=length(ccc[ccc%in%intersect(ccc,chongfugene)==F])
    missvalue=0
    missvalue=ndiff/(length(ccc)*507)
    tempp=data.frame()
    phit=0
    pmiss=0
    pscore=0
    score=c()
    
    for(k in 1:507){
      if(selectidpr[k,1]%in%ccc){
        phit=phit+selectidpr[k,3]
        
      }
      if(selectidpr[k,1]%in%ccc==F){
        pmiss=pmiss+missvalue
      }
      pscore=phit-pmiss
      p=c(phit,pmiss,pscore)
      tempp=rbind(tempp,p)
    }
    score=max(tempp[,3])
    scorelist=c(scorelist,score)
  }
  randomscore=rbind(randomscore,scorelist)
}


for(j in 42:82){
  pvalue=0
  for(i in 1:10000){
    if(as.numeric(randomscore[j-40,i])>=as.numeric(aascore[j,2])){
      pvalue=pvalue+1
    }
  }
  pvalue=pvalue/10000
  aascore[j,4]=pvalue
}
colnames(aascore)=c("mokuai","dafen","geshu","pvalue")
write.csv(aascore,file="C:/Users/mzy/Documents/学习/CHD病理/mokuaidafen.csv")

for(i in 1:82){
  print(i)
  ccc=c()
  
  for(j in 1:length(bb[[i]])){
    x=bb[[i]][j]
    x=sub('^.','',x)
    x=sub('^.','',x)
    x=sub('^.','',x)
    x=sub('^.','',x)
    ccc=c(ccc,x)}
  
  ndiff=length(ccc[ccc%in%intersect(ccc,chongfugene)==F])
  missvalue=0
  missvalue=ndiff/(length(ccc)*507)
  tempp=data.frame()
  phit=0
  pmiss=0
  pscore=0
  score=c()
  
  for(k in 1:507){
    if(selectidpr[k,1]%in%ccc){
      phit=phit+selectidpr[k,3]
      
    }
    if(selectidpr[k,1]%in%ccc==F){
      pmiss=pmiss+missvalue
    }
    pscore=phit-pmiss
    p=c(phit,pmiss,pscore)
    tempp=rbind(tempp,p)
  }
  tiff(filename = paste("plot_",i,".tiff", sep = ""),
       width = 3200, height = 3200, units = "px", pointsize = 12,
       compression = "lzw", res = 400)
  plot(tempp[,3],type='l',ylab="score", xlab="gene")
  dev.off()
  
}
#模块画图
nodes=c()
for(i in 1:48){
  if(i==5||i==11||i==45||i==48){
    for(j in (1:length(bb[[i]]))){
      nodeid=translateKEGGID2GeneID(bb[[i]][j]) 
      nodes=c(nodes,nodeid)
      nodeid=0
}
  }
}
nodeframe=as.data.frame(nodes)
for (i in 1:nrow(nodeframe)){
  if (i<167){
    nodeframe[i,2]=5
  }
  if(i>166&&i<335){
    nodeframe[i,2]=11}
if(i>334&&i<344){
  nodeframe[i,2]=45}
if(i>343){
  nodeframe[i,2]=48}
}

for (i in 1:354){
  for(j in 1:2313){
    if(nodeframe[i,1]==row.names(tempdf)[j]){
      nodeframe[i,3]=tempdf[j,2]
      }
  }
}
for (i in 1:354){
  if(is.na(nodeframe[i,3])==TRUE){
    nodeframe[i,3]=0
  }
}
for(i in 1:354){
  nodeframe[i,4]=bitr(nodeframe[i,1], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Rn.eg.db", drop = FALSE)[2]
}
write.csv(nodeframe,file="C:/Users/mzy/Documents/学习/CHD病理/nodes.csv")
