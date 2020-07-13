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

proteingene=gene[which(gene[,3]==TRUE),]
proteingene=proteingene[,1]
mrnagene=gene[which(gene[,2]==TRUE),]
mrnagene=mrnagene[,1]
for(i in 1:63){proteingene[i]=sub('^.','',proteingene[i])}
for(i in 1:2264){mrnagene[i]=sub('^.','',mrnagene[i])}
geneallname=gene[,1]
for(i in 1:2316){geneallname[i]=sub('^.','',geneallname[i])}
unique(geneallname) 

# x[duplicated(x)] 
# geneallname=c(proteingene,mrnagene)
# geneallname=geneallname[!duplicated(geneallname)]

 #修改名称去掉x
write.table(geneallname,file="C:/Users/mzy/Documents/学习/CHD病理/genealllist.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(proteingene,file="C:/Users/mzy/Documents/学习/CHD病理/proteingenelist.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(mrnagene,file="C:/Users/mzy/Documents/学习/CHD病理/mrnagenelist.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)

bothgene=c()
for (i in 1:52){if(proteingene[i]%in%mrnagene){bothgene=c(bothgene,proteingene[i])}}
for(i in 1:14){bothgene[i]=sub('^.','',bothgene[i])}
#找出了重复的基因14个

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

#取出来是那些基因的信息

# proteinbothrowname=list()
# for(i in 1:11){
#   proteinbothrowname[i]=paste("PM", row.names(proteinbothg)[i], sep="") 
# }
# mrnabothrowname=list()
# for(i in 1:11){
#   mrnabothrowname[i]=paste("MP", row.names(mrnabothg)[i], sep="") 
# }
# proteinrowname=list()
# for(i in 1:52){
#   proteinrowname[i]=paste("P", row.names(proteing)[i], sep="") 
# }
# mrnarowname=list()
# for(i in 1:2253){
#   mrnarowname[i]=paste("M", row.names(mrnag)[i], sep="") 
# }
# #修改基因名 转录的加上m 蛋白的加上p 

# proteinboth=data.frame(row.names =proteinbothrowname,F_1=proteinbothg$F_1,F_2=proteinbothg$F_2,F_3=proteinbothg$F_3,E_1=proteinbothg$E_1,E_2=proteinbothg$E_2 )
# mrnaboth=data.frame(row.names =mrnabothrowname,F_1=mrnabothg$J52,F_2=mrnabothg$J54,F_3=mrnabothg$J57,E_1=mrnabothg$M26,E_2=mrnabothg$M30)

# proteingg=data.frame(row.names =proteingene,F_1=proteing$F_1,F_2=proteing$F_2,F_3=proteing$F_3,E_1=proteing$E_1,E_2=proteing$E_2 )
# mrnagg=data.frame(row.names =mrnagene,F_1=mrnag$J52,F_2=mrnag$J54,F_3=mrnag$J57,E_1=mrnag$M26,E_2=mrnag$M30)
#定义新的dataframe 
finalgene=rbind(proteinboth,mrnaboth,proteingg,mrnagg)
# genefinal=rbind(proteingg,mrnagg)
#最后的记忆表达量数据
proteinall=rbind(proteinboth,proteingg)
mrnaall=rbind(mrnaboth,mrnagg)
nametemp1=row.names(proteinall)
nametemp2=row.names(mrnaall)
for(i in 1:63){nametemp1[i]=sub('^.','',nametemp1[i])}
for(i in 1:2264){nametemp2[i]=sub('^.','',nametemp2[i])}
for(i in 1:11){nametemp1[i]=sub('^.','',nametemp1[i])}
for(i in 1:11){nametemp2[i]=sub('^.','',nametemp2[i])}
row.names(proteinall)=nametemp1
row.names(mrnaall)=nametemp2

# 
# p.kendallcor=cor(t(proteinall),use="p",method=("kendall"))
# p.kendall=p.kendallcor[lower.tri(p.kendallcor)]
# hist(p.kendall)
# p.spearmancor=cor(t(proteinall),use="p",method=("spearman"))
# p.spearman=p.spearmancor[lower.tri(p.spearmancor)]
# hist(p.spearman)
# 
# m.kendallcor=cor(t(mrnaall),use="p",method=("kendall"))
# m.kendall=m.kendallcor[lower.tri(m.kendallcor)]
# hist(m.kendall)
# m.spearmancor=cor(t(mrnaall),use="p",method=("spearman"))
# m.spearman=m.spearmancor[lower.tri(m.spearmancor)]
# hist(m.spearman)
# 
# proteinalldisease=proteinall[,1:3]
# mrnaalldisease=mrnaall[,1:3]
# 
# pd.kendallcor=cor(t(proteinalldisease),use="p",method=("kendall"))
# pd.kendall=pd.kendallcor[lower.tri(pd.kendallcor)]
# hist(pd.kendall)
# pd.spearmancor=cor(t(proteinalldisease),use="p",method=("spearman"))
# pd.spearman=pd.spearmancor[lower.tri(pd.spearmancor)]
# hist(pd.spearman)
# 
# md.kendallcor=cor(t(mrnaalldisease),use="p",method=("kendall"))
# md.kendall=md.kendallcor[lower.tri(md.kendallcor)]
# hist(md.kendall)
# md.spearmancor=cor(t(mrnaalldisease),use="p",method=("spearman"))
# md.spearman=md.spearmancor[lower.tri(md.spearmancor)]
# hist(md.spearman)
# #疾病组 相关系数图
# n=0
# for(i in 1:63){
# c<-c(proteinall[i,1],proteinall[i,2],proteinall[i,3],proteinall[i,4],proteinall[i,5])
# if(shapiro.test(c)$p.value < 0.05)
#   {n=n+1}
# c<-c()
# }
# m=0
# for(i in 1:2264){
#   c<-c(mrnaall[i,1],mrnaall[i,2],mrnaall[i,3],mrnaall[i,4],mrnaall[i,5])
#   if(shapiro.test(c)$p.value < 0.05)
#   {m=m+1}
#   c<-c()
# }
#  
# kendallcor=cor(t(finalgene),use="p",method=("kendall"))
# kendallshuzi=kendallcor[lower.tri(kendallcor)]
# hist(kendallshuzi)
# spearmancor=cor(t(finalgene),use="p",method=("spearman"))
# spearmanshuzi=spearmancor[lower.tri(spearmancor)]
# hist(spearmanshuzi)

# library(WGCNA)
# kendallpvalue=corPvalueStudent(cor(t(genefinal),method=("kendall")),n=5)
# #spearman方式求pvalue矩阵
# names=row.names(genefinal)
# 
# holmpvalue=p.adjust(kendallpvalue,method="holm")
# holmpvalue=matrix(holmpvalue,2305,2305)
# row.names(holmpvalue)=names
# colnames(holmpvalue)=names 
# #p.adjust第一种方式holm
# hochbergpvalue=p.adjust(kendallpvalue,method="hochberg")
# hochbergpvalue=matrix(hochbergpvalue,2305,2305)
# row.names(hochbergpvalue)=names
# colnames(hochbergpvalue)=names 
# #p.adjust第二种方式hochberg
# bonferronipvalue=p.adjust(kendallpvalue,method="bonferroni")
# bonferronipvalue=matrix(bonferronipvalue,2305,2305)
# row.names(bonferronipvalue)=names
# colnames(bonferronipvalue)=names 
# #p.adjust第三种方式bonferroni
# BHpvalue=p.adjust(kendallpvalue,method="BH")
# BHpvalue=matrix(BHpvalue,2305,2305)
# row.names(BHpvalue)=names
# colnames(BHpvalue)=names 
# #p.adjust第四种方式BH
# BYpvalue=p.adjust(kendallpvalue,method="BY")
# BYpvalue=matrix(BYpvalue,2305,2305)
# row.names(BYpvalue)=names
# colnames(BYpvalue)=names 
# #p.adjust第五种方式BY
# fdrpvalue=p.adjust(kendallpvalue,method="fdr")
# fdrpvalue=matrix(fdrpvalue,2305,2305)
# row.names(fdrpvalue)=names
# colnames(fdrpvalue)=names 
# #p.adjust第五种方式fdr
# 
# holmpvalue[as.numeric(holmpvalue)>=0.05]=2
# holmpvalue[as.numeric(holmpvalue)<0.05]=1
# holmpvalue[holmpvalue==2]=0
# hochbergpvalue[as.numeric(hochbergpvalue)>=0.05]=2
# hochbergpvalue[as.numeric(hochbergpvalue)<0.05]=1
# hochbergpvalue[hochbergpvalue==2]=0
# bonferronipvalue[as.numeric(bonferronipvalue)>=0.05]=2
# bonferronipvalue[as.numeric(bonferronipvalue)<0.05]=1
# bonferronipvalue[bonferronipvalue==2]=0
# BHpvalue[as.numeric(BHpvalue)>=0.05]=2
# BHpvalue[as.numeric(BHpvalue)<0.05]=1
# BHpvalue[BHpvalue==2]=0
# BYpvalue[as.numeric(BYpvalue)>=0.05]=2
# BYpvalue[as.numeric(BYpvalue)<0.05]=1
# BYpvalue[BYpvalue==2]=0
# fdrpvalue[as.numeric(fdrpvalue)>=0.05]=2
# fdrpvalue[as.numeric(fdrpvalue)<0.05]=1
# fdrpvalue[fdrpvalue==2]=0
# #换成0-1矩阵
# for(i in 1:2305){fdrpvalue[i,i]=0}
# 
# library(igraph)
# g1=graph_from_adjacency_matrix(holmpvalue,mode =  "undirected")
# g2=graph_from_adjacency_matrix(hochbergpvalue,mode =  "undirected")
# g3=graph_from_adjacency_matrix(BHpvalue,mode =  "undirected")
# g4=graph_from_adjacency_matrix(BYpvalue,mode =  "undirected")
# g5=graph_from_adjacency_matrix(fdrpvalue,mode =  "undirected")
# 
# g1degree=degree(g1)
# g2degree=degree(g2)
# g3degree=degree(g3)
# g4degree=degree(g4)
# g5degree=degree(g1)
# 
# g1degreetable=table(g1degree)
# g2degreetable=table(g2degree)
# g3degreetable=table(g3degree)
# g4degreetable=table(g4degree)
# g5degreetable=table(g5degree)
# 
# y=as.numeric(rownames(g1degreetable))
# y=log(y)
# a=lm(y~log(g1degreetable))
# summary(a)
# 
# edge_density(g1)
# diameter(g1)
# max(components(g1)$csize)
# 
# pagerankvalue2=data.frame()
# names=row.names(genefinal)
# row.names(pagerankvalue2)=names
# 
# for(i in 1:2305){
#   x=fdrpvalue[-i,-i]
#   library(igraph)
#   gz=graph_from_adjacency_matrix(x,mode =  "undirected")
#   pagerankvalue2[i,1]=edge_density(gz)
#   pagerankvalue2[i,2]=diameter(gz)
#   pagerankvalue2[i,3]=max(components(gz)$csize)
#   gz=0
#   x=0
# }
# 
# pagerankvalue3=data.frame()
# names=row.names(genefinal)
# row.names(pagerankvalue3)=names
# 
# 
# pr=page.rank(g5)$vector
# df <- data.frame(Object = 1:2305, PageRank = pr)
# pr1=page.rank(g1)$vector
# df1 <- data.frame(Object = 1:2305, PageRank = pr1)
# names=row.names(df)
# gene.order=data.frame(F1=genefinal[,1][order(df[,2],decreasing=TRUE )],F2=genefinal[,2][order(df[,2],decreasing=TRUE )],F3=genefinal[,3][order(df[,2],decreasing=TRUE )],E1=genefinal[,4][order(df[,2],decreasing=TRUE )],E2=genefinal[,5][order(df[,2],decreasing=TRUE )])
# library(WGCNA)
# orderpvalue=corPvalueStudent(cor(t(gene.order),method=("kendall")),n=5)
# fdrorderpvalue=p.adjust(orderpvalue,method="fdr")
# fdrorderpvalue=matrix(fdrorderpvalue,2305,2305)
# row.names(fdrorderpvalue)=names
# fdrorderpvalue[as.numeric(fdrorderpvalue)>=0.05]=2
# fdrorderpvalue[as.numeric(fdrorderpvalue)<0.05]=1
# fdrorderpvalue[fdrorderpvalue==2]=0
# for(i in 1:2305){fdrorderpvalue[i,i]=0}
# 
# pagerankvalue7=data.frame()
# library(igraph)
# x=fdrorderpvalue
# for(i in 1:2305){
#   x=x[-1,-1]
#   g7=graph_from_adjacency_matrix(x,mode =  "undirected")
#   pagerankvalue7[i,1]=edge_density(g7)
#   pagerankvalue7[i,2]=diameter(g7)
#   pagerankvalue7[i,3]=max(components(g7)$csize)
#   g7=0
#  
# }
# 
# finalname=names[1:2160]
# for(i in 1:2160){finalname[i]=sub('^.','',finalname[i])}
# write.table(finalname,file="C:/Users/mzy/Documents/学习/CHD病理/genelist.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)



#以下是分别构建网络 合在一起 再去求pagerank 20161104
library(WGCNA)
mrnaallpvalue=corPvalueStudent(cor(t(mrnaall),method=("kendall")),n=5)
proteinallpvalue=corPvalueStudent(cor(t(proteinall),method=("kendall")),n=5)

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
uniongraph=union(graphmrna,graphprotein)
allpagerank=page.rank(uniongraph)$vector
tempdf<- data.frame(Object = 1:2313, PageRank =allpagerank)

#标准化
tempdf[,2]=as.numeric(tempdf[,2])
pagerankmax=max(tempdf[,2])
pagerankmin=min(tempdf[,2])
for(i in 1:2313){
  tempdf[i,3]=(tempdf[i,2]-pagerankmin)/(pagerankmax-pagerankmin)
}
tempdf=  tempdf[order(  tempdf[,3],decreasing=T),]


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


# for(i in 1:272){
#   selectidpr[i,1]=paste("rno:",selectidpr[i,1],sep="")
# }
# # 

#通路富集结果20161103
library(clusterProfiler)
keggresult=enrichKEGG(geneallname,organism="rno")
write.csv(summary(keggresult),file="C:/Users/mzy/Documents/学习/CHD病理/keggresult20161103.csv")
write.table(geneallname,file="C:/Users/mzy/Documents/学习/CHD病理/geneallname.txt",,quote = FALSE,col.names = FALSE,row.names = FALSE)

for (i in 1:11){mrnaboth[i,6]=((mrnaboth[i,4]+mrnaboth[i,5])/2)/((mrnaboth[i,1]+mrnaboth[i,2]+mrnaboth[i,3])/3)
mrnaboth[i,7]=abs(log2(mrnaboth[i,6]))}
for (i in 1:11){proteinboth[i,6]=((proteinboth[i,4]+proteinboth[i,5])/2)/((proteinboth[i,1]+proteinboth[i,2]+proteinboth[i,3])/3)
proteinboth[i,7]=abs(log2(proteinboth[i,6]))}
names(mrnaboth)[6]='foldchange'
names(proteinboth)[6]='foldchange'
names(mrnaboth)[7]='abs'
names(proteinboth)[7]='abs'
mrnaboth=mrnaboth[order(mrnaboth[,7],decreasing =TRUE ),]
proteinboth=proteinboth[order(proteinboth[,7],decreasing =TRUE ),]



 

library(KEGGgraph)
# tmp <- "rno05210.xml"
# retrieveKGML("05210", organism="rno", destfile=tmp,method = "internal")
# tmp1 <- "rno00380.xml"
# retrieveKGML("00380", organism="rno", destfile=tmp1,method = "internal")
# tmp2 <- "rno05212.xml"
# retrieveKGML("05212", organism="rno", destfile=tmp2,method = "internal")
filelist=list.files(path="C:/Users/mzy/Documents/xml")
for(i in 1:20){
  filelist[i]=paste("C:/Users/mzy/Documents/xml/",filelist[i],sep="")
}
for( i in 1:20){
  kegggraphg[[i]]=parseKGML2Graph(filelist[i],expandGenes=TRUE)
}

mergedKEGG=mergeKEGGgraphs(kegggraphg)
geneid=translateKEGGID2GeneID(nodes(mergedKEGG)) 

# selectgenename=row.names(finalgene)
# for(i in 1:2327){selectgenename[i]=sub('^.','',selectgenename[i])}
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

selectkeggid=translateGeneID2KEGGID(selectid,organism="rno")
finalmergedgraph=subGraph(selectkeggid,mergedKEGG)

fgraph=igraph.from.graphNEL(finalmergedgraph)
finalgraph=fgraph-names(degree(fgraph)[degree(fgraph)==0])
aa=edge.betweenness.community(finalgraph)
bb=walktrap.community(finalgraph)
cc=leading.eigenvector.community(finalgraph)

#模块分析 打分
aascore=data.frame()

for(i in 1:28){
  print(i)
  ccc=c()
  
  for(j in 1:length(cc[[i]])){
    x=cc[[i]][j]
    x=sub('^.','',x)
    x=sub('^.','',x)
    x=sub('^.','',x)
    x=sub('^.','',x)
    ccc=c(ccc,x)}
  tempallname=c(ccc,chongfugene)
  tempfre=as.data.frame(table(tempallname))
  ndiff=0
  for(h in 1:length(ccc)){
    if(tempfre[h,2]==1){
      ndiff=ndiff+1
    }
  }
  missvalue=0
  missvalue=ndiff/(length(ccc)*278)
  tempp=data.frame()
  phit=0
  pmiss=0
  pscore=0
  score=c()

  for(k in 1:278){
    if(selectidpr[k,1]%in%ccc){
      phit=phit+selectidpr[i,3]
      
    }
    if(selectidpr[k,1]%in%ccc==F){
      pmiss=pmiss+missvalue
    }
    pscore=phit-pmiss
    p=c(phit,pmiss,pscore)
    tempp=rbind(tempp,p)
  }
  
  score=c(paste("aa", i, sep = ""),max(tempp[,3]))
  aascore=rbind(aascore,score,stringsAsFactors = F )
  
}
#aascore即为各模块打分


# for(i in 1:29){
#   ccc=c()
#   for(j in 1:length(aa[[i]]))
#   {x=aa[[i]][j]
#    x=sub('^.','',x)
#    x=sub('^.','',x)
#    x=sub('^.','',x)
#    x=sub('^.','',x)
#    ccc=c(ccc,x)}
#   ccresult=enrichKEGG(ccc,organism="rno")
#   write.csv(summary(ccresult),file=paste("C:/Users/mzy/Documents/",i,".csv"),quote = TRUE)
#   ccresult=0
# }