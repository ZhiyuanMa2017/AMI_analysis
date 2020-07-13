library(ggplot2)
library(readxl)
GO =  read.csv("juleifenxi.csv",header=TRUE)  
GO = GO[0:30,]
p=ggplot(GO,aes(Fold.Enrichment,Term))
pbubble = p+ geom_point(aes(shape = Category,size=Count,color=PValue))
pr = pbubble+scale_color_gradient(low="red",high = "green")
pr = pr+labs(color=expression(PValue),size="Count",x="Fold Enrichment",y='',title="Top 30 of GO Enrichment")
pr + theme_bw()



kegg = read.csv("kegg.csv",header=TRUE) 
pk=ggplot(kegg,aes(Fold.Enrichment,Term))
pbubblek = pk + geom_point(aes(size=Count,color=PValue))
prk = pbubblek+labs(color=expression(PValue),size="Count",x="Fold Enrichment",y='',title="Enriched KEGG Pathway")
prk + theme_bw()