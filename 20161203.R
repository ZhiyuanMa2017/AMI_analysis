library(NOISeq)
data(Marioni)
myfactors = data.frame(Tissue = c("Kidney", "Liver", "Kidney", "Liver", "Liver", "Kidney", "Liver", "Kidney", "Liver", "Kidney"), TissueRun = c("Kidney_1","Liver_1", "Kidney_1", "Liver_1", "Liver_1", "Kidney_1", "Liver_1","Kidney_2", "Liver_2", "Kidney_2"), Run = c(rep("R1", 7), rep("R2",3)))
setwd("C:/Users/mzy/Documents/学习/CHD病理/CHD病理/新建文件夹")
a = list.files()
mrna = read.table(a[1],header = T,row.names=1)
mfactors = data.frame(Tissue = c("jiashoushu", "jiashoushu", "jiashoushu", "model", "model"))
mdata <- readData(data = mrna, factors = mfactors)
mnoiseq = noiseq(mdata, k = 0.5, norm = "rpkm", factor = "Tissue", pnr = 0.2,nss = 5, v = 0.02, lc = 1, replicates = "technical")
mynoiseq.deg = degenes(mnoiseq, q = 0.8, M = NULL)
for (i in 1:2327){
  for (j in 1:624){
    if (featuregene[i,1]==row.names(mynoiseq.deg)[j]){
    featuregene[i,15]=mynoiseq.deg[j,5]
    featuregene[i,16]=1-mynoiseq.deg[j,5]
  }
  }
}
write.xlsx(featuregene, "chongfu12032158.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
# chongfu=c()
# for (i in 1:2327){
#   if (featuregene[i,1]%in%row.names(mynoiseq.deg)){
#     chongfu=c(chongfu,featuregene[i,1])
#   }
# }
# chongfugene=data.frame(geneid=chongfu)
# for (i in 1:380){
#   for (j in 1:624){
#   if (chongfugene[i,1]==row.names(mynoiseq.deg)[j]){
#     chongfugene[i,2]=mynoiseq.deg[j,5]
#     chongfugene[i,3]=1-mynoiseq.deg[j,5]
#   }
#   }
# }

