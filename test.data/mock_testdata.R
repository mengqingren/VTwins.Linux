setwd("K:/CRC-Pair/sc2meta-master/sc2meta-master")



####### mock dataset
set.seed(12345)
dataset <- data.frame(matrix(runif(1200, min = 1e-5, max = 1),nrow = 120,ncol = 10))
colnames(dataset) <- paste("Feature",1:10,sep = '')
rownames(dataset) <- paste("Sample",1:120,sep = '')
dataset.normalized <- decostand(dataset,method = "total",1)
write.table(dataset.normalized,file = "test.data.txt",sep = '\t',quote = F)
####### mock phenodata
phe_data <- data.frame(id = paste("Sample",1:120,sep = ''),grp=rep(c("grp1","grp2"),c(60,60)))
write.table(phe_data,file = "test.phenodata.txt",sep = '\t',quote = F,row.names = F)



data <- read.table("test.data.txt",header = T,row.names = 1,sep = '\t')
phe_data <- read.table("test.phenodata.txt",header = T,sep = "\t")


pair_find(data=data,
          phenodata=phe_data)







