setwd("K:/CRC-Pair/")

library(tidyverse)
library(readxl)
library(ggpubr)
library(ggstar)
library(fdrtool)
library(qvalue)
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(ROCR)
library("caTools")
library(sampling)
library(scales)
library(ggpmisc)
library(ggrepel)
library(UpSetR)
library(gridExtra)
library(gridGraphics)
library(pheatmap)
library(ggthemes)
library(ggsci)
library(ComplexHeatmap)
library(ggrepel)
library(randomcoloR)


##################################  Unique Samples  ->  permutation 0.2-0.8 #############################################
# Update function #### Shuffle => Pair => Recall the shuffle Pair =>  more one interation
pair_find<-function(data=data,phenodata=data.frame(),k="euclidean",SavePath = NULL,ShuffleWstat = NULL, BoundarySample = NULL,BoundaryPair=NULL,ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8){ # colnames(phenodata) = c("id","grp"); grp1 = "Ctrl", grp2 = "Disease"
  suppressMessages(library(tidyverse))
  #suppressMessages(library(fdrtool))
  #suppressMessages(library(qvalue))
  phenodata <- phenodata %>% dplyr::arrange(grp) %>% data.frame()
  #data must be ordered as grp order
  groupnum <- table(phenodata$grp) %>% data.frame()
  RAN_num = groupnum[1,2]
  RAP_num = groupnum[2,2]
  RAN.Samples<-phenodata %>% dplyr::filter(grp=="grp1")
  RAP.Samples<-phenodata %>% dplyr::filter(grp=="grp2")
  RAN <- data[rownames(data) %in% RAN.Samples$id,] %>% data.frame()
  RAP <- data[rownames(data) %in% RAP.Samples$id,] %>% data.frame()
  
  n=dim(data)[1]
  num=floor(sqrt(RAN_num+RAP_num))
  results_whole=matrix(nrow=0,ncol=(num+1))
  for (i in 1:n) {
    new<-rbind(data[i,],data)
    new_dis<-dist(new,method = k)
    new_dis<-as.matrix(new_dis)
    wa_d<-sort(new_dis[,1])[3:(num+2)]
    new_row<-c(rownames(data)[i],wa_d)
    results_whole<-rbind(results_whole,new_row)
  }
  mean_whole=mean(as.numeric(as.vector(results_whole[,2:(num+1)])))
  sd_whole=sd(as.numeric(as.vector(results_whole[,2:(num+1)])))
  
  RAP_knn_num=floor(sqrt(RAP_num))
  RAN_knn_num=floor(sqrt(RAN_num))
  
  results_sample_pair<-matrix(nrow=0,ncol=(2*RAN_knn_num+2*RAP_knn_num+2))
  
  for (j in 1:(RAN_num+RAP_num)){
    sample_mean<-mean(as.numeric(results_whole[j,2:(num+1)]))
    if (sample_mean <= (mean_whole + sd_whole)){
      #print(results_whole[j,1])
      if (results_whole[j,1] %in% rownames(RAN)){
        new_RANN<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAN)
        #print(results_whole[j,1])
        new_RANN_dis<-dist(new_RANN,method = k)
        new_RANN_dis<-as.matrix(new_RANN_dis)
        wa_RANN_dis<-sort(new_RANN_dis[,1])[3:(RAN_knn_num+2)]
        knn_RANN_mean<-mean(as.numeric(wa_RANN_dis))
        knn_RANN_sd <- sd(as.numeric(wa_RANN_dis))
        
        new_RANP<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAP)
        new_RANP_dis<-dist(new_RANP,method = k)
        new_RANP_dis<-as.matrix(new_RANP_dis)
        wa_RANP_dis<-sort(new_RANP_dis[,1])[2:(RAP_knn_num+1)]
        knn_RANP_mean<-mean(as.numeric(wa_RANP_dis))
        knn_RANP_sd<-sd(as.numeric(wa_RANP_dis))
        
        name_NN<-c(rep(0,RAN_knn_num))
        name_NP<-c(rep(0,RAP_knn_num))
        
        if(knn_RANP_mean <= (knn_RANN_mean+knn_RANN_sd)){
          for (m in 1:RAN_knn_num){
            RANN_sample<-rownames(new_RANN_dis)[which(new_RANN_dis[,1] == wa_RANN_dis[m])]
            name_NN[m]<-RANN_sample
          }
          for (h in 1:RAP_knn_num){
            RANP_sample<-rownames(new_RANP_dis)[which(new_RANP_dis[,1] == wa_RANP_dis[h])]
            name_NP[h]<-RANP_sample
          }
          new.row<-c(as.character(results_whole[j,1]),rownames(data[which(rownames(data) == results_whole[j,1]),]),name_NN, as.vector(wa_RANN_dis), name_NP, as.vector(wa_RANP_dis))
          #print(new.row)
          results_sample_pair<-rbind(results_sample_pair,new.row)
        }
      } else {
        new_RAPN<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAN)
        new_RAPN_dis<-dist(new_RAPN,method = k)
        new_RAPN_dis<-as.matrix(new_RAPN_dis)
        wa_RAPN_dis<-sort(new_RAPN_dis[,1])[2:(RAN_knn_num+1)]
        knn_RAPN_mean<-mean(wa_RAPN_dis)
        
        new_RAPP<-rbind(data[which(rownames(data) == results_whole[j,1]),],RAP)
        new_RAPP_dis<-dist(new_RAPP,method = k)
        new_RAPP_dis<-as.matrix(new_RAPP_dis)
        wa_RAPP_dis<-sort(new_RAPP_dis[,1])[3:(RAP_knn_num+2)]
        knn_RAPP_mean<-mean(wa_RAPP_dis)
        knn_RAPP_sd<-sd(wa_RAPP_dis)
        
        name_PN<-c(rep(0,RAN_knn_num))
        name_PP<-c(rep(0,RAP_knn_num))
        
        if (knn_RAPN_mean <= (knn_RAPP_mean+knn_RAPP_sd)){
          for (m in 1:RAN_knn_num){
            RAPN_sample<-rownames(new_RAPN_dis)[which(new_RAPN_dis[,1] == wa_RAPN_dis[m])]
            name_PN[m]<-RAPN_sample
          }
          for (h in 1:RAP_knn_num){
            RAPP_sample<-rownames(new_RAPP_dis)[which(new_RAPP_dis[,1] == wa_RAPP_dis[h])]
            name_PP[h]<-RAPP_sample
          }
          new.row<-c(as.character(results_whole[j,1]),rownames(data[which(rownames(data) == results_whole[j,1]),]), name_PN, as.vector(wa_RAPN_dis), name_PP, as.vector(wa_RAPP_dis))
          #print(new.row)
          results_sample_pair<-rbind(results_sample_pair,new.row)
        }
      }
    }
  }
  cat("KNN FINISHED\n")
  
  if (dim(results_sample_pair)[1] > 1) {
    res<-results_sample_pair[,1:(2*RAN_knn_num+2*RAP_knn_num+2)]
  }else{
    cat("Nothing Done \n")
    return(NULL)
  }
  
  if (is.null(SavePath)) {
    SavePath="./"
  }else{
    if(!dir.exists(SavePath)){
      dir.create(SavePath)
    }
  }
  
  #res<-res[,c(2:(2+RAN_knn_num),(2+1+2*RAN_knn_num):(2+2*RAN_knn_num+RAP_knn_num))] %>% data.frame() %>% remove_rownames()
  res<-res[,c(2:(2+2*RAN_knn_num+2*RAP_knn_num))] %>% data.frame() %>% remove_rownames()
  if (!is.null(BoundarySample)) {
    write.csv(res,paste(SavePath,"/",BoundarySample,".csv",sep = ''),row.names = F)
  }
  #View(res)
  #return(res)
  #sum(str_detect(res$X1,groupnum$Var1[1]%>% as.vector()))
  #View(res)
  #write.csv(res,file = filename)
  pairinfor = matrix(ncol =3,nrow = 0)
  for (i in 1:dim(res)[1]) {
    indexpo1 = which(phenodata[,1] == res[i,1])
    indexpo2 = which(phenodata[,1] == res[i,2])
    if (phenodata[indexpo1,2] == phenodata[indexpo2,2]) {
      for (mkl in 1:RAP_knn_num) {
        pairinfor <- rbind(pairinfor,c(res[i,1],res[i,1+2*RAN_knn_num+mkl],res[i,1+2*RAN_knn_num+mkl+RAP_knn_num]))
      }
    }else{
      for (pkl in 1:RAN_knn_num) {
        pairinfor <- rbind(pairinfor,c(res[i,pkl+1],res[i,1],res[i,pkl+1+RAN_knn_num]))
      }
    }
  }
  
  Newmatrix <<- matrix(ncol = 3,nrow = 0)
  Extract_Dist <- function(PairData,...){
    if (dim(PairData)[1] > 0) {
      Findline <- PairData[which(PairData[,3] == min(PairData[,3])),]
      Newmatrix <<- rbind(Newmatrix,Findline)
      #View(Newmatrix)
      NewMidData <- PairData %>% filter(!(Ctl %in% Findline$Ctl) & !(Disease %in% Findline$Disease))
      #cat(dim(NewMidData)[1])
      #cat("\n")
      return(Extract_Dist(NewMidData))
    }else{
      return(Newmatrix)
    }
  }
  
  #View(pairinfor)
  pairinfor <- pairinfor %>% data.frame() %>% dplyr::rename(Ctl=1,Disease=2,Distance=3) %>% 
    mutate(Distance = as.numeric(as.character(Distance))) %>% dplyr::arrange(Ctl,Disease) %>% unique()
  #View(pairinfor)
  pairinfor <- Extract_Dist(pairinfor)
  pairinfor <- pairinfor %>% data.frame() %>% dplyr::select(-Distance)
  
  if (!is.null(BoundaryPair)) {
    write.csv(pairinfor,paste(SavePath,"/",BoundaryPair,".csv",sep = ''),row.names = F)
  }
  #write.csv(pairinfor,file = "test.csv",row.names = F)
  cat(paste("the redundant pair number is ",dim(pairinfor)[1],"\n",sep = ''))
  cat("PAIR FINISHED\n")
  #return(pairinfor)
  
  ###### Pair shuffle
  
  ShufflePair <- lapply(1:ShuffleTime, function(Shuffle){
    Random <- runif(1,DownPercent,Uppercent)
    n <- round(dim(pairinfor)[1]*Random,0)
    RandomIndex <- sample(1:dim(pairinfor)[1],n,replace = F)
    
    pairinforMiddata.mid1.1 <- pairinfor[-RandomIndex,] %>% data.frame()
    pairinforMiddata.mid2 <- pairinfor[RandomIndex,] %>% data.frame()
    pairinforMiddata.mid2.1 <- data.frame(Ctl = pairinforMiddata.mid2$Disease, Disease = pairinforMiddata.mid2$Ctl)
    pairinforMiddata.mid <- rbind(pairinforMiddata.mid1.1,pairinforMiddata.mid2.1) %>% data.frame()
  })
  
  cat("SHUFFLE PAIR FINISHED\n")
  
  ####### calculate W stat and rank
  FinalMatrix <- matrix(ncol = 2+ShuffleTime,nrow = 0) #Species, Original W, 
  MeanData <- matrix(ncol = 3,nrow = 0) # record mean value of pair
  for (i in 1:dim(data)[2]) {
    ### original cohort
    former<-rep(NA,dim(pairinfor)[1])
    latter<-rep(NA,dim(pairinfor)[1])
    for (j in 1:dim(pairinfor)[1]) {
      index_former<-which(rownames(data) == pairinfor$Ctl[j])
      former[j]=as.numeric(as.character(data[index_former,i]))
      index_latter<-which(rownames(data) == pairinfor$Disease[j])
      latter[j]=as.numeric(as.character(data[index_latter,i]))
    }
    Middata <- data.frame(Ctrl = former, CRC = latter)
    test<-wilcox.test(Middata$Ctrl,Middata$CRC,paired = TRUE)
    AllNW <- (dim(pairinfor)[1]*(dim(pairinfor)[1]+1))/2
    OriginalStat = AllNW - test$statistic
    
    Ctrlmean <- mean(Middata$Ctrl)
    CRCmean <- mean(Middata$CRC)
    MeanData <- rbind(MeanData,c(as.character(colnames(data)[i]),Ctrlmean,CRCmean))
    
    cat(as.character(colnames(data)[i]))
    cat("\n")
    #print(AllNW)
    
    ### calculate shuffle W stat
    #ShuffleStat <- rep(NA,ShuffleTime)
    
    ShuffleStat <- lapply(ShufflePair, function(x){
      former<-rep(NA,dim(x)[1])
      latter<-rep(NA,dim(x)[1])
      for (j in 1:dim(x)[1]) {
        index_former<-which(rownames(data) == x$Ctl[j])
        former[j]=as.numeric(as.character(data[index_former,i]))
        index_latter<-which(rownames(data) == x$Disease[j])
        latter[j]=as.numeric(as.character(data[index_latter,i]))
      }
      
      Middata.mid <- data.frame(Ctrl = former, CRC = latter)
      test1<-wilcox.test(Middata.mid$Ctrl,Middata.mid$CRC,paired = TRUE)
      #ShuffleStat[Iter] = AllNW - test1$statistic
      return(AllNW - test1$statistic)
    })
    ShuffleStat <- unlist(ShuffleStat)
    FinalMatrix <- rbind(FinalMatrix,c(as.character(colnames(data)[i]),OriginalStat,ShuffleStat))
  }
  colnames(FinalMatrix) <- c("Feature","OriginalSata.W",paste("ShuffleStat.W.",1:ShuffleTime,sep = ""))
  
  if (!is.null(ShuffleWstat)) {
    write.csv(pairinfor,paste(SavePath,"/",ShuffleWstat,".csv",sep = ''),row.names = F)
  }
  
  Mid.Matrix <- matrix(ncol = 6,nrow = 0)
  for (I.index in 1:dim(FinalMatrix)[1]) {
    Increasing.Rank.Min <- rank(as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "min")[1]
    Decreasing.Rank.Min <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "min")[1]
    Increasing.Rank.Max <- rank(as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "max")[1]
    Decreasing.Rank.Max <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "max")[1]
    Increasing.Rank.Average <- rank(as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "average")[1]
    Decreasing.Rank.Average <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:(ShuffleTime+2)])),ties.method= "average")[1]
    Mid.Matrix <- rbind(Mid.Matrix,c(Increasing.Rank.Min,Decreasing.Rank.Min,Increasing.Rank.Max,Decreasing.Rank.Max,Increasing.Rank.Average,Decreasing.Rank.Average))
  }
  Mid.Matrix <- Mid.Matrix %>% data.frame() %>% remove_rownames() %>% 
    dplyr::rename(Increasing.Rank.Min=1,Decreasing.Rank.Min=2,Increasing.Rank.Max=3,Decreasing.Rank.Max=4,Increasing.Rank.Average=5,Decreasing.Rank.Average=6)
  
  Mid.Matrix$Decre.maxRank.P <- Mid.Matrix$Decreasing.Rank.Max/(ShuffleTime+1)
  Mid.Matrix$Decre.aveRank.P <- Mid.Matrix$Decreasing.Rank.Average/(ShuffleTime+1)
  Mid.Matrix$Decre.minRank.P <- Mid.Matrix$Decreasing.Rank.Min/(ShuffleTime+1)
  Mid.Matrix$Incre.maxRank.P <- Mid.Matrix$Increasing.Rank.Max/(ShuffleTime+1)
  Mid.Matrix$Incre.aveRank.P <- Mid.Matrix$Increasing.Rank.Average/(ShuffleTime+1)
  Mid.Matrix$Incre.minRank.P <- Mid.Matrix$Increasing.Rank.Min/(ShuffleTime+1)
  
  Mid.Matrix$Species <- FinalMatrix[,1] #Feature
  
  Mid.Matrix <- Mid.Matrix %>% data.frame() %>%
    dplyr::arrange(Decre.minRank.P) %>% 
    dplyr::mutate(Decre.minRank.P.FDR = p.adjust(.$Decre.minRank.P,method = "BH",n=length(.$Decre.minRank.P))) %>%
    dplyr::arrange(Decre.maxRank.P) %>% 
    dplyr::mutate(Decre.maxRank.P.FDR = p.adjust(.$Decre.maxRank.P,method = "BH",n=length(.$Decre.maxRank.P))) %>%
    dplyr::arrange(Decre.aveRank.P) %>% 
    dplyr::mutate(Decre.aveRank.P.FDR = p.adjust(.$Decre.aveRank.P,method = "BH",n=length(.$Decre.aveRank.P)))
  
  MeanData <- MeanData %>% data.frame() %>% dplyr::rename(Species=1,Ctlmean=2,Dismean=3)
  Mid.Matrix <- merge(Mid.Matrix,MeanData,by="Species") %>% data.frame()
  
  cat("All done\n")
  return(Mid.Matrix)
}# END - function: Pair_Find
####### Merge Data with Speices Abundance and Metadata ########

#### @0 Download Data ####
Study <- c("FengQ_2015","HanniganGD_2017","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014")
suffix<-c(".metaphlan_bugs_list.stool",".pathabundance_relab.stool") #gene_family Data
#metadata <- c("subjectID","study_condition","disease","disease_subtype","country","number_bases","body_site","age","age_category","gender","BMI","")
dir.create("Data")
for (i in Study) {
  for (j in suffix) {
    datanew <- eval(parse(text = paste(i,j,"()",sep = "")))
    datanew@phenoData@data %>% write.csv(file = paste("Data/",i,".metadata.csv",sep = ""),quote = F,row.names = TRUE,col.names = TRUE)
    datanew@assayData$exprs %>% write.table(file = paste("Data/",i,j,".txt",sep = ""),quote = F,row.names = TRUE,col.names = TRUE,sep = '\t')
  }
}

#### @1 Merge MetaData ####
#list.files("Old/",pattern = "[.txt]")
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

MetaData <- data.frame()
for (name in Study) {
  metadata<-read.csv(paste("Data/",name,".metadata.csv",sep = ''),row.names = 1,check.names = F)
  MetaData<-metadata %>%  rownames_to_column() %>% filter(study_condition == "control" | study_condition == "CRC") %>% 
    mutate(study_condition = factor(study_condition,levels = c("control","CRC"))) %>%
    arrange(study_condition) %>% select(rowname,study_condition) %>% mutate(Study=name) %>% 
    rbind(MetaData)
}

#MetaData[(MetaData$Study != "PRJEB27928" & !str_detect(MetaData$rowname,"\\.21\\.0")),] %>% dim()
dir.create("./Species-8Study-20201010")
dir.create("./Pathway-8Study-20201010")

## Metadata have duplicated samples in ZellerG_2014 study and PRJEB27928. So we first redundant - exclude the redundant samples in ZellerG_2014
## ZellerG_2014 contain two datasets with two countries 
redundantsample <- MetaData %>% arrange(rowname) %>% select(rowname) %>% group_by(rowname) %>% summarise(count=n()) %>%
  filter(count > 1) %>% select(rowname)
MetaData <- MetaData %>% filter(! ((rowname %in% redundantsample$rowname) & (Study == "ZellerG_2014")))


write.csv(MetaData,"Species-8Study-20201010/EightStudies-MetaData.csv",row.names = F)
write.csv(MetaData,"Pathway-8Study-20201010/EightStudies-MetaData.csv",row.names = F)

#test1 <- read.csv(paste("Data/","ZellerG_2014",".metadata.csv",sep = ''),row.names = 1,check.names = F)
#test1 %>% rownames_to_column() %>% mutate(rowname3 = str_replace_all(rowname,"-",".")) %>% filter()


#### @2 Merge Species Abundance ####
SpeciesData <- read_excel("RawListSpecies.xlsx")
SpeciesData <- SpeciesData[str_detect(SpeciesData$Species,"k__Bacteria") & str_detect(SpeciesData$Species,"s__") & !str_detect(SpeciesData$Species,"t__") ,] %>% data.frame()
SpeciesData$Species <- SpeciesData$Species %>% as.character() %>% str_replace_all(".*s__","")

for (name in Study) {
  print(name)
  Speciesdata <- read.table(paste("Data/",name,".metaphlan_bugs_list.stool.txt",sep = ''),row.names = 1,check.names = F,sep = '\t',header = T)
  print(dim(Speciesdata))
  Speciesdata <- Speciesdata[str_detect(rownames(Speciesdata),"k__Bacteria") & str_detect(rownames(Speciesdata),"s__") & !str_detect(rownames(Speciesdata),"t__") ,] %>% data.frame()
  rownames(Speciesdata) = rownames(Speciesdata) %>% as.character() %>% str_replace_all(".*s__","")
  Speciesdata <- Speciesdata %>% rownames_to_column()
  SpeciesData <-  merge(SpeciesData,Speciesdata,by.x = "Species",by.y = "rowname",all = T)
  print(dim(SpeciesData))
}

SpeciesData[is.na(SpeciesData)] = 0

write.table(SpeciesData,"Species-8Study-20201010/EightStudies-SpeciesAbundance.txt",row.names = F,sep = '\t',quote = F)

#### @3 Merge Pathway Abundance ####
PathwayData <- read_excel("RawListPathway.xlsx")
PathwayData$Pathway <- PathwayData$Pathway %>% str_replace_all("-",".") #%>% str_replace_all(" ",".") %>% str_replace_all(":",".") %>% str_replace_all("\\(",".") %>% str_replace_all("\\)",".")
for (name in Study) {
  print(name)
  if (name == "PRJEB27928") {
    Pathwaydata <- read.csv(paste("Data/",name,".pathabundance_relab.csv",sep = ''),check.names = F,row.names = 1)
  }else{
    Pathwaydata <- read.table(paste("Data/",name,".pathabundance_relab.stool.txt",sep = ''),check.names = F,sep = '\t',header = T,row.names = 1)
  }
  
  Pathwaydata <- Pathwaydata[!str_detect(rownames(Pathwaydata),"\\|"),] %>% data.frame()
  Pathwaydata <- Pathwaydata %>% data.frame() %>% rownames_to_column()
  Pathwaydata$rowname <- Pathwaydata$rowname %>% str_replace_all("-",".")
  PathwayData <- merge(PathwayData,Pathwaydata,by.x = "Pathway",by.y = "rowname",all = T)
}

PathwayData[is.na(PathwayData)] = 0

write.table(PathwayData,"Pathway-8Study-20201010/EightStudies-PathwayAbundance.txt",row.names = F,sep = '\t')

################################## 
######### All ###############

SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',stringsAsFactors = F,check.names = F)
SpeciesData2 <- SpeciesData2 %>% dplyr::arrange(study_condition)

groupnum <- table(SpeciesData2$study_condition) %>% data.frame()
middata <- SpeciesData2 %>% dplyr::arrange(study_condition) %>% dplyr::select(-study_condition,-Study)
phe_data <- SpeciesData2 %>% dplyr::arrange(study_condition) %>% dplyr::select(study_condition) %>% rownames_to_column() %>%
  dplyr::rename(id=1,grp=2) %>%
  dplyr::arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
filename = paste("Unique.BoundarySmaple.Distance","All",".csv",sep="")
res<-pair_find(data=middata,phenodata = phe_data,k="euclidean",BoundarySample="All.BoundarySample.Species",BoundaryPair="All.BoundaryPair.Species")
write.csv(res,paste("Unique.Pair.Sample.","All",".Pair.Distance.csv",sep = ''),row.names = F)

pairinfor <- res
FinalMatrix <- matrix(ncol = 20002,nrow = 0) #Species, Original W, 
for (i in 1:dim(middata)[2]) {
  former<-rep(NA,dim(pairinfor)[1])
  latter<-rep(NA,dim(pairinfor)[1])
  for (j in 1:dim(pairinfor)[1]) {
    index_former<-which(rownames(middata) == pairinfor$Ctl[j])
    former[j]=as.numeric(as.character(middata[index_former,i]))
    index_latter<-which(rownames(middata) == pairinfor$Disease[j])
    latter[j]=as.numeric(as.character(middata[index_latter,i]))
  }
  Middata <- data.frame(Ctrl = former, CRC = latter)
  test<-wilcox.test(Middata$Ctrl,Middata$CRC,paired = TRUE)
  AllNW <- (dim(res)[1]*(dim(res)[1]+1))/2
  OriginalStat = AllNW - test$statistic
  #Output = list()
  #Output[[as.character(colnames(middata1)[i])]][["Original"]] = AllNW - test$statistic
  
  #Output[[as.character(colnames(middata1)[i])]][["Shuffle"]] = rep(NA,1000)
  #Output[[as.character(colnames(middata1)[i])]][["Percent"]] = rep(NA,1000)
  ShuffleStat <- rep(NA,10000)
  PercentCount <- rep(NA,10000)
  cat(as.character(colnames(middata)[i]))
  cat("\n")
  print(AllNW)
  for (Iter in 1:10000) {
    Random <- runif(1,0.2,0.8)
    n <- round(dim(res)[1]*Random,0)
    #Output[[as.character(colnames(middata1)[i])]][["Percent"]][Iter] = n
    PercentCount[Iter] = n
    
    RandomIndex <- sample(1:dim(res)[1],n,replace = F)
    
    Middata.mid1.1 <- Middata[-RandomIndex,] %>% data.frame()
    Middata.mid2 <- Middata[RandomIndex,] %>% data.frame()
    Middata.mid2.1 <- data.frame(Ctrl = Middata.mid2$CRC, CRC = Middata.mid2$Ctrl)
    Middata.mid <- rbind(Middata.mid1.1,Middata.mid2.1) %>% data.frame()
    
    test1<-wilcox.test(Middata.mid$Ctrl,Middata.mid$CRC,paired = TRUE)
    #Output[[as.character(colnames(middata1)[i])]][["Shuffle"]][Iter] = AllNW - test1$statistic
    ShuffleStat[Iter] = AllNW - test1$statistic
    cat(test1$statistic)
    cat("\n")
  }
  FinalMatrix <- rbind(FinalMatrix,c(as.character(colnames(middata)[i]),OriginalStat,ShuffleStat,PercentCount))
  cat(as.character(colnames(middata)[i]))
  cat("\n")
}

colnames(FinalMatrix) <- c("Feature","OriginalSata.W",paste("ShuffleStat.W.",1:10000,sep = ""),paste("PercentCount.Overturn.",1:10000,sep = ""))
Mid.Matrix <- matrix(ncol = 2,nrow = 0)
for (I.index in 1:dim(FinalMatrix)[1]) {
  Increasing.Rank <- rank(as.numeric(as.character(FinalMatrix[I.index,2:10002])),ties.method= "min")[1]
  Decreasing.Rank <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:10002])),ties.method= "min")[1]
  Mid.Matrix <- rbind(Mid.Matrix,c(Increasing.Rank,Decreasing.Rank))
}
colnames(Mid.Matrix) <- c("Increasing.Rank","Decreasing.Rank")
Final.Matrix <- cbind(Mid.Matrix,FinalMatrix)
write.csv(Final.Matrix,paste("Unique.Pair.Sample","All",".0.2.0.8.10000.Feature.Shuffle.Rank.csv",sep = ''),row.names = F)

RankData <- read.csv(paste("Unique.Pair.Sample","All",".0.2.0.8.10000.Feature.Shuffle.Rank.csv",sep = ''))
RandData.Pvalue <- data.frame(Increasing.Rank.Min = RankData$Increasing.Rank,Decreasing.Rank.Min = RankData$Decreasing.Rank)
RandData.Pvalue$Decre.minRank.P <- RandData.Pvalue$Decreasing.Rank.Min/10001

Decreasing.Rank.Average <- rep(NA,dim(RankData)[1])
Increasing.Rank.Average <- rep(NA,dim(RankData)[1])
Decreasing.Rank.Max <- rep(NA,dim(RankData)[1])
Increasing.Rank.Max <- rep(NA,dim(RankData)[1])
for (I.index in 1:dim(RankData)[1]) {
  Decreasing.Rank.Average[I.index] <- rank(-as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "average")[1]
  Increasing.Rank.Average[I.index] <- rank(as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "average")[1]
  Decreasing.Rank.Max[I.index] <- rank(-as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "max")[1]
  Increasing.Rank.Max[I.index] <- rank(as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "max")[1]
}
RandData.Pvalue$Decreasing.Rank.Average <- Decreasing.Rank.Average
RandData.Pvalue$Increasing.Rank.Average <- Increasing.Rank.Average
RandData.Pvalue$Decreasing.Rank.Max <- Decreasing.Rank.Max
RandData.Pvalue$Increasing.Rank.Max <- Increasing.Rank.Max
RandData.Pvalue$Decre.maxRank.P <- RandData.Pvalue$Decreasing.Rank.Max/10001
RandData.Pvalue$Decre.aveRank.P <- Decreasing.Rank.Average/10001
RandData.Pvalue$OriginalStat <- RankData$OriginalSata.W
RandData.Pvalue$Sig <- if_else(RandData.Pvalue$Decreasing.Rank.Average == RandData.Pvalue$Increasing.Rank.Average,"FALSE",if_else(RandData.Pvalue$Decre.aveRank.P <= 0.5,"TRUE","FALSE"))
RandData.Pvalue$Species <- RankData$Feature
write.csv(RandData.Pvalue,paste("Unique.Pair.Sample","All",".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''),row.names = F)

RandData <- read.csv(paste("Unique.Pair.Sample","All",".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''))
RandData <- RandData %>% 
  arrange(Decre.minRank.P) %>% 
  mutate(Decre.minRank.P.FDR = p.adjust(.$Decre.minRank.P,method = "BH",n=length(.$Decre.minRank.P))) %>%
  arrange(Decre.maxRank.P) %>% 
  mutate(Decre.maxRank.P.FDR = p.adjust(.$Decre.maxRank.P,method = "BH",n=length(.$Decre.maxRank.P))) %>%
  arrange(Decre.aveRank.P) %>% 
  mutate(Decre.aveRank.P.FDR = p.adjust(.$Decre.aveRank.P,method = "BH",n=length(.$Decre.aveRank.P)))

write.csv(RandData,paste("Unique.Pair.Sample","All",".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''),row.names = F)


######### Every Study #########
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
Study <- c("ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928") #"FengQ_2015","ThomasAM_2018a",
for (name in Study) {
  middata2 <- SpeciesData2 %>% filter(Study == name) %>% arrange(study_condition)
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  phe_data <- middata2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
    dplyr::rename(id=1,grp=2) %>%
    arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
  groupnum <- table(phe_data$grp) %>% data.frame()
  filename = paste("Unique.BoundarySmaple.Distance",name,".csv",sep="")
  res<-pair_find(data=middata,phe_data,k="euclidean")
  
  write.csv(res,paste("Unique.Pair.Sample.",name,".Pair.Distance.csv",sep = ''),row.names = F)
  
  pairinfor <- res
  FinalMatrix <- matrix(ncol = 20002,nrow = 0) #Species, Original W, 
  for (i in 1:dim(middata)[2]) {
    former<-rep(NA,dim(pairinfor)[1])
    latter<-rep(NA,dim(pairinfor)[1])
    for (j in 1:dim(pairinfor)[1]) {
      index_former<-which(rownames(middata) == pairinfor$Ctl[j])
      former[j]=as.numeric(as.character(middata[index_former,i]))
      index_latter<-which(rownames(middata) == pairinfor$Disease[j])
      latter[j]=as.numeric(as.character(middata[index_latter,i]))
    }
    Middata <- data.frame(Ctrl = former, CRC = latter)
    test<-wilcox.test(Middata$Ctrl,Middata$CRC,paired = TRUE)
    AllNW <- (dim(res)[1]*(dim(res)[1]+1))/2
    OriginalStat = AllNW - test$statistic
    #Output = list()
    #Output[[as.character(colnames(middata1)[i])]][["Original"]] = AllNW - test$statistic
    
    #Output[[as.character(colnames(middata1)[i])]][["Shuffle"]] = rep(NA,1000)
    #Output[[as.character(colnames(middata1)[i])]][["Percent"]] = rep(NA,1000)
    ShuffleStat <- rep(NA,10000)
    PercentCount <- rep(NA,10000)
    cat(as.character(colnames(middata)[i]))
    cat("\n")
    print(AllNW)
    for (Iter in 1:10000) {
      Random <- runif(1,0.2,0.8)
      n <- round(dim(res)[1]*Random,0)
      #Output[[as.character(colnames(middata1)[i])]][["Percent"]][Iter] = n
      PercentCount[Iter] = n
      
      RandomIndex <- sample(1:dim(res)[1],n,replace = F)
      
      Middata.mid1.1 <- Middata[-RandomIndex,] %>% data.frame()
      Middata.mid2 <- Middata[RandomIndex,] %>% data.frame()
      Middata.mid2.1 <- data.frame(Ctrl = Middata.mid2$CRC, CRC = Middata.mid2$Ctrl)
      Middata.mid <- rbind(Middata.mid1.1,Middata.mid2.1) %>% data.frame()
      
      test1<-wilcox.test(Middata.mid$Ctrl,Middata.mid$CRC,paired = TRUE)
      #Output[[as.character(colnames(middata1)[i])]][["Shuffle"]][Iter] = AllNW - test1$statistic
      ShuffleStat[Iter] = AllNW - test1$statistic
      cat(test1$statistic)
      cat("\n")
    }
    FinalMatrix <- rbind(FinalMatrix,c(as.character(colnames(middata)[i]),OriginalStat,ShuffleStat,PercentCount))
    cat(as.character(colnames(middata)[i]))
    cat("\n")
  }
  
  colnames(FinalMatrix) <- c("Feature","OriginalSata.W",paste("ShuffleStat.W.",1:10000,sep = ""),paste("PercentCount.Overturn.",1:10000,sep = ""))
  Mid.Matrix <- matrix(ncol = 2,nrow = 0)
  for (I.index in 1:dim(FinalMatrix)[1]) {
    Increasing.Rank <- rank(as.numeric(as.character(FinalMatrix[I.index,2:10002])),ties.method= "min")[1]
    Decreasing.Rank <- rank(-as.numeric(as.character(FinalMatrix[I.index,2:10002])),ties.method= "min")[1]
    Mid.Matrix <- rbind(Mid.Matrix,c(Increasing.Rank,Decreasing.Rank))
  }
  colnames(Mid.Matrix) <- c("Increasing.Rank","Decreasing.Rank")
  Final.Matrix <- cbind(Mid.Matrix,FinalMatrix)
  write.csv(Final.Matrix,paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.csv",sep = ''),row.names = F)
}

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928") #
RandData.Pvalue <- data.frame()
for (name in Study){
  RankData <- read.csv(paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.csv",sep = ''))
  RandData.Pvalue <- data.frame(Increasing.Rank.Min = RankData$Increasing.Rank,Decreasing.Rank.Min = RankData$Decreasing.Rank)
  RandData.Pvalue$Decre.minRank.P <- RandData.Pvalue$Decreasing.Rank.Min/10001
  
  Decreasing.Rank.Average <- rep(NA,dim(RankData)[1])
  Increasing.Rank.Average <- rep(NA,dim(RankData)[1])
  Decreasing.Rank.Max <- rep(NA,dim(RankData)[1])
  Increasing.Rank.Max <- rep(NA,dim(RankData)[1])
  for (I.index in 1:dim(RankData)[1]) {
    Decreasing.Rank.Average[I.index] <- rank(-as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "average")[1]
    Increasing.Rank.Average[I.index] <- rank(as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "average")[1]
    Decreasing.Rank.Max[I.index] <- rank(-as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "max")[1]
    Increasing.Rank.Max[I.index] <- rank(as.numeric(as.character(RankData[I.index,4:10004])),ties.method= "max")[1]
  }
  RandData.Pvalue$Decreasing.Rank.Average <- Decreasing.Rank.Average
  RandData.Pvalue$Increasing.Rank.Average <- Increasing.Rank.Average
  RandData.Pvalue$Decreasing.Rank.Max <- Decreasing.Rank.Max
  RandData.Pvalue$Increasing.Rank.Max <- Increasing.Rank.Max
  RandData.Pvalue$Decre.maxRank.P <- RandData.Pvalue$Decreasing.Rank.Max/10001
  RandData.Pvalue$Decre.aveRank.P <- Decreasing.Rank.Average/10001
  RandData.Pvalue$OriginalStat <- RankData$OriginalSata.W
  RandData.Pvalue$Sig <- if_else(RandData.Pvalue$Decreasing.Rank.Average == RandData.Pvalue$Increasing.Rank.Average,"FALSE",if_else(RandData.Pvalue$Decre.aveRank.P <= 0.5,"TRUE","FALSE"))
  RandData.Pvalue$Species <- RankData$Feature
  write.csv(RandData.Pvalue,paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''),row.names = F)
  
  RandData <- read.csv(paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''))
  
  RandData <- RandData %>% 
    arrange(Decre.minRank.P) %>% 
    mutate(Decre.minRank.P.FDR = p.adjust(.$Decre.minRank.P,method = "BH",n=length(.$Decre.minRank.P))) %>%
    arrange(Decre.maxRank.P) %>% 
    mutate(Decre.maxRank.P.FDR = p.adjust(.$Decre.maxRank.P,method = "BH",n=length(.$Decre.maxRank.P))) %>%
    arrange(Decre.aveRank.P) %>% 
    mutate(Decre.aveRank.P.FDR = p.adjust(.$Decre.aveRank.P,method = "BH",n=length(.$Decre.aveRank.P)))
  
  write.csv(RandData,paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''),row.names = F)
}


################################## Study Transfer Study  => Top 10,20,30,40,50 ##################################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
dir.create("StudyToStudy")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

### Average Rank -> RF model
for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    middata <- read.csv(paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''))
    middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
    mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
    middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-kcount,Pvalue)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## train test data
    TrainData <- SpeciesSelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self 0.6 => All
    set.seed(123)
    split = sample.split(TrainData$Label,SplitRatio = .6)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=2,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, TrainData)
    cm<-confusionMatrix(rf.pred,TrainData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,BaseModel = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,TrainData,type = "prob")
    pred<-prediction(predob[,2],TrainData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = name,BaseModel = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name,BaseModel = name) %>% rbind(ModelImportance)
    
    rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
    
    ## model for Other All Studies
    for (name2 in Study) {
      if (name != name2) {
        cat(paste(name,"\t",name2,"\n"))
        
        test.data <- SpeciesSelectData %>% filter(Study == name2) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
          select(-study_condition,-Study)
        
        control <- trainControl(method="repeatedcv",number=3,repeats=5)
        
        test.data$Label <- factor(test.data$Label,levels = c(0,1))
        
        rf.pred <- predict(rf.fit, test.data)
        cm<-confusionMatrix(rf.pred,test.data$Label)
        confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name2,BaseModel = name) %>%rbind(confusionMatrixdata)
        
        predob = predict(rf.fit,test.data,type = "prob")
        pred<-prediction(predob[,2],test.data$Label)
        perf<-performance(pred,'tpr','fpr')
        
        #extrac plot ROC data
        modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name2,BaseModel = name) %>%
          rbind(modelROC)
        
        auc<-performance(pred,"auc")
        auc<-unlist(slot(auc,"y.values"))
        model.AUC <- data.frame(Predict=name2,AUC=auc,BaseModel = name) %>% rbind(model.AUC)
        
        ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
          mutate(Rank = floor(rank(-MeanDecreaseGini)),StudyModel = name2,BaseModel = name) %>% rbind(ModelImportance)
      }
    }
  }
  
  write.csv(confusionMatrixdata,file = paste("StudyToStudy/Study2Study.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("StudyToStudy/Study2Study.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("StudyToStudy/Study2Study.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("StudyToStudy/Study2Study.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}

################################## Intra/Within Study CrossValidation => Top 10,20,30,40,50 #########################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
dir.create("IntraStudy")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")

ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()
for (kcount in seq(10,50,10)) {
  for (name in Study) {
    middata <- read.csv(paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''))
    #middataTop30 <- middata %>% arrange(Decre.aveRank.P) %>% top_n(-kcount,Decre.aveRank.P)
    middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
    mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
    middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-kcount,Pvalue)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## Inter Study CV5 
    Traindata <- SpeciesSelectData %>% filter(Study == name) %>% select(c(Study,study_condition,middataTop30$Species)) %>% 
      arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1))
    Traindata$Label <- factor(Traindata$Label,levels = c(0,1))
    
    for (repeatn in 1:10) {
      for (k in 1:5) {
        Folds <- createFolds(y = Traindata$study_condition,k = 5)
        
        TrainData <- Traindata[-Folds[[k]],] %>% data.frame() %>% 
          select(-study_condition,-Study)
        
        TestData <- Traindata[Folds[[k]],] %>% data.frame() %>% 
          select(-study_condition,-Study)
        
        TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
        TestData$Label <- factor(TestData$Label,levels = c(0,1))
        # cross training model
        control <- trainControl(method="repeatedcv",repeats=5)
        
        rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
        
        rf.pred <- predict(rf.fit, TestData)
        cm<-confusionMatrix(rf.pred,TestData$Label)
        confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>%rbind(confusionMatrixdata)
        
        predob = predict(rf.fit,TestData,type = "prob")
        pred<-prediction(predob[,2],TestData$Label)
        perf<-performance(pred,'tpr','fpr')
        
        #extrac plot ROC data
        modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>%
          rbind(modelROC)
        
        auc<-performance(pred,"auc")
        auc<-unlist(slot(auc,"y.values"))
        model.AUC <- data.frame(AUC=auc,Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>% rbind(model.AUC)
        
        ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
          mutate(Rank = floor(rank(-MeanDecreaseGini)),Study=name,RepeatTimes = repeatn,Fold=k,FeatureCount=kcount) %>% rbind(ModelImportance)
      }
    }
  }
}

write.csv(confusionMatrixdata,file = paste("IntraStudy/IntraStudy.CV5Top10-50",".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("IntraStudy/IntraStudy.CV5Top10-50",".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("IntraStudy/IntraStudy.CV5Top10-50",".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("IntraStudy/IntraStudy.CV5Top10-50",".RF.model.Importance.csv",sep = ''),row.names = F)


################################## CrossValidatation RF Model for Pooled Study => Top 10,20,30,40,50 #################################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% dplyr::arrange(study_condition)

ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()
for (kcount in seq(10,50,10)){
  middata <- read.csv(paste("All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",sep = ''))
  middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
  mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
  middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-kcount,Pvalue)
  #middataTop30 <- middata %>% arrange(Decre.aveRank.P) %>% top_n(-kcount,Decre.aveRank.P)
  for (i in 1:20) {
    print(i)
    for (k in 1:10) {
      Folds <- createFolds(y = SpeciesData2$study_condition,k = 10)
      #sample index for train and test data
      #Index <- sample(1:10,5,replace = F)
      #two for two
      TrainData <- SpeciesData2[-Folds[[k]],] %>% data.frame() %>% select(c(Study,study_condition,middataTop30$Species)) %>% 
        mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      TestData <- SpeciesData2[Folds[[k]],] %>% data.frame() %>% select(c(Study,study_condition,middataTop30$Species)) %>% 
        mutate(Label=if_else(study_condition =="control",0,1)) %>%
        select(-study_condition,-Study)
      
      TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
      TestData$Label <- factor(TestData$Label,levels = c(0,1))
      # cross training model
      control <- trainControl(method="repeatedcv",number=3,repeats=5)
      
      rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(rf.fit, TestData)
      cm<-confusionMatrix(rf.pred,TestData$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(RepeatTimes = i,Fold=k,FatureCount=kcount) %>%rbind(confusionMatrixdata)
      
      predob = predict(rf.fit,TestData,type = "prob")
      pred<-prediction(predob[,2],TestData$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(RepeatTimes = i,Fold=k,FatureCount=kcount) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(AUC=auc,RepeatTimes = i,Fold=k,FatureCount=kcount) %>% rbind(model.AUC)
      
      ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini)),RepeatTimes = i,Fold=k,FatureCount=kcount) %>% rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("AllStudy.CV10Top10-50",".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("AllStudy.CV10Top10-50",".RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("AllStudy.CV10Top10-50",".RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("AllStudy.CV10Top10-50",".RF.model.Importance.csv",sep = ''),row.names = F)


################################## LODO Analysis #####################################################
######## Differential Species ################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
dir.create("LODO")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("PRJEB27928") #"FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176",

for (name in Study) {
  cat(paste("LODO  ",name,"  Start\n",sep = ""))
  middata2 <- SpeciesData2 %>% filter(Study != name) %>% arrange(study_condition)
  
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  phe_data <- middata2 %>% select(study_condition) %>% rownames_to_column() %>%
    dplyr::rename(id=1,grp=2) %>%
    arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))

  groupnum <- table(phe_data$grp) %>% data.frame()
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential species
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  write.csv(wilcox_res,file = paste("LODO/LODO-Exclude.",name,"-Metaphlan2-wilcoxonTest.csv",sep = ''),row.names = F)
  
  # pair analysi wilcoxon pair test
  res<-pair_find(data=middata,phe_data,k="euclidean",SavePath = "LODO",BoundarySample=paste("LODO.BoundarySample.",name,sep = ""),BoundaryPair=paste("LODO.BoundaryPair.",name,sep = ""),ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
  write.csv(res,file = paste("LODO/LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),row.names = F)
  cat(paste("LODO  ",name,"  End\n",sep = ""))
}



######## LODO RF -> Ave Top 10-50 #############
setwd("K:/CRC-Pair/Unique.Pair.Permutation/LODO.Test")
SpeciesData2 <- read.table("../../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928") #

for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    ## Process Data
    #middata <- get(paste(name,".LODO.wilcoxonsign",sep = ''))
    middata <- read.csv(paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
    middataTop30 <- middata %>% arrange(Decre.aveRank.P) %>% top_n(-kcount,Decre.aveRank.P)
    #write.csv(middataTop50,file = paste(name,"-PairTop50.Species.csv",sep = ''),row.names = F)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## train test data
    TrainData <- SpeciesSelectData %>% filter(Study != name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self
    set.seed(123)
    split = sample.split(TrainData$Label,SplitRatio = .7)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=3,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, test_data)
    cm<-confusionMatrix(rf.pred,test_data$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,test_data,type = "prob")
    pred<-prediction(predob[,2],test_data$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict="Self",AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = "Self",ModelExcludeStudy = name) %>% rbind(ModelImportance)
    
    ## model for exclulded study
    TestData <- SpeciesSelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TestData$Label <- factor(TestData$Label,levels = c(0,1))
    
    control <- trainControl(method="repeatedcv",number=5,repeats=5)
    rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(rf.fit, TestData)
    cm<-confusionMatrix(rf.pred,TestData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(rf.fit,TestData,type = "prob")
    pred<-prediction(predob[,2],TestData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name,ModelExcludeStudy = name) %>%
      rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = name,ModelExcludeStudy = name) %>% rbind(ModelImportance)
  }
  
  write.csv(confusionMatrixdata,file = paste("LODO.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("LODO.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("LODO.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("LODO.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}


################################## Species Differential Analysis ###########
######## All -> Pair + Wolcoxon #########

setwd("K:/CRC-Pair/Unique.Pair.Permutation")

SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)

groupnum <- table(SpeciesData2$study_condition) %>% data.frame()
middata1 <- SpeciesData2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
phe_data <- SpeciesData2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
  dplyr::rename(id=1,grp=2) %>%
  arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))

#filename=paste("Species-8Study-20201010/All-8Study-Contained.Speices.pair.csv",sep = '')
res<-pair_find(data=middata1,phe_data,k="euclidean",ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
#ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
write.csv(res,file = paste("All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",sep = ''),row.names = F)

# all wilcoxon
rownames(middata1)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
wilcox_res<-matrix(ncol = 7,nrow = 0)
## wilcox test differential species
for (j in 1:dim(middata1)[2]) {
  test<-wilcox.test(middata1[,j][1:groupnum$Freq[1]],middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]])
  wilcox_res <- rbind(wilcox_res,c(colnames(middata1)[j],test$p.value,mean(middata1[,j][1:groupnum$Freq[1]]),median(middata1[,j][1:groupnum$Freq[1]]),mean(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),median(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),j))
}
wilcox_res <- as.data.frame(wilcox_res)
colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
wilcox_res <- na.omit(wilcox_res)
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
  mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()

wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
write.csv(wilcox_res,file = paste("All-8Study-Contained-Species-wilcoxon-test.csv",sep = ''),row.names = F)

####
#data=middata1;phenodata=phe_data;k="euclidean";ShuffleTime=10000;DownPercent = 0.2;Uppercent=0.8

######## Every Study -> Pair + Wolcoxon ########
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (name in Study) {
  middata2 <- SpeciesData2 %>% filter(Study == name) %>% arrange(study_condition)
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  phe_data <- middata2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
    dplyr::rename(id=1,grp=2) %>%
    arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
  groupnum <- table(phe_data$grp) %>% data.frame()
  #write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste(name,"-metaphlan2.SubjectID-Newname.csv",sep = ''),row.names = F)
  #rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential species
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  write.csv(wilcox_res,file = paste(name,"-Metaphlan2-wilcoxonTest.csv",sep = ''),row.names = F)
  assign(paste(name,".wilcoxon",sep = ''),wilcox_res)
  
  # pair analysi wilcoxon pair test
  filename=paste(name,".metaphaln2.pair.csv",sep = '')
  #res<-pair_find(data=middata,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
  res<-pair_find(data=middata,phe_data,k="euclidean")
  
  #assign(paste(name,".wilcoxonsign",sep = ''),res)
  write.csv(RandData,paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''),row.names = F)
}



#################################### LODO -> Species Permutation #############################
############# Differential Species #############
dir.create("LODO.Test")
setwd("K:/CRC-Pair/Unique.Pair.Permutation/LODO.Test")
SpeciesData2 <- read.table("../../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928") #

for (name in Study) {
  cat(paste("LODO  ",name,"  Start\n",sep = ""))
  middata2 <- SpeciesData2 %>% filter(Study != name) %>% arrange(study_condition)
  
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  phe_data <- middata2 %>% select(study_condition) %>% rownames_to_column() %>%
    dplyr::rename(id=1,grp=2) %>%
    arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
  
  groupnum <- table(phe_data$grp) %>% data.frame()
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential species
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  write.csv(wilcox_res,file = paste("LODO-Exclude.",name,"-Metaphlan2-wilcoxonTest.csv",sep = ''),row.names = F)
  
  # pair analysi wilcoxon pair test
  res<-pair_find(data=middata,phe_data,k="euclidean",ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
  write.csv(res,file = paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''),row.names = F)
  cat(paste("LODO  ",name,"  End\n",sep = ""))
}

############# LODO RF -> Ave Top 10-50 #############
setwd("K:/CRC-Pair/Unique.Pair.Permutation/LODO.Test")
SpeciesData2 <- read.table("../../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928") #

for (kcount in seq(10,50,10)) {
  ModelImportance <- data.frame()
  confusionMatrixdata <- data.frame()
  modelROC <- data.frame()
  model.AUC <- data.frame()
  
  for (name in Study) {
    ## Process Data
    #middata <- get(paste(name,".LODO.wilcoxonsign",sep = ''))
    middata <- read.csv(paste("LODO-Exclude.",name,"-Metaphlan2-PairWilcoxonSign.csv",sep = ''))
    middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
    mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
    middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-kcount,Pvalue)
    #middataTop30 <- middata %>% arrange(Decre.aveRank.P) %>% top_n(-kcount,Decre.aveRank.P)
    #write.csv(middataTop50,file = paste(name,"-PairTop50.Species.csv",sep = ''),row.names = F)
    
    SpeciesSelectData <- SpeciesData2 %>% select(c(Study,study_condition,middataTop30$Species)) 
    ## train test data
    TrainData <- SpeciesSelectData %>% filter(Study != name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
    
    ## model for self
    set.seed(123)
    split = sample.split(TrainData$Label,SplitRatio = .7)
    train_data = subset(TrainData,split == TRUE)
    test_data  = subset(TrainData,split == FALSE)
    
    TrainData$Label<-factor(TrainData$Label,levels = c(0,1))
    control <- trainControl(method="repeatedcv",number=3,repeats=5)
    
    train_data$Label <- factor(train_data$Label,levels = c(0,1))
    test_data$Label <- factor(test_data$Label,levels = c(0,1))
    
    fit.rf <- train(Label ~ .,data = train_data, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(fit.rf, test_data)
    cm<-confusionMatrix(rf.pred,test_data$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(fit.rf,test_data,type = "prob")
    pred<-prediction(predob[,2],test_data$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict = "Self",ModelExcludeStudy = name) %>%rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict="Self",AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = "Self",ModelExcludeStudy = name) %>% rbind(ModelImportance)
    
    ## model for exclulded study
    TestData <- SpeciesSelectData %>% filter(Study == name) %>% arrange(study_condition) %>% mutate(Label=if_else(study_condition =="control",0,1)) %>%
      select(-study_condition,-Study)
    TestData$Label <- factor(TestData$Label,levels = c(0,1))
    
    control <- trainControl(method="repeatedcv",number=5,repeats=5)
    rf.fit <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
    
    rf.pred <- predict(rf.fit, TestData)
    cm<-confusionMatrix(rf.pred,TestData$Label)
    confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% mutate(Predict = name,ModelExcludeStudy = name) %>%rbind(confusionMatrixdata)
    
    predob = predict(rf.fit,TestData,type = "prob")
    pred<-prediction(predob[,2],TestData$Label)
    perf<-performance(pred,'tpr','fpr')
    
    #extrac plot ROC data
    modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% mutate(Predict  = name,ModelExcludeStudy = name) %>%
      rbind(modelROC)
    
    auc<-performance(pred,"auc")
    auc<-unlist(slot(auc,"y.values"))
    model.AUC <- data.frame(Predict=name,AUC=auc,ModelExcludeStudy = name) %>% rbind(model.AUC)
    
    ModelImportance<-rf.fit$finalModel$importance %>% data.frame() %>% rownames_to_column() %>% 
      mutate(Rank = floor(rank(-MeanDecreaseGini)),Predict = name,ModelExcludeStudy = name) %>% rbind(ModelImportance)
  }
  
  write.csv(confusionMatrixdata,file = paste("LODO.Top",kcount,".RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
  write.csv(model.AUC,file = paste("LODO.Top",kcount,".RF.model.AUC.csv",sep = ''),row.names = F)
  write.csv(modelROC,file = paste("LODO.Top",kcount,".RF.model.ROC.csv",sep = ''),row.names = F)
  write.csv(ModelImportance,file = paste("LODO.Top",kcount,".RF.model.Importance.csv",sep = ''),row.names = F)
}


#################################### Species Differential Analysis ###########
######## All -> Pair + Wolcoxon #########

setwd("K:/CRC-Pair/Unique.Pair.Permutation")

SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)

groupnum <- table(SpeciesData2$study_condition) %>% data.frame()
middata1 <- SpeciesData2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
phe_data <- SpeciesData2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
  dplyr::rename(id=1,grp=2) %>%
  arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))

#filename=paste("Species-8Study-20201010/All-8Study-Contained.Speices.pair.csv",sep = '')
res<-pair_find(data=middata1,phe_data,k="euclidean",ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
#ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
write.csv(res,file = paste("All-8Study-Contained-Species-pair-wilcoxonsign-res.csv",sep = ''),row.names = F)

# all wilcoxon
rownames(middata1)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
wilcox_res<-matrix(ncol = 7,nrow = 0)
## wilcox test differential species
for (j in 1:dim(middata1)[2]) {
  test<-wilcox.test(middata1[,j][1:groupnum$Freq[1]],middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]])
  wilcox_res <- rbind(wilcox_res,c(colnames(middata1)[j],test$p.value,mean(middata1[,j][1:groupnum$Freq[1]]),median(middata1[,j][1:groupnum$Freq[1]]),mean(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),median(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),j))
}
wilcox_res <- as.data.frame(wilcox_res)
colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
wilcox_res <- na.omit(wilcox_res)
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
  mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()

wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
write.csv(wilcox_res,file = paste("All-8Study-Contained-Species-wilcoxon-test.csv",sep = ''),row.names = F)

####
data=middata1;phenodata=phe_data;k="euclidean";ShuffleTime=10000;DownPercent = 0.2;Uppercent=0.8

######## Every Study -> Pair + Wolcoxon ########
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
for (name in Study) {
  middata2 <- SpeciesData2 %>% filter(Study == name) %>% arrange(study_condition)
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  phe_data <- middata2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
    dplyr::rename(id=1,grp=2) %>%
    arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
  groupnum <- table(phe_data$grp) %>% data.frame()
  #write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste(name,"-metaphlan2.SubjectID-Newname.csv",sep = ''),row.names = F)
  #rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential species
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  write.csv(wilcox_res,file = paste(name,"-Metaphlan2-wilcoxonTest.csv",sep = ''),row.names = F)
  assign(paste(name,".wilcoxon",sep = ''),wilcox_res)
  
  # pair analysi wilcoxon pair test
  filename=paste(name,".metaphaln2.pair.csv",sep = '')
  #res<-pair_find(data=middata,RAN_num=groupnum$Freq[1],RAP_num=groupnum$Freq[2],k="euclidean")
  res<-pair_find(data=middata,phe_data,k="euclidean")
  
  #assign(paste(name,".wilcoxonsign",sep = ''),res)
  write.csv(RandData,paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''),row.names = F)
}


#################################### Sampling ###############################
######## Sampling Repeat 10 Times #########
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
dir.create("Sampling")
SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
SpeciesData2 <- SpeciesData2 %>% arrange(study_condition)
#Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200) #
for (k in c(1:10)) {
  for (i in SampleNumber) {
    print(k)
    print(i)
    library(sampling)
    SamplingIndex <- strata(SpeciesData2,stratanames = 'study_condition',size=c(i,i),method = 'srswor')
    SpeciesSelectData2 <- getdata(SpeciesData2,SamplingIndex) %>% arrange(study_condition)
    write.csv(SpeciesSelectData2,file = paste("Sampling/Sampling-",i,".",k,"-Metaphlan2.Species-Abundance.csv",sep = ""))
    ##
    #groupnum <- table(SpeciesSelectData$study_condition) %>% data.frame()
    SpeciesSelectData <- SpeciesSelectData2 %>% arrange(study_condition) %>% select(-study_condition,-Study,-ID_unit, -Prob,-Stratum)
    phe_data <- SpeciesSelectData2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
      dplyr::rename(id=1,grp=2) %>%
      arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
    groupnum <- table(phe_data$grp) %>% data.frame()
    #rownames(SpeciesSelectData)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
    
    wilcox_res<-matrix(ncol = 7,nrow = 0)
    ## wilcox test differential species
    for (j in 1:dim(SpeciesSelectData)[2]) {
      test<-wilcox.test(SpeciesSelectData[,j][1:groupnum$Freq[1]],SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]])
      wilcox_res <- rbind(wilcox_res,c(colnames(SpeciesSelectData)[j],test$p.value,mean(SpeciesSelectData[,j][1:groupnum$Freq[1]]),median(SpeciesSelectData[,j][1:groupnum$Freq[1]]),mean(SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]]),median(SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]]),j))
    }
    wilcox_res <- as.data.frame(wilcox_res)
    colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
    wilcox_res <- na.omit(wilcox_res)
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
      mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
    wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
    
    wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
    wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.05, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
    write.csv2(wilcox_res,file = paste("Sampling/Sampling-",i,".",k,"-Metaphlan2.Species-wilcoxonTest.csv",sep = ''),row.names = F)
    
    #filename=paste("Sampling-",i,".",k,".Speices.pair.csv",sep = '')
    res<-pair_find(data=SpeciesSelectData,phenodata = phe_data,k="euclidean",
                   SavePath = "Sampling",
                   BoundarySample=paste("Sampling.BoundarySample.",i,".",k,sep = ""),
                   BoundaryPair=paste("Sampling.BoundaryPair.",i,".",k,sep = ""),
                   ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
    #ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
    write.csv(res,file = paste("Sampling/Sampling-",i,".",k,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''),row.names = F)
  }
}

lapply(1:10, function(k) {
  lapply(SampleNumber, function(i){
    print(k)
    print(i)
    SamplingIndex <- strata(SpeciesData2,stratanames = 'study_condition',size=c(i,i),method = 'srswor')
    SpeciesSelectData2 <- getdata(SpeciesData2,SamplingIndex) %>% arrange(study_condition)
    write.csv(SpeciesSelectData2,file = paste("Sampling/Sampling-",i,".",k,"-Metaphlan2.Species-Abundance.csv",sep = ""))
    ##
    #groupnum <- table(SpeciesSelectData$study_condition) %>% data.frame()
    SpeciesSelectData <- SpeciesSelectData2 %>% arrange(study_condition) %>% select(-study_condition,-Study,-ID_unit, -Prob,-Stratum)
    phe_data <- SpeciesSelectData2 %>% arrange(study_condition) %>% select(study_condition) %>% rownames_to_column() %>%
      dplyr::rename(id=1,grp=2) %>%
      arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
    groupnum <- table(phe_data$grp) %>% data.frame()
    #rownames(SpeciesSelectData)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
    
    wilcox_res<-matrix(ncol = 7,nrow = 0)
    ## wilcox test differential species
    for (j in 1:dim(SpeciesSelectData)[2]) {
      test<-wilcox.test(SpeciesSelectData[,j][1:groupnum$Freq[1]],SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]])
      wilcox_res <- rbind(wilcox_res,c(colnames(SpeciesSelectData)[j],test$p.value,mean(SpeciesSelectData[,j][1:groupnum$Freq[1]]),median(SpeciesSelectData[,j][1:groupnum$Freq[1]]),mean(SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]]),median(SpeciesSelectData[,j][(1+groupnum$Freq[1]):dim(SpeciesSelectData)[1]]),j))
    }
    wilcox_res <- as.data.frame(wilcox_res)
    colnames(wilcox_res)=c("Species","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
    wilcox_res <- na.omit(wilcox_res)
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
      mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
    wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
    wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
    
    wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
    wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.05, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
    write.csv2(wilcox_res,file = paste("Sampling/Sampling-",i,".",k,"-Metaphlan2.Species-wilcoxonTest.csv",sep = ''),row.names = F)
    
    #filename=paste("Sampling-",i,".",k,".Speices.pair.csv",sep = '')
    res<-pair_find(data=SpeciesSelectData,phe_data,k="euclidean",
                   SavePath = "Sampling",
                   BoundarySample=paste("Sampling.BoundarySample.",i,".",k,sep = ""),
                   BoundaryPair=paste("Sampling.BoundaryPair.",i,".",k,sep = ""),
                   ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
    #ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
    write.csv(res,file = paste("Sampling/Sampling-",i,".",k,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''),row.names = F)
  })
})

######## Sampling RF model ##########
setwd("K:/CRC-Pair/Unique.Pair.Permutation/Sampling")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
SpeciesData2 <- read.table("../../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)

ModelImportance <- data.frame()
confusionMatrixdata <- data.frame()
modelROC <- data.frame()
model.AUC <- data.frame()

for (feature in seq(10,50,10)) {
  for (SampleN in SampleNumber) {
    for (i in 1:10) {
      middata <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      SampleData <- read.csv(paste("Sampling-",SampleN,".",i,"-Metaphlan2.Species-Abundance.csv",sep = ''),row.names = 1)
      
      middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
      mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
      mid2 <- mid %>% arrange(Pvalue) %>% top_n(-feature,Pvalue)
      
      #mid2 <- mid %>% top_n(-feature,Decre.aveRank.P) %>% arrange(Decre.aveRank.P) #%>% arrange(Enrichment)
      
      TrainData <- SpeciesData2[rownames(SpeciesData2) %in% rownames(SampleData),] %>% data.frame() %>% 
        mutate(Label = if_else(study_condition == "control",0,1)) %>% select(c(mid2$Species,Label))
      TrainData$Label <- factor(TrainData$Label,levels = c(0,1))
      
      
      TestData <- SpeciesData2[!rownames(SpeciesData2) %in% rownames(SampleData),] %>% data.frame()%>% 
        mutate(Label = if_else(study_condition == "control",0,1)) %>% select(c(mid2$Species,Label))
      TestData$Label <- factor(TestData$Label,levels = c(0,1))
      
      control <- trainControl(method="repeatedcv",repeats=5)
      
      fit.rf <- train(Label ~ .,data = TrainData, method = "rf", metric="Accuracy", trControl = control)
      
      rf.pred <- predict(fit.rf, TestData)
      cm<-confusionMatrix(rf.pred,TestData$Label)
      confusionMatrixdata<-data.frame(cbind(t(cm$overall),t(cm$byClass))) %>% 
        mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>% 
        rbind(confusionMatrixdata)
      
      predob = predict(fit.rf,TestData,type = "prob")
      pred<-prediction(predob[,2],TestData$Label)
      perf<-performance(pred,'tpr','fpr')
      
      #extrac plot ROC data
      modelROC<-data.frame(FPR = unlist(perf@x.values),TPR = unlist(perf@y.values)) %>% 
        mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>%
        rbind(modelROC)
      
      auc<-performance(pred,"auc")
      auc<-unlist(slot(auc,"y.values"))
      model.AUC <- data.frame(AUC=auc) %>% mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>% rbind(model.AUC)
      
      ModelImportance<-fit.rf$finalModel$importance %>% data.frame() %>% rownames_to_column()  %>% 
        mutate(Rank = floor(rank(-MeanDecreaseGini))) %>% mutate(FeatureCount = feature,Sampling = SampleN*2,RepeatTimes = i) %>%
        rbind(ModelImportance)
    }
  }
}

write.csv(confusionMatrixdata,file = paste("Sampling.RF.model.ConfusionMatrix.csv",sep = ''),row.names = F)
write.csv(model.AUC,file = paste("Sampling.RF.model.AUC.csv",sep = ''),row.names = F)
write.csv(modelROC,file = paste("Sampling.RF.model.ROC.csv",sep = ''),row.names = F)
write.csv(ModelImportance,file = paste("Sampling.RF.model.Importance.csv",sep = ''),row.names = F)


################################## PATHWAY #################################
dir.create("PATHWAY.PERMUTATION")
setwd("K:/CRC-Pair/Unique.Pair.Permutation/PATHWAY.PERMUTATION")

#### All Paired + Wilcoxon ####
PathwayData2 <- read.table("../../Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
#PathwayData2 <- PathwayData2 %>% arrange(study_condition)
NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2$study_condition <- factor(PathwayData2$study_condition,levels = c("control","CRC"))
phe_data <- PathwayData2 %>% select(study_condition) %>% rownames_to_column() %>%
  dplyr::rename(id=1,grp=2) %>%
  arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))

PathwayData2 <- PathwayData2[order(PathwayData2$study_condition),] %>% data.frame()

groupnum <- table(PathwayData2$study_condition) %>% data.frame()
middata1 <- PathwayData2[,-c(1:2)] %>% data.frame() #%>% select(-study_condition,-Study)
filename=paste("All-8Study-Contained.Pathway.pair.csv",sep = '')
res<-pair_find(data=middata1,phenodata=phe_data,k="euclidean",BoundarySample="All.BoundarySample",BoundaryPair = "All.BoundaryPair")
#ExcludeData<-res2 %>% filter(adj.fdr.p <= 0.01)  %>% mutate(ExcludeStudy = name) %>% rbind(ExcludeData)
res<- merge(res,NameData,by.x = "Species",by.y = "Rename",all.x =T)
write.csv(res,file = paste("All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv",sep = ''),row.names = F)


rownames(middata1)<-c(paste("Ctl",1:groupnum$Freq[1],sep = ""),paste("CRC",1:groupnum$Freq[2],sep = ""))
wilcox_res<-matrix(ncol = 7,nrow = 0)
## wilcox test differential Pathway
for (j in 1:dim(middata1)[2]) {
  test<-wilcox.test(middata1[,j][1:groupnum$Freq[1]],middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]])
  wilcox_res <- rbind(wilcox_res,c(colnames(middata1)[j],test$p.value,mean(middata1[,j][1:groupnum$Freq[1]]),median(middata1[,j][1:groupnum$Freq[1]]),mean(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),median(middata1[,j][(1+groupnum$Freq[1]):dim(middata1)[1]]),j))
}
wilcox_res <- as.data.frame(wilcox_res)
colnames(wilcox_res)=c("Pathway","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","SamplingCountForEachGroup")
wilcox_res <- na.omit(wilcox_res)
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
  mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()

wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
wilcox_res<- merge(wilcox_res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
write.csv(wilcox_res,file = paste("All-8Study-Contained-Pathway-wilcoxon-test.csv",sep = ''),row.names = F)

#### Every Study Pair + Wilcoxon ####
setwd("K:/CRC-Pair/Unique.Pair.Permutation/PATHWAY.PERMUTATION")
PathwayData2 <- read.table("../../Pathway-8Study-20201010/EightStudies-PathwayAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)
NameData <- data.frame(Origin = colnames(PathwayData2)[3:dim(PathwayData2)[2]],Rename=paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = ''))
write.csv(NameData,file = "Pathway.Rename.csv",row.names = F)

colnames(PathwayData2)[3:dim(PathwayData2)[2]] = paste("Pathway",1:(dim(PathwayData2)[2]-2),sep = '')

PathwayData2$study_condition <- factor(PathwayData2$study_condition,levels = c("control","CRC"))
PathwayData2 <- PathwayData2[order(PathwayData2$study_condition),] %>% data.frame()

Study <- c("FengQ_2015","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016") #
for (name in Study) {
  middata2 <- PathwayData2 %>% filter(Study == name) %>% arrange(study_condition)
  groupnum <- table(middata2$study_condition) %>% data.frame()
  
  phe_data <- middata2 %>% select(study_condition) %>% rownames_to_column() %>%
    dplyr::rename(id=1,grp=2) %>%
    arrange(grp) %>% mutate(grp = dplyr::case_when(grp == "control" ~ "grp1",grp == "CRC" ~ "grp2"))
  
  middata <- middata2 %>% arrange(study_condition) %>% select(-study_condition,-Study)
  
  # pair analysi wilcoxon pair test
  #filename=paste(name,".metaphaln2.pair.csv",sep = '')
  res<-pair_find(data=middata,phenodata = phe_data,k="euclidean",BoundarySample=paste("Pathway.",name,".BoundarySample.",name,sep = ""),BoundaryPair=paste("Pathway.",name,".BoundaryPair.",name,sep = ""),ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8)
  #assign(paste(name,".wilcoxonsign",sep = ''),res)
  #data=middata;phenodata=phe_data;k="euclidean";BoundarySample=paste("Pathway.",name,".BoundarySample.",name,sep = "");BoundaryPair=paste("Pathway.",name,".BoundaryPair.",name,sep = "");ShuffleTime=10000;DownPercent = 0.2;Uppercent=0.8
  # merge final Pathway Data
  res<- merge(res,NameData,by.x = "Species",by.y = "Rename",all.x =T)
  write.csv(res,file = paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),row.names = F)
  
  
  rownames(middata) <- c(paste(groupnum$Var1[1]%>% as.vector(),1:groupnum$Freq[1]),paste(groupnum$Var1[2]%>% as.vector(),1:groupnum$Freq[2])) %>% str_remove_all(" ")
  #write.csv(data.frame(NewName=rownames(middata),subjectid=rownames(middata2)),file = paste(name,"-Humann2.SubjectID-Newname.csv",sep = ''),row.names = F)
  wilcox_res<-matrix(ncol = 7,nrow = 0)
  ## wilcox test differential Pathway
  for (i in 1:dim(middata)[2]) {
    test<-wilcox.test(middata[,i][1:groupnum$Freq[1]],middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]])
    wilcox_res <- rbind(wilcox_res,c(colnames(middata)[i],test$p.value,mean(middata[,i][1:groupnum$Freq[1]]),median(middata[,i][1:groupnum$Freq[1]]),mean(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),median(middata[,i][(1+groupnum$Freq[1]):dim(middata)[1]]),name))
  }
  wilcox_res <- as.data.frame(wilcox_res)
  colnames(wilcox_res)=c("Pathway","Pvalue","Ctlmean","Ctlmedian","CRCmean","CRCmedian","Study")
  wilcox_res <- na.omit(wilcox_res)
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit() %>% mutate(pvalue = as.numeric(as.character(Pvalue)))  %>%
    mutate(Ctlmean = as.numeric(as.character(Ctlmean)),CRCmean=as.numeric(as.character(CRCmean)))
  wilcox_res$Pvalue<-as.numeric(as.character(wilcox_res$Pvalue))
  wilcox_res <- wilcox_res %>% arrange(Pvalue) %>% na.omit()
  
  wilcox_res$adj.fdr.p<-p.adjust(wilcox_res$Pvalue,method="fdr",n=length(wilcox_res$Pvalue))
  wilcox_res$"Enrichment" = if_else(wilcox_res$adj.fdr.p <= 0.01, ifelse(wilcox_res$Ctlmean > wilcox_res$CRCmean,"Ctl","CRC"),"Nodiff")
  # merge final Pathway Data
  wilcox_res<- merge(wilcox_res,NameData,by.x = "Pathway",by.y = "Rename",all.x =T)
  
  write.csv(wilcox_res,file = paste(name,"-Humann2-wilcoxonTest.csv",sep = ''),row.names = F)
  #assign(paste(name,".wilcoxon",sep = ''),wilcox_res)
  
}

###########################################################################################################################
########### Fig 1 - MIN #############
#### dataset ######
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
Fig1Data <- data.frame()
IntersectSpecies <- read_excel("MetaIntersectSpecies.xlsx")
IntersectSpecies$Species <- IntersectSpecies$Species %>% str_replace_all(" ","_")

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (name in Study) {
  index <- which(Study == name)
  middata <- read.csv(paste(name,"-Metaphlan2-wilcoxonTest.csv",sep = ''))
  #Wilcon Test
  Fig1Data  <- na.omit(middata) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0)) %>% 
    filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".1",sep = ''),Test="Wilcoxon") %>% rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(Fig1Data)
  #Pair
  middata <- read.csv(paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''))
  Fig1Data.mid  <- na.omit(middata) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min == 1)))
  Fig1Data.mid$Incre.minRank.P <- Fig1Data.mid$Increasing.Rank.Min/10001
  Fig1Data.mid2 <- data.frame(Species = c(Fig1Data.mid$Species,Fig1Data.mid$Species),Pvalue = c(Fig1Data.mid$Incre.minRank.P,Fig1Data.mid$Decre.minRank.P))
  Fig1Data.mid3 <- Fig1Data.mid2 %>% data.frame() %>% dplyr::arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()
  
  
  Fig1Data <- Fig1Data.mid3 %>% data.frame() %>% dplyr::arrange(Pvalue) %>% data.frame() %>% mutate(Rank = rank(Pvalue,ties.method = "min")) %>% 
    filter(Species %in% IntersectSpecies$Species) %>% mutate(Study = name) %>%
    select(Species,Study,Pvalue,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".2",sep = ''),Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(Fig1Data)
}

#Meta Wilcoxon
MetaData <- read.csv("../Species-8Study-20201010/All-8Study-Contained-Species-wilcoxon-test.csv")
Fig1Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0),Study="Meta") %>% 
  filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.1",Test="Wilcoxon") %>% rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind(Fig1Data )

#Meta Pair
MetaData <- read.csv("Unique.Pair.SampleAll.0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv")
Fig1Data.mid  <- na.omit(MetaData) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min== 1)))
Fig1Data.mid$Incre.minRank.P <- Fig1Data.mid$Increasing.Rank.Min/10001
Fig1Data.mid2 <- data.frame(Species = c(Fig1Data.mid$Species,Fig1Data.mid$Species),Pvalue = c(Fig1Data.mid$Incre.minRank.P,Fig1Data.mid$Decre.minRank.P))
Fig1Data.mid3 <- Fig1Data.mid2 %>% data.frame() %>% arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()

Fig1Data <- Fig1Data.mid3 %>% data.frame() %>% arrange(Pvalue) %>% mutate(Rank = round(rank(Pvalue,ties.method = "min"),0)) %>% 
  filter(Species %in% IntersectSpecies$Species) %>% mutate(Study="Meta") %>%
  dplyr::select(Species,Study,Pvalue,Rank) %>% 
  mutate(Country = "Meta.2",Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind(Fig1Data)

write.csv(Fig1Data,"Figure/Figure1.Data.csv")

Pdata <- Fig1Data %>% select(Species,Country,FDR,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata$Species <- factor(Pdata$Species,levels = IntersectSpecies$Species)
Pdata <- Pdata[order(Pdata$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,FDR) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata,"Figure/Fig1a.csv")

#### Fig 1 ####
Pdata.1 <- read.csv("Figure/Fig1a.csv",row.names = 1,stringsAsFactors = F)
Pdata.1 <- Pdata.1 %>% arrange(Meta.2,Meta.1)
#Pdata <- t(Pdata)
#rownames(Pdata.1) <- rownames(Pdata.1) %>% str_replace_all("_"," ")
Pdata <- read.csv("Figure/Fig1a.csv",row.names = 1,stringsAsFactors = F)
Pdata <- Pdata %>% rownames_to_column() %>%
  gather(key = "Condition",value = "FDR",-rowname)

labels = c("<0.001","0.001-0.01","0.01-0.05",">0.05")
breaks <- c(-1,0.001,0.01,0.05,1)
breaks2 <- c(0,001,0.01,0.05)
#mid <- data.frame(labels=labels,breaks=breaks2)
Pdata$labels = cut(Pdata$FDR,breaks,labels,ordered_result = T)
#Pdata<-merge(Pdata,mid,by= "labels")
#Pdata$alpha <- Pdata$FDR/Pdata$breaks
Pdata$rowname <-factor(Pdata$rowname,levels = rownames(Pdata.1))
interval.cols <- c("#8B1A1A","#FF6A6A","#FFC1C1","#EEE9E9")#brewer.pal(6,"Set2")
names(interval.cols) <- levels(Pdata$labels)

Pdata.2 <- Pdata %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  mutate(rowname = factor(rowname,levels = rownames(Pdata.1))) %>% arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))

## text note => False
TextFunc <- function(dat, col = "black", fontsize = 8, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}
###
col_cat <- c("<0.001"="#8B1A1A","0.001-0.01"="#FF6A6A","0.01-0.05"="#FFC1C1",">0.05"="#EEE9E9")

### Sperate the data for twins cohort and original cohort
### Fig1 Top
Country.Heatmap <- c("AUS","CHI","FRA","GER","ITA1","ITA2","JAP","USA","Meta")
Pdata.2.mid1 <- Pdata.2 %>% select(ends_with(".1")) %>% select(paste(Country.Heatmap,".1",sep = ''))
rownames(Pdata.2.mid1) <- rownames(Pdata.2.mid1) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid1) = Country.Heatmap
Pdata.2.mid1 <- t(Pdata.2.mid1)

Pdata.1.mid1 <- Pdata.1 %>% select(ends_with(".1")) %>% select(paste(Country.Heatmap,".1",sep = ''))
rownames(Pdata.1.mid1) <- rownames(Pdata.1.mid1) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid1) = Country.Heatmap
Pdata.1.mid1 <- t(Pdata.1.mid1)


pdf("Figure/Fig1a.Top.pdf",width = 7,height = 5)
Heatmap(Pdata.2.mid1, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid1)
)
dev.off()

### Fig1 Middle
Pdata.2.mid2 <- Pdata.2 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.2.mid2) <- rownames(Pdata.2.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid2) = Country.Heatmap
Pdata.2.mid2 <- t(Pdata.2.mid2)

Pdata.1.mid2 <- Pdata.1 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.1.mid2) <- rownames(Pdata.1.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid2) = Country.Heatmap
Pdata.1.mid2 <- t(Pdata.1.mid2)


pdf("Figure/Fig1a.Middle.pdf",width = 7,height = 5)
Heatmap(Pdata.2.mid2, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid2)
)
dev.off()

########### Fig 1 - AVE #############
#### dataset ######
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
Fig1Data <- data.frame()
IntersectSpecies <- read_excel("MetaIntersectSpecies.xlsx")
IntersectSpecies$Species <- IntersectSpecies$Species %>% str_replace_all(" ","_")

Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (name in Study) {
  index <- which(Study == name)
  middata <- read.csv(paste(name,"-Metaphlan2-wilcoxonTest.csv",sep = ''))
  #Wilcon Test
  Fig1Data  <- na.omit(middata) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0)) %>% 
    filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".1",sep = ''),Test="Wilcoxon") %>% rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(Fig1Data)
  #Pair
  middata <- read.csv(paste("Unique.Pair.Sample",name,".0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv",sep = ''))
  Fig1Data.mid  <- na.omit(middata) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min == 1)))
  Fig1Data.mid$Incre.aveRank.P <- Fig1Data.mid$Increasing.Rank.Average/10001
  Fig1Data.mid2 <- data.frame(Species = c(Fig1Data.mid$Species,Fig1Data.mid$Species),Pvalue = c(Fig1Data.mid$Incre.aveRank.P,Fig1Data.mid$Decre.aveRank.P))
  Fig1Data.mid3 <- Fig1Data.mid2 %>% data.frame() %>% dplyr::arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()
  
  
  Fig1Data <- Fig1Data.mid3 %>% data.frame() %>% dplyr::arrange(Pvalue) %>% data.frame() %>% mutate(Rank = rank(Pvalue,ties.method = "min")) %>% 
    filter(Species %in% IntersectSpecies$Species) %>% mutate(Study = name) %>%
    select(Species,Study,Pvalue,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".2",sep = ''),Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(Fig1Data)
}

#Meta Wilcoxon
MetaData <- read.csv("../Species-8Study-20201010/All-8Study-Contained-Species-wilcoxon-test.csv")
Fig1Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0),Study="Meta") %>% 
  filter(Species %in% IntersectSpecies$Species) %>% select(Species,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.1",Test="Wilcoxon") %>% rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind(Fig1Data )

#Meta Pair
MetaData <- read.csv("Unique.Pair.SampleAll.0.2.0.8.10000.Feature.Shuffle.Rank.Pvalue.csv")
Fig1Data.mid  <- na.omit(MetaData) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min== 1)))
Fig1Data.mid$Incre.aveRank.P <- Fig1Data.mid$Increasing.Rank.Average/10001
Fig1Data.mid2 <- data.frame(Species = c(Fig1Data.mid$Species,Fig1Data.mid$Species),Pvalue = c(Fig1Data.mid$Incre.aveRank.P,Fig1Data.mid$Decre.aveRank.P))
Fig1Data.mid3 <- Fig1Data.mid2 %>% data.frame() %>% arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()

Fig1Data <- Fig1Data.mid3 %>% data.frame() %>% arrange(Pvalue) %>% mutate(Rank = round(rank(Pvalue,ties.method = "min"),0)) %>% 
  filter(Species %in% IntersectSpecies$Species) %>% mutate(Study="Meta") %>%
  dplyr::select(Species,Study,Pvalue,Rank) %>% 
  mutate(Country = "Meta.2",Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind(Fig1Data)

write.csv(Fig1Data,"Figure/Figure1.AVE.Data.csv")

Pdata <- Fig1Data %>% select(Species,Country,FDR,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata$Species <- factor(Pdata$Species,levels = IntersectSpecies$Species)
Pdata <- Pdata[order(Pdata$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,FDR) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata,"Figure/Fig1a.AVE.csv")

#### Fig 1 ####
Pdata.1 <- read.csv("Figure/Fig1a.AVE.csv",row.names = 1,stringsAsFactors = F)
Pdata.1 <- Pdata.1 %>% arrange(Meta.2,Meta.1)
#Pdata <- t(Pdata)
#rownames(Pdata.1) <- rownames(Pdata.1) %>% str_replace_all("_"," ")
Pdata <- read.csv("Figure/Fig1a.AVE.csv",row.names = 1,stringsAsFactors = F)
Pdata <- Pdata %>% rownames_to_column() %>%
  gather(key = "Condition",value = "FDR",-rowname)

labels = c("<0.001","0.001-0.01","0.01-0.05",">0.05")
breaks <- c(-1,0.001,0.01,0.05,1)
breaks2 <- c(0,001,0.01,0.05)
#mid <- data.frame(labels=labels,breaks=breaks2)
Pdata$labels = cut(Pdata$FDR,breaks,labels,ordered_result = T)
#Pdata<-merge(Pdata,mid,by= "labels")
#Pdata$alpha <- Pdata$FDR/Pdata$breaks
Pdata$rowname <-factor(Pdata$rowname,levels = rownames(Pdata.1))
interval.cols <- c("#8B1A1A","#FF6A6A","#FFC1C1","#EEE9E9")#brewer.pal(6,"Set2")
names(interval.cols) <- levels(Pdata$labels)

Pdata.2 <- Pdata %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  mutate(rowname = factor(rowname,levels = rownames(Pdata.1))) %>% arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))

## text note => False
TextFunc <- function(dat, col = "black", fontsize = 8, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}
###
col_cat <- c("<0.001"="#8B1A1A","0.001-0.01"="#FF6A6A","0.01-0.05"="#FFC1C1",">0.05"="#EEE9E9")

### Sperate the data for twins cohort and original cohort

### Fig1 Middle
Pdata.2.mid2 <- Pdata.2 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.2.mid2) <- rownames(Pdata.2.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid2) = Country.Heatmap
Pdata.2.mid2 <- t(Pdata.2.mid2)

Pdata.1.mid2 <- Pdata.1 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.1.mid2) <- rownames(Pdata.1.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid2) = Country.Heatmap
Pdata.1.mid2 <- t(Pdata.1.mid2)


pdf("Figure/Fig1a.Middle.AVE.pdf",width = 7,height = 5)
Heatmap(Pdata.2.mid2, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid2)
)
dev.off()

########### Fig S1 - AVE #########
Fig1Data <- read.csv("Figure/Figure1.AVE.Data.csv",row.names = 1)
Pdata2 <- Fig1Data %>% select(Species,Country,Rank,Test) %>% arrange(Country) %>% arrange(desc(Test))
#Pdata2$Species <- factor(Pdata2$Species,levels = IntersectSpecies$Species)
Pdata2 <- Pdata2[order(Pdata2$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,Rank) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata2,"Figure/Fig1b.csv")

dPdata2 <- read.csv("Figure/Fig1b.csv",row.names = 1)
dPdata3 <- dPdata2 %>% rownames_to_column() %>%
  gather(key = "Condition",value = "Rank",-rowname)

labels = c("1-10","10-20","20-30","30-40","40-50",">50")
breaks <- c(0,10,20,30,40,50,500)
dPdata3$labels = cut(dPdata3$Rank,breaks,labels,ordered_result = T)

TextFunc <- function(dat, col = "black", fontsize = 8, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(round(dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(round(dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}

Country.Heatmap <- c("AUS","CHI","FRA","GER","ITA1","ITA2","JAP","USA","Meta")

dPdata4 <- dPdata3 %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  select(rowname,ends_with(".2")) %>% 
  mutate(rowname = factor(.$rowname,levels = rownames(dPdata2))) %>%
  arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>% 
  select(paste(Country.Heatmap,".2",sep = ''))
  
colnames(dPdata4) = Country.Heatmap

dPdata5 <- dPdata3 %>% select(Rank,rowname,Condition) %>% spread(Condition,Rank) %>% 
  select(rowname,ends_with(".2")) %>% 
  mutate(rowname = factor(.$rowname,levels = rownames(dPdata2))) %>%
  arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>% 
  select(paste(Country.Heatmap,".2",sep = ''))
colnames(dPdata5) = Country.Heatmap

col_cat <- c("1-10"="#8B1A1A","10-20"="#EE0000","20-30"="#FF6A6A","30-40"="#FFC1C1","40-50"="#CDC9C9",">50"="#EEE9E9")
pdf("Figure/FigS1.Right.AVE.pdf",width = 6,height = 3)
Heatmap(dPdata4, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(dPdata5)
)
dev.off()


########### Fig 2 #############
#### Figure 2a => Intra Study CV 5 Repeat 10 Times ####
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")
NewStudy <- data.frame(Study,StudyPos)

IntraData <- read.csv("IntraStudy/IntraStudy.CV5Top10-50.RF.model.AUC.csv")
IntraData$FeatureCount <- factor(IntraData$FeatureCount,levels = seq(10,50,10))
IntraData2 <- IntraData %>% dplyr::group_by(Study,FeatureCount) %>% dplyr::summarise(MeanAUC = mean(AUC))
IntraData2 <- merge(IntraData2,NewStudy,by = "Study",all.x = T)

p1 <- ggplot(IntraData2,aes(FeatureCount,MeanAUC))+
  geom_point(aes(group=StudyPos,color=StudyPos),size=3)+
  geom_line(aes(group=StudyPos,color=StudyPos),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  labs(x="No. of features used",y="AUC")+ylim(0.6,1)+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))
ggsave(p1,filename = "Figure/Fig2a.pdf",width = 5,height = 4)

#### Figure 2b => Study to Study Transfer + LODO ####
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
###
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # 计算长度
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # 重命名  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  # 计算标准偏差
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
###
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")
NewStudy <- data.frame(Study,StudyPos)
Fig2bData <- data.frame()
for (feature in seq(10,50,10)) {
  data <- read.csv(paste("StudyToStudy/Study2Study.Top",feature,".RF.model.AUC.csv",sep = ""))
  Fig2bData <- data %>% filter(Predict != BaseModel) %>% dplyr::select(-Predict) %>% dplyr::rename(AUC=1,Study=2) %>%
    dplyr::mutate(Condition="Study2study",FeatureCount = feature) %>% rbind(Fig2bData)
  
  data <- read.csv(paste("LODO.Test/LODO.Top",feature,".RF.model.AUC.csv",sep = ''))
  Fig2bData <- data %>% filter(Predict != "Self") %>% dplyr::select(-Predict) %>% dplyr::rename(AUC=1,Study=2) %>%
    dplyr::mutate(Condition="LODOValidation", FeatureCount = feature) %>% rbind(Fig2bData)
}

Fig2bData <- merge(Fig2bData,NewStudy,by="Study",all.x = T)
write.csv(Fig2bData,"Figure/Fig2b.data.csv",row.names = F)

Fig2bData <- read.csv("Figure/Fig2b.data.csv")
for (feature in seq(10,50,10)) {
  mid <- Fig2bData %>% filter(FeatureCount == feature)
  mid$Condition <- factor(mid$Condition,levels = c("Study2study","LODOValidation"))
  mid2 <- summarySE(mid, measurevar="AUC", groupvars=c("StudyPos","Condition"))
  mid2$Condition <- factor(mid2$Condition,levels = c("Study2study","LODOValidation"))
  p3<-ggplot(mid2, aes(x=StudyPos, y=AUC, fill=Condition)) + 
    geom_bar(position=position_dodge(),stat='identity',colour='black') +
    geom_jitter(data=mid,aes(StudyPos,AUC),color="grey60",position=position_dodge(width=0.85))+
    geom_errorbar(aes(x=StudyPos,ymin=AUC-se, ymax=AUC+se),
                  width=.3,size=1, # 设置误差线的宽度 
                  position=position_dodge(.8))+
    theme_few()+
    labs(y="AUC",x="")+scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))+
    theme(legend.position = "top",legend.title = element_blank())+
    scale_fill_manual(values=c( 'white','grey'))
  
  ggsave(p3,filename = paste("Figure/Fig2b.Features.",feature,".2.pdf",sep = ''),width = 6,height = 5)
}

#### Figure 2c => LODO Importance #####
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

Studyselect <- c("AUS","CHI","FRA","GER","ITA1","ITA2","JAP","USA")
  
for (feature in c(10,20,30,40,50)) {
  mid <- read.csv(paste("LODO.Test/LODO.Top",feature,".RF.model.Importance.csv",sep = ''))
  mid2 <- mid %>% filter(Predict != "Self") %>% dplyr::select(-MeanDecreaseGini,-Predict) %>% spread(ModelExcludeStudy,Rank) %>%
    remove_rownames() %>% column_to_rownames("rowname")
  
  mid2 <- mid2[apply(mid2, 1, function(x){sum(is.na(x))}) <= 1,] %>% data.frame()
  
  DiffExplore <- data.frame()
  #confirm the species difference
  for (name in Study) {
    data<-read.csv(paste("LODO.Test/LODO-Exclude.",name,"-Metaphlan2-wilcoxonTest.csv",sep = ''))
    DiffExplore <- data[data$Species %in% rownames(mid2),] %>% data.frame() %>% dplyr::select(Species,Enrichment) %>% 
      #mutate(Study=name) %>% 
      rbind(DiffExplore)
  }
  
  mid2 <- mid2 %>% select(Study)
  colnames(mid2) = StudyPos
  mid2 <- mid2 %>% dplyr::select(Studyselect)
  #sort factor
  mid3 <- apply(mid2,1,function(x) sum(na.omit(x))) %>% sort()
  mid4 <- mid2 %>% data.frame() %>% rownames_to_column() %>% mutate(rowname = factor(.$rowname,levels = names(mid3))) %>% 
    arrange(rowname) %>% mutate(rowname = str_replace_all(rowname,"_"," ")) %>% remove_rownames() %>% column_to_rownames("rowname")
  
  pdf(paste("Figure/Fig2C.LODO-ExcludeStudy.Importance.",feature,".Re.pdf",sep = ''),height = 10,width = 8)
  pheatmap::pheatmap(mid4, 
                     display_numbers = round(mid4,0), #matrix(ifelse(abs(tdata) > 0.4, "*", ""), nrow(tdata)), 
                     #annotation_row = annotation_row, 
                     #annotation_colors = anno_color,
                     color = colorRampPalette(c("#8B1A1A","Red","RosyBrown", "white"))(100),cluster_rows = F,cluster_cols = F,
                     border_color = "black",
                     number_color = "#DCDCDC",
                     cellwidth = 24,
                     cellheight = 12,
                     fontsize_row = 12,
                     fontsize_col = 12,
                     na_col = "white",
                     legend = T,angle_col = 45)
  dev.off()
}

########### Fig 3  ################
#### Fig 3a  => Sampling Intersect Meta  core 7 Species ####
Species7 <- c("Parvimonas_unclassified","Gemella_morbillorum","Peptostreptococcus_stomatis",
              "Fusobacterium_nucleatum","Parvimonas_micra","Porphyromonas_asaccharolytica",
              "Clostridium_symbiosum")
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)

for (feature in seq(10,50,10)) {
  DataList=list()
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      middata <- read.csv(paste("Sampling/Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
      mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
      mid2 <- mid %>% arrange(Pvalue) %>% top_n(-feature,Pvalue)
      #mid2 <- mid %>% top_n(-feature,Decre.aveRank.P) %>% arrange(Decre.aveRank.P) #%>% arrange(Enrichment)
      Mid[i] <- length(intersect(mid2$Species,Species7))
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  Data <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="IntersectSpecies",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X"))
  
  p<-ggboxplot(Data,x="Sampling",y="IntersectSpecies",color = "Sampling",add = "jitter")+
    labs(x="Sampling",y="Intersect Species Count")+
    theme_few() + theme(legend.position = "none")
  
  ggsave(p,filename = paste("Figure/Sampling.Top.",feature,".IntersectMeta7.pdf",sep = ''),width = 3,height = 3)
  
}
#### Fig S4  => Sampling Intersect Meta  core 13 Species ####
setwd("K:/CRC-Pair/Unique.Pair.Permutation")

MetaSpecies <- read_excel("MetaIntersectSpecies.xlsx")
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)

for (feature in seq(10,50,10)) {
  DataList=list()
  for (SampleN in SampleNumber) {
    Mid <- rep(NA,10)
    for (i in 1:10) {
      middata <- read.csv(paste("Sampling/Sampling-",SampleN,".",i,"-Metaphlan2.Species-PairwilcoxonSign-res.csv",sep = ''))
      #mid2 <- mid %>% top_n(-feature,Decre.aveRank.P) %>% arrange(Decre.aveRank.P) #%>% arrange(Enrichment)
      middata$Incre.aveRank.P <- middata$Increasing.Rank.Average/10001
      mid <- data.frame(Pvalue = c(middata$Incre.aveRank.P,middata$Decre.aveRank.P),Species = c(middata$Species,middata$Species))
      mid2 <- mid %>% arrange(Pvalue) %>% top_n(-feature,Pvalue)
      Mid[i] <- length(intersect(mid2$Species,MetaSpecies$Species))
    }
    DataList[[as.character(SampleN*2)]] = Mid
  }
  Data <- data.frame(DataList) %>% mutate(RepeatTimes = 1:10)%>% gather(key="Sampling",value="IntersectSpecies",-RepeatTimes) %>%
    mutate(Sampling = str_remove_all(Sampling,"X"))
  
  p<-ggboxplot(Data,x="Sampling",y="IntersectSpecies",color = "Sampling",add = "jitter")+
    labs(x="Sampling",y="Intersect Species Count")+
    theme_few() + theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45))
  
  ggsave(p,filename = paste("Figure/Figure3a.Sampling.Top.",feature,".Species.IntersectMeta.pdf",sep = ''),width = 3,height = 3)
  
}

#### Fig 3b  => Sampling AUC ####
SampleNumber <- c(15,20,30,40,50,60,70,80,90,100,200)
AUCdata <- read.csv("Sampling/Sampling.RF.model.AUC.csv")

AUCdata$Sampling <- factor(AUCdata$Sampling,levels = SampleNumber*2)
AUCdata$FeatureCount <- factor(AUCdata$FeatureCount,levels = seq(10,50,10))

for (feature in seq(10,50,10)) {
  AUCdata2 <- AUCdata %>% filter(FeatureCount == feature)
  p<-ggplot(data=AUCdata2,aes(Sampling,AUC,color=Sampling)) + geom_boxplot()+geom_jitter(width = 0.1,height = 0.05,size=1)+ 
    theme_few() + 
    theme(legend.position = "none")+ylim(0,1)+
    theme(axis.text = element_text(size=13),axis.title = element_text(size=15),axis.text.x = element_text(angle = 45))
  
  ggsave(p,filename = paste("Figure/Figure3b.Sampling.RF.Top.",feature,".model0-1.AUC.pdf",sep = ""),height = 3,width = 3)
}

########### Fig S2 ##################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
AllSpecies <- read.csv("All-8Study-Contained-Species-pair-wilcoxonsign-res.csv")
AllSpecies$Incre.aveRank.P <- AllSpecies$Increasing.Rank.Average/10001
AllSpecies2 <- AllSpecies %>% mutate(Enriched = if_else(Incre.aveRank.P <= 0.05,"Ctrl",if_else(Decre.aveRank.P <= 0.05,"Disease","N.S.")))

## calculate log2FC
PairSample <- read.csv("All.BoundaryPair.Species.csv")

SpeciesData2 <- read.table("../Species-8Study-20201010/EightStudies-SpeciesAbundance-Group.txt",sep = '\t',row.names = 1,stringsAsFactors = F)

Res.Data <- data.frame()
for (name in AllSpecies2$Species) {
  Ctrlmean <- rep(NA,dim(PairSample)[1])
  CRCmean <- rep(NA,dim(PairSample)[1])
  for (i in 1:dim(PairSample)[1]) {
    Ctrlmean[i] <- SpeciesData2[PairSample$Ctl[i],name]
    CRCmean[i] <- SpeciesData2[PairSample$Disease[i],name]
  }
  Res.Data <- rbind(Res.Data,c(name,mean(Ctrlmean),mean(CRCmean)))
}
Res.Data2 <- Res.Data %>% dplyr::rename(Species=1,Ctrlmean=2,CRCmean=3) %>% mutate(log2FC = log2(as.numeric(CRCmean)/as.numeric(Ctrlmean)))

AllSpecies3 <- merge(AllSpecies2,Res.Data2,by="Species")

library(readxl)
Meta13 <- read_xlsx("MetaIntersectSpecies.xlsx") %>% mutate(id = Species) 
AllSpecies4 <- AllSpecies3 %>% mutate(ID = if_else(Species %in% Meta13$Species,Species,""))

mid <- data.frame(Pvalue = c(AllSpecies4$Incre.aveRank.P,AllSpecies4$Decre.aveRank.P),Species = c(AllSpecies4$Species,AllSpecies4$Species))
middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-30,Pvalue)
AllSpecies5 <- AllSpecies4 %>% mutate(Label = if_else((Species %in% middataTop30$Species) & (ID == ""),Species,"")) %>% 
  mutate(Shape = if_else(ID == "","B","A"))

write.csv(AllSpecies5,"Figure/FigS2.data.csv")
library(ggthemes)
library(ggrepel)
p1<-ggplot(AllSpecies5,aes(x=log2FC,y=-log10(Decre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'CRC-enriched P'))
ggsave(p1,filename = "Figure/FigS2.Decre.pdf",height = 10,width = 10)


p2<-ggplot(AllSpecies5,aes(x=log2FC,y=-log10(Incre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'Ctrl-enriched P'))
ggsave(p2,filename = "Figure/FigS2.Incre.pdf",height = 10,width = 10)


p1<-ggplot(AllSpecies5,aes(x=log2FC,y=-log10(Decre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'CRC-enriched P'))
ggsave(p1,filename = "Figure/FigS2.Decre.Rev.pdf",height = 5,width = 5)


p2<-ggplot(AllSpecies5,aes(x=log2FC,y=-log10(Incre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=2,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=3,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'Ctrl-enriched P'))
ggsave(p2,filename = "Figure/FigS2.Incre.Rev.pdf",height = 5,width = 5)

p1<-ggplot(AllSpecies5,aes(x=log2FC,y=-log10(Decre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=4.5,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=7,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black",max.overlaps = 20)+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'CRC-enriched P'))
ggsave(p1,filename = "Figure/FigS2.Decre.Test.pdf",height = 10,width = 10)


p2<-ggplot(AllSpecies5,aes(x=log2FC,y=-log10(Incre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=4.5,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label),size=7,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black",max.overlaps = 20)+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'Ctrl-enriched P'))
ggsave(p2,filename = "Figure/FigS2.Incre.Test.pdf",height = 10,width = 10)


########### Fig 4 ###########################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
AllPathway <- read.csv("PATHWAY.PERMUTATION/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")

AllPathway$Incre.aveRank.P <- AllPathway$Increasing.Rank.Average/10001
AllPathway2 <- AllPathway %>% mutate(Enriched = if_else(Incre.aveRank.P <= 0.01,"Ctrl",if_else(Decre.aveRank.P <= 0.01,"Disease","N.S.")))
AllPathway3 <- AllPathway2 %>% mutate(log2FC = log2(as.numeric(Dismean)/as.numeric(Ctlmean)))

#mid <- data.frame(Pvalue = c(AllPathway3$Incre.aveRank.P,AllPathway3$Decre.aveRank.P),Species = c(AllPathway3$Species,AllPathway3$Species))
#middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-30,Pvalue)

p1<-ggplot(AllPathway3,aes(x=log2FC,y=-log10(Decre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=4.5,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label2),size=7,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.01)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'CRC-enriched P'))
ggsave(p1,filename = "Figure/Fig4a.Decre.pdf",height = 10,width = 10)

p2<-ggplot(AllPathway3,aes(x=log2FC,y=-log10(Incre.aveRank.P),shape=Shape,color=Enriched))+geom_point(size=4.5,alpha=0.6)+theme_few()+
  scale_color_manual(values = c("#0072B5","#BC3C28","grey60"))+
  #scale_size_continuous(range = c(1,10))+
  geom_text_repel(aes(label = Label2),size=7,box.padding=unit(0.5, "lines"),arrow = arrow(length=unit(0.008, "npc")),
                  alpha=1,force = 3,direction = "both",#angle=ifelse(data$Color=="Disease",5,if_else(data$Color=="Ctl",-5,0)),
                  #nudge_x = ifelse(data$Color == "Ctl", -1, 0), nudge_y = ifelse(data$Color=="Disease", 1, 0)
                  max.iter = 3e3,fontface="bold",point.padding=unit(0.01, "lines"),
                  segment.color="black")+
  #scale_x_continuous(limits = c(-2,4))+
  geom_hline(aes(yintercept=-log10(0.01)),linetype=2,color="black")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=24),
        legend.position = "none")+
  labs(x=expression(Log[2]*'Fold Change'),y=expression(-Log[10]*'Ctrl-enriched P'))
ggsave(p2,filename = "Figure/Fig4a.Incre.pdf",height = 10,width = 10)



########### Fig S3 ############
setwd("K:/CRC-Pair/Unique.Pair.Permutation/LODO.Test")
AllAUC <- data.frame()
for (feature in seq(10,50,10)) {
  mid <- read.csv(paste("LODO.Top",feature,".RF.model.AUC.csv",sep = ''))
  AllAUC <- mid %>% mutate(FeatureCount = feature) %>% rbind(AllAUC)
}
AllAUC$FeatureCount <- factor(AllAUC$FeatureCount,levels = seq(10,50,10))
AllAUC$"Condition" = if_else(AllAUC$Predict == "Self","Self","ExcludeStudy")
p<-ggplot(AllAUC%>%filter(Condition == "Self"),aes(FeatureCount,AUC))+
  geom_point(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),size=3)+
  geom_line(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  xlab("No. of features used")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))+ylim(0.5,1)

ggsave(p,filename = "LODO-ExcludeStudy.Self.AUC.pdf",height = 4,width = 6)

p<-ggplot(AllAUC%>%filter(Condition == "ExcludeStudy"),aes(FeatureCount,AUC))+
  geom_point(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),size=3)+
  geom_line(aes(group=ModelExcludeStudy,color=ModelExcludeStudy),linetype="dashed",size=1)+
  theme_few() + theme(legend.position = "top",legend.title = element_blank())+
  xlab("No. of features used")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size=15))+ylim(0.5,1)

ggsave(p,filename = "../Figure/FigS3.pdf",height = 4,width = 4)

########### Fig S5 - Min ##########
setwd("K:/CRC-Pair/Unique.Pair.Permutation/PATHWAY.PERMUTATION")
AllPathway <- read.csv("All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")
AllPathway$Incre.aveRank.P <- AllPathway$Increasing.Rank.Average/10001
AllPathway2 <- AllPathway %>% mutate(Enriched = if_else(Incre.aveRank.P <= 0.05,"Ctrl",if_else(Decre.aveRank.P <= 0.05,"Disease","N.S.")))
AllPathway3 <- AllPathway2 %>% mutate(log2FC = log2(as.numeric(Dismean)/as.numeric(Ctlmean)))
AllPathway4 <- AllPathway3 %>% filter(Shape == "A")
## Dataset
FigS5Data <- data.frame()
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (name in Study) {
  index <- which(Study == name)
  middata <- read.csv(paste(name,"-Humann2-wilcoxonTest.csv",sep = ''))
  #Wilcon Test
  FigS5Data  <- na.omit(middata) %>% data.frame() %>% dplyr::arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0)) %>% 
    filter(Pathway %in% AllPathway4$Species) %>% dplyr::select(Pathway,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".1",sep = ''),Test="Wilcoxon") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(FigS5Data)
  #Pair
  middata <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
  FigS5Data.mid  <- na.omit(middata) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min == 1)))
  FigS5Data.mid$Incre.minRank.P <- FigS5Data.mid$Increasing.Rank.Min/10001
  FigS5Data.mid2 <- data.frame(Species = c(FigS5Data.mid$Species,FigS5Data.mid$Species),Pvalue = c(FigS5Data.mid$Incre.minRank.P,FigS5Data.mid$Decre.minRank.P))
  FigS5Data.mid3 <- FigS5Data.mid2%>% data.frame() %>% arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()
  
  
  FigS5Data <- FigS5Data.mid3%>% data.frame() %>% arrange(Pvalue) %>% mutate(Rank = round(rank(Pvalue,ties.method = "min"),0)) %>% 
    filter(Species %in% AllPathway4$Species) %>% mutate(Study = name) %>%
    select(Species,Study,Pvalue,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".2",sep = ''),Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(FigS5Data)
}

#Meta Wilcoxon
MetaData <- read.csv("All-8Study-Contained-Pathway-wilcoxon-test.csv")
FigS5Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0),Study="Meta") %>% 
  filter(Pathway %in% AllPathway4$Species) %>% select(Pathway,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.1",Test="Wilcoxon") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind( FigS5Data)

#Meta Pair
MetaData <- read.csv("All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")
FigS5Data.mid  <- na.omit(MetaData) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min == 1)))
FigS5Data.mid$Incre.minRank.P <- FigS5Data.mid$Increasing.Rank.Min/10001
FigS5Data.mid2 <- data.frame(Species = c(FigS5Data.mid$Species,FigS5Data.mid$Species),Pvalue = c(FigS5Data.mid$Incre.minRank.P,FigS5Data.mid$Decre.minRank.P))
FigS5Data.mid3 <- FigS5Data.mid2%>% data.frame() %>% arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()

FigS5Data <- FigS5Data.mid3%>% data.frame() %>% arrange(Pvalue) %>% mutate(Rank = round(rank(Pvalue,ties.method = "min"),0)) %>% 
  filter(Species %in% AllPathway4$Species) %>% mutate(Study="Meta") %>%
  dplyr::select(Species,Study,Pvalue,Rank) %>% 
  mutate(Country = "Meta.2",Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind(FigS5Data)

write.csv(FigS5Data,"../Figure/FigureS5.Data.csv")

### Figure
FigS5Data <- read.csv("../Figure/FigureS5.Data.csv",row.names = 1,stringsAsFactors = F)
Pdata <- FigS5Data %>% select(Species,Country,FDR,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata$Species <- factor(Pdata$Species,levels = AllPathway4$Species)
Pdata <- Pdata[order(Pdata$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,FDR) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata,"../Figure/FigS5.Left.csv")

Pdata.1 <- read.csv("../Figure/FigS5.Left.csv",row.names = 1,stringsAsFactors = F)
#Pdata.1 <- Pdata.1 %>% dplyr::arrange(Meta.2,Meta.1)
#Pdata <- t(Pdata)
#rownames(Pdata.1) <- rownames(Pdata.1) %>% str_replace_all("_"," ")
Pdata <- read.csv("../Figure/FigS5.Left.csv",row.names = 1,stringsAsFactors = F)
Pdata <- Pdata %>% rownames_to_column() %>%
  gather(key = "Condition",value = "FDR",-rowname)

labels = c("<0.001","0.001-0.01","0.01-0.05",">0.05")
breaks <- c(-1,0.001,0.01,0.05,1)
breaks2 <- c(0,001,0.01,0.05)
#mid <- data.frame(labels=labels,breaks=breaks2)
Pdata$labels = cut(Pdata$FDR,breaks,labels,ordered_result = T)
#Pdata<-merge(Pdata,mid,by= "labels")
#Pdata$alpha <- Pdata$FDR/Pdata$breaks
Pdata$rowname <-factor(Pdata$rowname,levels = AllPathway4$Species)
interval.cols <- c("#8B1A1A","#FF6A6A","#FFC1C1","#EEE9E9")#brewer.pal(6,"Set2")
names(interval.cols) <- levels(Pdata$labels)

Pdata.2 <- Pdata %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  mutate(rowname = factor(rowname,levels = AllPathway4$Species)) %>% dplyr::arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))

## text note => False
TextFunc <- function(dat, col = "black", fontsize = 9, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}
###
col_cat <- c("<0.001"="#8B1A1A","0.001-0.01"="#FF6A6A","0.01-0.05"="#FFC1C1",">0.05"="#EEE9E9")

### FigS5 left
Country.Heatmap <- c("AUS","CHI","FRA","GER","ITA1","ITA2","JAP","USA","Meta")
Pdata.2.mid1 <- Pdata.2 %>% select(ends_with(".1")) %>% select(paste(Country.Heatmap,".1",sep = ''))
rownames(Pdata.2.mid1) <- rownames(Pdata.2.mid1) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid1) = Country.Heatmap
#Pdata.2.mid1 <- t(Pdata.2.mid1)

Pdata.1.mid1 <- Pdata.1 %>% select(ends_with(".1")) %>% select(paste(Country.Heatmap,".1",sep = ''))
rownames(Pdata.1.mid1) <- rownames(Pdata.1.mid1) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid1) = Country.Heatmap
#Pdata.1.mid1 <- t(Pdata.1.mid1)


pdf("../Figure/FigS5.Left.pdf",width = 7,height = 4)
Heatmap(Pdata.2.mid1, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid1)
)
dev.off()

### FigS5 Right
Pdata.2.mid2 <- Pdata.2 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.2.mid2) <- rownames(Pdata.2.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid2) = Country.Heatmap
#Pdata.2.mid2 <- t(Pdata.2.mid2)

Pdata.1.mid2 <- Pdata.1 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.1.mid2) <- rownames(Pdata.1.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid2) = Country.Heatmap
#Pdata.1.mid2 <- t(Pdata.1.mid2)


pdf("../Figure/FigS5.Right.pdf",width = 7,height = 4)
Heatmap(Pdata.2.mid2, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid2)
)
dev.off()


########### Fig S5 - AVE #####################
setwd("K:/CRC-Pair/Unique.Pair.Permutation/PATHWAY.PERMUTATION")
AllPathway <- read.csv("All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")
AllPathway$Incre.aveRank.P <- AllPathway$Increasing.Rank.Average/10001
AllPathway2 <- AllPathway %>% mutate(Enriched = if_else(Incre.aveRank.P <= 0.05,"Ctrl",if_else(Decre.aveRank.P <= 0.05,"Disease","N.S.")))
AllPathway3 <- AllPathway2 %>% mutate(log2FC = log2(as.numeric(Dismean)/as.numeric(Ctlmean)))
AllPathway4 <- AllPathway3 %>% filter(Shape == "A")
## Dataset
FigS5Data <- data.frame()
Study <- c("FengQ_2015","ThomasAM_2018a","ThomasAM_2018b","VogtmannE_2016","YuJ_2015","ZellerG_2014","PRJDB4176","PRJEB27928")
StudyPos <- c("AUS","ITA1","ITA2","USA","CHI","FRA","JAP","GER")

for (name in Study) {
  index <- which(Study == name)
  middata <- read.csv(paste(name,"-Humann2-wilcoxonTest.csv",sep = ''))
  #Wilcon Test
  FigS5Data  <- na.omit(middata) %>% data.frame() %>% dplyr::arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0)) %>% 
    filter(Pathway %in% AllPathway4$Species) %>% dplyr::select(Pathway,Study,adj.fdr.p,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".1",sep = ''),Test="Wilcoxon") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(FigS5Data)
  #Pair
  middata <- read.csv(paste(name,"-Humann2-PairWilcoxonSign.csv",sep = ''),stringsAsFactors = F)
  FigS5Data.mid  <- na.omit(middata) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min == 1)))
  FigS5Data.mid$Incre.aveRank.P <- FigS5Data.mid$Increasing.Rank.Average/10001
  FigS5Data.mid2 <- data.frame(Species = c(FigS5Data.mid$Species,FigS5Data.mid$Species),Pvalue = c(FigS5Data.mid$Incre.aveRank.P,FigS5Data.mid$Decre.aveRank.P))
  FigS5Data.mid3 <- FigS5Data.mid2%>% data.frame() %>% arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()
  
  
  FigS5Data <- FigS5Data.mid3%>% data.frame() %>% arrange(Pvalue) %>% mutate(Rank = round(rank(Pvalue,ties.method = "min"),0)) %>% 
    filter(Species %in% AllPathway4$Species) %>% mutate(Study = name) %>%
    select(Species,Study,Pvalue,Rank) %>% 
    mutate(Country = paste(StudyPos[index],".2",sep = ''),Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
    rbind(FigS5Data)
}

#Meta Wilcoxon
MetaData <- read.csv("All-8Study-Contained-Pathway-wilcoxon-test.csv")
FigS5Data  <- na.omit(MetaData) %>% data.frame() %>% arrange(adj.fdr.p) %>% mutate(Rank = round(rank(adj.fdr.p,ties.method = "min"),0),Study="Meta") %>% 
  filter(Pathway %in% AllPathway4$Species) %>% select(Pathway,Study,adj.fdr.p,Rank) %>% 
  mutate(Country = "Meta.1",Test="Wilcoxon") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind( FigS5Data)

#Meta Pair
MetaData <- read.csv("All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")
FigS5Data.mid  <- na.omit(MetaData) %>% data.frame() %>% filter(!((Increasing.Rank.Min == 1) & (Decreasing.Rank.Min == 1)))
FigS5Data.mid$Incre.aveRank.P <- FigS5Data.mid$Increasing.Rank.Average/10001
FigS5Data.mid2 <- data.frame(Species = c(FigS5Data.mid$Species,FigS5Data.mid$Species),Pvalue = c(FigS5Data.mid$Incre.aveRank.P,FigS5Data.mid$Decre.aveRank.P))
FigS5Data.mid3 <- FigS5Data.mid2%>% data.frame() %>% arrange(Species) %>% group_by(Species) %>% top_n(-1,Pvalue) %>% unique()

FigS5Data <- FigS5Data.mid3%>% data.frame() %>% arrange(Pvalue) %>% mutate(Rank = round(rank(Pvalue,ties.method = "min"),0)) %>% 
  filter(Species %in% AllPathway4$Species) %>% mutate(Study="Meta") %>%
  dplyr::select(Species,Study,Pvalue,Rank) %>% 
  mutate(Country = "Meta.2",Test="Pair") %>% dplyr::rename(Species=1,Study=2,FDR=3,Rank=4,Country=5,Test=6) %>%
  rbind(FigS5Data)

write.csv(FigS5Data,"../Figure/FigureS5.AVE.Data.csv")

### Figure
FigS5Data <- read.csv("../Figure/FigureS5.AVE.Data.csv",row.names = 1,stringsAsFactors = F)
Pdata <- FigS5Data %>% select(Species,Country,FDR,Test) %>% arrange(Country) %>% arrange(desc(Test))
Pdata$Species <- factor(Pdata$Species,levels = AllPathway4$Species)
Pdata <- Pdata[order(Pdata$Species),] %>% data.frame() %>% select(-Test) %>% spread(Country,FDR) %>% remove_rownames() %>%
  column_to_rownames("Species")
write.csv(Pdata,"../Figure/FigS5.Left.AVE.csv")

Pdata.1 <- read.csv("../Figure/FigS5.Left.AVE.csv",row.names = 1,stringsAsFactors = F)
#Pdata.1 <- Pdata.1 %>% dplyr::arrange(Meta.2,Meta.1)
#Pdata <- t(Pdata)
#rownames(Pdata.1) <- rownames(Pdata.1) %>% str_replace_all("_"," ")
Pdata <- read.csv("../Figure/FigS5.Left.AVE.csv",row.names = 1,stringsAsFactors = F)
Pdata <- Pdata %>% rownames_to_column() %>%
  gather(key = "Condition",value = "FDR",-rowname)

labels = c("<0.001","0.001-0.01","0.01-0.05",">0.05")
breaks <- c(-1,0.001,0.01,0.05,1)
breaks2 <- c(0,001,0.01,0.05)
#mid <- data.frame(labels=labels,breaks=breaks2)
Pdata$labels = cut(Pdata$FDR,breaks,labels,ordered_result = T)
#Pdata<-merge(Pdata,mid,by= "labels")
#Pdata$alpha <- Pdata$FDR/Pdata$breaks
Pdata$rowname <-factor(Pdata$rowname,levels = AllPathway4$Species)
interval.cols <- c("#8B1A1A","#FF6A6A","#FFC1C1","#EEE9E9")#brewer.pal(6,"Set2")
names(interval.cols) <- levels(Pdata$labels)

Pdata.2 <- Pdata %>% select(labels,rowname,Condition) %>% spread(Condition,labels) %>% 
  mutate(rowname = factor(rowname,levels = AllPathway4$Species)) %>% dplyr::arrange(rowname) %>%
  remove_rownames() %>% column_to_rownames("rowname") %>%
  select(colnames(Pdata.1))

## text note => False
TextFunc <- function(dat, col = "black", fontsize = 9, numdat = TRUE,digit = 2){
  if(numdat == TRUE){
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }else{
    function(j, i, x, y, width, height, fill){
      grid.text(sprintf("%.0e", dat[i, j]), x, y, gp = gpar(fontsize = fontsize, col  = col))
    }
  }}
###
col_cat <- c("<0.001"="#8B1A1A","0.001-0.01"="#FF6A6A","0.01-0.05"="#FFC1C1",">0.05"="#EEE9E9")

### FigS5 left
Country.Heatmap <- c("AUS","CHI","FRA","GER","ITA1","ITA2","JAP","USA","Meta")
Pdata.2.mid1 <- Pdata.2 %>% select(ends_with(".1")) %>% select(paste(Country.Heatmap,".1",sep = ''))
rownames(Pdata.2.mid1) <- rownames(Pdata.2.mid1) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid1) = Country.Heatmap
#Pdata.2.mid1 <- t(Pdata.2.mid1)

Pdata.1.mid1 <- Pdata.1 %>% select(ends_with(".1")) %>% select(paste(Country.Heatmap,".1",sep = ''))
rownames(Pdata.1.mid1) <- rownames(Pdata.1.mid1) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid1) = Country.Heatmap
#Pdata.1.mid1 <- t(Pdata.1.mid1)


pdf("../Figure/FigS5.Left.AVE.pdf",width = 7,height = 4)
Heatmap(Pdata.2.mid1, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid1)
)
dev.off()

### FigS5 Right
Pdata.2.mid2 <- Pdata.2 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.2.mid2) <- rownames(Pdata.2.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.2.mid2) = Country.Heatmap
#Pdata.2.mid2 <- t(Pdata.2.mid2)

Pdata.1.mid2 <- Pdata.1 %>% select(ends_with(".2")) %>% select(paste(Country.Heatmap,".2",sep = ''))
rownames(Pdata.1.mid2) <- rownames(Pdata.1.mid2) %>% str_replace_all("_"," ")
colnames(Pdata.1.mid2) = Country.Heatmap
#Pdata.1.mid2 <- t(Pdata.1.mid2)


pdf("../Figure/FigS5.Right.AVE.pdf",width = 7,height = 4)
Heatmap(Pdata.2.mid2, rect_gp = gpar(lwd = 1, col = "black"), 
        name = "FDR",
        col = col_cat,
        show_row_names = T,
        show_column_names = T,na_col="white",
        cell_fun = TextFunc(Pdata.1.mid2)
)
dev.off()

#### Figure 4 => boxplot ####
setwd("D:/CRC-Pair/Pathway-8Study-20201010/Res")
PairData <- read.csv("../All-8Study-Contained.Pathway.pair.csv")
AbunData <- read.table("../EightStudies-PathwayAbundance-Group.txt",row.names = 1,header = T,sep = '\t')

#PWY.5667..CDP.diacylglycerol.biosynthesis.I #Pathway223
#PWY0.1319..CDP.diacylglycerol.biosynthesis.II #Pathway467
#PWY.4984..urea.cycle #Pathway162
#PWY.6123..inosine.5..phosphate.biosynthesis.I #Pathway275
#PWY.6124..inosine.5..phosphate.biosynthesis.II #Pathway276
#PWY.7234..inosine.5..phosphate.biosynthesis.III #Pathway403
middata <- AbunData[rownames(AbunData) %in% c(unique(PairData$Ctl),unique(PairData$Disease)),]
middata2 <- middata[,c(1,2,225,164,277,278,405,469)]

p1<-ggplot(middata2,aes(x=study_condition,y=PWY.5667..CDP.diacylglycerol.biosynthesis.I,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "CDP-diacylglycerol biosynthesis I")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p1,filename = "../../Figure/Figure4b.CDP.pdf",width = 3.5,height = 3)

p6<-ggplot(middata2,aes(x=study_condition,y=PWY0.1319..CDP.diacylglycerol.biosynthesis.II,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "CDP-diacylglycerol biosynthesis II")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p6,filename = "../../Figure/Figure4b.CDPII.pdf",width = 3.5,height = 3)

p2<-ggplot(middata2,aes(x=study_condition,y=PWY.4984..urea.cycle,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "urea cycle")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p2,filename = "../../Figure/Figure4b.ureacycle.pdf",width = 3.5,height = 3)

p3<-ggplot(middata2,aes(x=study_condition,y=PWY.6123..inosine.5..phosphate.biosynthesis.I,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis I")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p3,filename = "../../Figure/Figure4b.IMP1.pdf",width = 3.5,height = 3)

p4<-ggplot(middata2,aes(x=study_condition,y=PWY.6124..inosine.5..phosphate.biosynthesis.II,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis II")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p4,filename = "../../Figure/Figure4b.IMP2.pdf",width = 3.5,height = 3)

p5<-ggplot(middata2,aes(x=study_condition,y=PWY.7234..inosine.5..phosphate.biosynthesis.III,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis III")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p5,filename = "../../Figure/Figure4b.IMP3.pdf",width = 3.5,height = 3)

## Sum the same pathway
middata2$"IMP biosynthesis I+II+III" <- middata2$PWY.6123..inosine.5..phosphate.biosynthesis.I+
  middata2$PWY.6124..inosine.5..phosphate.biosynthesis.II+
  middata2$PWY.7234..inosine.5..phosphate.biosynthesis.III
  
p8<-ggplot(middata2,aes(x=study_condition,y=`IMP biosynthesis I+II+III`,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "IMP biosynthesis I+II+III")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))#+stat_compare_means(comparisons=c("control","CRC"),label = "p.value",hide.ns = F)

ggsave(p8,filename = "../../Figure/Figure4b.IMP1+2+3.pdf",width = 3.5,height = 3)

middata2$"CDP-diacylglycerol biosynthesis I+II" = middata2$PWY.5667..CDP.diacylglycerol.biosynthesis.I + 
  middata2$PWY0.1319..CDP.diacylglycerol.biosynthesis.II

p7<-ggplot(middata2,aes(x=study_condition,y=`CDP-diacylglycerol biosynthesis I+II`,color=study_condition)) +
  geom_boxplot()+
  geom_jitter(width = 0.1,height = 0,size=0.3)+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Relative Abundance",title = "CDP-diacylglycerol biosynthesis I+II")+
  theme(axis.text = element_text(size=13),axis.title = element_text(size=15))+
  scale_color_manual(values = c("#1E90FF","#DC143C"))

ggsave(p7,filename = "../../Figure/Figure4b.CDP1+2.pdf",width = 3.5,height = 3)


###########################################################################################################################


###########################################################################################################################

###########################################################################################################################



cal_w <- function(ctrl,crc){
  Difference <- ctrl - crc
  ranks <- rank(abs(Difference),ties.method = "average")
  signs <- sign(Difference)
  sign.rank <- signs*ranks
  W.P <- sum(sign.rank[sign.rank>0])
  return(W.P)
}




wilcox.test(c(85,70,40,65,80,75,55,20),c(75,50,50,40,20,65,40,25),paired = T)


sign(c(85,70,40,65,80,75,55,20)-c(75,50,50,40,20,65,40,25))


cal_w(c(85,70,40,65,80,75,55,20),c(75,50,50,40,20,65,40,25))
cal_w(c(118,134,130,124,105,130,130,132,123,128,126,140,135,126,132),c(125,132,138,120,125,127,136,139,131,132,135,136,128,127,130))
cal_w(Middata.mid$Ctrl,Middata.mid$CRC)



for (i in 1:10) {
  for (j in 1:8) {
    cat(i)
    if (j==4) {
      break()
    }
  }
}













