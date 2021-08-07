###############################################################################################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")


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

#### Pair - Shuffle Relative Abundance -> less one interation
pair_find<-function(data=data,phenodata=data.frame(),k="euclidean",SavePath = NULL,ShuffleWstat = NULL, BoundarySample = NULL,BoundaryPair=NULL,ShuffleTime=10000,DownPercent = 0.2,Uppercent=0.8){ # colnames(phenodata) = c("id","grp"); grp1 = "Ctrl", grp2 = "Disease"
  suppressMessages(library(tidyverse))
  #suppressMessages(library(fdrtool))
  #suppressMessages(library(qvalue))
  phenodata <- phenodata %>% dplyr::arrange(grp) %>% data.frame()
  #data must be ordered as grp order
  groupnum <- table(phenodata$grp) %>% data.frame()
  RAN_num = groupnum[1,2]
  RAP_num = groupnum[2,2]
  RAN<-data[1:RAN_num,]
  RAP<-data[(RAN_num+1):(RAN_num+RAP_num),]
  n=dim(data)[1]
  num=floor(sqrt(RAN_num+RAP_num))
  results_whole=matrix(nrow=0,ncol=(num+1))
  for (i in 1:n) {
    #cat(i)
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
  pairinfor <- pairinfor %>% data.frame() %>% dplyr::rename(Ctl=1,Disease=2,Distance=3) %>% mutate(Distance = as.numeric(as.character(Distance))) %>% dplyr::arrange(Ctl,Disease) %>% unique()
  #View(pairinfor)
  pairinfor <- Extract_Dist(pairinfor)
  pairinfor <- pairinfor %>% data.frame() %>% dplyr::select(-Distance)
  if (!is.null(BoundaryPair)) {
    write.csv(pairinfor,paste(SavePath,"/",BoundaryPair,".csv",sep = ''),row.names = F)
  }
  #write.csv(pairinfor,file = "test.csv",row.names = F)
  cat(paste("the redundant pair number is ",dim(pairinfor)[1],"\n",sep = ''))
  cat("PAIR FINISHED\n")
  #View(pairinfor)
  #return(pairinfor)
  
  FinalMatrix <- matrix(ncol = 20002,nrow = 0) #Species, Original W, 
  for (i in 1:dim(data)[2]) {
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
    
    ShuffleStat <- rep(NA,10000)
    cat(as.character(colnames(data)[i]))
    #cat("\n")
    #print(AllNW)
    for (Iter in 1:ShuffleTime) {
      Random <- runif(1,DownPercent,Uppercent)
      n <- round(dim(pairinfor)[1]*Random,0)
      RandomIndex <- sample(1:dim(pairinfor)[1],n,replace = F)
      
      Middata.mid1.1 <- Middata[-RandomIndex,] %>% data.frame()
      Middata.mid2 <- Middata[RandomIndex,] %>% data.frame()
      Middata.mid2.1 <- data.frame(Ctrl = Middata.mid2$CRC, CRC = Middata.mid2$Ctrl)
      Middata.mid <- rbind(Middata.mid1.1,Middata.mid2.1) %>% data.frame()
      
      test1<-wilcox.test(Middata.mid$Ctrl,Middata.mid$CRC,paired = TRUE)
      #Output[[as.character(colnames(middata1)[i])]][["Shuffle"]][Iter] = AllNW - test1$statistic
      ShuffleStat[Iter] = AllNW - test1$statistic
      #cat(test1$statistic)
      #cat("\n")
    }
    FinalMatrix <- rbind(FinalMatrix,c(as.character(colnames(data)[i]),OriginalStat,ShuffleStat))
    #cat(as.character(colnames(data)[i]))
    
    #cat("\n")
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
    dplyr::rename(Increasing.Rank.Min=1,Decreasing.Rank.Min=2,
                  Increasing.Rank.Max=3,Decreasing.Rank.Max=4,
                  Increasing.Rank.Average=5,Decreasing.Rank.Average=6)
  
  Mid.Matrix$Decre.maxRank.P <- Mid.Matrix$Decreasing.Rank.Max/(ShuffleTime+1)
  Mid.Matrix$Decre.aveRank.P <- Mid.Matrix$Decreasing.Rank.Average/(ShuffleTime+1)
  Mid.Matrix$Decre.minRank.P <- Mid.Matrix$Decreasing.Rank.Min/(ShuffleTime+1)
  Mid.Matrix$Species <- FinalMatrix[,1] #Feature
  
  Mid.Matrix <- Mid.Matrix %>% data.frame() %>%
    dplyr::arrange(Decre.minRank.P) %>% 
    dplyr::mutate(Decre.minRank.P.FDR = p.adjust(.$Decre.minRank.P,method = "BH",n=length(.$Decre.minRank.P))) %>%
    dplyr::arrange(Decre.maxRank.P) %>% 
    dplyr::mutate(Decre.maxRank.P.FDR = p.adjust(.$Decre.maxRank.P,method = "BH",n=length(.$Decre.maxRank.P))) %>%
    dplyr::arrange(Decre.aveRank.P) %>% 
    dplyr::mutate(Decre.aveRank.P.FDR = p.adjust(.$Decre.aveRank.P,method = "BH",n=length(.$Decre.aveRank.P)))
  
  cat("All done\n")
  return(Mid.Matrix)
}# END - function: Pair_Find

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


#################################### 



#################################### 


#################################### 


#################################### 


#################################### 



####################### Fig 4 ###########################
setwd("K:/CRC-Pair/Unique.Pair.Permutation")
AllPathway <- read.csv("PATHWAY.PERMUTATION/All-8Study-Contained-Pathway-pair-wilcoxonsign-res.csv")

AllPathway$Incre.aveRank.P <- AllPathway$Increasing.Rank.Average/10001
AllPathway2 <- AllPathway %>% mutate(Enriched = if_else(Incre.aveRank.P <= 0.01,"Ctrl",if_else(Decre.aveRank.P <= 0.01,"Disease","N.S.")))
AllPathway3 <- AllPathway2 %>% mutate(log2FC = log2(as.numeric(Dismean)/as.numeric(Ctlmean)))

#mid <- data.frame(Pvalue = c(AllPathway3$Incre.aveRank.P,AllPathway3$Decre.aveRank.P),Species = c(AllPathway3$Species,AllPathway3$Species))
#middataTop30 <- mid %>% arrange(Pvalue) %>% top_n(-30,Pvalue)











