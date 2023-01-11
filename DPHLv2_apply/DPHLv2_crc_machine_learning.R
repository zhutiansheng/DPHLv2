###detdiff prot
rf<-read.csv("reviewedfull_diff_prot.csv",header = T)
lib3<-read.csv("3lib_diff_prot.csv",header = T)
useprot<-lib3[lib3$X%in%setdiff(lib3$X,rf$prot),]
write.csv(useprot,"dphl2_3lib_spec_diff_prot.csv",row.names = F)



#####################################################################
#四个库over的蛋白做机器学习

file_list_up<-list.files(path = "E:/Qsync/work/DPHL2/0415_DPHL2_apply_new/crcd/analysis4_90/analysis_90_no_gene_name/4lib_diff_prot/",pattern = "protUP")
prot_name_up<-c()
for (i in file_list_up) {
  indsel<-read.csv(i,header = T)
  prot_name_up<-c(prot_name_up,indsel$X)
  
}       
prot_num_up<-data.frame(table(prot_name_up))

prot_ovlap_up<-prot_num_up[prot_num_up$Freq==4,]


file_list_down<-list.files(path = "E:/Qsync/work/DPHL2/0415_DPHL2_apply_new/crcd/analysis4_90/analysis_90_no_gene_name/4lib_diff_prot/",pattern = "protDown_")
prot_name_down<-c()
for (i in file_list_down) {
  indsel<-read.csv(i,header = T)
  prot_name_down<-c(prot_name_down,indsel$X)
  
}       
prot_num_down<-data.frame(table(prot_name_down))

prot_ovlap_down<-prot_num_down[prot_num_down$Freq==4,]

###########################################################################
#用reviewedfull的样本做预测
refull<-read.csv("reviewedfull_uni_matrix.csv",header = T)

protmat_up<-refull[,colnames(refull)%in%prot_ovlap_up$prot_name_up]
protmat_down<-refull[,colnames(refull)%in%prot_ovlap_down$prot_name_down]
protmat<-cbind(protmat_up,protmat_down)
row.names(protmat)<-refull$X
protmat$label<-refull$pat
protmat$label<-gsub("N",0,protmat$label)
protmat$label<-gsub("P",1,protmat$label)
protmat$label<-as.factor(protmat$label)
#####################################################################################
#机器学习

#随机森林

library(readr)
library(readxl)
library(stringr)
library(magrittr)
library("randomForest")
library(pROC)
set.seed(123)
nn<-sample(241, 200, replace = FALSE, prob = NULL)
train_samp<-protmat[nn,]
test_samp<-protmat[-nn,]
Data_ML <- as.data.frame(train_samp[,c(1427,1:1426)])
train_set<-Data_ML



source("datamining_library_ge20200306.R")

###################################################


#-------------------------------------
#验证模型
library(caret)
accu <- c()
for (i in seq(4,6,0.5)) {
  set.seed(2022.7)
  tmpRF2 <- randomForest(label ~ . ,data=train_set,importance=T,ntree=1000,nodesize=5)
  result <- data.frame(importance(tmpRF2,type=1))
  result1 <- row.names(result)[result$MeanDecreaseAccuracy>i]
  set.seed(2022.7)
  for (seed in runif(10,1,1000)) {
    # for (fd in 3:6) {
    set.seed(seed)
    folds <- createFolds(train_set$label,5)
    n=0
    for(fold in folds){
      n=n+1
      #fold=folds[[8]]
      valids <- train_set[fold,result1]
      valids$label <- train_set$label[fold]
      trains <- train_set[setdiff(1:dim(train_set)[1],fold),result1]
      trains$label <- train_set$label[setdiff(1:dim(train_set)[1],fold)]
      trains$label <- as.factor(trains$label)
      # for (ntree in seq(600,1000,200)) {
      set.seed(2022.7)
      tmpRF <- randomForest(as.factor(label) ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
      fea <- data.frame(importance(tmpRF,type=1))
      for (dec in seq(4,6,0.5)) {
        feature <- row.names(fea)[fea$MeanDecreaseAccuracy>dec]
        if(length(feature)>1){
          train2 <- trains[,feature]
          train2$label <- trains$label
          tmpRF3 <- randomForest(as.factor(label) ~ . ,data=train2,importance=T,ntree=1000,nodesize=5)
          predicted <- predict(tmpRF3,valids,type='prob')
          predicts <- t(apply(predicted,1,function(v){v/sum(v)}))
          colnames(predicts) <- colnames(predicted)
          predicts <- data.frame(predicts,check.names=F)
          predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
          predicts$observed <- valids$label
          ROC <- roc(predicts$observed, as.numeric(predicts$`1`))
          auc <- as.numeric(ROC[["auc"]])
          acc <- sum(predicts$predicted==predicts$observed)
          accu <- rbind(accu,c(i,seed,5,n,1000,dec,acc/length(fold),auc))
          
          
        }
      }
    }
  }
}
write.csv(accu,paste0("DPHL2_crcd_revfull_4lib_overlap_diffprot_rf_grid_search_ntree1000_20221219.csv"),row.names = F)
# accu<-read.csv("DPHL2_crcd_revfull_4lib_overlap_diffprot_rf_grid_search_ntree1000.csv",header = T)
accu<-data.frame(accu)
# accu_fd5<-accu[accu$n==2,]
colnames(accu)<-c("i","seed","fd","n","ntree","dec","acc/length(fold)","auc")
mean <- aggregate(accu,by=list(accu$i,accu$seed,accu$fd,accu$ntree,accu$dec),mean)
tmp3 <- mean[which(mean[,12]>0.98 & mean[,13]>0.98), ]
tmp3$label<-paste(tmp3$i,tmp3$seed,tmp3$dec,sep = "_")
accu$label<-paste(accu$i,accu$seed,accu$dec,sep = "_")
for (jj in 1:nrow(tmp3)) {
  indsel<-accu[accu$label==tmp3$label[jj],]
  indsel1<-dplyr::arrange(indsel,dplyr::desc(`acc/length(fold)`),dplyr::desc(auc))
  tmp3$n[jj]<-indsel1$n[1]
}


accu3 <- c()
for (j in 1:nrow(tmp3)) {
  #j=1
  set.seed(2022.7)
  tmpRF2 <- randomForest(as.factor(label) ~ . ,data=train_set,importance=T,ntree=1000,nodesize=5)
  result <- data.frame(importance(tmpRF2,type=1))
  result1 <- row.names(result)[result$MeanDecreaseAccuracy>tmp3[j,6]]
  
  
  set.seed(tmp3[j,7])
  folds <- createFolds(train_set$label,tmp3[j,8])
  fold=folds[[tmp3[j,9]]]
  valids <- train_set[fold,result1]
  valids$label <- train_set$label[fold]
  trains <- train_set[setdiff(1:dim(train_set)[1],fold),result1]
  trains$label <- train_set$label[setdiff(1:dim(train_set)[1],fold)]
  trains$label <- as.factor(trains$label)
  set.seed(2022.7)
  tmpRF <- randomForest(as.factor(label) ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
  
  fea <- data.frame(importance(tmpRF,type=1))
  feature <- row.names(fea)[fea$MeanDecreaseAccuracy>tmp3[j,11]]
  train2 <- trains[,feature]
  train2$label <- trains$label
  tmpRF3 <- randomForest(as.factor(label) ~ . ,data=train2,importance=T,ntree=1000,nodesize=5)
  
  # for (time in 1:100) {
  # sm <- sample(1:42,30)
  # valids2 <- train_set[sm,]
  # valids2_label <- train_set$label[sm]
  valids2<-valids[,feature]
  valids2_label <- valids$label
  predicted <- predict(tmpRF3,valids2,type='prob')
  predicts <- t(apply(predicted,1,function(v){v/sum(v)}))
  colnames(predicts) <- colnames(predicted)
  predicts <- data.frame(predicts,check.names=F)
  predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
  predicts$observed <- valids2_label
  # write.csv(predicts,"ovc_prot_69_valid_pred.csv",row.names = T)
  train3<-train2[,1:(ncol(train2)-1)]
  predicted_train <- predict(tmpRF3,train3,type='prob')
  predicts_train <- t(apply(predicted_train,1,function(v){v/sum(v)}))
  colnames(predicts_train) <- colnames(predicted_train)
  predicts_train <- data.frame(predicts_train,check.names=F)
  predicts_train$predicted <- apply(predicts_train,1,function(v){names(v)[max(v)==v]})
  predicts_train$observed <- train2$label
  # write.csv(predicts_train,"ovc_prot_69_train_pred.csv",row.names = T)
  acc <- sum(predicts$predicted==predicts$observed)
  accu3 <- rbind(accu3,c(tmp3[j,6],tmp3[j,7],tmp3[j,8],tmp3[j,9],tmp3[j,10],tmp3[j,11],acc/length(valids2_label)))
  
  
}
accu3<-data.frame(accu3)
names(accu3)<-c(names(tmp3[6:11]),"acc")
write.csv(accu3,paste0("DPHL2_crcd_revfull_4lib_overlap_diffprot_rf_grid_search_select_ntree1000_20221219.csv"),row.names = F)

# tmp5 <- read.csv("test_120severeallresultv4.csv")
# accu4 <- accu3[,-1]
# mean <- aggregate(accu3,by=list(accu3[,1],accu3[,2],accu3[,3],accu3[,4],accu3[,5],accu3[,6]),mean)
# sd <- aggregate(accu3,by=list(accu3[,1],accu3[,2],accu3[,3],accu3[,4],accu3[,5],accu3[,6]),sd)

# best <- mean[which(mean$V6==1 & sd$V6 ==0),]


###########################


test_set1<-test_samp[,c(1427,1:1426)]
# colnames(test_set1)[1]<-"label"





sum <-roc<-accu4<-sum_train<-roc_train<- c()

for(nn in 1:nrow(accu3)){
  # nn=3
  # nn=17
  # nn=40
  set.seed(2022.7)
  tmpRF2 <- randomForest(as.factor(label) ~ . ,data=train_set,importance=T,ntree=1000,nodesize=5)
  result <- data.frame(importance(tmpRF2,type=1))
  result1 <- row.names(result)[result$MeanDecreaseAccuracy>accu3[nn,1]]
  	     write.csv(result1,"20221219_dphl2_severe_step1_feature_ntree1000_40.csv")
  	     pdf("20221219_dphl2_RF_step1_severeimportant_select_ntree1000_40.pdf")
          varImpPlot(tmpRF2,n.var=min(length(result1), nrow(tmpRF2$importance)))
         dev.off()
  
  set.seed(accu3[nn,2])
  folds <- createFolds(train_set$label,accu3[nn,3])
  fold=folds[[accu3[nn,4]]]
  valids <- train_set[fold,result1]
  valids$label <- train_set$label[fold]
  trains <- train_set[setdiff(1:dim(train_set)[1],fold),result1]
  trains$label <- train_set$label[setdiff(1:dim(train_set)[1],fold)]
  trains$label <- as.factor(trains$label)
  set.seed(2022.7)
  tmpRF <- randomForest(as.factor(label) ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
  fea <- data.frame(importance(tmpRF,type=1))
  feature <- row.names(fea)[fea$MeanDecreaseAccuracy>accu3[nn,6]]
  	     write.csv(feature,"20221219_dphl2_severe_step2_feature_ntree1000_40.csv")
  pdf("20221912_dphl2_RF_step2_severeimportant_select_ntree1000_40.pdf")
  varImpPlot(tmpRF,n.var=min(length(feature), nrow(tmpRF$importance)))
  dev.off()
  dmc1<-fea[fea$MeanDecreaseAccuracy>accu3[nn,6],]
  dmc1<-cbind(feature,dmc1)
  # write.csv(dmc1,"dphl2_MeanDecreaseAccuracy_dia_feature_ntree1000_24.csv",row.names = F)
  train2 <- trains[,feature]
  train2$label<-trains$label
  # train2_2<-train_mat[,feature]
  # sum(is.na(train2_2))/(ncol(train2_2)*nrow(train2_2))
  # train2_2$label <- train_mat$grouping
  # write.csv(train2_2,"zzovc_train_feature.csv",row.names = F)
  test_set2 <-test_set1[,feature]
  sum(is.na(test_set2))/(ncol(test_set2)*nrow(test_set2))
  # test_set2$label<-test_set1$label
  # write.csv(test_set,"zzovc_test_feature.csv",row.names = F)
  test_label<-as.factor(test_set1$label)
  # test_set2[is.na(test_set2)]<-9.6
  tmpRF3 <- randomForest(as.factor(label) ~ . ,data=train2,importance=T,ntree=1000,nodesize=5)
  ########################
  train3<-trains[,feature]
  train_label<-trains$label
  predicted_train <- predict(tmpRF3,train3,type='prob') 
  predicts_train <- t(apply(predicted_train,1,function(v){v/sum(v)}))
  colnames(predicts_train) <- colnames(predicted_train)
  predicts_train <- data.frame(predicts_train,check.names=F)
  predicts_train$predicted <- apply(predicts_train,1,function(v){names(v)[max(v)==v]})
  predicts_train$observed <- train_label
  # write.csv(predicts,"dphl2_prot_21_test_pred_ntree1000_24.csv",row.names = T)
  predicts_train$batch<-row.names(train3)
  # write.csv(predicts_train,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","dphl2_acc_tabel_24_ntree1000.csv"))
  acc_train <- sum(predicts_train$predicted==predicts_train$observed)
  acc_train/length(train_label)
  sum_train <- c(sum_train,acc/length(train_label))
  
  ROC_train <- roc(predicts_train$observed,predicts_train$`1`)
  roc_train<-c(roc_train,ROC_train[["auc"]])
  pdf(paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","ROC_severe_dphl2_40_ntree1000_20221219_auc1_train_acc0.9875776.pdf"))
  plot.roc(ROC_train,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
           main="RF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
  dev.off()
  # write.csv(feature,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","severe_feature_dphl2_24_ntree1000.csv"))
  ########################
  predicted <- predict(tmpRF3,test_set2,type='prob') 
  predicts <- t(apply(predicted,1,function(v){v/sum(v)}))
  colnames(predicts) <- colnames(predicted)
  predicts <- data.frame(predicts,check.names=F)
  predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
  predicts$observed <- test_label
  # write.csv(predicts,"dphl2_prot_21_test_pred_ntree1000_24.csv",row.names = T)
  predicts$batch<-row.names(test_samp)
  # write.csv(predicts,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","dphl2_acc_tabel_24_ntree1000.csv"))
  acc <- sum(predicts$predicted==predicts$observed)
  acc/length(test_label)
  sum <- c(sum,acc/length(test_label))
  
  ROC <- roc(predicts$observed,predicts$`1`)
  roc<-c(roc,ROC[["auc"]])
  pdf(paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","ROC_severe_dphl2_40_ntree1000_auc0.9034_test_acc_0.9268293.pdf"))
  plot.roc(ROC,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
           main="RF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
  dev.off()
  # write.csv(feature,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","severe_feature_dphl2_24_ntree1000.csv"))
  
  accu4 <- rbind(accu4,c(accu3[nn,1],accu3[nn,2],accu3[nn,3],accu3[nn,4],accu3[nn,6],accu3[nn,7],acc_train/length(train_label),ROC_train[["auc"]],acc/length(valids2_label),ROC[["auc"]]))
  
}
# acc_auc<-cbind(sum,roc)
accu4<-data.frame(accu4)
names(accu4)<-c(names(accu3[c(1:4,6:7)]),"acc_train","ROC_train","acc_test","ROC_test")
write.csv(accu4,"test_acc_auc_new_group_42_tissue_70miss_0.8min_1_5fc_zscore_20221219.csv") 	   
####################################################################################
#train predict 散点图

predicted1<-cbind(train3,predicts_train)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
row.names(predicted1)<-predicted1$batch
pdf("20221219_dphl2_rf_train_pred_40_ntree1000.pdf",width = 5.5,height = 5)
plot(predicted1$`1`,predicted1$mean,xlim = c(0,1))
points(predicted1$`1`[predicted1$observed==0],predicted1$mean[predicted1$observed==0], col="red", pch=19, cex=1)
# text(predicted1$`1`[predicted1$observed==0],predicted1$mean[predicted1$observed==0], labels =predicted1$batch[predicted1$observed==0],col=1, bg = "grey", cex=0.8,font = 3)
points(predicted1$`1`[predicted1$observed==1],predicted1$mean[predicted1$observed==1], col = "blue", pch=19,cex=1)
# text(predicted1$`1`[predicted1$observed==1],predicted1$mean[predicted1$observed==1], labels =predicted1$batch[predicted1$observed==1],col=1, bg = "grey", cex=0.8,font = 3)
abline(v=0.5,lty=2,lwd=1)
dev.off()
####################################################################################
#test predict 散点图

predicted1<-cbind(test_set2,predicts)
predicted1$mean<-apply(predicted1[,1:14],1,mean)
row.names(predicted1)<-predicted1$batch
pdf("20221219_dphl2_rf_test_pred_40_ntree1000.pdf",width = 5.5,height = 5)
plot(predicted1$`1`,predicted1$mean,xlim = c(0,1))
points(predicted1$`1`[predicted1$observed==0],predicted1$mean[predicted1$observed==0], col="red", pch=19, cex=1)
# text(predicted1$`1`[predicted1$observed==0],predicted1$mean[predicted1$observed==0], labels =predicted1$batch[predicted1$observed==0],col=1, bg = "grey", cex=0.8,font = 3)
points(predicted1$`1`[predicted1$observed==1],predicted1$mean[predicted1$observed==1], col = "blue", pch=19,cex=1)
# text(predicted1$`1`[predicted1$observed==1],predicted1$mean[predicted1$observed==1], labels =predicted1$batch[predicted1$observed==1],col=1, bg = "grey", cex=0.8,font = 3)
abline(v=0.5,lty=2,lwd=1)
dev.off()

##################################################################################
#train sample 的acc roc
train3<-train2[,-ncol(train2)]
predicted_train <- predict(tmpRF3,train3,type='prob') 
predicts_train <- t(apply(predicted_train,1,function(v){v/sum(v)}))
colnames(predicts_train) <- colnames(predicted_train)
predicts_train <- data.frame(predicts_train,check.names=F)
predicts_train$predicted <- apply(predicts_train,1,function(v){names(v)[max(v)==v]})
predicts_train$observed <- train2$label
batch_label<-as.numeric(row.names(predicts_train))
predicts_train$batch<-train_mat$batchId[batch_label]
# write.csv(predicts,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","acc_tabel_70miss_0.8min_1_5fc_69_zscore.csv"))
acc <- sum(predicts_train$predicted==predicts_train$observed)
acc/length(train2$label)
sum <- c(sum,acc/length(train2$label))

ROC <- roc(predicts_train$observed,predicts_train$`1`)
roc<-c(roc,ROC[["auc"]])
pdf(paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","ROC_severe_70miss_0.8min_1_5fc_69_zscore_train.pdf"))
plot.roc(ROC,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="RF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()

####################################################################################
#validation sample 的acc roc
valids_3<-valids[,feature]
predicted_vali <- predict(tmpRF3,valids_3,type='prob') 
predicts_vali <- t(apply(predicted_vali,1,function(v){v/sum(v)}))
colnames(predicts_vali) <- colnames(predicted_vali)
predicts_vali <- data.frame(predicts_vali,check.names=F)
predicts_vali$predicted <- apply(predicts_vali,1,function(v){names(v)[max(v)==v]})
predicts_vali$observed <- valids$label
batch_label<-as.numeric(row.names(predicts_vali))
predicts_vali$batch<-train_mat$batchId[batch_label]
# write.csv(predicts,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","acc_tabel_70miss_0.8min_1_5fc_69_zscore.csv"))
acc <- sum(predicts_vali$predicted==predicts_vali$observed)
acc/length(valids$label)
sum <- c(sum,acc/length(valids$label))

ROC <- roc(predicts_vali$observed,predicts_vali$`1`)
roc<-c(roc,ROC[["auc"]])
pdf(paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","ROC_severe_70miss_0.8min_1_5fc_69_zscore_vali.pdf"))
plot.roc(ROC,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="RF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()






################################################################################
#用新的数据input模型中
vali_mat<-read.csv("SJOVA_65files.csv",header = T)
vali_samp_inf<-read.csv("validation_sample_inf.csv",header = T)
vali_pat_inf<-read.csv("validation_patient_inf.csv",header = T)
vali_mat<-merge(vali_mat,vali_samp_inf,by = "batch_id",all.x = T)
vali_mat<-merge(vali_mat,vali_pat_inf,by= "ID",all.x = T)
###############################################################
#技术重复取平均，在这里patientid一样的算技术重复
#暂时未做
######################################################
feature_prot<-read.csv("102_1_364.434086317429_5_1_1000_1.5_severe_feature_70miss_1_5fc_102_zscore.csv",header=T)
vali_set<-vali_mat[,colnames(vali_mat)%in%row.names(ipainput)]
sum(is.na(vali_set))/(ncol(vali_set)*nrow(vali_set))
vali_set[is.na(vali_set)]<-na_to
vali_set1<-data.frame(matrix(NA,nrow = nrow(vali_set),ncol = ncol(vali_set)))
names(vali_set1)<-names(vali_set)
for (i in 1:length(temp_mean)) {
  aa<-vali_set[,i]
  vali_set1[,i]<-sapply(aa,function(e){(e-temp_mean[i])/temp_sd[i]} )
  
}
vali_set2<-vali_set1[,colnames(vali_set1)%in%feature_prot$x]

vali_set_label<-vali_mat$grouping
# vali_set$grouping<-vali_mat$grouping
# write.csv(vali_set,"88_validation_feature.csv",row.names = F)
vali_set_label<-gsub("res",0,vali_set_label)
vali_set_label<-gsub("sen",1,vali_set_label)
vali_set_label<- as.factor(vali_set_label)

predicted1 <- predict(tmpRF3,vali_set2,type='prob') 
predicts1 <- t(apply(predicted1,1,function(v){v/sum(v)}))
colnames(predicts1) <- colnames(predicted1)
predicts1 <- data.frame(predicts1,check.names=F)
predicts1$predicted <- apply(predicts1,1,function(v){names(v)[max(v)==v]})
predicts1$observed <- vali_set_label
predicts1$batch<-vali_mat$batch_id
write.csv(predicts1,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","acc_tabel_validation_102.csv"))
acc <- sum(predicts1$predicted==predicts1$observed)
acc/length(vali_set_label)
# sum <- c(sum,acc/length(vali_set_label))

ROC <- roc(predicts1$observed,predicts1$`1`)
roc<-c(roc,ROC[["auc"]])
pdf(paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","ROC_severe_validation_102.pdf"))
plot.roc(ROC,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="RF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
# write.csv(feature,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","severe_feature_102.csv"))


################################################################################
#用新的数据input模型中_使用去批次效应版
vali_mat<-read.csv("SJOVA_65files_corrected.csv",header = T)
vali_samp_inf<-read.csv("validation_sample_inf.csv",header = T)
vali_pat_inf<-read.csv("validation_patient_inf.csv",header = T)
vali_mat<-merge(vali_mat,vali_samp_inf,by = "batch_id",all.x = T)
vali_mat<-merge(vali_mat,vali_pat_inf,by= "ID",all.x = T)
###############################################################
#技术重复取平均，在这里patientid一样的算技术重复
#暂时未做
######################################################
feature_prot<-read.csv("69_1.5_816.161676998017_5_4_1000_0.5_severe_feature_70miss_1_5fc_102_zscore.csv",header=T)
vali_set<-vali_mat[,colnames(vali_mat)%in%row.names(ipainput)]
sum(is.na(vali_set))/(ncol(vali_set)*nrow(vali_set))
vali_set[is.na(vali_set)]<-na_to
vali_set1<-data.frame(matrix(NA,nrow = nrow(vali_set),ncol = ncol(vali_set)))
names(vali_set1)<-names(vali_set)
for (i in 1:length(temp_mean)) {
  aa<-vali_set[,i]
  vali_set1[,i]<-sapply(aa,function(e){(e-temp_mean[i])/temp_sd[i]} )
  
}
vali_set2<-vali_set1[,colnames(vali_set1)%in%feature_prot$x]

vali_set_label<-vali_mat$grouping

vali_set_label<-gsub("res",0,vali_set_label)
vali_set_label<-gsub("sen",1,vali_set_label)
vali_set_label<- as.factor(vali_set_label)



predicted1 <- predict(tmpRF3,vali_set2,type='prob') 
predicts1 <- t(apply(predicted1,1,function(v){v/sum(v)}))
colnames(predicts1) <- colnames(predicted1)
predicts1 <- data.frame(predicts1,check.names=F)
predicts1$predicted <- apply(predicts1,1,function(v){names(v)[max(v)==v]})
predicts1$observed <- vali_set_label
predicts1$batch<-vali_mat$AD
write.csv(predicts1,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","acc_tabel_validation_corr_19.csv"))
acc <- sum(predicts1$predicted==predicts1$observed)
acc/length(vali_set_label)
# sum <- c(sum,acc/length(vali_set_label))

ROC <- roc(predicts1$observed,predicts1$`1`)
roc<-c(roc,ROC[["auc"]])
pdf(paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","ROC_severe_validation_corr_102.pdf"))
plot.roc(ROC,col = "blue3",ylim=c(0,1), print.auc=TRUE, print.thres="best",
         main="RF ROC",legacy.axes = TRUE,print.auc.cex=1.2)
dev.off()
# write.csv(feature,paste0(nn,"_",accu3[nn,1],"_",accu3[nn,2],"_",accu3[nn,3],"_",accu3[nn,4],"_",accu3[nn,5],"_",accu3[nn,6],"_","severe_feature.csv"))






















