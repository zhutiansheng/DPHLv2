rm(list = ls())
options(stringsAsFactors = F)
# first_category_name = list.files() 
crcdinf<-read.csv("20190522CRCD_sampleinf.csv",stringsAsFactors = F,header = T)
# genename<-read.csv("uniprot-human-gene_name.csv",stringsAsFactors = F,header = T)
# first_category_name<-c("reviewedfull","reviewedsemi","isoformfull","isoformsemi","dphl1","library_free")
filename<-c("RF_crcd_sublib_diann_prot_90%.csv","RS_crcd_sublib_diann_prot_90%.csv","IF_crcd_sublib_diann_prot_90%.csv","IS_crcd_sublib_diann_prot_90%.csv","DPHL1_crcd_sublib_diann_prot_90%.csv","library_free")
for(i in filename ){
  protmat<-read.csv(i,header = T)
protmat[is.na(protmat)]<-min(protmat[,2:ncol(protmat)],na.rm = T)
protmat[,2:ncol(protmat)]<-log2(protmat[,2:ncol(protmat)])
  # protmat<-data.frame(t(readin))
  # protmat<-readin[,apply(!is.na(readin),2,sum)>nrow(readin)*0.2]
  # for (j in 2:ncol(protmat)) {
  #   indsel<-which(genename$Entry==colnames(protmat)[j])
  #   if(length(indsel)>0){
  #     colnames(protmat)[j]<-paste(colnames(protmat)[j],genename$Gene[indsel],sep = "_")
  #   }
  # }
  
  protmat2<-merge(protmat,crcdinf,by.x = "prot",by.y = "batch_number",all.x = T)
  protmat2<-protmat2[,2:(ncol(protmat2)-7)]
  delmat<-protmat2[duplicated(protmat2$PatientId),]
  #选择蛋白平均值
  MeanNA<-data.frame(matrix(NA,nrow=length(unique(protmat2$PatientId)),ncol=(ncol(protmat2)-1)))
  colnames(MeanNA)<-colnames(protmat2)[1:(ncol(protmat2)-1)]
  row.names(MeanNA)<-unique(protmat2$PatientId)
  for (j in 1:nrow(MeanNA)) {
    indsel<-which(protmat2$PatientId==row.names(MeanNA)[j])
    if(length(indsel)>1){
      MeanNA[j,]<-log2(apply(2^protmat2[indsel,1:(ncol(protmat2)-1)],2,mean,na.rm=T)) 
      
    }
    else
      MeanNA[j,]<-protmat2[indsel,1:(ncol(protmat2)-1)]
  }
  MeanNA$patID<-row.names(MeanNA)
  crcdinf1<-crcdinf[!duplicated(crcdinf$PatientId),]
  unimat<-merge(MeanNA,crcdinf1,by.x = "patID",by.y = "PatientId",all.x = T)
  unimat$pat<-sapply(strsplit(unimat$patID,""),function(e){e[1]})
  protmat1<-unimat[,c(2:(ncol(unimat)-9),ncol(unimat))]
  row.names(protmat1)<-unimat$patID
  readCont<-protmat1[grepl("N",protmat1$pat),]
  readCont<-readCont[,-ncol(readCont)]
  readPat<-protmat1[grepl("P",protmat1$pat),]
  readPat<-readPat[,-ncol(readPat)]
  
  FC<-vector();
  pvalue<-vector();
  for(x in 1:ncol(readCont)){
    if(sum(!is.na(readCont[,x]))>1&sum(!is.na(readPat[,x]))>1){
      
      FC[x]<-mean(2^readPat[,x],na.rm=T)/mean(2^readCont[,x],na.rm = T)
      pvalue[x]<-t.test(readPat[,x],readCont[,x], paxred = F,  var.equal = F)$p.value
    }
  }
  pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
  FC_PAd<-cbind(FC,pvalueAd)
  row.names(FC_PAd)<-colnames(readCont)
  write.csv(FC_PAd,paste0("FC_PAd_",i))
  hist(pvalueAd)
  hist(log2(FC),breaks = 1000)
  dev.off()
  library(ggplot2)
  pdf(paste0(i,"_vocano.pdf"),width =8,height = 6)
  plot(log2(FC),-log10(pvalueAd),
       pch=19, col="gray",
       xlim=c(-8,12),
       main="patient vs control"
       ,xlab="log2 FC",ylab="-log10 pvalueAd",
       cex.axis = 1.5,cex.lab=1.5)
  up <- which(pvalueAd < 0.01 & FC > 2)
  down <- which(pvalueAd < 0.01 & FC <  0.5)
  points(log2(FC[up]), -log10(pvalueAd[up]), col=1, bg = "red", pch=21, cex=1)
  points(log2(FC[down]), -log10(pvalueAd[down]), col = "blue", pch=19,cex=1)
  abline(h=-log10(0.01),v=c(-log2(2),log2(2)),lty=2,lwd=1)
  dev.off()
  protUP<-FC_PAd[up,]
  protdown<-FC_PAd[down,]
  diff_prot<-rbind(protUP,protdown)
  write.csv(diff_prot,paste0("diff_prot_",i))
  # write.csv(protUP,paste0("protUP_",i))
  # write.csv(protdown,paste0("protDown_",i))
  upprot<-cbind(protmat1[,colnames(protmat1)%in%row.names(protUP)],protmat1$pat)
  downprot<-cbind(protmat1[,colnames(protmat1)%in%row.names(protdown)],protmat1$pat)
  write.csv(upprot,paste0("protUP_matrix_",i))
  write.csv(downprot,paste0("protDown_matrix_",i))
  }

#6个库dia search结果差异蛋白的overlap
#down
first_category_name<-c("reviewedfull","reviewedsemi","isoformfull","isoformsemi","dphl1","library_free")
filename<-paste0("protDown_",first_category_name,"_crcd_diann_prot_90.csv")
c<-list()
m<-1
for (i in filename) {
  indsel<-read.csv(i,header = T)
  c[[m]]<-indsel$X
  m<-m+1
}
uniP<-unique(c(unlist(c[[1]]),unlist(c[[2]]),unlist(c[[3]]),unlist(c[[4]]),unlist(c[[5]]),unlist(c[[6]])))
unimat<-data.frame(matrix(0,nrow = length(uniP),ncol = 6))
row.names(unimat)<-uniP
colnames(unimat)<-first_category_name
for (j in 1:ncol(unimat)) {
  for (i in 1:nrow(unimat)) {
    pr_this <- row.names(unimat)[i]
    if(pr_this%in%c[[j]]){
      unimat[i,j]<-1
    } 
}

}

###heatmap

#圆型heatmap
library(circlize)
library(grid)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(0,1), c( "#8F82BC","#F19E9E"))

pdf("20220103_DPHL2_diann_crcd_375Down_prot_circheatmap.pdf",width=8,height=8)
circos.heatmap(unimat,col =col_fun1,cluster =T,rownames.side = "outside",dend.side = "inside",bg.border = T)
circos.clear()
# lgd = Legend(title = "heatmap", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()
dev.off()


#upset图
library(UpSetR)

listinput <- list(reviewedFull=c[[1]],
                  reviwedSemi=c[[2]],
                  isoformFull=c[[3]],
                  isoformSemi=c[[4]],
                  DPHL1=c[[5]],
                library_free=c[[6]])
pdf('20220103_DPHL2_diann_crcd_down_prot.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()

###########################
#5个库的应用结果及差异蛋白的数量情况比较
temp<-read.csv("5_sublib_diff_summary_2.csv",header = T)

library(RColorBrewer)
library(randomcoloR)
library(reshape2)
library(ggplot2)
library(plyr)
protsum<-melt(temp)

a<-ggplot(protsum,aes(x=lib,y=value,fill=variable))+geom_bar(stat="identity",position = 'stack',colour="black")+
  geom_text(aes(label=value),size=4,position = position_stack())+guides(fill=guide_legend(reverse=F))+
  scale_fill_manual(values=c("burlywood1","gold2","gray74"))+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave("5lib_prot_up_down.pdf", plot=a,device = NULL, width = 12, height = 8)

#######################################################################
#missing value
library(ggplot2)
missing<-read.csv("5_sublib_diff_summary.csv",header = T)
missing$lib<-factor(missing$lib,levels=c("DPHL1","IS","IF","RS","RF"))
fig<-ggplot(missing,aes(x=lib,y=missing))+geom_bar(aes(y=missing,fill=lib),stat = "identity",width=0.6)+coord_flip()+
  ylim(0,0.5)+theme_minimal()+theme(legend.position = "none")
ggsave("dphl2_crcd__sublib_missing_value.pdf",fig,width = 8,height = 8)

RF<-read.csv("diff_prot_RF_crcd_sublib_diann_prot_90%.csv",header = T)
RS<-read.csv("diff_prot_RS_crcd_sublib_diann_prot_90%.csv",header = T)
IF<-read.csv("diff_prot_IF_crcd_sublib_diann_prot_90%.csv",header = T)
IS<-read.csv("diff_prot_IS_crcd_sublib_diann_prot_90%.csv",header = T)
DPHL1<-read.csv("diff_prot_DPHL1_crcd_sublib_diann_prot_90%.csv",header = T)
# Libfree<-read.csv("Libfree_CRCD_diff_prot.csv",header = T)

#upset图
library(UpSetR)

listinput <- list(reviewedFull=RF$X,
                  reviwedSemi=RS$X,
                  isoformFull=IF$X,
                  isoformSemi=IS$X,
                  DPHL1=DPHL1$X)
pdf('DPHL2_diann_crcd_sublib_diff_prot.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 5, order.by = "freq",text.scale =1.5)
dev.off()




#venn图
library(gplots)
library(VennDiagram)
venn1<-venn.diagram(listinput,resolution = 300, imagetype = "png",  cat.fontface=1,fontfamily=1,cat.cex=2,cex=2,
                    main="Protein down overlap of 6 librarys",main.cex = 2, main.fontface = 2, main.fontfamily = 3,filename = NULL)

# pdf(file="20220103_dphl_protein_down_overlap_lib.pdf")
grid.draw(venn1)#代码保存数字丢失，故需要手动保存
dev.off()



unimat$type<-paste(unimat$reviewedfull,unimat$reviewedsemi,unimat$isoformfull,unimat$isoformsemi,unimat$dphl1,sep = "_")
unimat1<-unimat[unimat$type=="1_1_1_1_1_1",]
write.csv(unimat1,"20220103_DPHL2_6lib_overlap_prot_down.csv",row.names = T)



#UP
filename<-paste0("protUP_",first_category_name,"_crcd_diann_prot_90.csv")
c<-list()
m<-1
for (i in filename) {
  indsel<-read.csv(i,header = T)
  c[[m]]<-indsel$X
  m<-m+1
}
uniP<-unique(c(unlist(c[[1]]),unlist(c[[2]]),unlist(c[[3]]),unlist(c[[4]]),unlist(c[[5]]),unlist(c[[6]])))
unimat<-data.frame(matrix(0,nrow = length(uniP),ncol = 6))
row.names(unimat)<-uniP
colnames(unimat)<-first_category_name
for (j in 1:ncol(unimat)) {
  for (i in 1:nrow(unimat)) {
    pr_this <- row.names(unimat)[i]
    if(pr_this%in%c[[j]]){
      unimat[i,j]<-1
    } 
  }
  
}

###heatmap

#圆型heatmap
library(circlize)
library(grid)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(0,1), c(  "#8F82BC","#F19E9E"))

pdf("20220103_DPHL2_diann_crcd_763UP_prot_circheatmap.pdf",width=8,height=8)
circos.heatmap(unimat,col =col_fun1,cluster =T,rownames.side = "outside",dend.side = "inside",bg.border = T)
circos.clear()
# lgd = Legend(title = "heatmap", col_fun = col_fun1)
# grid.draw(lgd)
# circos.clear()
dev.off()



#upset图
library(UpSetR)

listinput <- list(reviewedFull=c[[1]],
                  reviwedSemi=c[[2]],
                  isoformFull=c[[3]],
                  isoformSemi=c[[4]],
                  DPHL1=c[[5]],
                  library_free=c[[6]])
pdf('20220103_DPHL2_diann_crcd_up_prot.pdf',height = 8,width = 10)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()






#venn图
library(gplots)
library(VennDiagram)
venn1<-venn.diagram(listinput,resolution = 300, imagetype = "png",  cat.fontface=1,fontfamily=1,cat.cex=2,cex=2,
                    main="Protein up overlap of 6 librarys",main.cex = 2, main.fontface = 2, main.fontfamily = 3,filename =NULL)

grid.draw(venn1)#代码保存数字丢失，故需要手动保存
dev.off()

unimat$type<-paste(unimat$reviewedfull,unimat$reviewedsemi,unimat$isoformfull,unimat$isoformsemi,unimat$dphl1,sep = "_")
unimat1<-unimat[unimat$type=="1_1_1_1_1_1",]
write.csv(unimat1,"20220103_DPHL2_6lib_overlap_prot_up.csv",row.names = T)



sum<-read.csv("20220103_6lib_analysis_result_summary.csv",header = T)

library(ggplot2)
library(reshape2)
library(dplyr)
data_m <- melt(sum[,1:3], id.vars=c("X"))
p <- ggplot(data_m, aes(x=X, y=value))+ geom_bar(stat="identity", position="dodge", aes(fill=variable))+geom_text(aes(label=value))
p

ggsave("20220103_DPHL2_6lib_result_summary.pdf", plot=p,device = NULL, width = 12, height = 8)



# ###up中的蛋白太多了，各取前300，做overlap
# filename<-paste0("protUP_ 20210617_DPHL2_",first_category_name,"_DIANN_crcd_24raw_result.pr_matrix.tsv_prot.csv")
# c<-list()
# m<-1
# for (i in filename) {
#   indsel<-read.csv(i,header = T)
#   indsel1<-indsel[order(indsel$pvalueAd,decreasing = F),]
#   c[[m]]<-indsel1$X[1:300]
#   m<-m+1
# }
# uniP<-unique(c(unlist(c[[1]]),unlist(c[[2]]),unlist(c[[3]]),unlist(c[[4]])))
# unimat<-data.frame(matrix(0,nrow = length(uniP),ncol = 4))
# row.names(unimat)<-uniP
# colnames(unimat)<-first_category_name
# for (j in 1:ncol(unimat)) {
#   for (i in 1:nrow(unimat)) {
#     pr_this <- row.names(unimat)[i]
#     if(pr_this%in%c[[j]]){
#       unimat[i,j]<-1
#     } 
#   }
#   
# }
# 
# unimat$type<-paste(unimat$reviewedfull,unimat$reviewedsemi,unimat$isoformfull,unimat$isoformsemi,sep = "_")
# unimat1<-unimat[unimat$type=="1_1_1_1",]
# write.csv(unimat1,"20210923_DPHL2_4lib_top300_prot_overlap_up.csv",row.names = T)






















######################################################old################################################################
refucrcd<-read.csv("20210510_DIANN_DPHL2_prot_crcd_reviewedFull_uniprot.csv",stringsAsFactors = F,header = T)
crcdinf<-read.csv("crcd_sample_information.csv",stringsAsFactors = F,header = T)
genename<-read.csv("uniprot-human-gene_name.csv",stringsAsFactors = F,header = T)
protmat<-refucrcd[,apply(!is.na(refucrcd),2,sum)>nrow(refucrcd)*0.2]
for (i in 1:ncol(protmat)) {
  indsel<-which(genename$Entry==colnames(protmat)[i])
  if(length(indsel)>0){
    colnames(protmat)[i]<-paste(colnames(protmat)[i],genename$Gene[indsel],sep = "_")
  }
}
# protmat1<-log2(protmat[,2:4671])
# protmat[is.na(protmat)]<-1
#先把定量的值log2一下,NA填充为log2后最小值的0.5倍
protmat1<-log2(protmat[,2:5325])
protmat1[is.na(protmat1)]<-min(protmat1,na.rm=T)*0.1
protmat1$prot<-protmat$prot
protmat2<-merge(protmat1,crcdinf,by.x = "prot",by.y = "Sample.name",all.x = T)
row.names(protmat2)<-protmat2$prot
readCont<-protmat2[grepl("Benign prostatic hyperplasia tissue",protmat2$Tissue.type),]
readCont<-readCont[,2:5325]
readPat<-protmat2[grepl("Prostate cancer",protmat2$Tissue.type),]
readPat<-readPat[,2:5325]
FC<-vector();
pvalue<-vector();
for(i in 1:ncol(readCont)){
  if(sum(!is.na(readCont[,i]))>1&sum(!is.na(readPat[,i]))>1){
    FC[i]<-mean(readPat[,i])/mean(readCont[,i])
    pvalue[i]<-t.test(readPat[,i],readCont[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}
pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(readCont)
write.csv(FC_PAd,"20210512_DPHL2_crcd_dainn_reviewfull_FC_PAd.csv")
hist(pvalueAd)
hist(log2(FC),breaks = 1000)


library(ggplot2)
pdf("20210512_DPHL2_crcd_diann_reviewfull_volcano_FC_pvalueAd.pdf",width =8,height = 6)
plot(log2(FC),-log10(pvalueAd),
     pch=19, col="gray",
     xlim=c(-10,10),
     main=paste("patient vs control")
     ,xlab="log2 FC",ylab="-log10 pvalueAd",
     cex.axis = 1.5,cex.lab=1.5)
up <- which(pvalueAd < 0.05 & log2(FC) > 2)
down <- which(pvalueAd < 0.05 & log2(FC) <  -2)
points(log2(FC[up]), -log10(pvalueAd[up]), col=1, bg = "red", pch=21, cex=1)
points(log2(FC[down]), -log10(pvalueAd[down]), col = "blue", pch=19,cex=1)
abline(h=-log10(0.05),v=c(-log2(4),log2(4)),lty=2,lwd=1)
dev.off()
protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]
write.csv(protUP,"20210512_DPHL2_crcd_diann_reviewfull_protUP.csv")
write.csv(protdown,"20210512_DPHL2_crcd_diann_reviewfull_protdown.csv")

# pdf("0622_BoxplotUp_pvalueAd.pdf",width =8,height = 6)
# for(i in 1:length(up)){
#   
#   sel1<-as.numeric(readCont[,up[i]])
#   sel2<-as.numeric(readPat[,up[i]])
#   boxplot(list(sel1,sel2),na.rm=T,main=colnames(readCont)[up[i]],names=c("control","patient"),
#           pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.8),cex.axis=1.5)
#   stripchart(list(sel1,sel2),vertical = T,add=T,pch=19,jitter = 0.2,method = "jitter",cex=1
#              ,col=c("#A5C0DF","#EA8C8C"))
#   legend("bottom",bty="n",legend = c(paste("log2FC: ",round(log2(FC[up[i]]),6)),paste("pvalueAd:",round(pvalueAd[up[i]],6))),cex=1.2)
#   
# }
# dev.off()
# 
# pdf("0622_BoxplotDown_pvalueAd.pdf",width =8,height = 6)
# for(i in 1:length(down)){
#   sel1<-as.numeric(readCont[,down[i]])
#   sel2<-as.numeric(readPat[,down[i]])
#   boxplot(list(sel1,sel2),na.rm=T,main=colnames(readCont)[down[i]],names=c("control","patient"),
#           pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.8),cex.axis=1.5)
#   stripchart(list(sel1,sel2),vertical = T,add=T,pch=19,jitter = 0.2,method = "jitter",cex=1
#              ,col=c("#A5C0DF","#EA8C8C"))
#   legend("bottom",bty="n",legend = c(paste("log2FC: ",round(log2(FC[up[i]]),6)),paste("pvalueAd:",round(pvalueAd[up[i]],6))),cex=1.2)
#   
# }
# dev.off()

library(Biostrings)
fasta<-readBStringSet("2021-05-06-nodecoys-UP000005640_reviewed_human_20200401_no_irt_fragpipe.fasta",
                      format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

a<-data.frame(fasta@ranges@NAMES)

b<-data.frame(sapply(a$fasta.ranges.NAMES,function(e){unlist(strsplit(e," "))[1]}))#fasta中的|符号代码识别不了，只能导出在excel里面改

write.csv(b,"reviewed_fasta_prot_list.csv",row.names = F)


fasta<-readBStringSet("2021-05-06-nodecoys-UP000005640_isoform_reviewed_human_20200401_no_irt_fragpipe.fasta",
                      format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

a<-data.frame(fasta@ranges@NAMES)

b<-data.frame(sapply(a$fasta.ranges.NAMES,function(e){unlist(strsplit(e," "))[1]}))#fasta中的|符号代码识别不了，只能导出在excel里面改

write.csv(b,"isoform_fasta_prot_list.csv",row.names = F)


#isofull
#up
revfst<-read.csv("reviewed_fasta_prot_list.csv",header = T)
isofup<-read.csv("protUP_isoformfull_crcd_diann_prot_90.csv",header = T)
isofup$X<-gsub("\\.","-",isofup$X)
upprot<-setdiff(isofup$X,revfst$prot)
write.csv(upprot,"isofull_up_isoform.csv",row.names = F)

#down
isofdown<-read.csv("protDown_isoformfull_crcd_diann_prot_90.csv",header = T)
isofdown$X<-gsub("\\.","-",isofdown$X)
downprot<-setdiff(isofdown$X,revfst$prot)
write.csv(downprot,"isofull_down_isoform.csv",row.names = F)


#isosemi
#up
revfst<-read.csv("reviewed_fasta_prot_list.csv",header = T)
isofup<-read.csv("protUP_isoformsemi_crcd_diann_prot_90.csv",header = T)
isofup$X<-gsub("\\.","-",isofup$X)
upprot<-setdiff(isofup$X,revfst$prot)
write.csv(upprot,"isosemi_up_isoform.csv",row.names = F)

#down
isofdown<-read.csv("protDown_isoformsemi_crcd_diann_prot_90.csv",header = T)
isofdown$X<-gsub("\\.","-",isofdown$X)
downprot<-setdiff(isofdown$X,revfst$prot)
write.csv(downprot,"isosemi_down_isoform.csv",row.names = F)


#semi peptide protein 
#reviewedsemi
#up
revsemi<-read.csv("crcd_diann_reviewedsemi_pep_prot.csv",header = T)
revpepfull<-read.csv("20210526_DPHL2_reviewed_fasta_pep_nodecoy_noirt_miss1&no_miss.csv",header = T)
revsup<-read.csv("protUP_reviewedsemi_crcd_diann_prot_90.csv",header = T)
revsup$X<-gsub("\\.","-",revsup$X)
uppep<-revsemi[revsemi$prot%in%revsup$X,]
semipep<-uppep[uppep$pep%in%revpepfull$x,]
a<-unique(semipep$prot)
onlsemi<-uppep[!uppep$prot%in%a,]
b<-unique(onlsemi$prot)
write.csv(b,"reviewedsemi_up_only_semi_pep_prot.csv",row.names = F)

#down
revsdown<-read.csv("protDown_reviewedsemi_crcd_diann_prot_90.csv",header = T)
revsdown$X<-gsub("\\.","-",revsdown$X)
downpep<-revsemi[revsemi$prot%in%revsdown$X,]
semipep<-downpep[downpep$pep%in%revpepfull$x,]
a<-unique(semipep$prot)
onlsemi<-downpep[!downpep$prot%in%a,]
b<-unique(onlsemi$prot)
write.csv(b,"reviewedsemi_down_only_semi_pep_prot.csv",row.names = F)




#isoformsemi
#up
revsemi<-read.csv("crcd_diann_isoformsemi_pep_prot.csv",header = T)
revpepfull<-read.csv("20210526_DPHL2_isoform_fasta_pep_nodecoy_noirt_miss1&no_miss.csv",header = T)
revsup<-read.csv("protUP_isoformsemi_crcd_diann_prot_90.csv",header = T)
revsup$X<-gsub("\\.","-",revsup$X)
uppep<-revsemi[revsemi$prot%in%revsup$X,]
semipep<-uppep[uppep$pep%in%revpepfull$x,]
a<-unique(semipep$prot)
onlsemi<-uppep[!uppep$prot%in%a,]
b<-unique(onlsemi$prot)
write.csv(b,"isoformsemi_up_only_semi_pep_prot.csv",row.names = F)

#down
revsdown<-read.csv("protDown_isoformsemi_crcd_diann_prot_90.csv",header = T)
revsdown$X<-gsub("\\.","-",revsdown$X)
downpep<-revsemi[revsemi$prot%in%revsdown$X,]
semipep<-downpep[downpep$pep%in%revpepfull$x,]
a<-unique(semipep$prot)
onlsemi<-downpep[!downpep$prot%in%a,]
b<-unique(onlsemi$prot)
write.csv(b,"isoformsemi_down_only_semi_pep_prot.csv",row.names = F)

##################################################################################################
#已将差异蛋白转成resrsq格式，与下载的关于结直肠癌fasta中的蛋白匹配，找出与结直肠癌密切相关的isoform差异蛋白
fast<-read.csv("crcd_bioprot_list_isoform.csv",header = T)
isofup<-read.csv("crcd_diann_isofull_up_isoform_list.csv",header = T)
uselist<-isofup[isofup$To%in%fast$res_name,]
write.csv(uselist,"crcd_diann_isofull_up_isoformlist_in_Colorectal_cancer.csv",row.names = F)


isofup<-read.csv("crcd_diann_isosemi_down_isoform_list.csv",header = T)
uselist<-isofup[isofup$To%in%fast$res_name,]
write.csv(uselist,"crcd_diann_isosemi_down_isoformlist_in_Colorectal_cancer.csv",row.names = F)


isofup<-read.csv("crcd_diann_isosemi_up_isoform_list.csv",header = T)
uselist<-isofup[isofup$To%in%fast$res_name,]
write.csv(uselist,"crcd_diann_isosemi_up_isoformlist_in_Colorectal_cancer.csv",row.names = F)


isofup<-read.csv("crcd_diann_isofull_down_isoform_list.csv",header = T)
uselist<-isofup[isofup$To%in%fast$res_name,]
write.csv(uselist,"crcd_diann_isofull_down_isoformlist_in_Colorectal_cancer.csv",row.names = F)



