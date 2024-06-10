# partial correlation 
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(Hmisc)
require(ppcor)
require(readxl)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/pathway/")
setwd(resdir)
# omics data 
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
load("omics_data_clean.RData")
subjects=rownames(omicsmat)
nsubjects=length(subjects)
# metadata
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
subphenogrp=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","modified_DI_heyjun","ie_heyjun","hepatic_IR_heyjun","HOMA_IR_heyjun","ffa_avg_heyjun")
nclasses=length(subphenogrp)
metadatashow=as.data.frame(matrix(NA,nrow=nsubjects,ncol=length(subphenogrp)+1))
rownames(metadatashow)=subjects
colnames(metadatashow)=c("study_id",subphenogrp)
for(subj in subjects){
    match_id=paste0("STUDYID-",str_pad(subj,3,pad="0"))
    metadatashow[subj,]=metadata[metadata[,"study_id"]==match_id,c("study_id",subphenogrp)]
}
metadatashow=metadatashow[,-1]
metadatashow=metadatashow[complete.cases(metadatashow),]
# select significant features
nomics=dim(omicsmat)[2]
pcor_mat=matrix(NA,nrow=nomics,ncol=nclasses)
colnames(pcor_mat)=subphenogrp
rownames(pcor_mat)=featvec
p_mat=pcor_mat
tabstat_cor=c()
for(omicsi in seq(nomics)){
    xmat=metadatashow
    yvec=omicsmat[rownames(metadatashow),omicsi]
    nonaind=(!is.na(yvec))
    if(length(unique(yvec[nonaind]))<=1){
        next
    }
    datadf=cbind(xmat[nonaind,],yvec[nonaind])
    lastind=dim(datadf)[2]
    corres=pcor(datadf,method="spearman")
    pcor_mat[omicsi,]=corres$estimate[lastind,-lastind]
    p_mat[omicsi,]=corres$p.value[lastind,-lastind]
    tabstat_cor=rbind(tabstat_cor,data.frame(type=rep(omicsvec[omicsi],times=nclasses),clinic=subphenogrp,feature=rep(featvec[omicsi],times=nclasses),correlation=pcor_mat[omicsi,],pvalue=p_mat[omicsi,]))
}
omicslist=unique(omicsvec)
tabstat_cor$p.adj=rep(NA,times=dim(tabstat_cor)[1])
for(classi in seq(nclasses)){
    for(omicsele in omicslist){
        ind=which(tabstat_cor[,"type"]==omicsele&tabstat_cor[,"clinic"]==subphenogrp[classi])
        tabstat_cor[ind,"p.adj"]=p.adjust(tabstat_cor[ind,"pvalue"],"fdr")
    }
}
save(tabstat_cor,file="pcor_subtypes.RData")