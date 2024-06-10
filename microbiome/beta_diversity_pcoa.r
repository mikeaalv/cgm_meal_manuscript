# beta diversity scatter plot
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(Hmisc)
require(readxl)
require(dplyr)
require(data.table)
require(ggfortify)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/")
predatadir=paste0(pardir,"result/microbiome/eda/")
resdir=paste0(pardir,"result/microbiome/diversity/")
setwd(resdir)
# 
load(paste0(predatadir,"microbiome_processed.RData"))
div_tab=read.table("diversity_table.tsv",header=TRUE)
colnames(div_tab)[1]="sample"
div_tab$sample=str_remove(string=div_tab$sample,pattern="\\_trimmed\\.fq\\.gz\\_metaphlan\\_bugs\\_list$")
# meta data
metadata_clinic=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
grouplist=c("sex.factor","ethnicity","a1c_t2d_status_bl","sspg_status_heyjun","di_3classes_heyjun","ie_3_classes_heyjun","FFA_3classes_heyjun","hepatic_ir_3classes_heyjun")
metadata_clinic=metadata_clinic[,c("study_id",grouplist)]
metadata_clinic$transformedid<-metadata_clinic$study_id%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()
for(groupfeat in grouplist){
    metadata_clinic[which(metadata_clinic[,groupfeat]=="Unknown"),groupfeat]=NA
}
metacomb=merge(sampfilematchmeta[,c("subject_id","files")],metadata_clinic,by.x="subject_id",by.y="transformedid")
metacomb$sample=str_remove(string=metacomb$files,pattern="\\_trimmed\\.fq\\.gz$")
types=c("bc","jaccard")
selecpc=c("PC1","PC2")
# 
for(type in types){
    seleccol=paste0(type,selecpc)
    pc_tab=div_tab[,c("sample",seleccol)]
    plotdf=merge(metacomb[,c("sample",grouplist)],pc_tab,by="sample")
    colnames(plotdf)=str_remove(string=colnames(plotdf),pattern=paste0("^",type))
    for(grp in grouplist){
        p<-ggplot(plotdf,aes_string(x="PC1",y="PC2",color=grp))+geom_point()
        ggsave(paste0(type,"_",grp,"_scatter.pdf"),plot=p) 
    }
}