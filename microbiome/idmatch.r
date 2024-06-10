# sample and filename matching
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
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/")
qcdir=paste0(pardir,"result/microbiome/fastqc/")
resdir=paste0(pardir,"result/microbiome/eda/")
setwd(resdir)
# 
file_tab_trim=read.table(paste0(qcdir,"trim_all_list.txt"),header=FALSE)
files=file_tab_trim[,9]
# id matching table
expr_tab=c()
temptab=as.data.frame(read_excel(paste0(pardir,"data/microbiome_code/CGM_Samples.xlsx"),sheet="Pool1"))
temptab=temptab[,c("subject_id","weight","reserve_well")]
colnames(temptab)=c("subject_id","weight","vial_id")
temptab$plate=rep("pool_1",times=dim(temptab)[1])
filesvec=list()
for(rowi in seq(nrow(temptab))){
    matind=str_which(string=files,pattern=fixed(paste0("-CGM",temptab[rowi,"vial_id"],"-")))
    if(length(matind)==0){
        filesvec=c(filesvec,NA)
    }else{
        filesvec=c(filesvec,list(files[matind]))
    }
}
temptab$files=filesvec
expr_tab=temptab
temptab=as.data.frame(read_excel(paste0(pardir,"data/microbiome_code/CGM_Samples.xlsx"),sheet="CGM_E23"))
temptab=temptab[,c("subject_id","weight","file_id","plate...9","filename")]
colnames(temptab)=c("subject_id","weight","vial_id","plate","match_str")
filesvec=list()
for(rowi in seq(nrow(temptab))){
    # mask1=str_detect(string=files,pattern=fixed(paste0("-",temptab[rowi,"vial_id"],"-")))
    # mask2=str_detect(string=files,pattern=fixed(paste0("-",temptab[rowi,"plate"],"-")))
    # matind=which(mask1&mask2)
    mask=str_detect(string=files,pattern=fixed(temptab[rowi,"match_str"]))
    matind=which(mask)
    if(length(matind)==0){
        filesvec=c(filesvec,NA)
    }else{
        filesvec=c(filesvec,list(files[matind]))
    }
}
temptab$files=filesvec
expr_tab=rbind(expr_tab,temptab[,c("subject_id","weight","vial_id","plate","files")])
expr_tab[,"subject_id"]=str_remove(string=expr_tab[,"subject_id"],pattern="\\.0$")
save(expr_tab,file="sample_matching.RData")