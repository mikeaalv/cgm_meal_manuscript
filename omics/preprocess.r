# reformatting of omics features
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(readxl)
require(dplyr)
require(data.table)
require(reshape2)
require(masscleaner)
require(ggfortify)
require(Rcpm)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/omics/")
data_pre=paste0(pardir,"data/")
setwd(resdir)
omicsset=c("metabolomics","lipidomics")
# batches of samples
batchtab_list=list()
batchtab_list[["lipidomics"]]=as.data.frame(read_excel(paste0(data_pre,"metabolomics_lipidomics/CGM batches.xlsx"),sheet="lipidomics"))
batchtab_list[["metabolomics"]]=as.data.frame(read_excel(paste0(data_pre,"metabolomics_lipidomics/CGM batches.xlsx"),sheet="metabolomics"))
samp_match_tab=as.data.frame(read_excel(paste0(data_pre,"metabolomics_lipidomics/Sample list.xlsx"),sheet="Randomized sample list"))
date_form<-samp_match_tab[,"Visit Date"] %>% format(.,"%m.%d.%Y") %>% str_replace_all(string=.,pattern="^0+",replacement="") %>% str_replace_all(string=.,pattern="\\.0+",replacement=".")
sample_id=paste0("X",samp_match_tab[,"SubjectID"],"_",date_form)
id_conv=data.frame(intern_id=samp_match_tab[,"sample ID"],sample_id=sample_id)
batchconv_list=list()
for(omics in omicsset){
    loctab=merge(batchtab_list[[omics]],id_conv,by.x="sample",by.y="intern_id")
    loctab=loctab[,c("sample_id","prep batch","Acquisition batch")]
    colnames(loctab)=c("sample_id","prep_batch","acquisition_batch")
    batchconv_list[[omics]]=loctab
}
# 
omicslist_arra=list()
omicslist_arra_unscale=list()
annolist=list()
############# metabolomics dataset #############
metatab=read.table(paste0(data_pre,"metabolomics_lipidomics/CGM_plasma_metabolomics_all_Anno_metID(1).csv"),sep=",",header=TRUE)
colnameall=colnames(metatab)
col_sele=colnameall[str_detect(string=colnameall,pattern="X\\d+\\_\\d+\\.\\d+\\.\\d+")]
datamat=as.matrix(metatab[,col_sele])
# PQN normalization
datamat=t(pqn(t(datamat)))
# keep level 1 and 2
level_mask=metatab[,"Level"]==1 | metatab[,"Level"]==2
name_mask=!is.na(metatab[,"Compound.name"])
metatab=metatab[level_mask&name_mask,]
datamat=datamat[level_mask&name_mask,]
# if multiple peaks match to the same compound, select the closest ones by level and merge 
countab=table(metatab[,"Compound.name"])
dupmask=countab>1
# the single match ones
name_sin=names(countab[!dupmask])
sin_mask=metatab[,"Compound.name"]%in%name_sin
datamat_nrp=datamat[sin_mask,]
annotations=metatab[sin_mask,"Compound.name"]
# merge adducts
for(ind_dup in which(dupmask)){
    compd=names(countab)[ind_dup]
    matchind=which(metatab[,"Compound.name"]==compd)
    levels=metatab[matchind,"Level"];
    bestlevel=min(levels);
    feat_bestlev_ind=matchind[which(levels==bestlevel)];
    if(length(feat_bestlev_ind)>1){
        tempmat=colSums(datamat[feat_bestlev_ind,])
    }else{
        tempmat=datamat[feat_bestlev_ind,]
    }
    datamat_nrp=rbind(datamat_nrp,tempmat)
    annotations=c(annotations,compd)
}
# merge mulitple batches by scaling median
sampid=colnames(datamat_nrp)
batchtab=batchconv_list[["metabolomics"]]
batchtab=batchtab[batchtab[,"sample_id"]%in%sampid,]
batchtab=batchtab[match(sampid,batchtab[,"sample_id"]),]
split_dataarray<-purrr::map(unique(batchtab$acquisition_batch),
    .f=function(thebatch){
        temp_sample_info<-batchtab%>%dplyr::filter(acquisition_batch==thebatch)
        datamat_nrp[,temp_sample_info$sample_id,drop=FALSE]
    })
#calculating the scale factor
correct_factor<-purrr::map(split_dataarray,function(x) {return(apply(x,1,median))})
correct_factor<-lapply(correct_factor,function(x) {correct_factor[[1]]/x})
correct_factor<-correct_factor%>%
    lapply(function(x){
        x[is.na(x)]=1
        x[is.infinite(x)]=1
        x
    })
split_dataarray<-purrr::map2(split_dataarray,correct_factor,function(x,y){x*y})
split_dataarray_merg=dplyr::bind_cols(split_dataarray)
# center and scale transformation for each feature
split_dataarray_merg2=t(scale(log(t(split_dataarray_merg))))
# PCA plot 
# split_dataarray_merg_t=t(split_dataarray_merg)
# pca_res=prcomp(split_dataarray_merg_t,center=FALSE)
# autoplot(pca_res)
omicslist_arra[["metabolomics"]]=split_dataarray_merg2
omicslist_arra_unscale[["metabolomics"]]=t(log(t(split_dataarray_merg)))
anno_tab=metatab[match(annotations,metatab[,"Compound.name"]),]
anno_tab=anno_tab[,!colnames(anno_tab)%in%col_sele]
annolist[["metabolomics"]]=anno_tab
############# lipidomics dataset #############
lipidtab=read.table(paste0(data_pre,"metabolomics_lipidomics/CGM_Plasma_Lipid_Species(1).csv"),sep=",",header=TRUE)
colnameall=colnames(lipidtab)
col_sele=colnameall[str_detect(string=colnameall,pattern="X\\d+\\_\\d+\\.\\d+\\.\\d+")]
datamat=as.matrix(lipidtab[,col_sele])
# PQN normalization
datamat=t(pqn(t(datamat)))
annotations=lipidtab[,"X"]
datamat_nrp=datamat
# merge mulitple batches by scaling median
sampid=colnames(datamat_nrp)
batchtab=batchconv_list[["lipidomics"]]
batchtab=batchtab[batchtab[,"sample_id"]%in%sampid,]
batchtab=batchtab[match(sampid,batchtab[,"sample_id"]),]
split_dataarray<-purrr::map(unique(batchtab$acquisition_batch),
    .f=function(thebatch){
        temp_sample_info<-batchtab%>%dplyr::filter(acquisition_batch==thebatch)
        datamat_nrp[,temp_sample_info$sample_id,drop=FALSE]
    })
#calculating the scale factor
correct_factor<-purrr::map(split_dataarray,function(x) {return(apply(x,1,median))})
correct_factor<-lapply(correct_factor,function(x) {correct_factor[[1]]/x})
correct_factor<-correct_factor%>%
    lapply(function(x){
        x[is.na(x)]=1
        x[is.infinite(x)]=1
        x
    })
split_dataarray<-purrr::map2(split_dataarray,correct_factor,function(x,y){x*y})
split_dataarray_merg=dplyr::bind_cols(split_dataarray)
# center and scale transformation for each feature
split_dataarray_merg2=t(scale(log(t(split_dataarray_merg))))
# PCA plot 
# pca_res=prcomp(t(split_dataarray_merg),center=FALSE)
# autoplot(pca_res,data=batchtab,colour="prep_batch")
# autoplot(pca_res,data=batchtab,colour="acquisition_batch")
omicslist_arra[["lipidomics"]]=split_dataarray_merg2
omicslist_arra_unscale[["lipidomics"]]=t(log(t(split_dataarray_merg)))
annolist[["lipidomics"]]=data.frame("Compound.name"=annotations)
############# proteomics dataset #############
proteomicstab=read.table(paste0(data_pre,"proteomics/StanfordProject_FahimAbbasi_NPX_2022-08-24_CGM1study.csv"),sep=";",header=TRUE)
proteomicstab[,"subjectid"]<-proteomicstab[,"SampleID"]%>%str_extract(string=.,pattern="^\\d+\\-\\d+")%>%str_extract(string=.,pattern="-\\d+")%>%str_remove(string=.,pattern="-0+")%>%paste0("X",.)
subjects=unique(proteomicstab[,"subjectid"])
proteins=unique(proteomicstab[,"OlinkID"])
dataarray_prot=matrix(NA,nrow=length(proteins),ncol=length(subjects))
colnames(dataarray_prot)=subjects
anno_tab=c()
for(proti in seq(length(proteins))){
    prot=proteins[proti]
    loctab=proteomicstab[proteomicstab[,"OlinkID"]==prot,]
    loctab=loctab[match(subjects,loctab[,"subjectid"]),]
    dataarray_prot[proti,]=loctab[,"NPX"]
    anno_tab=rbind(anno_tab,loctab[1,c("OlinkID","UniProt","Assay","MissingFreq","Panel","Panel_Lot_Nr","PlateID")])
}
omicslist_arra[["proteomics"]]=t(scale(t(dataarray_prot)))
omicslist_arra_unscale[["proteomics"]]=dataarray_prot
annolist[["proteomics"]]=anno_tab
save(omicslist_arra,annolist,omicslist_arra_unscale,file="metabolomics_lipidomics_processed.RData")