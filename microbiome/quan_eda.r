# taxonomy quatnfication and plot
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
resdir=paste0(pardir,"result/microbiome/eda/")
matdir=paste0(pardir,"result/microbiome/matrix_res/")
convdir=paste0(pardir,"result/microbiome/feature_id_conv/")
setwd(resdir)
# 
load("sample_matching.RData")
tab_taxon=read.table(paste0(matdir,"species_full.tsv"),header=TRUE,sep="\t",check.names=FALSE)
allclades=tab_taxon[,"clade_name"]
ref_ind=str_which(string=allclades,pattern="halotolerans")
sele_ind=c(seq(length(allclades)),ref_ind)
taxon_g=allclades[sele_ind]
# 
mat_rel=as.matrix(tab_taxon[,seq(from=3,to=dim(tab_taxon)[2])])
mat_rel=mat_rel[sele_ind,]
rownames(mat_rel)=taxon_g
# remove samples with no measured reference speicies
ref_zero_mask=mat_rel[nrow(mat_rel),]==0
mat_rel=mat_rel[,!ref_zero_mask]
# meta data
expr_tab[,"weight"]=as.numeric(expr_tab[,"weight"])
expr_tab=expr_tab[!is.na(expr_tab[,"weight"]),]
files=colnames(mat_rel)
expr_tab_long=as.data.frame(tidyr::unnest(expr_tab,cols=files))
expr_tab_long[,"files"]=str_remove(string=expr_tab_long[,"files"],pattern=fixed("_clean.fastq.gz"))
# refilter both meta data and data table for the intersection part
files<-files%>%str_remove(string=.,pattern="\\_metaphlan\\_bugs\\_list$")
matchind=match(files,expr_tab_long[,"files"])
nonnamaks=!is.na(matchind)
files=files[nonnamaks]
mat_rel=mat_rel[,nonnamaks]
metadata=expr_tab_long[matchind[nonnamaks],]
# 
sampvec=rep("cgm",times=dim(metadata)[1])
matchlist=list("ipop"="-","blank"="blank","mock"="mock","generous_donor"="generous_donor")
for(samptype in names(matchlist)){
    ind=str_which(string=metadata[,"subject_id"],pattern=fixed(matchlist[[samptype]]))
    sampvec[ind]=samptype
}
metadata$samptypes=sampvec
# divided by the reference species
mat_rel_ref=t(t(mat_rel)/mat_rel[nrow(mat_rel),])
#remove reference
mat_rel_ref=mat_rel_ref[-nrow(mat_rel),]
taxon_g=taxon_g[-length(taxon_g)]
# normalize by weight
mat_rel_ref_norm=t(t(mat_rel_ref)/metadata[,"weight"])
# select the groups
gmask=str_detect(string=allclades,pattern="\\|g\\_\\_")
smask=str_detect(string=allclades,pattern="\\|s\\_\\_")
amask=str_detect(string=allclades,pattern="^k\\_\\_")
list_ind=list(g=which(gmask&(!smask)),s=which(smask),all=which(amask))
mat_list=list()
for(grptype in names(list_ind)){
    mat_rel_ref_norm_here=mat_rel_ref_norm[list_ind[[grptype]],]
    # pca absolute quantification
    pca_res_ref=prcomp(t(mat_rel_ref_norm_here))
    p<-autoplot(pca_res_ref,data=metadata,colour='plate')
    ggsave(paste0("pca_abs_quant",grptype,".pdf"),plot=p)
    pdf(paste0("hist_abs_quant",grptype,".pdf"))
    hist(c(mat_rel_ref_norm_here))
    dev.off()
    lowthres=sort(unique(c(mat_rel_ref_norm_here)),decreasing=FALSE)[2]
    mat_rel_ref_norm_log=log(mat_rel_ref_norm_here+lowthres/2)
    pca_res_ref=prcomp(t(mat_rel_ref_norm_log))
    p<-autoplot(pca_res_ref,data=metadata,colour='plate')
    ggsave(paste0("pca_abs_quant_log",grptype,".pdf"),plot=p)
    pdf(paste0("hist_abs_quant_log",grptype,".pdf"))
    hist(c(mat_rel_ref_norm_log))
    dev.off()
    pdf(paste0("hist_abs_quant_log_no0",grptype,".pdf"))
    hist(c(log(mat_rel_ref_norm_here)))
    dev.off()
    colnames(mat_rel_ref_norm_log)=files
    mat_list[[paste0("species_abs_",grptype)]]=mat_rel_ref_norm_log
}
# 
list_sources=list("go"="go_full.tsv","kegg"="kegg_full.tsv","pathway"="pathway_full_cpm.tsv","pfam"="pfam_full.tsv","rxn"="rxn_full.tsv")
filesall=files
missingthres=0.9
summfeat_list=list()
missing_ratio_list=list()
avg_perc_list=list()
# clean each source of microbiome data
for(type in names(list_sources)){
    tab_temp=read.table(paste0(matdir,list_sources[[type]]),header=TRUE,sep="\t",check.names=FALSE,comment.char="",quote="")
    # filter samples
    cols=colnames(tab_temp)
    allsamp=str_replace_all(string=cols[-1],pattern="\\.gz.*$",replacement=".gz")
    if(length(setdiff(filesall,allsamp))>0){
        stop("missing samples")
    }
    tab_temp=tab_temp[,c(1,match(filesall,allsamp)+1)]
    # filter features 
    allfeats=tab_temp[,1]
    smask=str_detect(string=allfeats,pattern="\\.s\\_\\_")
    sumfeatvec=allfeats[!smask]
    sumfeatvec=sumfeatvec[!str_detect(string=sumfeatvec,pattern="(UNMAPPED)|(UNINTEGRATED)|(UNGROUPED)|(unclassified)")]
    splitfeatvec=allfeats[smask]
    collfeat=c()
    tab_mat_split=tab_temp[match(splitfeatvec,allfeats),-1]
    ratiovec=rowSums(tab_mat_split==0)/dim(tab_mat_split)[2]
    collfeat=splitfeatvec[ratiovec<missingthres]
    selefeat=c(sumfeatvec,collfeat)
    mat_rel=as.matrix(tab_temp[match(selefeat,allfeats),-1])
    rownames(mat_rel)=selefeat
    colnames(mat_rel)=filesall
    # 
    missing_ratio_list[[type]]=rowSums(mat_rel==0)/dim(mat_rel)[2]
    avg_perc_list[[type]]=rowMeans(apply(mat_rel,2,percent_rank))
    # 
    minvec=sort(unique(c(mat_rel)),decreasing=FALSE)
    if(minvec[1]==0){
        lowthres=minvec[2]
    }else{
        lowthres=minvec[1]
    }
    mat_rel_log=log(mat_rel+lowthres/2)
    # pca
    pca_res_ref=prcomp(t(mat_rel_log[sumfeatvec,]))
    p<-autoplot(pca_res_ref,data=metadata,colour='plate')
    ggsave(paste0("pca_abs_quant_log_",type,".pdf"),plot=p)
    # histogram
    pdf(paste0("hist_abs_quant_log_",type,".pdf"))
    hist(c(mat_rel_log[sumfeatvec,]))
    dev.off()
    mat_list[[type]]=mat_rel_log
    summfeat_list[[type]]=sumfeatvec
}
# pca plot with subtypes
metadata_clinic=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
grouplist=c("sex.factor","ethnicity","a1c_t2d_status_bl","sspg_status_heyjun","di_3classes_heyjun","ie_3_classes_heyjun","FFA_3classes_heyjun","hepatic_ir_3classes_heyjun")
metadata_clinic=metadata_clinic[,c("study_id",grouplist)]
metadata_clinic$transformedid<-metadata_clinic$study_id%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()
for(groupfeat in grouplist){
    metadata_clinic[which(metadata_clinic[,groupfeat]=="Unknown"),groupfeat]=NA
}
# 
summfeat_list[["species_abs_g"]]=rownames(mat_list[["species_abs_g"]])
summfeat_list[["species_abs_s"]]=rownames(mat_list[["species_abs_s"]])
summfeat_list[["species_abs_all"]]=rownames(mat_list[["species_abs_all"]])
# 
types=names(mat_list)
featvec_micro=c()
for(type in types){
    mat=as.data.frame(t(mat_list[[type]]))
    summfeat=summfeat_list[[type]]
    pca_res=prcomp(mat[,summfeat])
    # metadata reformat
    samples=rownames(mat)
    subjectids=metadata[match(samples,metadata[,"files"]),"subject_id"]
    meta_tranform_df=as.data.frame(metadata_clinic[match(subjectids,metadata_clinic[,"transformedid"]),])
    # PCA plot with metabolic subtypes
    for(grp in grouplist){
        p<-autoplot(pca_res,data=meta_tranform_df,colour=grp)
        ggsave(paste0(type,"_",grp,".pdf"),plot=p) 
    }
    featvec_micro=c(featvec_micro,rep(type,times=dim(mat)[2]))
}
subjects=unique(subjectids)
subjects=subjects[!str_detect(string=subjects,pattern="\\-")]
ncols=sum(sapply(mat_list,function(x) dim(x)[1]))
micromat=matrix(NA,nrow=length(subjects),ncol=ncols)
rownames(micromat)=subjects
colnames(micromat)=unlist(sapply(mat_list,rownames))
for(subj in subjects){
    tempmat=c()
    filelist=metadata[metadata[,"subject_id"]==subj,"files"]
    for(type in types){
        sourmat=mat_list[[type]]
        tempmat=c(tempmat,apply(sourmat[,filelist,drop=FALSE],1,mean))
    }
    micromat[subj,]=tempmat
}
weighttab=metadata[,c("subject_id","weight")]%>%group_by(subject_id)%>%summarize(meanw=mean(weight),.groups="keep")%>%as.data.frame()
# id conversion (name) table
needconvlist=c("go","kegg","pfam","rxn")
idname_match_list=list()
for(idtype in needconvlist){
    tab_temp=read.table(paste0(convdir,idtype,"_full_idconv.tsv"),header=TRUE,sep="\t",check.names=FALSE,comment.char="",quote="")
    namestr<-tab_temp[,1]%>%str_remove(string=.,pattern="\\|.*$")%>%unique(.)
    namestr=namestr[str_detect(string=namestr,pattern=": ")]
    tabsplit=str_split(string=namestr,pattern=": ",simplify=TRUE)
    idname_match_list[[idtype]]=tabsplit
}
sampfilematchmeta=metadata
save(summfeat_list,micromat,mat_list,featvec_micro,weighttab,idname_match_list,missing_ratio_list,avg_perc_list,sampfilematchmeta,file="microbiome_processed.RData")
