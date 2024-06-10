# compare finding related groups of mitigators 
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(Hmisc)
require(npreg)
require(readxl)
require(dplyr)
require(data.table)
require(ComplexHeatmap)
require(reshape2)
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/")
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
load("heatmap_plot.RData")
# only subject with mitigator data
mat_mitigator_effect=na.omit(mat_mitigator_effect)
subjects=rownames(mat_mitigator_effect)
mitigator_type_mat=matrix(0,nrow=length(subjects),ncol=4)
colnames(mitigator_type_mat)=c("fiber_spiker","fiber_mitigator","fat_mitigator","protein_mitigator")
rownames(mitigator_type_mat)=subjects
# 
fsp_mask=(!is.na(mat_mitigator_effect[,"Rice+Fiber"]))&mat_mitigator_effect[,"Rice+Fiber"]>0
mitigator_type_mat[,"fiber_spiker"]=ifelse(fsp_mask,1,0)
bes_supress_ind_list=apply(mat_mitigator_effect,1,function(x) which.min(x))
bes_supress_ind_list[sapply(bes_supress_ind_list,length)==0]=NA
bes_supress_ind=unlist(bes_supress_ind_list)
indmatch=seq(3)
names(indmatch)=c("fat","fiber","protein")
for(miti in names(indmatch)){
    seleind=which(bes_supress_ind==indmatch[miti]&mitigator_type_mat[,"fiber_spiker"]!=1)
    mitigator_type_mat[seleind,paste0(miti,"_mitigator")]=1
}
# compare between clinical data (set group comparing with others combined)
collists=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","bmi_avg_all","systolic_avg_all","diastolic_avg_all","fructosamine_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","modified_DI_heyjun","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","hepatic_IR_heyjun","ogtt_t_120_avg_all","hscrp_avg_all","creatinine_avg_all","albumin_avg_all","C_peptide_metabolic_testing_mean","GLP1_metabolic_testing_mean","GIP_metabolic_testing_mean","Glucagon_pmol_L_metabolic_testing_mean")
stat_coll_clinic=c()
subjects_sele=rownames(mitigator_type_mat[mitigator_type_mat[,"fiber_spiker"]!=1,])
metatab=metadatashow[,c("record_id",collists)]
metatab=metatab[metatab[,"record_id"]%in%subjects_sele,]
for(mittype in colnames(mitigator_type_mat)[-1]){
    tmptab=c()
    for(feature in collists){
        subvec=rownames(mitigator_type_mat[mitigator_type_mat[,mittype]==1,])
        gp_mask=metatab[,"record_id"]%in%subvec
        xvec=metatab[gp_mask,feature]
        yvec=metatab[!gp_mask,feature]
        xvec=xvec[!is.na(xvec)]
        yvec=yvec[!is.na(yvec)]
        if(length(xvec)<=3){
            next
        }
        tmptab=rbind(tmptab,data.frame(foods=mittype,feature=feature,pval=t.test(xvec,yvec)$p.value,delta=mean(xvec,na.rm=TRUE)-mean(yvec,na.rm=TRUE)))
    }
    tmptab$padj=p.adjust(tmptab$pval,method="fdr")
    stat_coll_clinic=rbind(stat_coll_clinic,tmptab)
}
write.table(stat_coll_clinic,paste0("mitig_type_clinic.tsv"),row.names=FALSE)
# compare between omics data: metaboloimcs, lipidomics, proteomics
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
load(paste0(pardir,"result/pathway/omics_data_clean.RData"))
nondupmask=!duplicated(featvec)
omicsmat=omicsmat[subjects_sele,nondupmask]
omicsvec=omicsvec[nondupmask]
featvec=featvec[nondupmask]
colnames(omicsmat)=featvec
# 
groupvec=colnames(mitigator_type_mat)[-1]
testtab=merge(omicsmat,mitigator_type_mat[subjects_sele,groupvec],by="row.names")
rownames(testtab)=testtab[,1]
testtab=testtab[,-1]
# 
nomics=dim(omicsmat)[2]
nmitigsig=length(groupvec)
# 
delta_mat=matrix(NA,nrow=nomics,ncol=nmitigsig)
colnames(delta_mat)=groupvec
p_mat=delta_mat
list_sig=vector(mode="list",length=nmitigsig*2)
names(list_sig)=paste0(rep(c("spikh","spikl"),each=nmitigsig),rep(groupvec,times=2))
stattab=c()
for(foodsigi in seq(nmitigsig)){
    pvec=c()
    detltavec=c()
    for(omicsi in seq(nomics)){
        locdf=testtab[,c(featvec[omicsi],groupvec[foodsigi])]
        colnames(locdf)=c("value","group")
        nonaind=(!is.na(locdf[,1]))&(!is.na(locdf[,2]))
        locdf=locdf[nonaind,]
        if(length(unique(locdf[,"value"]))<=1 | min(table(locdf[,"group"]))<=3){
            pvec=c(pvec,NA)
            detltavec=c(detltavec,NA)
        }else{
            pval=t.test(value~group,data=locdf)$p.value
            delta=mean(locdf[locdf[,"group"]==1,"value"],na.rm=TRUE)-mean(locdf[locdf[,"group"]==0,"value"],na.rm=TRUE)
            pvec=c(pvec,pval)
            detltavec=c(detltavec,delta)
        }
    }
    delta_mat[,foodsigi]=detltavec
    p_mat[,foodsigi]=pvec
    padj=rep(NA,times=length(pvec))
    namsk=is.na(pvec)
    for(theomics in unique(omicsvec)){
        themask=omicsvec==theomics
        locmask=themask&(!namsk)
        padj[locmask]=p.adjust(pvec[locmask],method="fdr")
    }
    sigmask=which(pvec<0.05)
    len=length(sigmask)
    if(len>0){
        compdnames=paste0(omicsvec[sigmask],featvec[sigmask])
        dirc=ifelse(detltavec[sigmask]>0,"spikh","spikl")
        listnames=paste0(dirc,groupvec[foodsigi])
        for(listname_i in seq(length(listnames))){
            listname=listnames[listname_i]
            list_sig[[listname]]=c(list_sig[[listname]],compdnames[listname_i])
        }
    }
    stattab=rbind(stattab,data.frame(food=rep(groupvec[foodsigi],times=len),feature=featvec[sigmask],omics=omicsvec[sigmask],delta=detltavec[sigmask],pval=pvec[sigmask],padj=padj[sigmask]))
}
write.table(stattab,paste0("mitig_type_omics.tsv"),row.names=FALSE)