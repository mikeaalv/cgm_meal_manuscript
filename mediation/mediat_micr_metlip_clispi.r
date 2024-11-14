# mediation analysis
# microbiome -> metabolomics,lipidomics -> clincal,spikes
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(igraph)
require(Hmisc)
require(npreg)
require(readxl)
require(circlize)
require(ComplexHeatmap)
require(mediation)
# 
require(foreach)
require(doSNOW)
require(Rodin)
cl<-makeSOCKcluster(8)
registerDoSNOW(cl)
# 
set.seed(1)
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/mediation/")
setwd(resdir)
##  microbiome->metabolomics
load(paste0(pardir,"result/microbiome/asso_other_omics/association_metabolomics_microbiome.RData"))
subtab1=tabstat_cor[,c("feature","metabolites","padj")]
colnames(subtab1)=c("microbiome","metabolites","padj1")
## metabolomics->clinical
load(paste0(pardir,"/result/pathway/association_subtype_omics.RData"))
subtab2=tabstat_cor[,c("feature","clinic","padj")]
colnames(subtab2)=c("metabolites","outcome","padj2")
## metabolomics->spikes
load(paste0(pardir,"result/pathway/association_omics_food_all.RData"))
subtab2_2=stat_cor[,c("feature","food","padj")]
colnames(subtab2_2)=c("metabolites","outcome","padj2")
subtab2=rbind(subtab2,subtab2_2)
## microbiome->clinical
load(paste0(pardir,"result/microbiome/metabolic_subtypes/association_subtype_microbiome.RData"))
subtab3=tabstat_cor[,c("feature","clinic","padj")]
colnames(subtab3)=c("microbiome","outcome","padj3")
# microbiome->spikes
load(paste0(pardir,"result/microbiome/spikes/association_spike_microbiome.RData"))
subtab3_2=tabstat_cor[,c("feature","spikes","padj")]
colnames(subtab3_2)=c("microbiome","outcome","padj3")
subtab3=rbind(subtab3,subtab3_2)
# 
mergtab=merge(merge(subtab1,subtab2,by="metabolites"),subtab3,by=c("outcome","microbiome"))
# 
mergtab_sele=mergtab[mergtab[,"padj3"]<0.1|mergtab[,"padj1"]<0.05,]#mergtab[,"padj3"]<0.1|mergtab[,"padj1"]<0.05
# 
load(paste0(pardir,"result/microbiome/eda/microbiome_processed.RData"))
mergtab_sele=mergtab_sele[mergtab_sele[,"microbiome"]%in%colnames(micromat),]
covar_cols=c("sex","bmi_avg_all","age_today_avg_all")
mergtab_sele=mergtab_sele[!mergtab_sele[,"outcome"]%in%covar_cols,]
save(mergtab_sele,file="start_mediation_list.RData")
# load data matrix
## microbiome  
micromat_decr=apply(micromat,2,function(x){ifelse(x>quantile(x,0.2),1,0)})
subjects=rownames(micromat)
nsubjects=length(subjects)
# omics
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
load(paste0(pardir,"result/pathway/omics_data_clean.RData"))
metab_ind=which(omicsvec=="metabolomics"|omicsvec=="lipidomics")
metab_mat=omicsmat[,metab_ind]
featvec_metab=featvec[metab_ind]
colnames(metab_mat)=featvec_metab
# clinical data
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
metadata$sex=ifelse(metadata[,"sex.factor"]=="Male",1,0)
class_num_cols=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","modified_DI_heyjun","ie_heyjun","hepatic_IR_heyjun","HOMA_IR_heyjun","ffa_avg_heyjun","HOMA_S_heyjun","HOMA_B_heyjun","albumin_avg_all","creatinine_avg_all","hscrp_avg_all","ldl_avg_all","hdl_avg_all","triglyceride_avg_all","fructosamine_avg_all","insulin_fasting_avg_all","alt_sgpt_avg_all","systolic_avg_all","diastolic_avg_all")
allmetacols=c(class_num_cols,covar_cols)
metadataclinic=as.data.frame(matrix(NA,nrow=nsubjects,ncol=length(allmetacols)+1))
rownames(metadataclinic)=subjects
colnames(metadataclinic)=c("study_id",allmetacols)
for(subj in subjects){
    match_id=paste0("STUDYID-",str_pad(subj,3,pad="0"))
    metadataclinic[subj,]=metadata[metadata[,"study_id"]==match_id,c("study_id",allmetacols)]
}
metadataclinic=as.matrix(metadataclinic[,-1])
# spike data 
load(paste0(pardir,"result/cgm_meal/heatmap_plot.RData"))
food_sig_mat=cbind(mat_quantile,mat_mitigator_effect)
# 
subejects=Reduce(intersect,list("micro"=rownames(micromat_decr),"metali"=rownames(metab_mat),"clic"=rownames(metadataclinic),"spike"=rownames(food_sig_mat)))
inputmat=micromat_decr[subejects,]
mediationmat=metab_mat[subejects,]
outcomemat=cbind(metadataclinic[subejects,class_num_cols],food_sig_mat[subejects,])
covarmat=metadataclinic[subejects,covar_cols]
wvec=weighttab[match(subejects,weighttab[,"subject_id"]),"meanw"]
covarmat=cbind(covarmat,wvec)
colnames(covarmat)[dim(covarmat)[2]]="samp_weight"
# 
progress<-function(n) cat(sprintf("task %d is complete\n",n))
opts<-list(progress=progress)
tabstat_list<-foreach(rowi=seq(nrow(mergtab_sele)),.packages=c("mediation"),.options.snow=opts)%dopar%{
    set.seed(1)
    locrecrd=mergtab_sele[rowi,]
    locdf=as.data.frame(cbind(data.frame(outcome=outcomemat[,locrecrd[,"outcome"]],metabolite=mediationmat[,locrecrd[,"metabolites"]],microbiome=inputmat[,locrecrd[,"microbiome"]]),covarmat))
    locdf=na.omit(locdf)
    mediator_lm=lm(metabolite ~ microbiome + sex + bmi_avg_all + samp_weight,data=locdf)
    outcome_lm=lm(outcome ~ microbiome + metabolite + sex + bmi_avg_all + samp_weight,data=locdf)
    md_model=mediate(mediator_lm,outcome_lm,treat="microbiome",mediator="metabolite")
    modelsumr=summary(md_model)
    # plot(md_model)
    data.frame(microbiomes=locrecrd[,"microbiome"],metabolites=locrecrd[,"metabolites"],outcome=locrecrd[,"outcome"],mediationeffect=modelsumr$d0,mediationeffect_p=modelsumr$d0.p,directeffect=modelsumr$z0,directeffect_p=modelsumr$z0.p,proportion=modelsumr$n0,proportion_p=modelsumr$n0.p,step1coef=summary(mediator_lm)$coefficients["microbiome","Estimate"],step2coef=summary(outcome_lm)$coefficients["metabolite","Estimate"],ci=paste0(formatC(modelsumr$d0.ci,format="E",digits=2),collapse=" "))
}
stattab=Reduce("rbind",tabstat_list)
# reverse run
tabstat_list_rev<-foreach(rowi=seq(nrow(mergtab_sele)),.packages=c("mediation"),.options.snow=opts)%dopar%{
    set.seed(1)
    locrecrd=mergtab_sele[rowi,]
    locdf=as.data.frame(cbind(data.frame(outcome=outcomemat[,locrecrd[,"outcome"]],metabolite=mediationmat[,locrecrd[,"metabolites"]],microbiome=inputmat[,locrecrd[,"microbiome"]]),covarmat))
    locdf=na.omit(locdf)
    mediator_lm=lm(outcome ~ microbiome + sex + bmi_avg_all + samp_weight,data=locdf)
    outcome_lm=lm(metabolite ~ microbiome + outcome + sex + bmi_avg_all + samp_weight,data=locdf)
    md_model=mediate(mediator_lm,outcome_lm,treat="microbiome",mediator="outcome")
    modelsumr=summary(md_model)
    # plot(md_model)
    data.frame(microbiomes=locrecrd[,"microbiome"],metabolites=locrecrd[,"metabolites"],outcome=locrecrd[,"outcome"],mediationeffect_rev=modelsumr$d0,mediationeffect_rev_p=modelsumr$d0.p)
}
stattab_rev=Reduce("rbind",tabstat_list_rev)
# pvalue adjust and clean up
outcomelist=unique(stattab[,"outcome"])
stattab_padj=c()
for(outcomele in outcomelist){
    loctab=stattab[stattab[,"outcome"]==outcomele,]
    stattab_rev_loc=stattab_rev[stattab_rev[,"outcome"]==outcomele,]
    remtab=stattab_rev_loc[p.adjust(stattab_rev_loc[,"mediationeffect_rev_p"],method="fdr")<0.05,]
    indvec=c()
    for(rowi in seq(nrow(remtab))){
        locrecrd=remtab[rowi,]
        matchind=which(loctab[,"microbiomes"]==locrecrd[,"microbiomes"]&loctab[,"metabolites"]==locrecrd[,"metabolites"]&loctab[,"outcome"]==locrecrd[,"outcome"])
        indvec=c(indvec,matchind)
    }
    if(length(indvec)>0){
        loctab=loctab[-indvec,]
    }
    loctab$padj=p.adjust(loctab[,"mediationeffect_p"],method="fdr")
    stattab_padj=rbind(stattab_padj,loctab)
}
save(stattab_padj,file="mediation_micro_metab_clincspike.RData")
# match functional descriptions
idmatchtab=Reduce(rbind,idname_match_list)
rownamevec=stattab_padj[,"microbiomes"]
matchind=match(rownamevec,idmatchtab[,1])
notnaind=which(!is.na(matchind))
rownamevec[notnaind]=idmatchtab[matchind[notnaind],2]
stattab_padj$descriptions=rownamevec
# dominate speceis 
tab_select=stattab_padj[stattab_padj[,"padj"]<0.2&stattab_padj[,"outcome"]%in%class_num_cols,]
# select first reactions for Citramalic acid and Threonic acid
tab_select=tab_select[order(tab_select[,"padj"]),]
rmind=c()
for(filtcmpd in c("Citramalic acid","Threonic acid")){
    nonfirind=tab_select[,"metabolites"]==filtcmpd
    rmind=c(rmind,which(nonfirind)[-1])
}
tab_select=tab_select[-rmind,]
thres_sele=0.1
allfeat=colnames(micromat)
tab_select_new=c()
for(feati in seq(nrow(tab_select))){
    feat=tab_select[feati,"microbiomes"]
    featselecspec=allfeat[str_detect(string=allfeat,pattern=fixed(paste0(feat,"|")))]
    quanvec=apply(micromat[,featselecspec,drop=FALSE],2,mean)
    nsele=ceiling(length(quanvec)*thres_sele)
    featdom=names(sort(quanvec,decreasing=TRUE)[seq(nsele)])
    tab_select_new=rbind(tab_select_new,data.frame(microbiome=featdom,metabolites=rep(tab_select[feati,"metabolites"],times=nsele),outcome=rep(tab_select[feati,"outcome"],times=nsele)))
}
# print(feature_spec_species)
# oberve and change the table manually accordingly
# tab_select_new=data.frame(microbiome=feature_spec_species,metabolites=c("DAG(20:0/20:0)","DAG(20:0/20:0)","C9H15NO3","C9H15NO3","C9H15NO3","C9H15NO3"),outcome=c("ffa_avg_heyjun","ffa_avg_heyjun","hepatic_IR_heyjun","hepatic_IR_heyjun","hepatic_IR_heyjun","insulin_fasting_avg_all"))
tab_specwis=c()
for(rowi in seq(nrow(tab_select_new))){
    locrecrd=tab_select_new[rowi,]
    locdf=as.data.frame(cbind(data.frame(outcome=outcomemat[,locrecrd[,"outcome"]],metabolite=mediationmat[,locrecrd[,"metabolites"]],microbiome=inputmat[,locrecrd[,"microbiome"]]),covarmat))
    locdf=na.omit(locdf)
    mediator_lm=lm(metabolite ~ microbiome + sex + bmi_avg_all + samp_weight,data=locdf)
    outcome_lm=lm(outcome ~ microbiome + metabolite + sex + bmi_avg_all + samp_weight,data=locdf)
    md_model=mediate(mediator_lm,outcome_lm,treat="microbiome",mediator="metabolite")
    modelsumr=summary(md_model)
    temptab=data.frame(microbiomes=locrecrd[,"microbiome"],metabolites=locrecrd[,"metabolites"],outcome=locrecrd[,"outcome"],mediationeffect=modelsumr$d0,mediationeffect_p=modelsumr$d0.p,directeffect=modelsumr$z0,directeffect_p=modelsumr$z0.p,proportion=modelsumr$n0,proportion_p=modelsumr$n0.p)
    tab_specwis=rbind(tab_specwis,temptab)
}
save(tab_specwis,tab_select,file="targeted_grp.RData")
write.table(tab_select,file="mediation_sig_clinic.tsv",row.names=FALSE)