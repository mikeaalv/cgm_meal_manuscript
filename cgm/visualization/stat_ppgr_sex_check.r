# Check sex as co-variate and interaction terms
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
require(plyr)
require(dplyr)
require(data.table)
require(ComplexHeatmap)
require(ggpubr)
require(circlize)
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
# 
feattab=read.table("cgm_foods_manual_features.csv",header=TRUE,sep=",")
feattab=feattab[,c("subject","foods","rep","mitigator","peak_value","baseline_glucose","AUC_above_baseline","time_to_peak")]
feattab$peak_relative=feattab$peak_value-feattab$baseline_glucose
metadata=as.data.frame(read.table(paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
metadata=metadata[,c("study_id","sex.factor")]
metadata[,"study_id"]<-metadata[,"study_id"]%>%str_replace_all(string=.,pattern="STUDYID",replacement="")%>%as.numeric()
colnames(metadata)=c("subject","sex")
metadata=metadata[!is.na(metadata[,"sex"]),]
# foods spikers (worst food)
food_list_carbs=c("Beans","Berries","Bread","Grapes","Pasta","Potatoes","Rice","Rice+Fat","Rice+Fiber","Rice+Protein")
feattab=feattab[feattab[,"foods"]%in%food_list_carbs,]
plotfeats=c("peak_relative","AUC_above_baseline","time_to_peak")
for(pfeat in plotfeats){
    feattab_loc=merge(feattab,metadata,by="subject")
    feattab_loc=feattab_loc[,c("subject","foods","sex",pfeat)]
    p<-ggbarplot(feattab_loc,x="sex",y=pfeat,color="sex",add=c("mean_se","jitter"),facet.by="foods",short.panel.labs=FALSE,ylim=c(0,max(feattab_loc[,pfeat])*1.2))+stat_compare_means(aes(label=paste0("p = ",after_stat(p.format))),hide.ns=FALSE)
    ggsave(plot=p,paste0(resdir,"bar_ppgr_sex_",pfeat,".pdf"))
}
# stat model with sex added
load("heatmap_plot.RData")
# data preparation 
temptab=summtab2[match(subjects,summtab2[,"subject"]),c("subject","ricespiker","potato_vs_grape_val")]
temptab$subject=subjects
stattab=temptab
temptab=mat_mitigator_effect
colnames(temptab)=c("Fat","Fiber","Protein")
stattab=cbind(stattab,temptab[subjects,])
metadatashow[,"subject"]<-metadatashow[,"study_id"]%>%str_remove(string=.,pattern="STUDYID")%>%as.numeric()
temptab=metadatashow[,c("subject","age_today_avg_all","bmi_avg_all","sspg_status_heyjun","DI_status_heyjun","ethnicity","a1c_avg_all","fbg_avg_all","ffa_avg_heyjun","sspg_avg_all","systolic_avg_all","diastolic_avg_all","ethnicity","modified_DI_heyjun","sex.factor")]
stattab=cbind(stattab,temptab[match(subjects,temptab[,"subject"]),-1])
temptab=as.data.frame(mat_worstfood[,"Bread",drop=FALSE])
colnames(temptab)="breadspiker"
stattab=cbind(stattab,temptab[subjects,,drop=FALSE])
temptab=mat_quantile[,c("Potatoes","Rice")]
stattab=cbind(stattab,temptab[subjects,])
# 
stattab$"IsAsian"=ifelse(str_detect(stattab[,"ethnicity"],pattern="Asian"),0.5,-0.5)
stattab$"IsBreadspiker"=ifelse(stattab[,"breadspiker"]==1,0.5,-0.5)
stattab$"IsRicespiker"=ifelse(stattab[,"ricespiker"]=="ricespiker",0.5,-0.5)
stattab$"IsMuscleIR"=ifelse(stattab[,"sspg_status_heyjun"]=="IR",0.5,-0.5)
stattab$"IsBetaDys"=ifelse(stattab[,"DI_status_heyjun"]=="BC_dys",0.5,-0.5)
stattab$sex=ifelse(stattab[,"sex.factor"]=="Female",-0.5,0.5)
# rename columns
colnames(stattab)=c("subject","ricespiker","potato_vs_grape_ratio","Fat_mitigation_effect","Fiber_mitigation_effect","Protein_mitigation_effect","age","bmi","sspg_status","DI_status","ethnicity","a1c","fbg","FFA","sspg","systolic_bp","diastolic_bp","ethnicity.1","DI","sex.factor","breadspiker","Potatoes","Rice","IsAsian","IsBreadspiker","IsRicespiker","IsMuscleIR","IsBetaDys","sex")
# model list
modelist=list("potato_IR"="Potatoes~IsMuscleIR","potato_beta"="Potatoes~IsBetaDys","Fat_IR"="Fat_mitigation_effect~IsMuscleIR","Fiber_IR"="Fiber_mitigation_effect~IsMuscleIR","Protein_IR"="Protein_mitigation_effect~IsMuscleIR","Fat_Bcell"="Fat_mitigation_effect~IsBetaDys","Fiber_Bcell"="Fiber_mitigation_effect~IsBetaDys","Protein_Bcell"="Protein_mitigation_effect~IsBetaDys","Rice_Asian"="Rice~IsAsian","Fat_Rice"="Fat_mitigation_effect~IsRicespiker","Fiber_Rice"="Fiber_mitigation_effect~IsRicespiker","Protein_Rice"="Protein_mitigation_effect~IsRicespiker","SSPG_PG"="sspg~potato_vs_grape_ratio","FFA_PG"="FFA~potato_vs_grape_ratio","A1C_PG"="a1c~potato_vs_grape_ratio","FBG_PG"="fbg~potato_vs_grape_ratio","Sys_Bread"="systolic_bp~IsBreadspiker","Dia_Bread"="diastolic_bp~IsBreadspiker","PGR_IR"="potato_vs_grape_ratio~IsMuscleIR","PGR_beta"="potato_vs_grape_ratio~IsBetaDys","comb_miti_protein"="Protein_mitigation_effect~DI+sspg","comb_miti_fat"="Fat_mitigation_effect~DI+sspg","comb_miti_fiber"="Fiber_mitigation_effect~DI+sspg")
addterms="+bmi+age+sex"
statcoll=c()
for(modelname in names(modelist)){
    formula_mlm_str=paste0(modelist[[modelname]],addterms)
    formula_mlm=as.formula(formula_mlm_str)
    termsall=attr(terms(lme4:::nobars(formula_mlm)),"term.labels")
    targetterm=termsall[1]
    interterm=paste0(targetterm,":sex")
    formula_mlm_str=paste0(formula_mlm_str," + ",interterm)
    formula_mlm=as.formula(formula_mlm_str)
    lmdol=lm(formula_mlm,stattab)
    coeftab=summary(lmdol)$coefficients
    statcoll=rbind(statcoll,data.frame(Y=all.vars(formula_mlm)[1],X=targetterm,coef=coeftab[targetterm,"Estimate"],pval=coeftab[targetterm,"Pr(>|t|)"],coef_sex=coeftab["sex","Estimate"],p_sex=coeftab["sex","Pr(>|t|)"],p_sex_inter=coeftab[interterm,"Pr(>|t|)"],coef_sex_inter=coeftab[interterm,"Estimate"],formula=formula_mlm_str,statistics=coeftab[targetterm,"t value"],ste=coeftab[targetterm,"Std. Error"],Df=summary(lmdol)$df[2]))
}
# reformat 
save(statcoll,file="comb_model_assocheck_sexadd.RData")
write.table(statcoll,file="comb_model_assocheck_sexadd.txt",row.names=FALSE)