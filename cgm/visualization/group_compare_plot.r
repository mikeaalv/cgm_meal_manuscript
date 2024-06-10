# compare finding related groups of spikers, etc
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
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
load("heatmap_plot.RData")
# mitigator effect for rice spikers
plot_df=as.data.frame(mat_mitigator_effect)
ricespikers=rownames(mat_worstfood)[mat_worstfood[,"Rice"]==1]
plot_df$rice_spiker=rownames(plot_df)%in%ricespikers
plot_df_long=melt(plot_df,id=c("rice_spiker"))
p<-ggplot(plot_df_long,aes(x=variable,y=value,fill=rice_spiker))+geom_boxplot()+xlab("mitigator tyeps")+ylab("mitigator effect")
ggsave(paste0("rice_spiker_mitigator_effect_box.pdf"),plot=p)
# t.test 
p_asso_vec=c()
for(col in colnames(mat_mitigator_effect)){
    loctab=plot_df_long[plot_df_long[,"variable"]==col,]
    p_asso_vec=c(p_asso_vec,t.test(value~rice_spiker,data=loctab)$p.value)
}
names(p_asso_vec)=colnames(mat_mitigator_effect)
p.adjust(p_asso_vec,method="fdr")
# all spiker groups
plot_df=as.data.frame(mat_mitigator_effect)
spikervec=rep(NA,times=dim(plot_df)[1])
names(spikervec)=rownames(plot_df)
for(spikfood in colnames(mat_worstfood)){
    spikertye=rownames(mat_worstfood)[mat_worstfood[,spikfood]==1]
    spikervec[spikertye]=spikfood
}
plot_df$spikertype=as.factor(spikervec)
plot_df_long2=melt(plot_df,id=c("spikertype"))
p<-ggplot(plot_df_long2,aes(x=variable,y=value,fill=spikertype))+geom_boxplot()+xlab("mitigator tyeps")+ylab("mitigator effect")
ggsave(paste0("spiker_group_mitigator_effect_box.pdf"),plot=p)
# 
pvec=c()
mitigators=unique(plot_df_long[,"variable"])
for(mitigator in mitigators){
    loctab=plot_df_long[plot_df_long[,"variable"]==mitigator,]
    loctab=loctab[!is.na(loctab[,"value"]),]
    aovtest=aov(value~rice_spiker,loctab)
    pvec=c(pvec,summary(aovtest)[[1]][["Pr(>F)"]][1])
}
names(pvec)=mitigators
p.adjust(pvec,method="fdr")
# blood test and metabolic testing results among different spiker groups
featlist=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","sspg_avg_all","systolic_avg_all","diastolic_avg_all","ffa_avg_heyjun","ie_heyjun")
plot_df=metadatashow_plot[,featlist]
vec_spiker=rep(NA,times=dim(mat_worstfood)[1])
for(spiker in colnames(mat_worstfood)){
    vec_spiker[mat_worstfood[,spiker]==1]=spiker
}
plot_df$spiker=vec_spiker
plot_df_long=melt(plot_df,id=c("spiker"),na.rm=TRUE)
p<-ggplot(plot_df_long,aes(x=variable,y=value,fill=spiker))+geom_boxplot()+xlab("blood test values")+ylab("values")
ggsave(paste0("spiker_bloodtestvalue_box.pdf"),plot=p)
pvec=c()
variables=unique(plot_df_long[,"variable"])
for(var in variables){
    loctab=plot_df_long[plot_df_long[,"variable"]==var,]
    loctab=loctab[!is.na(loctab[,"value"]),]
    aovtest=aov(value~spiker,loctab)
    pvec=c(pvec,summary(aovtest)[[1]][["Pr(>F)"]][1])
}
names(pvec)=variables
# Rice spiker vs Asian population
stat_df=metadatashow_plot[,c("study_id","ethnicity")]
stat_df$asian=ifelse(stat_df$ethnicity=="Asian","Asian","non-Asian")
stat_df$rice_spiker=rownames(stat_df)%in%ricespikers
print(fisher.test(x=as.factor(stat_df[,"asian"]),y=as.factor(stat_df[,"rice_spiker"])))
# bar plot
stat_df$rice_peaks=ifelse(stat_df$rice_spiker,"spiker","non-spiker")
n_df_ch=as.data.frame(table(stat_df[,c("asian","rice_peaks")]))
p<-ggplot(n_df_ch,aes(x=asian,y=Freq,fill=rice_peaks))+geom_bar(stat="identity",position="dodge")+ylab("number")
ggsave(paste0("rice_spiker_ethnicity_freq.pdf"),plot=p)
# mitigation effect bar(increase/decrease)
plotdf=c()
featureslist=colnames(mat_mitigator_effect)
for(feathere in featureslist){
    subtab=mat_mitigator_effect[,feathere]
    valuen=c(sum(subtab>0,na.rm=TRUE),-sum(subtab<0,na.rm=TRUE))
    plotdf=rbind(plotdf,data.frame(feature=c(feathere,feathere),value=valuen,direction=c("posi","nega")))
}
p=ggplot(plotdf,aes(x=feature,y=value,fill=direction))+geom_bar(stat="identity",position="identity")
ggsave(plot=p,filename=paste0(resdir,"mitigator_numb.pdf"),device="pdf")