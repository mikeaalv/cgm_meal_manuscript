# compare finding related groups of spikers carb-response-types, etc (based on the quantile result)
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
ricespikers=rownames(mat_worstfood_thre)[mat_worstfood_thre[,"Rice"]==1]
plot_df$rice_spiker=rownames(plot_df)%in%ricespikers
plot_df_long=melt(plot_df,id=c("rice_spiker"))
p<-ggplot(plot_df_long,aes(x=variable,y=value,fill=rice_spiker))+geom_boxplot()+xlab("mitigator types")+ylab("mitigator effect")
ggsave(paste0("rice_spiker_mitigator_effect_box_thres.pdf"),plot=p)
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
# Rice spiker vs Asian population
stat_df=metadatashow_plot[,c("study_id","ethnicity")]
stat_df$asian=ifelse(stat_df$ethnicity=="Asian","Asian","non-Asian")
stat_df$rice_spiker=rownames(stat_df)%in%ricespikers
print(fisher.test(x=as.factor(stat_df[,"asian"]),y=as.factor(stat_df[,"rice_spiker"])))
# bar plot
stat_df$rice_peaks=ifelse(stat_df$rice_spiker,"spiker","non-spiker")
n_df_ch=as.data.frame(table(stat_df[,c("asian","rice_peaks")]))
p<-ggplot(n_df_ch,aes(x=asian,y=Freq,fill=rice_peaks))+geom_bar(stat="identity",position="dodge")+ylab("number")
ggsave(paste0("rice_spiker_ethnicity_freq_thres.pdf"),plot=p)
# meta data comparison between groups
plotdf=c()
for(food in food_list_carbs){
    subtab1=mat_worstfood_thre[mat_worstfood_thre[,food]==1,food,drop=FALSE]
    subtab2=metadatashow_plot[rownames(subtab1),c(numcolssub,"ethnicity")]
    subtab2$worstfood=rep(food,times=dim(subtab2)[1])
    plotdf=rbind(plotdf,subtab2)
}
foodord=c("Pasta","Grapes","Rice","Bread","Potatoes","Berries","Beans")
plotdf[,"worstfood"]=factor(plotdf[,"worstfood"],level=foodord)
comp_feat_list=list("Grapes"=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","sspg_avg_all"),
                    "Potatoes"=c("a1c_avg_all","fbg_avg_all","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","bmi_avg_all"),
                    "Bread"=c("systolic_avg_all","diastolic_avg_all"))
for(center_gr in names(comp_feat_list)){
    compfeat=comp_feat_list[[center_gr]]
    comparisons=list()
    for(food in setdiff(food_list_carbs,center_gr)){
        comparisons[[length(comparisons)+1]]=c(food,center_gr)
    }
    for(colfeat in compfeat){
        p<-ggboxplot(plotdf,x="worstfood",y=colfeat)+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
        ggsave(plot=p,paste0(resdir,"boxplot_spiketype_meta",colfeat,"_compare_",center_gr,".thres.pdf"))
    }
}