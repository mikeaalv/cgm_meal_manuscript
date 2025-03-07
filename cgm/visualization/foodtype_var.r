# consistency of carb-response-type assignment
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
# foods spikers (worst food)
food_list_carbs=c("Beans","Berries","Bread","Grapes","Pasta","Potatoes","Rice")
summartab_mean <- feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(meanval=mean(peak_relative),.groups="keep") %>% group_by(subject) %>% summarize(worstfood=foods[which.max(meanval)],.groups="keep") %>% as.data.frame()
ressummtab=c()
for(therep in c(1,2)){
    feattab_loc=feattab[feattab[,"rep"]==therep,]
    summartab_1st <- feattab_loc %>% filter(foods %in% food_list_carbs) %>% group_by(subject) %>% summarize(worstfood=foods[which.max(peak_relative)],.groups="keep") %>% as.data.frame()
    colnames(summartab_1st)[2]="worst_1st"
    summartab_2nd <- feattab_loc %>% filter(foods %in% food_list_carbs) %>% group_by(subject) %>% summarize(worstfood=foods[order(peak_relative,decreasing=TRUE)[2]],.groups="keep") %>% as.data.frame()
    colnames(summartab_2nd)[2]="worst_2nd"
    rescomp_tab=merge(merge(summartab_mean,summartab_1st,by="subject"),summartab_2nd,by="subject")
    ressummtab=rbind(ressummtab,rescomp_tab)
}
statdf=c()
for(thefood in food_list_carbs){
    subtab=ressummtab[ressummtab[,"worstfood"]==thefood,]
    nrep=nrow(subtab)
    statdf=rbind(statdf,data.frame(food=thefood,"worstfood"=sum(subtab[,"worst_1st"]==thefood,na.rm=TRUE)/nrep,"worstfood_2nd"=sum(subtab[,"worst_2nd"]==thefood,na.rm=TRUE)/nrep))
}
statdf=statdf[!is.na(statdf[,"worstfood"]),]
statdf_long=melt(statdf,id=c("food"))
statdf_long$variable=factor(statdf_long$variable,levels=c("worstfood_2nd","worstfood"))
p<-ggplot(statdf_long,aes(x=food,y=value,fill=variable))+geom_bar(stat="identity",position="stack")+ylab("ratio")
ggsave(paste0("rep_carb_resp_type_repro_proport.pdf"),plot=p)
savdf=statdf_long
write.table(statdf_long,file=paste0("figext4.txt"))