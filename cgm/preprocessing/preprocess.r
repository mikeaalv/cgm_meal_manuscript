#  CGM data smoothing
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
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
data_pre=paste0(pardir,"result/cgm_meal_preprocessing/")
setwd(resdir)
# 
tab_cgm=read.table(paste0(data_pre,"cgm_foods.csv"),header=TRUE,sep=",")
tab_cgm=as.data.frame(tab_cgm)
# summartab<- tab_cgm %>% group_by(subject,foods,rep) %>% summarize(min=min(mins_since_start),max=max(mins_since_start),n=n(),.groups="keep")
# sort(summartab[["max"]])[1:60]
# sort(summartab[["min"]])
tab_cgm=na.omit(tab_cgm)
rangeselec=c(-25,170)
gap_h=30#minutes
ntime=40
timegrid=seq(from=rangeselec[1],to=rangeselec[2],length.out=ntime)
groupsplits=split(tab_cgm,list(tab_cgm$subject,tab_cgm$foods,tab_cgm$rep))
groupsplits=groupsplits[sapply(groupsplits,function(x) dim(x)[1])>1]
tab_cgm_smoothed_list=c()
groupname=names(groupsplits)
for(groupi in seq(length(groupsplits))){#
    locdf=groupsplits[[groupi]]
    timevec=locdf[,"mins_since_start"]
    start_flag=min(timevec)>rangeselec[1]
    end_flag=max(timevec)<rangeselec[2]
    gap_flag=any(abs(diff(timevec))>gap_h)
    if(start_flag | end_flag | gap_flag){
        print(paste0(groupi," ",groupname[[groupi]],start_flag,end_flag,gap_flag))
        next
    }
    yvec=locdf[,"glucose"]
    ssmod=ss(x=timevec,y=yvec,all.knots=TRUE,m=3,method="REML")
    yvec_new=predict(ssmod,x=timegrid)$y
    # pdf(paste0(groupi,".pdf"))
    # plot(timevec,yvec)
    # lines(ssmod)
    # dev.off()
    # 
    temptab=data.frame(glucose=yvec_new,subject=unique(locdf[,"subject"]),foods=unique(locdf[,"foods"]),mitigator=unique(locdf[,"mitigator"]),food=unique(locdf[,"food"]),rep=unique(locdf[,"rep"]),mins_since_start=timegrid,dates=unique(as.Date(locdf[,"datetime_local"])))
    tab_cgm_smoothed_list[[length(tab_cgm_smoothed_list)+1]]=temptab
}
tab_cgm_smoothed=bind_rows(tab_cgm_smoothed_list)
write.table(tab_cgm_smoothed,file="cgm_foods_smooth.txt",row.names=FALSE)