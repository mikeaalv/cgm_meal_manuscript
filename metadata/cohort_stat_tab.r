# number of samples and features for all data sets
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
# 
projdir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/"
setwd(paste0(projdir,"result/metadata/"))
# metadata
meta_data=read.table(file=paste0(projdir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE)
meta_data=meta_data[meta_data[,"study_id"]!="STUDYID-SOMESUBJECTID",]
meta_data[,"subject"]<-meta_data[,"study_id"]%>%str_remove(string=.,pattern="STUDYID-")%>%as.numeric()
# cgm
cgmdata=read.table(file=paste0(projdir,"result/cgm_meal/cgm_foods_smooth.txt"),header=TRUE)
id_cgm=unique(cgmdata[,"subject"])
meta_data=meta_data[meta_data[,"subject"]%in%id_cgm,]
# cohort statistics table
tabstat=c()
numlist=list("Age, y"="age_today_avg_all","BMI, kg/m2"="bmi_avg_all","Systolic Blood Pressure, mmHg"="systolic_avg_all","Diastolic Blood Pressure, mmHg"="diastolic_avg_all","HbA1c, %"="a1c_avg_all","Fasting Plasma Glucose, mg/dL"="fbg_avg_all","Fasting Insulin, mmol/L"="insulin_fasting_avg_all","Fructosamine umol/L"="fructosamine_avg_all","Total Cholesterol, mg/dL"="cholesterol_total_avg_all","Triglyceride, mg/dL"="triglyceride_avg_all","HDL, mg/dL"="hdl_avg_all","LDL, mg/dL"="ldl_avg_all","hs-CRP, mg/L"="hscrp_avg_all","ALT/SGPT, U/L"="alt_sgpt_avg_all","Creatine,mg/dL"="creatinine_avg_all")
classlist=list("Sex, F/M"="sex.factor","Ethnicity, Caucasian/Asian/Hispanic/Others"="ethnicity","Diabetes/PreDM/Euglycemia"="a1c_t2d_status_bl")
classlevels=list("sex.factor"=c("Female","Male"),"ethnicity"=c("White","Asian","Hispanic or Latino","Asian, White"),"a1c_t2d_status_bl"=c("T2D","preDM","Normal"))
for(col in names(numlist)){
    valvec=meta_data[,numlist[[col]]]
    meanval=mean(valvec,na.rm=TRUE)
    sdval=sd(valvec,na.rm=TRUE)
    tabstat=rbind(tabstat,data.frame(variable=col,output=paste0(signif(meanval,digits=3),"\u00b1",signif(sdval,digits=3))))
}
for(col in names(classlist)){
    valvec=meta_data[,classlist[[col]]]
    valvec=valvec[!is.na(valvec)]
    numvec=table(valvec)[classlevels[[classlist[[col]]]]]
    tabstat=rbind(tabstat,data.frame(variable=col,output=paste0(c(numvec),collapse="/")))
}
write.table(tabstat,file="summ_stat.txt",row.names=FALSE)

