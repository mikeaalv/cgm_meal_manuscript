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
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/meal_prop/")
setwd(resdir)
datadir=paste0(comp,"RAWDATAFOLDER")
# 
dirlist=list.files(path=datadir,pattern="\\.cron.csv",recursive=TRUE,ignore.case=TRUE)
subjectids_vec=str_extract(string=dirlist,pattern="^\\d+-\\d+")
subjectids=unique(subjectids_vec)
filetab=data.frame(subjectid=subjectids_vec,cronometerfile=dirlist)
navec=rep(NA,times=length(subjectids))
# 
tab_food_prop=data.frame("subject_id"=subjectids,"rice"=navec,"berries"=navec,"grapes"=navec,"bread"=navec,"pasta"=navec,"beans"=navec,"potatoes"=navec,"proteins"=navec,"FrenchFries"=navec,"chicken"=navec,"fish"=navec)
tab_food_energy_prop=tab_food_prop
matchpattern=list("rice"="rice","berries"="berr","grapes"="grape","bread"="bread","pasta"=c("pasta","Angel","macaroni","spaghett","penne","lasagna"),"beans"="bean","potatoes"=c("potato","hash brown","fries"),"FrenchFries"="fries","chicken"="chicken","fish"=c("fish","salmon"))
for(subj in subjectids){
    files=filetab[filetab[,"subjectid"]==subj,"cronometerfile"]
    tab_merg=c()
    for(thfile in files){
        loctab=read.table(paste0(datadir,"/",thfile),header=TRUE,sep=",",quote="\"")
        if(ncol(loctab)<10){
            next
        }
        loctab=loctab[,c("Food.Name","Energy..kcal.","Net.Carbs..g.","Protein..g.")]
        tab_merg=rbind(tab_merg,loctab)
    }
    subjind=which(tab_food_prop[,"subject_id"]==subj)
    # foods match
    for(food in names(matchpattern)){
        patternvec=matchpattern[[food]]
        inds=c()
        for(patt in patternvec){
            ind=str_which(string=tab_merg[,"Food.Name"],pattern=fixed(patt,ignore_case=TRUE))
            inds=c(inds,ind)
        }
        tab_food_prop[subjind,food]=length(inds)/dim(tab_merg)[1]
        tab_food_energy_prop[subjind,food]=sum(tab_merg[inds,"Energy..kcal."],na.rm=TRUE)/sum(tab_merg[,"Energy..kcal."],na.rm=TRUE)
    }
    # proteins
    tab_food_prop[subjind,"proteins"]=sum(tab_merg[,"Protein..g."],na.rm=TRUE)/sum(tab_merg[,"Net.Carbs..g."],na.rm=TRUE)
}
save(tab_food_prop,tab_food_energy_prop,file="meal_prop.RData")
