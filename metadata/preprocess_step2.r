# process the metadata from manual cleaning and merge with other data sources
# the 2nd step of metadata processing
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
require(data.table)
require(reshape2)
# read meta data raw file
projdir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/"
dirdata=paste0(projdir,"data/metadata/")
metatab=as.data.frame(read_excel(paste0(dirdata,"metadata_upd_042023_manuclean.xlsx"),na=c("","NA","na","N/A","NaN")))
tab_metabolictesting=read.table(paste0(projdir,"data/washu_combined_batches_09302021_metabolic_testing.csv"),header=TRUE,sep=",")
tab_metabolictesting=tab_metabolictesting[tab_metabolictesting[,"TimePoint"]==0,c("SubjectID","Insulin","C_peptide","GLP1","GIP","Glucagon_pmol_L","VitaminD")]
summartab<-tab_metabolictesting %>% group_by(SubjectID) %>% summarize(across(where(is.numeric),~ mean(.x,na.rm=TRUE),.names="{.col}_metabolic_testing_mean"))%>%as.data.frame()
summartab$study_id<-summartab[,"SubjectID"]%>%str_replace_all(string=.,pattern=fixed("S"),replacement="")%>%str_pad(.,3,pad="0")%>%paste0("STUDYID-",.)
metatab2=merge(metatab,summartab,by="study_id",all=TRUE)
# correct spelling
for(col in c("a1c_t2d_status_bl","fbg_t2d_status_bl")){
    vectemp=metatab2[,col]
    vectemp[vectemp=="PreDM"]="preDM"
    metatab2[,col]=vectemp
}
for(col in c("HBP_Med_bl","Statin_bl")){
    vectemp=metatab2[,col]
    vectemp[vectemp=="0"]="no"
    vectemp[vectemp=="1"]="yes"
    metatab2[,col]=vectemp
}
# # fill avg col with bl values if missing
# listfill=list("a1c_avg_all"="a1c_bl","fbg_avg_all"="fbg_bl")
# for(tocol in names(listfill)){
#     fromcol=listfill[[tocol]]
#     inds=which(is.na(metatab2[,tocol])&(!is.na(metatab2[,fromcol])))
#     metatab2[inds,tocol]=metatab2[inds,fromcol]
# }
write.table(metatab2,file=paste0(dirdata,"metadata_clean_all_ver2.tsv"),row.names=FALSE)
# QC of the metadata table 
metatab2_test=metatab2
clasvec=sapply(metatab2_test,class)
print(table(clasvec))
chind=which(clasvec=="character")
for(ind in chind){
    metatab2_test[,ind]=as.factor(metatab2_test[,ind])
}
print(summary(metatab2_test))