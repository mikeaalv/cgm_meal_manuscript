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
require(data.table)
require(reshape2)
require(ComplexUpset)
# 
projdir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/"
setwd(paste0(projdir,"result/"))
# metadata
meta_data=read.table(file=paste0(projdir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE)
targmis_tab=meta_data[,c("study_id","age_today_avg_all","bmi_avg_all","a1c_avg_all","sspg_avg_all","ogtt_t_120_avg_all","iigi_t_120_avg_all")]
misind=which(rowSums(is.na(targmis_tab))>1)
misids<-targmis_tab[misind,"study_id"]%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()
# cgm
cgmdata=read.table(file=paste0(projdir,"result/cgm_meal/cgm_foods_smooth.txt"),header=TRUE)
id_cgm=unique(cgmdata[,"subject"])
id_mitigators=unique(cgmdata[cgmdata[,"mitigator"]!="","subject"])
# omics data
lipidomics=read.table(file=paste0(projdir,"data/metabolomics_lipidomics/CGM_Plasma_Lipid_Species(1).csv"),header=TRUE,sep=",")
id_lipidomics=colnames(lipidomics)[-1]%>%str_extract(string=.,pattern="^X\\d+")%>%str_replace_all(string=.,pattern="^X",replacement="")%>%as.numeric()%>%unique()
metabomics=read.table(file=paste0(projdir,"data/metabolomics_lipidomics/CGM_plasma_metabolomics_all_Anno_metID(1).csv"),header=TRUE,sep=",")
id_metabomics=colnames(metabomics)[-1]%>%str_extract(string=.,pattern="^X\\d+")%>%na.omit()%>%str_replace_all(string=.,pattern="^X",replacement="")%>%as.numeric()%>%unique()
proteomics=read.table(file=paste0(projdir,"data/proteomics/StanfordProject_FahimAbbasi_NPX_2022-08-24_CGM1study.csv"),header=TRUE,sep=";")
id_proteomics<-proteomics[,"SampleID"]%>%str_extract(string=.,pattern="^STUDYID-\\d+")%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()%>%unique()
load(paste0(projdir,"result/microbiome/eda/microbiome_processed.RData"))
id_microbiome=unique(rownames(micromat))
allids_useful=Reduce(union,list(id_cgm,id_lipidomics,id_metabomics,id_proteomics,id_microbiome))
id_nometa=intersect(allids_useful,misids)
# missing data record
restab=targmis_tab[targmis_tab[,"study_id"]%in%paste0("STUDYID-",str_pad(id_nometa,3,pad="0")),]
write.table(restab,file="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/result/missingdata.csv",sep = ",",row.names=FALSE)
# fitbit data
load("/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/archive/Fitbit_data_IRpred/cgm1/rawdata/sleep_data.RData")
fitbitdata=as.data.frame(data)
id_fitbit<-fitbitdata[,"subject.id"]%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()%>%unique()
# upset plot of all data sources
metadata_term_list=list(demographic=c("age_today_avg_all","bmi_avg_all"),blood=c("a1c_avg_all"),ogtt=c("ogtt_t_120_avg_all"),sspg=c("sspg_avg_all"),iigi=c("iigi_t_120_avg_all"))
metadata_id_list=list()
metadata_id_list[["cgm"]]=id_cgm
metadata_id_list[["mitigator"]]=id_mitigators
metadata_id_list[["metab_lipidomics"]]=intersect(id_lipidomics,id_metabomics)
metadata_id_list[["proteomics"]]=id_proteomics
metadata_id_list[["microbiome"]]=id_microbiome
metadata_id_list[["fitbit"]]=id_fitbit
for(sourcetype in names(metadata_term_list)){
    ind=which(rowSums(is.na(targmis_tab[,metadata_term_list[[sourcetype]],drop=FALSE]))==0)
    ids<-targmis_tab[ind,"study_id"]%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()%>%unique()
    metadata_id_list[[sourcetype]]=ids
}
save(metadata_id_list,file=paste0(projdir,"data/metadata/data_avail.RData"))
allids=Reduce(union,metadata_id_list)
covars=names(metadata_id_list)
tabsig=matrix(NA,nrow=length(allids),ncol=length(metadata_id_list))
rownames(tabsig)=allids
colnames(tabsig)=covars
count_val=c()
for(sourcetype in covars){
    tabsig[,sourcetype]=allids%in%metadata_id_list[[sourcetype]]
    count_val=c(count_val,sum(tabsig[,sourcetype]))
}
# df bar
count_val=c(count_val,
            sum(apply(tabsig[,c("cgm","metab_lipidomics","microbiome","demographic","blood","ogtt","sspg","iigi")],1,all)),
            sum(apply(tabsig[,c("cgm","demographic","blood")],1,all)),
            sum(apply(tabsig[,c("cgm","metab_lipidomics","microbiome","demographic","blood")],1,all)),
            sum(apply(tabsig[,c("mitigator","metab_lipidomics","microbiome","demographic","blood")],1,all)))
bartypes=c(covars,"all_exp_fitbit_proteomics","cgm_demo_blood","cgm_demo_blood_metlip_micro","mitigator_demo_blood_metlip_micro")
bar_df=data.frame(type=factor(bartypes,levels=bartypes),count=count_val)
p<-ggplot(bar_df,aes(x=type,y=count))+geom_bar(stat="identity")+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
ggsave(plot=p,filename=paste0("datasource_bar_intersct.pdf"))
# mark other data availablity
tab_newval=as.data.frame(read_excel(paste0(projdir,"data/metadata/Data_Status_Summary_for_Metadata_052623.xlsx")))
sub_id_nece=tab_newval[,"Participants"]
metadata_id_list[["newvalidate"]]=sub_id_nece
needed_subj=Reduce(union,metadata_id_list[c("cgm","metab_lipidomics","proteomics","newvalidate")])
needed_subj_ids=paste0("STUDYID-",str_pad(needed_subj,3,pad="0"))
meta_data$needed=meta_data[,"study_id"]%in%needed_subj_ids
#update needed samples redo./....move code
sampneeds_tab=read.table(paste0(projdir,"data/metadata/sample_needs_full.txt"),sep="\t",header=TRUE)
updneeded_vec=sampneeds_tab[match(sampneeds_tab[,"study_id"],meta_data[,"study_id"]),"needed"]
meta_data$needed=meta_data$needed|updneeded_vec
write.table(meta_data,file=paste0(projdir,"data/metadata/metadata_clean_all.tmp.tsv"),sep="\t",row.names=FALSE)
