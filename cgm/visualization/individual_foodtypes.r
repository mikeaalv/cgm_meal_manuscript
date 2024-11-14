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
require(lsr)
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
summartab <- feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(meanval=mean(peak_relative),.groups="keep") %>% group_by(subject) %>% summarize(worstfood=foods[which.max(meanval)],.groups="keep") %>% as.data.frame()
mat_worstfood=matrix(0,nrow=nrow(summartab),ncol=length(food_list_carbs))
rownames(mat_worstfood)=paste0(summartab[,1])
colnames(mat_worstfood)=food_list_carbs
for(rowi in seq(dim(summartab)[1])){
    mat_worstfood[paste0(summartab[rowi,1]),summartab[rowi,2]]=1
}
subjects=rownames(mat_worstfood)
nsubjects=length(subjects)
save(summartab,file="worst.food.RData")
# targeted comparison on spikers
sumtab_quan<-feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(meanval=mean(peak_relative),.groups="keep") %>% group_by(subject)
summtab2<-sumtab_quan%>%summarize(ricespiker=ifelse(foods[which.max(meanval)]=="Rice","ricespiker","noricespiker"),potato_vs_grape_val=meanval[foods=="Potatoes"]/meanval[foods=="Grapes"],.groups="keep")%>% as.data.frame()
summtab2$potato_vs_grape=ifelse(summtab2$potato_vs_grape_val>1,"potatospiker","grapespiker")
save(summtab2,file="spiker_ben.RData")
# quantile matrix
feattab$peakquantile=percent_rank(feattab$peak_relative)
feattab_avg_q <- feattab%>%filter(foods %in% food_list_carbs)%>%group_by(subject,foods)%>%summarise(avg_peak_q=mean(peakquantile),.groups="keep")%>%as.data.frame()
mat_quantile=matrix(NA,nrow=nsubjects,ncol=length(food_list_carbs))
rownames(mat_quantile)=subjects
colnames(mat_quantile)=food_list_carbs
for(rowi in seq(dim(feattab_avg_q)[1])){
    mat_quantile[paste0(feattab_avg_q[rowi,1]),feattab_avg_q[rowi,2]]=feattab_avg_q[rowi,3]
}
save(feattab_avg_q,mat_quantile,file="worst.food.quantile.RData")
# quantile worst food
mat_worstfood_thre=matrix(0,nrow=nrow(summartab),ncol=length(food_list_carbs))
rownames(mat_worstfood_thre)=paste0(summartab[,1])
colnames(mat_worstfood_thre)=food_list_carbs
for(food in food_list_carbs){
    nselect=floor(sum(!is.na(mat_quantile[,food]))*0.2)
    indsele=order(mat_quantile[,food],decreasing=TRUE)[seq(nselect)]
    mat_worstfood_thre[indsele,food]=1
}
# matrix of three features mean [peak_relative, AUC, time_to_peak]
summartab <- feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(peakmean=mean(peak_relative),aucmean=mean(AUC_above_baseline),time2peakmean=mean(time_to_peak),.groups="keep") %>% group_by(subject) %>% as.data.frame()
featlist_intr=c("peakmean","aucmean","time2peakmean")
featfullmat=matrix(NA,nrow=nsubjects,ncol=length(featlist_intr)*length(food_list_carbs))
rownames(featfullmat)=subjects
colnames(featfullmat)=paste0(rep(featlist_intr,each=length(food_list_carbs)),rep(food_list_carbs,times=length(featlist_intr)))
for(subjecti in seq(nsubjects)){
    locatab=summartab[summartab[,"subject"]==subjects[subjecti],]
    locvec=c()
    for(feat in featlist_intr){
        matchind=sapply(food_list_carbs,function(x){
            ind=which(locatab[,"foods"]==x)
            ifelse(length(ind)==0,NA,ind)
        })
        locvec=c(locvec,scale(locatab[matchind,feat]))#scaled within each subject
    }
    featfullmat[subjecti,]=locvec
}
# matrix of three features mean [peak_relative, AUC, time_to_peak] with mitigators
fullfoodlist=c(food_list_carbs,"Rice+Fat","Rice+Fiber","Rice+Protein")
summartab <- feattab %>% filter(foods %in% fullfoodlist) %>% group_by(subject,foods) %>% summarize(peakmean=mean(peak_relative),aucmean=mean(AUC_above_baseline),time2peakmean=mean(time_to_peak),.groups="keep") %>% group_by(subject) %>% as.data.frame()
featlist_intr=c("peakmean","aucmean","time2peakmean")
featfullmat_full=matrix(NA,nrow=nsubjects,ncol=length(featlist_intr)*length(fullfoodlist))
rownames(featfullmat_full)=subjects
colnames(featfullmat_full)=paste0(rep(featlist_intr,each=length(fullfoodlist)),rep(fullfoodlist,times=length(featlist_intr)))
for(subjecti in seq(nsubjects)){
    locatab=summartab[summartab[,"subject"]==subjects[subjecti],]
    locvec=c()
    for(feat in featlist_intr){
        matchind=sapply(fullfoodlist,function(x){
            ind=which(locatab[,"foods"]==x)
            ifelse(length(ind)==0,NA,ind)
        })
        locvec=c(locvec,scale(locatab[matchind,feat]))#scaled within each subject
    }
    featfullmat_full[subjecti,]=locvec
}
neword=order(apply(featfullmat_full,1,function(x) sum(is.na(x))))
frontind=neword[which(rownames(featfullmat_full)[neword]=="SOMEID")]
neword_plot=c(frontind,setdiff(neword,frontind))
# featfullmat=scale(featfullmat)
# mitigators effect upon pure rice
food_list_mitig=c("Rice","Rice+Fat","Rice+Fiber","Rice+Protein")
# mitigator_list=c("Fat","Fiber","Protein")
meanvaltab <- feattab %>% filter(foods %in% food_list_mitig) %>% group_by(subject,foods) %>% summarize(meanval=mean(peak_relative),.groups="keep") %>% as.data.frame()
# threshold=1
mat_mitigator_effect=matrix(NA,nrow=length(subjects),ncol=3)
colnames(mat_mitigator_effect)=food_list_mitig[-1]
rownames(mat_mitigator_effect)=subjects
for(subject_i in seq(nsubjects)){
    subject=subjects[subject_i]
    loctab=meanvaltab[meanvaltab[,"subject"]==subject,]
    if(nrow(loctab)<=1){
        # print(subject)
        next
    }
    ricelevel=loctab[loctab[,"foods"]=="Rice","meanval"]
    # delta_t=ricelevel*threshold
    delta_t=ricelevel
    for(mitigator in loctab[loctab[,"foods"]%in%food_list_mitig[-1],"foods"]){
        delta=loctab[loctab[,"foods"]==mitigator,"meanval"]-ricelevel
        # if(delta>delta_t){
        #     mat_mitigator_effect[subject_i,mitigator]=1
        # }else if(delta < -delta_t){
        #     mat_mitigator_effect[subject_i,mitigator]=-1
        # }
        mat_mitigator_effect[subject,mitigator]=delta/delta_t
    }
}
# mitigators shift upon pure rice
meanvaltab <- feattab %>% filter(foods %in% food_list_mitig) %>% group_by(subject,foods) %>% summarize(meanshift=mean(time_to_peak),.groups="keep") %>% as.data.frame()
mat_mitigator_shift=matrix(NA,nrow=length(subjects),ncol=3)
colnames(mat_mitigator_shift)=food_list_mitig[-1]
rownames(mat_mitigator_shift)=subjects
for(subject_i in seq(nsubjects)){
    subject=subjects[subject_i]
    loctab=meanvaltab[meanvaltab[,"subject"]==subject,]
    if(nrow(loctab)<=1){
        next
    }
    ricelevel=loctab[loctab[,"foods"]=="Rice","meanshift"]
    for(mitigator in loctab[loctab[,"foods"]%in%food_list_mitig[-1],"foods"]){
        delta=loctab[loctab[,"foods"]==mitigator,"meanshift"]-ricelevel
        mat_mitigator_shift[subject,mitigator]=delta/ricelevel
    }
}
# heatmap clustered by mean features
h_3feat_anno=HeatmapAnnotation(feature=rep(featlist_intr,each=length(food_list_carbs)),annotation_name_side="left")
h_3feat=Heatmap(featfullmat,name="3featfood",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=FALSE,row_names_side="left",top_annotation=h_3feat_anno)#
h_worstfood_bar=HeatmapAnnotation(bar2=anno_barplot(colSums(mat_worstfood),height=unit(2,"cm")))
h_worstfood=Heatmap(mat_worstfood,name="worstfood",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE,top_annotation=h_worstfood_bar)#
col_fun_miti=colorRamp2(c(-1,0,1),c("blue","white","red"))
h_mitigator=Heatmap(mat_mitigator_effect,name="mitigator",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE,col=col_fun_miti)#
h_food_quantile=Heatmap(mat_quantile,name="food_quantile",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE)#
pdf(paste0("heatmap_food_all.pdf"))
draw(h_3feat + h_worstfood + h_food_quantile + h_mitigator,main_heatmap="3featfood")
dev.off()
# heatmap original order
h_3feat_anno_full=HeatmapAnnotation(feature=rep(featlist_intr,each=length(fullfoodlist)),annotation_name_side="left")
h_3feat_2=Heatmap(featfullmat_full,name="3featfood2",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,row_names_side="left",top_annotation=h_3feat_anno_full)#
pdf(paste0("heatmap_food_all_oriorder.pdf"))
draw(h_3feat_2 + h_worstfood + h_food_quantile + h_mitigator,main_heatmap="3featfood2",row_order=neword_plot)
dev.off()
# heatmap clustered by worst food
h_3feat2=Heatmap(featfullmat,name="3featfood",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,row_names_side="left",top_annotation=h_3feat_anno)#
h_worstfood2=Heatmap(mat_worstfood,name="worstfood",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE,top_annotation=h_worstfood_bar)#
pdf(paste0("heatmap_food_all2.pdf"))
draw(h_3feat2 + h_worstfood2 + h_food_quantile + h_mitigator,main_heatmap="worstfood")
dev.off()
# heatmap clustered by quantile of standardize food
h_food_quantile2=Heatmap(mat_quantile,name="food_quantile",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE)#
pdf(paste0("heatmap_food_all_sortquantile.pdf"))
draw(h_worstfood + h_food_quantile2 + h_mitigator,main_heatmap="food_quantile")
dev.off()
save(featfullmat,file="featuremean.RData")
# co-plot with whole metadata
metadata=as.data.frame(read.table(paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
cols=colnames(metadata)
coltypes=sapply(metadata,class)
numcols=cols[coltypes=="numeric"|coltypes=="integer"]
addcols=c("t2d_status_ada_bl","a1c_t2d_status_bl","fbg_t2d_status_bl","ogtt_t2d_status_heyjun","sspg_status_heyjun","DI_status_heyjun","IE_status_heyjun","HIR_status_heyjun","HOMA_IR_group_heyjun","ie_3_classes_heyjun","FFA_3classes_heyjun","di_3classes_heyjun","hepatic_ir_3classes_heyjun","sex.factor","ethnicity")
showcols=c(numcols,addcols)
metadatashow=as.data.frame(matrix(NA,nrow=nsubjects,ncol=length(showcols)+1))
rownames(metadatashow)=subjects
colnames(metadatashow)=c("study_id",showcols)
for(subj in subjects){
    match_id=paste0("STUDYID",str_pad(subj,3,pad="0"))
    metadatashow[subj,]=metadata[metadata[,"study_id"]==match_id,c("study_id",showcols)]
}
save(metadatashow,file="meta_idconv.RData")
metadatashow_plot=metadatashow
for(col2fac in addcols){#conert categories to numbers for easy heatmap vis
    metadatashow_plot[,col2fac]=as.numeric(as.factor(metadatashow_plot[,col2fac]))
}
for(col in showcols){
    metadatashow_plot[,col]=scale(metadatashow_plot[,col])[,1]
}
h_worstfood=Heatmap(mat_worstfood,name="worstfood",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE,top_annotation=h_worstfood_bar)#
h_metadata=Heatmap(metadatashow_plot[,-1],name="metadata",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE,column_names_gp=grid::gpar(fontsize=6))#
pdf(paste0("heatmap_metadata.pdf"))
draw(h_worstfood + h_mitigator + h_metadata,main_heatmap="worstfood")
dev.off()
# coplot with metadata selected subset
numcolssub=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","bmi_avg_all","systolic_avg_all","diastolic_avg_all","fructosamine_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","modified_DI_heyjun","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","hepatic_IR_heyjun")
#"alt_sgpt_avg_all","creatinine_avg_all","albumin_avg_all","Insulin_metabolic_testing_mean","C_peptide_metabolic_testing_mean","GLP1_metabolic_testing_mean","GIP_metabolic_testing_mean","Glucagon_pmol_L_metabolic_testing_mean"
chcolssub=c("hp_hypertension","hp_hyperlipidemia","sex.factor","ethnicity")
metadatashow_plot=metadatashow
metadatashow_plot$sspg_avg_all_unscale=metadatashow_plot$sspg_avg_all
metadatashow_plot$modified_DI_heyjun_unscale=metadatashow_plot$modified_DI_heyjun
metadatashow_plot$bmi_avg_all_unscale=metadatashow_plot$bmi_avg_all
for(col in numcolssub){
    metadatashow_plot[,col]=scale(metadatashow_plot[,col])[,1]
}
for(col in chcolssub){
    metadatashow_plot[,col]=as.character(metadatashow_plot[,col])
}
h_worstfood=Heatmap(mat_worstfood,name="worstfood",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE,top_annotation=h_worstfood_bar)#
h_metadata_num=Heatmap(metadatashow_plot[,numcolssub],name="metadata_num",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE)#
colors=structure(1:8,names=c("0","1","Male","Female","White","Asian","Asian, White","Hispanic or Latino")) 
h_metadata_ch=Heatmap(as.matrix(metadatashow_plot[,chcolssub]),name="metadata_ch",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE,col=colors)#
pdf(paste0("heatmap_metadata_subset.pdf"))
draw(h_worstfood + h_mitigator + h_metadata_num + h_metadata_ch,main_heatmap="worstfood")
dev.off()
pdf(paste0("heatmap_metadata_subset_sortquantile.pdf"))
draw(h_worstfood2 + h_food_quantile2 + h_mitigator + h_metadata_num + h_metadata_ch,main_heatmap="food_quantile")
dev.off()
# sorted by mitigator fat effect 
h_mitigator_fat=Heatmap(mat_mitigator_effect,name="mitigator",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_order=order(mat_mitigator_effect[,1]),row_names_side="left")#
pdf(paste0("heatmap_metadata_recordermitigator.pdf"))
draw(h_mitigator_fat + h_metadata_num + h_metadata_ch,main_heatmap="mitigator")
dev.off()
# sorted by mitigator fat shift 
h_mitigator_fat_shift=Heatmap(mat_mitigator_shift,name="mitigator_shift",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_order=order(mat_mitigator_shift[,1]),row_names_side="left")#
pdf(paste0("heatmap_metadata_recordermitigator_shift.pdf"))
draw(h_mitigator_fat_shift + h_metadata_num + h_metadata_ch,main_heatmap="mitigator_shift")
dev.off()
# sorted by mitigator fiber 
h_mitigator_fiber=Heatmap(mat_mitigator_effect,name="mitigator",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_order=order(mat_mitigator_effect[,2]),row_names_side="left")#
pdf(paste0("heatmap_metadata_recordermitigator_fiber.pdf"))
draw(h_mitigator_fiber + h_metadata_num + h_metadata_ch,main_heatmap="mitigator")
dev.off()
# coplot with omics data
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
# average for different samples of the same participant
ncolomics=sum(sapply(omicslist_arra,function(x) dim(x)[1]))
omicsmat=matrix(NA,nrow=nsubjects,ncol=ncolomics)
rownames(omicsmat)=subjects
omics=c("metabolomics","lipidomics","proteomics")
for(subj in subjects){
    tempmat=c()
    for(omic in omics){
        sourmat=omicslist_arra[[omic]]
        ind_data<-colnames(sourmat)%>%str_which(string=.,pattern=paste0("X",subj))
        tempmat=c(tempmat,apply(sourmat[,ind_data,drop=FALSE],1,mean))
    }
    omicsmat[subj,]=tempmat
}
h_worstfood=Heatmap(mat_worstfood,name="worstfood",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE,top_annotation=h_worstfood_bar,width=unit(20,"mm"))#
h_omics=Heatmap(omicsmat,name="omics",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=FALSE)#
pdf(paste0("heatmap_omics.pdf"))
draw(h_worstfood + h_omics,main_heatmap="worstfood")
dev.off()
save(list=ls(all.names=TRUE),file="heatmap_plot.RData")
# barplot and test of metadata between spiking groups
setwd(resdir)
load(paste0(resdir,"worst.food.RData"))
worstfood_group=summartab
plotdf=c()
foodsele=unique(worstfood_group[,"worstfood"])
for(food in foodsele){
    subtab1=worstfood_group[worstfood_group[,"worstfood"]==food,]
    subtab2=metadatashow_plot[as.character(subtab1[,"subject"]),c(numcolssub,"ethnicity")]
    plotdf=rbind(plotdf,cbind(subtab1,subtab2))
}
foodord=c("Pasta","Grapes","Rice","Bread","Potatoes")
plotdf[,"worstfood"]=factor(plotdf[,"worstfood"],level=foodord)
# stat test comparison
comp_feat_list=list("Grapes"=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","sspg_avg_all"),
                    "Potatoes"=c("a1c_avg_all","fbg_avg_all","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","bmi_avg_all","modified_DI_heyjun","hepatic_IR_heyjun"),
                    "Bread"=c("systolic_avg_all","diastolic_avg_all"))
for(center_gr in names(comp_feat_list)){
    compfeat=comp_feat_list[[center_gr]]
    comparisons=list()
    for(food in setdiff(foodsele,center_gr)){
        comparisons[[length(comparisons)+1]]=c(food,center_gr)
    }
    for(colfeat in compfeat){
        p<-ggboxplot(plotdf,x="worstfood",y=colfeat,add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
        ggsave(plot=p,paste0(resdir,"boxplot_spiketype_meta",colfeat,"_compare_",center_gr,".pdf"))
    }
}
# box plot without scale
plotdf=c()
foodsele=unique(worstfood_group[,"worstfood"])
for(food in foodsele){
    subtab1=worstfood_group[worstfood_group[,"worstfood"]==food,]
    subtab2=metadatashow[as.character(subtab1[,"subject"]),c(numcolssub,"ethnicity")]
    plotdf=rbind(plotdf,cbind(subtab1,subtab2))
}
foodord=c("Pasta","Grapes","Rice","Bread","Potatoes")
plotdf[,"worstfood"]=factor(plotdf[,"worstfood"],level=foodord)
# stat test comparison
comp_feat_list=list("Grapes"=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","sspg_avg_all"),
                    "Potatoes"=c("a1c_avg_all","fbg_avg_all","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","bmi_avg_all","modified_DI_heyjun","hepatic_IR_heyjun"),
                    "Bread"=c("systolic_avg_all","diastolic_avg_all"))
for(center_gr in names(comp_feat_list)){
    compfeat=comp_feat_list[[center_gr]]
    comparisons=list()
    for(food in setdiff(foodsele,center_gr)){
        comparisons[[length(comparisons)+1]]=c(food,center_gr)
    }
    for(colfeat in compfeat){
        p<-ggboxplot(plotdf,x="worstfood",y=colfeat,add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
        ggsave(plot=p,paste0(resdir,"boxplot_spiketype_meta",colfeat,"_compare_",center_gr,"unscale.pdf"))
    }
}
# test whether is asian vs rice spiker
testtab=plotdf[,c("worstfood","ethnicity")]
testtab$ricespiker=testtab[,"worstfood"]=="Rice"
testtab$asian=testtab[,"ethnicity"]=="Asian"
print(fisher.test(x=testtab$ricespiker,y=testtab$asian))
# correlation matrix of all spiking patterns
corr_mat_quantile=cor(mat_quantile,method="spearman",use="pairwise.complete.obs")
col_fun=colorRamp2(c(0,1),c("white","red"))
h_food_quantile_corr=Heatmap(corr_mat_quantile,name="food_quantile_corr",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,col=col_fun)#
pdf(paste0("heatmap_food_all_sortquantile_corr.pdf"))
draw(h_food_quantile_corr,main_heatmap="food_quantile_corr")
dev.off()
# correlation matrix of spiking and mitigator patterns
foodmat=merge(mat_quantile,mat_mitigator_effect,by="row.names")
corr_mat_quantile=cor(foodmat[,seq(from=2,to=dim(foodmat)[2])],method="spearman",use="pairwise.complete.obs")
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
h_food_quantile_corr=Heatmap(corr_mat_quantile,name="food_quantile_corr",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,col=col_fun)#
pdf(paste0("heatmap_corr_food_mitigator.pdf"))
draw(h_food_quantile_corr,main_heatmap="food_quantile_corr")
dev.off()
# spike comparison between metabolic subtypes
sumtab_quan=as.data.frame(sumtab_quan)
mat_peak=matrix(NA,nrow=length(subjects),ncol=length(food_list_carbs))
rownames(mat_peak)=subjects
colnames(mat_peak)=food_list_carbs
for(rowi in seq(nrow(sumtab_quan))){
    mat_peak[paste0(sumtab_quan[rowi,"subject"]),sumtab_quan[rowi,"foods"]]=sumtab_quan[rowi,"meanval"]
}
# add new def of IR and IE
sspg_status_new=rep(NA,times=dim(metadatashow_plot)[1])
sspg_status_new[metadatashow_plot$sspg_avg_all_unscale>120]="IR"
sspg_status_new[metadatashow_plot$sspg_avg_all_unscale<100]="IS"
metadatashow_plot$sspg_status_new=sspg_status_new
sspg_status_new2=rep(NA,times=dim(metadatashow_plot)[1])
sspg_status_new2[metadatashow_plot$sspg_avg_all_unscale>120]="IR"
sspg_status_new2[metadatashow_plot$sspg_avg_all_unscale<=120]="IS"
metadatashow_plot$sspg_status_new2=sspg_status_new2
IE_status_new=rep(NA,times=dim(metadatashow_plot)[1])
ievec=metadatashow_plot$modified_DI_heyjun_unscale
IE_status_new[ievec<quantile(ievec,0.4,na.rm=TRUE)]="Dysfunctional"
IE_status_new[ievec>quantile(ievec,0.6,na.rm=TRUE)]="Normal"
metadatashow_plot$IE_status_new=IE_status_new
a1c_t2d_status_bl_2=metadatashow_plot$a1c_t2d_status_bl
a1c_t2d_status_bl_2[a1c_t2d_status_bl_2=="T2D"]="preDM"
metadatashow_plot$a1c_t2d_status_bl_2=a1c_t2d_status_bl_2
# mixed group
t2d_ir_mix=rep(NA,times=dim(metadatashow_plot)[1])
for(typ1 in c("Normal","preDM")){
    for(typ2 in c("IS","IR")){
        ind=which(metadatashow_plot$a1c_t2d_status_bl_2==typ1&metadatashow_plot$sspg_status_heyjun==typ2)
        namstr=paste0(typ1,"_",typ2)
        t2d_ir_mix[ind]=namstr
    }
}
metadatashow_plot$t2d_ir_mix=t2d_ir_mix
# 
metadatashow_plot$bmi_group=ifelse(metadatashow_plot$bmi_avg_all_unscale>=25,"Overweight","Normal")
# 
subtyvec=c("t2d_status_ada_bl","a1c_t2d_status_bl","a1c_t2d_status_bl_2","sspg_status_heyjun","DI_status_heyjun","IE_status_heyjun","HOMA_IR_group_heyjun","FFA_3classes_heyjun","hepatic_ir_3classes_heyjun","sspg_status_new","IE_status_new","sspg_status_new2","t2d_ir_mix","sex.factor","ethnicity","bmi_group")
colnames(mat_mitigator_effect)=c("Fat","Fiber","Protein")
plotdf=cbind(metadatashow_plot[subjects,c(subtyvec,"modified_DI_heyjun","sspg_avg_all")],mat_peak)
plotdf=cbind(plotdf,mat_mitigator_effect)
foodslistcheck=c(colnames(mat_peak),colnames(mat_mitigator_effect))
orderlist=list("t2d_status_ada_bl"=c("Normal","preDM","T2D"),
               "a1c_t2d_status_bl"=c("Normal","preDM","T2D"),
               "a1c_t2d_status_bl_2"=c("Normal","preDM"),
               "sspg_status_heyjun"=c("IS","IR"),
               "DI_status_heyjun"=c("BC_normal","BC_inter","BC_dys"),
               "IE_status_heyjun"=c("IE_normal","IE_inter","IE_dys"),
               "HOMA_IR_group_heyjun"=c("IS","IR"),
               "FFA_3classes_heyjun"=c("Adipocyte_IS","Adipocyte_Intermediate","Adipocyte_IR"),
               "hepatic_ir_3classes_heyjun"=c("Hepatic_IS","Hepatic_Intermediate","Hepatic_IR"),
               "sspg_status_new"=c("IS","IR"),
               "IE_status_new"=c("Normal","Dysfunctional"),
               "sspg_status_new2"=c("IS","IR"),
               "t2d_ir_mix"=c("Normal_IS","Normal_IR","preDM_IS","preDM_IR"),
               "sex.factor"=c("Female","Male"),
               "ethnicity"=c("White","Asian"),
               "bmi_group"=c("Normal","Overweight"))
for(subtcol in subtyvec){
    plotdf2=plotdf[!is.na(plotdf[,subtcol]),]
    sublistv=plotdf2[,subtcol]
    sublistuniq=unique(sublistv)
    sublistuniq=intersect(sublistuniq,orderlist[[subtcol]])
    comparisons=combn(sublistuniq,2,simplify=FALSE)
    # reorder the plot dataframe
    plotdf2[,subtcol]=factor(plotdf2[,subtcol],levels=orderlist[[subtcol]])
    for(food in foodslistcheck){
        p<-ggboxplot(plotdf2,x=subtcol,y=food,add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
        ggsave(plot=p,paste0(resdir,"boxplot_metabotype_meta",food,"_compare_",subtcol,".pdf"))
    }
}
# potato vs grapes
plotdf$potatovsgrapes=plotdf$Potatoes/plotdf$Grapes
plotdf2=plotdf[!is.na(plotdf[,"sspg_status_new"]),]
sublistv=plotdf2[,"sspg_status_new"]
sublistuniq=unique(sublistv)
comparisons=combn(sublistuniq,2,simplify=FALSE)
# plot(plotdf[ind,"sspg_avg_all_unscale"],plotdf[ind,"Potatoes"])
p<-ggboxplot(plotdf2,x="sspg_status_new",y="potatovsgrapes",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
ggsave(plot=p,paste0(resdir,"boxplot_metabotype_meta_potatogrape_compare_sspg_all.pdf"))
# 
plotdf2=plotdf[!is.na(plotdf[,"sspg_status_new2"]),]
sublistv=plotdf2[,"sspg_status_new2"]
sublistuniq=unique(sublistv)
comparisons=combn(sublistuniq,2,simplify=FALSE)
# plot(plotdf[ind,"sspg_avg_all_unscale"],plotdf[ind,"Potatoes"])
p<-ggboxplot(plotdf2,x="sspg_status_new2",y="potatovsgrapes",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
ggsave(plot=p,paste0(resdir,"boxplot_metabotype_meta_potatogrape_compare_sspg_all_old_thres.pdf"))
# 
plotdf2=plotdf[!is.na(plotdf[,"sspg_status_heyjun"]),]
sublistv=plotdf2[,"sspg_status_heyjun"]
sublistuniq=unique(sublistv)
comparisons=combn(sublistuniq,2,simplify=FALSE)
# plot(plotdf[ind,"sspg_avg_all_unscale"],plotdf[ind,"Potatoes"])
p<-ggboxplot(plotdf2,x="sspg_status_heyjun",y="potatovsgrapes",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
ggsave(plot=p,paste0(resdir,"boxplot_metabotype_meta_potatogrape_compare_sspg_old.pdf"))
# t.test between groups of metabolic functions in food spikes
testtab=c()
subtyvec=c("a1c_t2d_status_bl","sspg_status_heyjun","DI_status_heyjun","IE_status_heyjun","HOMA_IR_group_heyjun","FFA_3classes_heyjun","hepatic_ir_3classes_heyjun")
for(subtcol in subtyvec){
    plotdf2=plotdf[!is.na(plotdf[,subtcol])&plotdf[,subtcol]!="Unknown",]
    sublistv=plotdf2[,subtcol]
    sublistuniq=unique(sublistv)
    comparisons=combn(sublistuniq,2,simplify=FALSE)
    # reorder the plot dataframe
    for(food in foodslistcheck){
        for(compi in seq(length(comparisons))){
            comppair=comparisons[[compi]]
            loctab=plotdf2[plotdf2[,subtcol]%in%comppair,c(food,subtcol)]
            colnames(loctab)=c("value","group")
            loctab=loctab[!is.na(loctab[,"value"]),]
            templist=list()
            if(min(table(loctab[,"group"]))>2){
                teststat=t.test(value~group,data=loctab)
                templist[["pval"]]=teststat$p.value
                templist[["statistics"]]=teststat$statistic
                templist[["ci"]]=paste0(formatC(teststat$conf.int,format="E",digits=2),collapse=" ")
                templist[["Df"]]=teststat$parameter
                templist[["effectsz"]]=cohensD(value~group,data=loctab)
            }else{
                templist[["pval"]]=NA
                templist[["statistics"]]=NA
                templist[["ci"]]=NA
                templist[["Df"]]=NA
                templist[["effectsz"]]=NA
            }
            gp1ind=loctab[,"group"]==comppair[1]
            delta=mean(loctab[!gp1ind,"value"],na.rm=TRUE)-mean(loctab[gp1ind,"value"],na.rm=TRUE)
            testtab=rbind(testtab,data.frame(food=food,category=subtcol,group1=comppair[1],group2=comppair[2],delta=delta,pval=templist[["pval"]],statistics=templist[["statistics"]],ci=templist[["ci"]],Df=templist[["Df"]],effectsz=templist[["effectsz"]]))
        }
    }
}
testtab_adj=c()
for(subtcol in subtyvec){
    subtab=testtab[testtab[,"category"]==subtcol,]
    subtab$padj=p.adjust(subtab$pval,method="fdr")
    testtab_adj=rbind(testtab_adj,subtab)
}
save(testtab_adj,file="stat_gp_test_food_subtype.RData")
write.table(testtab_adj,file="stat_gp_test_food_subtype.tsv",row.names=FALSE)
testtab_adj_sele=testtab_adj[testtab_adj[,"category"]%in%c("sspg_status_heyjun","DI_status_heyjun"),]
write.table(testtab_adj_sele,file="stat_gp_test_food_subtype_selec.tsv",row.names=FALSE)
# check food spike with both sspg and DI in one model
statcoll=c()
 for(food in foodslistcheck){
    formula_mlm_str=paste0(food," ~ modified_DI_heyjun+sspg_avg_all")
    formula_mlm=as.formula(formula_mlm_str)
    lmdol=lm(formula_mlm,plotdf)
    coeftab=summary(lmdol)$coefficients
    vars=all.vars(formula_mlm)
    statcoll=rbind(statcoll,data.frame(Y=rep(vars[1],times=2),X=vars[2:3],coef=coeftab[2:3,"Estimate"],pval=coeftab[2:3,"Pr(>|t|)"],formula=rep(formula_mlm_str,times=2)))
}
statcoll$padj=p.adjust(statcoll[,"pval"],method="fdr")
save(statcoll,file="sspg_and_di_foodspike.RData")
# paired test for mitigators for each groups
subtyvec=c("sspg_status_heyjun","DI_status_heyjun")
feat_to_compares=c("peak_relative","AUC_above_baseline","time_to_peak")
stat_coll_pair=c()
mitig_comb=c("Rice+Fiber","Rice+Fat","Rice+Protein")
for(subtcol in subtyvec){
    plotdf2=plotdf[!is.na(plotdf[,subtcol])&plotdf[,subtcol]!="Unknown",]
    sublistv=plotdf2[,subtcol]
    sublistuniq=unique(sublistv)
    ggplottablong=c()
    for(subg in sublistuniq){
        selecsub=rownames(plotdf2[plotdf2[,subtcol]==subg,])
        feattabsub=feattab[feattab[,"subject"]%in%selecsub,]
        for(feature in feat_to_compares){
            for(foodthe in mitig_comb){
                xtab <- feattabsub %>% filter(foods %in% c("Rice")) %>% group_by(subject) %>% summarize(meanval=mean(.data[[feature]]),.groups="keep") %>% group_by(subject) %>% as.data.frame()
                ytab <- feattabsub %>% filter(foods %in% c(foodthe)) %>% group_by(subject) %>% summarize(meanval=mean(.data[[feature]]),.groups="keep") %>% group_by(subject) %>% as.data.frame()
                subjs=ytab[,"subject"]
                xorder=match(subjs,xtab[,"subject"])
                teststat=t.test(xtab[xorder,"meanval"],ytab[,"meanval"],paired=TRUE)
                stat_coll_pair=rbind(stat_coll_pair,data.frame(category=subtcol,group=subg,foods=foodthe,feature=feature,pval=teststat$p.value,statistics=teststat$statistic,ci=paste0(formatC(teststat$conf.int,format="E",digits=2),collapse=" "),Df=teststat$parameter,effectsz=cohensD(xtab[,"meanval"],ytab[,"meanval"]),delta=mean(ytab[,"meanval"])-mean(xtab[xorder,"meanval"])))
            }
            temptab=feattabsub[feattabsub[,"foods"]%in%c("Rice",mitig_comb),c("subject",feature,"foods")]
            colnames(temptab)[2]="value"
            temptab$cgmfeature=feature
            temptab$subtype=subg
            ggplottablong=rbind(ggplottablong,temptab)
        }
    }
    # bar plot 
    for(feature in feat_to_compares){
        loctab=ggplottablong[ggplottablong[,"cgmfeature"]==feature,]
        loctab$foods=factor(loctab$foods,levels=c("Rice","Rice+Fiber","Rice+Protein","Rice+Fat"))
        cdata<-ddply(loctab,c("cgmfeature","subtype","foods"),summarise,
                N=length(value),
                mean=mean(value),
                sd=sd(value),
                se=sd/sqrt(N)
        )
        p<-ggplot(cdata,aes(x=foods,y=mean,fill=subtype))+geom_bar(stat="identity",position="dodge",width=0.2)+geom_errorbar(aes(ymin=mean-2*se,ymax=mean+2*se),position=position_dodge(.9),width=0.2)+geom_point(data=loctab,aes(x=foods,y=value),position=position_jitter(seed=1,width=0.1))+facet_wrap(~subtype,ncol=1,scales="free")+ylab(feature)+theme_bw()+theme(axis.line=element_line(colour="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background=element_blank()) 
        ggsave(paste0(subtcol,feature,"mitigation.pdf"),plot=p)
    }
}
stat_coll_pair=stat_coll_pair[stat_coll_pair[,"group"]!="BC_inter",]
stat_coll_pair_adj=c()
for(subtcol in subtyvec){
    plotdf2=plotdf[!is.na(plotdf[,subtcol])&plotdf[,subtcol]!="Unknown",]
    sublistv=plotdf2[,subtcol]
    sublistuniq=unique(sublistv)
    for(subg in sublistuniq){
        for(foodthe in mitig_comb){
            subtab=stat_coll_pair[stat_coll_pair[,"foods"]==foodthe&stat_coll_pair[,"category"]==subtcol&stat_coll_pair[,"group"]==subg,]
            subtab$padj=p.adjust(subtab$pval,method="fdr")
            stat_coll_pair_adj=rbind(stat_coll_pair_adj,subtab)
        }
    }
}
write.table(stat_coll_pair_adj,file="foodvs_ttest.mitig_pair_subt.txt",row.names=FALSE)