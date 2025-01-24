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
require(ggpubr)
require(circlize)
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
# 
feattab=read.table("cgm_foods_manual_features.csv",header=TRUE,sep=",")
feattab=feattab[,c("subject","foods","rep","mitigator","peak_value","baseline_glucose","AUC_above_baseline","time_to_peak","glucose_at_120_mins")]
checkfeats=c("AUC_above_baseline","glucose_at_120_mins","time_to_peak","peak_abs")
for(chfeat in checkfeats){
    if(chfeat=="AUC_above_baseline"){
        feattab$quan_relative=feattab$AUC_above_baseline
    }else if(chfeat=="time_to_peak"){
        feattab$quan_relative=feattab$time_to_peak
    }else if(chfeat=="peak_abs"){
        feattab$quan_relative=feattab$peak_value
    }else{
        feattab$quan_relative=feattab$glucose_at_120_mins-feattab$baseline_glucose
    }
    # foods spikers (worst food)
    food_list_carbs=c("Beans","Berries","Bread","Grapes","Pasta","Potatoes","Rice")
    summartab <- feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(meanval=mean(quan_relative),.groups="keep") %>% group_by(subject) %>% summarize(worstfood=foods[which.max(meanval)],.groups="keep") %>% as.data.frame()
    mat_worstfood=matrix(0,nrow=nrow(summartab),ncol=length(food_list_carbs))
    rownames(mat_worstfood)=paste0(summartab[,1])
    colnames(mat_worstfood)=food_list_carbs
    for(rowi in seq(dim(summartab)[1])){
        mat_worstfood[paste0(summartab[rowi,1]),summartab[rowi,2]]=1
    }
    subjects=rownames(mat_worstfood)
    nsubjects=length(subjects)
    save(summartab,file=paste0("worst.food.",chfeat,".RData"))
    # quantile matrix
    feattab$peakquantile=percent_rank(feattab$quan_relative)
    feattab_avg_q <- feattab%>%filter(foods %in% food_list_carbs)%>%group_by(subject,foods)%>%summarise(avg_peak_q=mean(peakquantile),.groups="keep")%>%as.data.frame()
    mat_quantile=matrix(NA,nrow=nsubjects,ncol=length(food_list_carbs))
    rownames(mat_quantile)=subjects
    colnames(mat_quantile)=food_list_carbs
    for(rowi in seq(dim(feattab_avg_q)[1])){
        mat_quantile[paste0(feattab_avg_q[rowi,1]),feattab_avg_q[rowi,2]]=feattab_avg_q[rowi,3]
    }
    save(feattab_avg_q,mat_quantile,file=paste0("worst.food.quantile.",chfeat,".RData"))
    # matrix of three features mean [quan_relative, AUC, time_to_peak]
    summartab <- feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(peakmean=mean(quan_relative),aucmean=mean(AUC_above_baseline),time2peakmean=mean(time_to_peak),.groups="keep") %>% group_by(subject) %>% as.data.frame()
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
    # matrix of three features mean [quan_relative, AUC, time_to_peak] with mitigators
    fullfoodlist=c(food_list_carbs,"Rice+Fat","Rice+Fiber","Rice+Protein")
    summartab <- feattab %>% filter(foods %in% fullfoodlist) %>% group_by(subject,foods) %>% summarize(peakmean=mean(quan_relative),aucmean=mean(AUC_above_baseline),time2peakmean=mean(time_to_peak),.groups="keep") %>% group_by(subject) %>% as.data.frame()
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
    frontind=neword[which(rownames(featfullmat_full)[neword]=="56")]
    neword_plot=c(frontind,setdiff(neword,frontind))
    # featfullmat=scale(featfullmat)
    # mitigators effect upon pure rice
    food_list_mitig=c("Rice","Rice+Fat","Rice+Fiber","Rice+Protein")
    # mitigator_list=c("Fat","Fiber","Protein")
    meanvaltab <- feattab %>% filter(foods %in% food_list_mitig) %>% group_by(subject,foods) %>% summarize(meanval=mean(quan_relative),.groups="keep") %>% as.data.frame()
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
    h_mitigator=Heatmap(mat_mitigator_effect,name="mitigator",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE)#
    h_food_quantile=Heatmap(mat_quantile,name="food_quantile",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE)#
    pdf(paste0("heatmap_food_all.",chfeat,".pdf"))
    draw(h_3feat + h_worstfood + h_food_quantile + h_mitigator,main_heatmap="3featfood")
    dev.off()
    datavec=colSums(mat_worstfood)
    savedf=data.frame(food=names(datavec),count=datavec)
    write.table(savedf,file=paste0("fig2c_",chfeat,".txt"),row.names=FALSE)
    # heatmap original order
    h_3feat_anno_full=HeatmapAnnotation(feature=rep(featlist_intr,each=length(fullfoodlist)),annotation_name_side="left")
    h_3feat_2=Heatmap(featfullmat_full,name="3featfood2",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,row_names_side="left",top_annotation=h_3feat_anno_full)#
    pdf(paste0("heatmap_food_all_oriorder.",chfeat,".pdf"))
    draw(h_3feat_2 + h_worstfood + h_food_quantile + h_mitigator,main_heatmap="3featfood2",row_order=neword_plot)
    dev.off()
    # heatmap clustered by worst food
    h_3feat2=Heatmap(featfullmat,name="3featfood",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,row_names_side="left",top_annotation=h_3feat_anno)#
    h_worstfood2=Heatmap(mat_worstfood,name="worstfood",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE,top_annotation=h_worstfood_bar)#
    pdf(paste0("heatmap_food_all2.",chfeat,".pdf"))
    draw(h_3feat2 + h_worstfood2 + h_food_quantile + h_mitigator,main_heatmap="worstfood")
    dev.off()
    # heatmap clustered by quantile of standardize food
    h_food_quantile2=Heatmap(mat_quantile,name="food_quantile",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=FALSE,show_column_names=TRUE)#
    pdf(paste0("heatmap_food_all_sortquantile.",chfeat,".pdf"))
    draw(h_worstfood + h_food_quantile2 + h_mitigator,main_heatmap="food_quantile")
    dev.off()
    save(featfullmat,file=paste0("featuremean.",chfeat,".RData"))
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
        match_id=paste0("STUDYID-",str_pad(subj,3,pad="0"))
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
    pdf(paste0("heatmap_metadata.",chfeat,".pdf"))
    draw(h_worstfood + h_mitigator + h_metadata,main_heatmap="worstfood")
    dev.off()
    # coplot with metadata selected subset
    numcolssub=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","bmi_avg_all","systolic_avg_all","diastolic_avg_all","fructosamine_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","modified_DI_heyjun","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","hepatic_IR_heyjun")
    #"alt_sgpt_avg_all","creatinine_avg_all","albumin_avg_all"
    # "Insulin_metabolic_testing_mean","C_peptide_metabolic_testing_mean","GLP1_metabolic_testing_mean","GIP_metabolic_testing_mean","Glucagon_pmol_L_metabolic_testing_mean"
    chcolssub=c("hp_hypertension","hp_hyperlipidemia","sex.factor","ethnicity")
    metadatashow_plot=metadatashow
    metadatashow_plot$sspg_avg_all_unscale=metadatashow_plot$sspg_avg_all
    metadatashow_plot$modified_DI_heyjun_unscale=metadatashow_plot$modified_DI_heyjun
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
    pdf(paste0("heatmap_metadata_subset.",chfeat,".pdf"))
    draw(h_worstfood + h_mitigator + h_metadata_num + h_metadata_ch,main_heatmap="worstfood")
    dev.off()
    pdf(paste0("heatmap_metadata_subset_sortquantile.",chfeat,".pdf"))
    draw(h_worstfood2 + h_food_quantile2 + h_mitigator + h_metadata_num + h_metadata_ch,main_heatmap="food_quantile")
    dev.off()
    # sorted by mitigator fat effect 
    h_mitigator_fat=Heatmap(mat_mitigator_effect,name="mitigator",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_order=order(mat_mitigator_effect[,1]),row_names_side="left")#
    pdf(paste0("heatmap_metadata_recordermitigator.",chfeat,".pdf"))
    draw(h_mitigator_fat + h_metadata_num + h_metadata_ch,main_heatmap="mitigator")
    dev.off()
    # sorted by mitigator fat shift 
    h_mitigator_fat_shift=Heatmap(mat_mitigator_shift,name="mitigator_shift",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_order=order(mat_mitigator_shift[,1]),row_names_side="left")#
    pdf(paste0("heatmap_metadata_recordermitigator_shift.",chfeat,".pdf"))
    draw(h_mitigator_fat_shift + h_metadata_num + h_metadata_ch,main_heatmap="mitigator_shift")
    dev.off()
    # sorted by mitigator fiber 
    h_mitigator_fiber=Heatmap(mat_mitigator_effect,name="mitigator",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_order=order(mat_mitigator_effect[,2]),row_names_side="left")#
    pdf(paste0("heatmap_metadata_recordermitigator_fiber.",chfeat,".pdf"))
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
    pdf(paste0("heatmap_omics.",chfeat,".pdf"))
    draw(h_worstfood + h_omics,main_heatmap="worstfood")
    dev.off()
    save(list=ls(all.names=TRUE),file=paste0("heatmap_plot.",chfeat,".RData"))
    # barplot and test of metadata between spiking groups
    setwd(resdir)
    load(paste0(resdir,"worst.food.",chfeat,".RData"))
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
    plotdf=plotdf[!is.na(plotdf[,"worstfood"]),]
    # stat test comparison
    comp_feat_list=list("Grapes"=c("a1c_avg_all","fbg_avg_all","insulin_fasting_avg_all","cholesterol_total_avg_all","ldl_avg_all","hdl_avg_all","sspg_avg_all"),
                        "Potatoes"=c("a1c_avg_all","fbg_avg_all","ffa_avg_heyjun","sspg_avg_all","ie_heyjun","bmi_avg_all"),
                        "Bread"=c("systolic_avg_all","diastolic_avg_all"))
    for(center_gr in names(comp_feat_list)){
        compfeat=comp_feat_list[[center_gr]]
        comparisons=list()
        for(food in setdiff(foodord,center_gr)){
            comparisons[[length(comparisons)+1]]=c(food,center_gr)
        }
        for(colfeat in compfeat){
            p<-ggboxplot(plotdf,x="worstfood",y=colfeat,add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
            ggsave(plot=p,paste0(resdir,"boxplot_spiketype_meta",colfeat,"_compare_",center_gr,chfeat,".pdf"))
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
    pdf(paste0("heatmap_food_all_sortquantile_corr.",chfeat,".pdf"))
    draw(h_food_quantile_corr,main_heatmap="food_quantile_corr")
    dev.off()
    # correlation matrix of spiking and mitigator patterns
    foodmat=merge(mat_quantile,mat_mitigator_effect,by="row.names")
    corr_mat_quantile=cor(foodmat[,seq(from=2,to=dim(foodmat)[2])],method="spearman",use="pairwise.complete.obs")
    col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
    h_food_quantile_corr=Heatmap(corr_mat_quantile,name="food_quantile_corr",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=FALSE,col=col_fun)#
    pdf(paste0("heatmap_corr_food_mitigator.",chfeat,".pdf"))
    draw(h_food_quantile_corr,main_heatmap="food_quantile_corr")
    dev.off()
    # quantification comparison between metabolic subtypes
    sumtab_quan<-feattab %>% filter(foods %in% food_list_carbs) %>% group_by(subject,foods) %>% summarize(meanval=mean(quan_relative),.groups="keep") %>% group_by(subject)
    sumtab_quan=as.data.frame(sumtab_quan)
    mat_peak=matrix(NA,nrow=length(subjects),ncol=length(food_list_carbs))
    rownames(mat_peak)=subjects
    colnames(mat_peak)=food_list_carbs
    for(rowi in seq(nrow(sumtab_quan))){
        mat_peak[paste0(sumtab_quan[rowi,"subject"]),sumtab_quan[rowi,"foods"]]=sumtab_quan[rowi,"meanval"]
    }
    sspg_status_new=rep(NA,times=dim(metadatashow_plot)[1])
    sspg_status_new[metadatashow_plot$sspg_avg_all_unscale>120]="IR"
    sspg_status_new[metadatashow_plot$sspg_avg_all_unscale<100]="IS"
    metadatashow_plot$sspg_status_new=sspg_status_new
    IE_status_new=rep(NA,times=dim(metadatashow_plot)[1])
    ievec=metadatashow_plot$modified_DI_heyjun_unscale
    IE_status_new[ievec<quantile(ievec,0.4,na.rm=TRUE)]="Dysfunctional"
    IE_status_new[ievec>quantile(ievec,0.6,na.rm=TRUE)]="Normal"
    metadatashow_plot$IE_status_new=IE_status_new
    subtyvec=c("t2d_status_ada_bl","a1c_t2d_status_bl","sspg_status_heyjun","DI_status_heyjun","IE_status_heyjun","HOMA_IR_group_heyjun","FFA_3classes_heyjun","hepatic_ir_3classes_heyjun","sspg_status_new","IE_status_new")
    colnames(mat_mitigator_effect)=c("Fat","Fiber","Protein")
    plotdf=cbind(metadatashow_plot[subjects,c(subtyvec,"modified_DI_heyjun","sspg_avg_all")],mat_peak)
    plotdf=cbind(plotdf,mat_mitigator_effect)
    foodslistcheck=c(colnames(mat_peak),colnames(mat_mitigator_effect))
    orderlist=list("t2d_status_ada_bl"=c("Normal","preDM","T2D"),
                "a1c_t2d_status_bl"=c("Normal","preDM","T2D"),
                "sspg_status_heyjun"=c("IS","IR"),
                "DI_status_heyjun"=c("BC_normal","BC_inter","BC_dys"),
                "IE_status_heyjun"=c("IE_normal","IE_inter","IE_dys"),
                "HOMA_IR_group_heyjun"=c("IS","IR"),
                "FFA_3classes_heyjun"=c("Adipocyte_IS","Adipocyte_Intermediate","Adipocyte_IR"),
                "hepatic_ir_3classes_heyjun"=c("Hepatic_IS","Hepatic_Intermediate","Hepatic_IR"),
                "sspg_status_new"=c("IS","IR"),
                "IE_status_new"=c("Normal","Dysfunctional"))
    for(subtcol in subtyvec){
        plotdf2=plotdf[!is.na(plotdf[,subtcol]),]
        sublistv=plotdf2[,subtcol]
        sublistuniq=unique(sublistv)
        comparisons=combn(sublistuniq,2,simplify=FALSE)
        # reorder the plot dataframe
        plotdf2[,subtcol]=factor(plotdf2[,subtcol],levels=orderlist[[subtcol]])
        for(food in foodslistcheck){
            p<-ggboxplot(plotdf2,x=subtcol,y=food,add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
            ggsave(plot=p,paste0(resdir,"boxplot_metabotype_meta",food,"_compare_",subtcol,".",chfeat,".pdf"))
        }
    }
}
