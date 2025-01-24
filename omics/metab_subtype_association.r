# annova test between subtype for each omics features
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(igraph)
require(Hmisc)
require(npreg)
require(ppcor)
require(readxl)
# 
require(clusterProfiler)
require(org.Hs.eg.db)
require(ReactomePA)
require(metpath)
require(ComplexHeatmap)
require(pathview)
require(circlize)
require(KEGG.db)
# 
require(foreach)
require(doSNOW)
require(Rodin)
cl<-makeSOCKcluster(8)
registerDoSNOW(cl)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/pathway/")
setwd(resdir)
source("/Users/yuewu/Documents/GitHub/meal_cgm/omics/func_lipid_enrich.r",local=TRUE)
# omics data 
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
load("omics_data_clean.RData")
subjects=rownames(omicsmat)
nsubjects=length(subjects)
# metadata
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
class_num_cols=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","modified_DI_heyjun","ie_heyjun","hepatic_IR_heyjun","HOMA_IR_heyjun","ffa_avg_heyjun","HOMA_S_heyjun","HOMA_B_heyjun","albumin_avg_all","creatinine_avg_all","hscrp_avg_all","ldl_avg_all","hdl_avg_all","triglyceride_avg_all","fructosamine_avg_all","insulin_fasting_avg_all","alt_sgpt_avg_all","systolic_avg_all","diastolic_avg_all","bmi_avg_all")
subphenogrp=class_num_cols[1:9]
metadatashow=as.data.frame(matrix(NA,nrow=nsubjects,ncol=length(class_num_cols)+1))
rownames(metadatashow)=subjects
colnames(metadatashow)=c("study_id",class_num_cols)
for(subj in subjects){
    match_id=paste0("STUDYID-",str_pad(subj,3,pad="0"))
    metadatashow[subj,]=metadata[metadata[,"study_id"]==match_id,c("study_id",class_num_cols)]
}
metadatashow=metadatashow[,-1]
# select significant features
nclasses=length(class_num_cols)
nomics=dim(omicsmat)[2]
cor_mat=matrix(NA,nrow=nomics,ncol=nclasses)
colnames(cor_mat)=class_num_cols
rownames(cor_mat)=featvec
p_mat=cor_mat
list_sig=vector(mode="list",length=nclasses*2)
names(list_sig)=paste0(rep(c("posicorr","negacorr"),each=nclasses),rep(class_num_cols,times=2))
for(omicsi in seq(nomics)){
    for(classi in seq(nclasses)){
        xvec=metadatashow[,classi]
        yvec=omicsmat[,omicsi]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        cor_mat[omicsi,classi]=corres$estimate
        p_mat[omicsi,classi]=corres$p.value
        if(corres$p.value<0.05){
            compdname=paste0(omicsvec[omicsi],featvec[omicsi])
            dirc=ifelse(corres$estimate>0,"posicorr","negacorr")
            listname=paste0(dirc,class_num_cols[classi])
            list_sig[[listname]]=c(list_sig[[listname]],compdname)
        }
    }
}
# significant table 
tabstat_cor=c()
omicslist=unique(omicsvec)
for(classi in seq(nclasses)){
    for(omicsele in omicslist){
        xvec=metadatashow[,classi]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        omicsind=which(omicsvec==omicsele)
        corvec=rep(NA,times=length(omicsind))
        features=featvec[omicsind]
        names(corvec)=featvec[omicsind]
        pval=corvec
        omicsmatloc=omicsmat[nonaind,omicsind]
        ind_dochange=apply(omicsmatloc,2,function(x) var(x,na.rm=TRUE))!=0
        for(coli in which(ind_dochange)){
            corres=cor.test(x=xvec[nonaind],y=omicsmatloc[,coli],method="spearman")
            pval[coli]=corres$p.value
            corvec[coli]=corres$estimate
        }
        p.adj=p.adjust(pval,"fdr")
        seleind=which(pval<0.05)
        seleind=seleind[!is.na(seleind)]
        len=length(seleind)
        tabstat_cor=rbind(tabstat_cor,data.frame(type=rep(omicsele,times=len),clinic=rep(class_num_cols[classi],times=len),feature=features[seleind],correlation=corvec[seleind],pvalue=pval[seleind],padj=p.adj[seleind]))
    }
}
save(tabstat_cor,file="association_subtype_omics.RData")
# metabolomics reformat kegg and hmdb ids
loctab=annolist[["metabolomics"]]
matchlist=list("HMDB.ID"="HMDB_New_DB","KEGG.ID"="KEGG_New_DB")
for(idcol in names(matchlist)){
    namask=is.na(loctab[,idcol])
    loctab[namask,idcol]=loctab[namask,matchlist[[idcol]]]
}
## correct the mismatched ids
qcinds=which(str_detect(string=loctab[,"HMDB.ID"],pattern="^C") | str_detect(string=loctab[,"KEGG.ID"],pattern="^HMDB"))
matchpatlist=list("HMDB.ID"="^HMDB","KEGG.ID"="^C")
for(qcind in qcinds){
    idlist=loctab[qcind,c("HMDB.ID","KEGG.ID")]
    for(matid in names(matchpatlist)){
        idvec=idlist[str_which(string=idlist,pattern=matchpatlist[[matid]])]
        loctab[qcind,matid]=ifelse(length(idvec)==0,NA,idvec)
    }
}
annolist[["metabolomics"]]=loctab
# ensemble id and entrez for proteomics
loctab=annolist[["proteomics"]]
loctab$uniprot=loctab$UniProt
idtrans=list("ENTREZID"="entrezid","ENSEMBL"="ensembl")
for(idname in names(idtrans)){
    locvec=mapIds(org.Hs.eg.db,keys=loctab$uniprot,column=idname,keytype="UNIPROT",multiVals="first")
    locvec[sapply(locvec,is.null)]<-NA
    loctab[[idtrans[[idname]]]]=unlist(locvec)
}
# keep only metabolomics related protein
loctab=loctab[loctab[,"Panel"]=="Cardiometabolic",]
annolist[["proteomics"]]=loctab
prot_sele_list=paste0("proteomics",loctab[,"UniProt"])
for(groupi in seq(length(list_sig))){
    locvec=list_sig[[groupi]]
    ind_prot=str_detect(string=locvec,pattern="proteomics")
    ind_prot_select=locvec %in% prot_sele_list
    list_sig[[groupi]]=locvec[(!ind_prot)|(ind_prot_select)]
}
# lipidomics 
loctab=annolist[["lipidomics"]]
lipidconvtab=read.table(paste0(pardir,"data/metabolomics_lipidomics/Lipomatdb.txt"),header=TRUE,sep="\t")
lipidconvtab[,"Lipid_Name"]=str_replace_all(string=lipidconvtab[,"Lipid_Name"],pattern="\\_",replacement="/")
addids=lipidconvtab[match(loctab[,1],lipidconvtab[,"Lipid_Name"]),c("HMDB_ID","KEGG_ID")]
addids[,"KEGG_ID"]=str_remove(addids[,"KEGG_ID"],pattern="\xa0")
colnames(addids)=c("HMDB.ID","KEGG.ID")
loctab=cbind(loctab,addids)
annolist[["lipidomics"]]=loctab
save(list_sig,annolist,file="correlated_features_subclass.RData")
# 
feat_select=list("metabolomics"=c("Compound.name","Compound.name"),"lipidomics"=c("Compound.name","Compound.name"),"proteomics"=c("UniProt","UniProt"))
lipid_universe=unique(annolist[["lipidomics"]][,1])
progress<-function(n) cat(sprintf("task %d is complete\n", n))
opts<-list(progress=progress)
groups=names(list_sig)
enrich_df_list_comb<-foreach(ind=seq(length(list_sig)),.packages=c("clusterProfiler","org.Hs.eg.db","ReactomePA","metpath","stringr","Rodin"),.options.snow=opts)%dopar%{
    cluster_ele=list_sig[[ind]]
    clust_len=length(cluster_ele)
    navec=rep(list(NA),times=clust_len)
    id_df=tibble(featid=cluster_ele,clutnames=navec,ensembl=navec,uniprot=navec,entrezid=navec,keggid=navec,hmdbid=navec)
    # formaute id transform as a dataframe
    for(featname_i in seq_len(clust_len)){
        featname=cluster_ele[featname_i]
        matchomics=str_extract(featname,names(feat_select))
        matchomics=matchomics[!is.na(matchomics)]
        featname=str_replace_all(string=featname,pattern=fixed(matchomics),replacement="")
        ftset=feat_select[[matchomics]]
        temptab=annolist[[matchomics]]
        tabind=which(temptab[,ftset[1]]==featname)
        if(length(tabind)>1){
            tabind=tabind[1]
        }
        id_df[["clutnames"]][featname_i]=temptab[[ftset[2]]][tabind]
        if(str_detect(string=matchomics,pattern=fixed("metabolomics")) | str_detect(string=matchomics,pattern=fixed("lipidomics"))){
            id_df[["keggid"]][featname_i]=temptab[["KEGG.ID"]][tabind]
            id_df[["hmdbid"]][featname_i]=temptab[["HMDB.ID"]][tabind]
        }else if(str_detect(string=matchomics,pattern=fixed("proteomics"))){
            for(geneidtype in c("ensembl","uniprot","entrezid")){
                geneidvec=temptab[[geneidtype]][tabind]
                protid=unique(geneidvec[!is.na(geneidvec)])
                id_df[[geneidtype]][featname_i]=ifelse(length(protid)==0,NA,protid)
            }
        }
    }
    functionlist=list("go"={X<-function(x) enrichGO(gene=idvec,OrgDb="org.Hs.eg.db",keyType="ENSEMBL",ont=x,pvalueCutoff=0.05,pAdjustMethod="fdr")},
                        "kegg"={X<-function(x) enrichKEGG(gene=idvec,organism="hsa",keyType="entrezid",pvalueCutoff=0.05,pAdjustMethod="fdr",use_internal_data=T)},
                        "reactome"={X<-function(x) ReactomePA::enrichPathway(gene=idvec,organism="human",pvalueCutoff=0.05,pAdjustMethod="fdr")})
    genematch_set=data.frame("names"=c("BP","MF","CC","KEGG","reactome"),
                            "functionid"=c("go","go","go","kegg","reactome"),
                            "idtype"=c("ensembl","ensembl","ensembl","entrezid","entrezid"))
    reslist=vector(mode="list")
    for(rowi in seq(dim(genematch_set)[1])){
        idtype=genematch_set[rowi,"idtype"]
        idvec=unlist(id_df[[idtype]])
        idvec=idvec[!is.na(idvec)]
        name=genematch_set[rowi,"names"]
        if(length(idvec)==0){
            matchres=NULL
        }else{
            matchres=functionlist[[genematch_set[rowi,"functionid"]]](name)
        }
        if(!is.null(matchres)){
            tab=matchres@result
            matchind=list()
            for(pathi in seq(dim(tab)[1])){
                pathgenes=str_split(string=tab[pathi,"geneID"],pattern="\\/")[[1]]
                if(length(pathgenes)==1 && pathgenes==""){
                    matchind[[length(matchind)+1]]=""
                }else{
                    locind=which(sapply(id_df[[idtype]],function(x) any(str_detect(x,pattern=pathgenes))))
                    matchind[[length(matchind)+1]]=id_df[["featid"]][locind]
                }
            }
            matchres@result$matchind=matchind
        }
        reslist[[name]]=matchres
    }
    if(length(reslist)>0){
        rescolname=colnames(reslist[[1]]@result)
    }else{
        rescolname=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","matchind")
    }
    if(length(reslist)>0){
        loc_enr_tab=merge_result(reslist)@compareClusterResult
    }else{
        loc_enr_tab=data.frame(matrix(ncol=length(rescolname),nrow=0))
        colnames(loc_enr_tab)=rescolname
    }
    # enrichment metabolites
    functionlist=list("kegg"={X<-function(x) enrich_kegg(query_id=idvec,query_type="compound",id_type="KEGG",pathway_database=kegg_hsa_pathway,p_cutoff=0.05,p_adjust_method="fdr",threads=1)},
                    "hmdb"={X<-function(x) enrich_hmdb(query_id=idvec,query_type="compound",id_type="HMDB",pathway_database=hmdb_pathway,p_cutoff=0.05,p_adjust_method="fdr",threads=1)})
    genematch_set=data.frame("names"=c("kegg","hmdb"),"functionid"=c("kegg","hmdb"),"idtype"=c("keggid","hmdbid"))
    reslist=vector(mode="list")
    for(rowi in seq(dim(genematch_set)[1])){
        idtype=genematch_set[rowi,"idtype"]
        idvec=unlist(id_df[[idtype]])
        idvec=idvec[!is.na(idvec)]
        idvec %>% str_split(string=.,pattern="\\|") %>% unlist() -> idvec
        idvec=unique(idvec[idvec!="0" & !str_detect(idvec,pattern="[:punct:]") & idvec!=""])
        name=genematch_set[rowi,"names"]
        matchres=functionlist[[genematch_set[rowi,"functionid"]]](name)
        if(!is.null(matchres)){
            tab=matchres@result
            if(length(which(tab[,"p_value_adjust"]<0.05))==0){
                matchres=NULL
            }else{
                matchind=list()
                for(pathi in seq(dim(tab)[1])){
                    pathgenes=str_split(string=tab[pathi,"mapped_id"],pattern=";")[[1]]
                    if(length(pathgenes)==1 && pathgenes==""){
                        matchind[[length(matchind)+1]]=""
                    }else{
                        locind=which(sapply(id_df[[idtype]],function(x) any(str_detect(x,pattern=pathgenes))))
                        matchind[[length(matchind)+1]]=id_df[["featid"]][locind]
                    }
                }
                matchres@result$matchind=matchind
            }

        }
        reslist[[name]]=matchres
    }
    tablist=list()
    if(length(reslist)>0){
        rescolname=colnames(reslist[[1]]@result)
    }else{
        rescolname=c("pathway_id","pathway_name","describtion","pathway_class","p_value","p_value_adjust","all_id","all_number","mapped_id","mapped_number","mapped_percentage","matchind")
    }
    for(fieldname in c("kegg","hmdb")){
        if(is.null(reslist[[fieldname]])){
            df=data.frame(matrix(ncol=length(rescolname),nrow=0))
            colnames(df)=rescolname
            tablist[[fieldname]]=df
        }else{
            tablist[[fieldname]]=reslist[[fieldname]]@result
        }
    }
    meta_enr_all=rbind(tablist[[1]],tablist[[2]])
    meta_enr_all=meta_enr_all[meta_enr_all[,"p_value_adjust"]<=0.05,]
    # enrichment lipidomics
    omics="lipidomics"
    idvec=unlist(id_df[["featid"]])
    lipidind=str_which(string=idvec,pattern=omics)
    idvec=str_replace_all(string=idvec[lipidind],pattern=omics,replacement="")
    if(length(idvec)==0){
        lipid_enr_all=NULL
    }else{
        lipid_enr_all=lipid_enrich(query=data.frame(lipid=idvec),universe=data.frame(lipid=lipid_universe))
        lipid_enr_all$"p.adjust"=p.adjust(lipid_enr_all[,"p-value"],method="fdr")
        lipid_enr_all=lipid_enr_all[lipid_enr_all[,"p.adjust"]<=0.05,]
    }
    if(is.null(lipid_enr_all)){
        rescolname=c("Type","Classifier","Count.query","Count.universe","%.query","%.universe","p-value","FDR.q-value","Fold.change","p.adjust")
        lipid_enr_all=data.frame(matrix(ncol=length(rescolname),nrow=0))
        colnames(lipid_enr_all)=rescolname
    }
    # 
    meta_enr_all$cluster_id=rep(groups[ind],times=dim(meta_enr_all)[1])
    loc_enr_tab$cluster_id=rep(groups[ind],times=dim(loc_enr_tab)[1])
    lipid_enr_all$cluster_id=rep(groups[ind],times=dim(lipid_enr_all)[1])
    enrichcomb=list(tempdf=data.frame(id=cluster_ele,name=unlist(id_df[["clutnames"]])),nonmeta=loc_enr_tab,meta=meta_enr_all,lipid=lipid_enr_all)
    enrichcomb
}
enrichment_df_nometa=c()
enrichment_df_meta=c()
enrichment_df_lipid=c()
templist=list()
for(ind in seq(length(enrich_df_list_comb))){
    enrichment_df_nometa=rbind(enrichment_df_nometa,enrich_df_list_comb[[ind]][["nonmeta"]])
    enrichment_df_meta=rbind(enrichment_df_meta,enrich_df_list_comb[[ind]][["meta"]])
    enrichment_df_lipid=rbind(enrichment_df_lipid,enrich_df_list_comb[[ind]][["lipid"]])
    templist[[length(templist)+1]]=enrich_df_list_comb[[ind]][["tempdf"]]
}
namechange=colnames(enrichment_df_meta)
namechange[namechange=="pathway_name"]="Description"
namechange[namechange=="p_value_adjust"]="p.adjust"
colnames(enrichment_df_meta)=namechange
namechange=colnames(enrichment_df_nometa)
namechange[namechange=="geneID"]="mapped_id"
colnames(enrichment_df_nometa)=namechange
namechange=colnames(enrichment_df_lipid)
namechange[namechange=="Classifier"]="Description"
colnames(enrichment_df_lipid)=namechange
clust_res_anno=list(anno=templist,enrich=list(nonmeta=enrichment_df_nometa,meta=enrichment_df_meta,lipid=enrichment_df_lipid))
# 
save(clust_res_anno,file="enrichment_subclass.RData")
# summary assocition through 
baset2d=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all")
subt2d=c("sspg_avg_all","modified_DI_heyjun","ie_heyjun","hepatic_IR_heyjun","ffa_avg_heyjun")
fevec=list(posi=c(),nega=c())
featlist_overlap=list("allt2d"=fevec,"baset2d"=fevec,"sspg_avg_all"=fevec,"modified_DI_heyjun"=fevec,"ie_heyjun"=fevec,"hepatic_IR_heyjun"=fevec,"ffa_avg_heyjun"=fevec)
for(dirc in names(fevec)){
    selegp=paste0(dirc,"corr",baset2d)
    allfeat=unlist(list_sig[selegp])
    allfeat_tb=table(allfeat)
    featlist_overlap[["baset2d"]][[dirc]]=names(allfeat_tb[allfeat_tb==3])
}
for(dirc in names(fevec)){
    selegp=paste0(dirc,"corr",class_num_cols)
    allfeat=unlist(list_sig[selegp])
    allfeat_tb=table(allfeat)
    featlist_overlap[["allt2d"]][[dirc]]=names(allfeat_tb[allfeat_tb>=6])
}
for(dirc in names(fevec)){
    for(subt in subt2d){
        selegp=paste0(dirc,"corr",subt)
        allfeat=unlist(list_sig[selegp])
        otherfeat=unlist(list_sig[names(list_sig)!=selegp])
        featlist_overlap[[subt]][[dirc]]=setdiff(allfeat,otherfeat)
    }
}
# heatmap correlation metabolomics and lipidomics, padj<0.05, features appear >=2
tabstat_cor_sub=tabstat_cor[tabstat_cor[,"clinic"]%in%subphenogrp,]
for(omicch in c("lipidomics","metabolomics")){
    temptab=tabstat_cor_sub[tabstat_cor_sub[,"type"]==omicch&tabstat_cor_sub[,"padj"]<0.05,]
    valvec=table(temptab[,"feature"])
    selefeat=names(valvec)[valvec>1]
    p_mat_sele=p_mat[selefeat,]
    cor_mat_sele=cor_mat[selefeat,]
    # the original fdr value 
    p_mat_sele_fdr=matrix(1,nrow=length(selefeat),ncol=dim(cor_mat_sele)[2])
    rownames(p_mat_sele_fdr)=selefeat
    colnames(p_mat_sele_fdr)=colnames(cor_mat_sele)
    for(subtyp in colnames(cor_mat_sele)){
        for(feat in selefeat){
            padj=tabstat_cor[tabstat_cor[,"feature"]==feat&tabstat_cor[,"clinic"]==subtyp,"padj"]
            if(length(padj)>0){
                p_mat_sele_fdr[feat,subtyp]=padj
            }
        }
    }
    h=Heatmap(cor_mat_sele,name="corrmat",cluster_columns=TRUE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,
            cell_fun=function(j,i,x,y,w,h,fill){
                if(p_mat_sele[i,j]<0.05){
                    if(p_mat_sele_fdr[i,j]<0.05){
                        grid.text("**",x,y)
                    }else{
                        grid.text("*",x,y)
                    }
                }
            },
            col=circlize::colorRamp2(c(-1,0,1),c("blue","white","red")))#
    pdf(paste0("omics_corr_heatmap_metalipi_clinic.",omicch,".pdf"),width=7,height=15)#,width=7,height=30
    draw(h,main_heatmap="corrmat")
    dev.off()
    savdf=cor_mat_sele
    colnames(savdf)=str_remove(string=colnames(savdf),pattern="(\\_avg\\_all)|(\\_heyjun)")
    rownames(savdf)=str_remove(string=rownames(savdf),pattern="^L\\-")
    if(omicch=="metabolomics"){
        rownames(savdf)[rownames(savdf)=="C20:4,DC FA"]="Unknown1"
    }
    write.table(savdf,file=paste0("figext8_9",omicch,".txt"))
}