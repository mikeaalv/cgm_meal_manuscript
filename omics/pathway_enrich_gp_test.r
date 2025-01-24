# groupwise test 
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
require(ggpubr)
require(lsr)
# 
require(clusterProfiler)
require(org.Hs.eg.db)
require(ReactomePA)
require(metpath)
require(ComplexHeatmap)
require(pathview)
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
load(paste0(pardir,"result/cgm_meal/heatmap_plot.RData"))
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
source("/Users/yuewu/Documents/GitHub/meal_cgm/omics/func_lipid_enrich.r",local=TRUE)
load("omics_data_clean.RData")
# remove proteomics duplications
nondupmask=!duplicated(featvec)
omicsmat=omicsmat[,nondupmask]
omicsvec=omicsvec[nondupmask]
featvec=featvec[nondupmask]
# select sig change features
colnames(omicsmat)=featvec
groupvec=c("Bread","Grapes","Potatoes","Rice","Pasta","Mitigator")
testtab=merge(omicsmat,mat_worstfood[,groupvec[seq(5)]],by="row.names")
rownames(testtab)=testtab[,1]
testtab=testtab[,-1]
mitig_tab=as.data.frame(ifelse(rowMeans(mat_mitigator_effect)>0,1,0))
colnames(mitig_tab)="Mitigator"
testtab=merge(testtab,mitig_tab,by="row.names")
rownames(testtab)=testtab[,1]
testtab=testtab[,-1]
# 
nomics=dim(omicsmat)[2]
nfoodsig=length(groupvec)
foodtypes=groupvec
# 
delta_mat=matrix(NA,nrow=nomics,ncol=nfoodsig)
colnames(delta_mat)=foodtypes
p_mat=delta_mat
list_sig=vector(mode="list",length=nfoodsig*2)
names(list_sig)=paste0(rep(c("spikh","spikl"),each=nfoodsig),rep(foodtypes,times=2))
stattab=c()
for(foodsigi in seq(nfoodsig)){
    pvec=c()
    statvec=c()
    civec=c()
    Dfvec=c()
    effectvec=c()
    detltavec=c()
    for(omicsi in seq(nomics)){
        locdf=testtab[,c(featvec[omicsi],foodtypes[foodsigi])]
        colnames(locdf)=c("value","group")
        nonaind=(!is.na(locdf[,1]))&(!is.na(locdf[,2]))
        locdf=locdf[nonaind,]
        if(length(unique(locdf[,"value"]))<=1 | min(table(locdf[,"group"]))<=3){
            pvec=c(pvec,NA)
            detltavec=c(detltavec,NA)
            statvec=c(statvec,NA)
            civec=c(civec,NA)
            Dfvec=c(Dfvec,NA)
            effectvec=c(effectvec,NA)
        }else{
            teststat=t.test(value~group,data=locdf)
            # pval=wilcox.test(value~group,data=locdf)$p.value
            delta=mean(locdf[locdf[,"group"]==1,"value"],na.rm=TRUE)-mean(locdf[locdf[,"group"]==0,"value"],na.rm=TRUE)
            pvec=c(pvec,teststat$p.value)
            detltavec=c(detltavec,delta)
            statvec=c(statvec,teststat$statistic)
            civec=c(civec,paste0(formatC(teststat$conf.int,format="E",digits=2),collapse=" "))
            Dfvec=c(Dfvec,teststat$parameter)
            effectvec=c(effectvec,cohensD(value~group,data=locdf))
        }
    }
    delta_mat[,foodsigi]=detltavec
    p_mat[,foodsigi]=pvec
    padj=rep(NA,times=length(pvec))
    namsk=is.na(pvec)
    for(theomics in unique(omicsvec)){
        themask=omicsvec==theomics
        locmask=themask&(!namsk)
        padj[locmask]=p.adjust(pvec[locmask],method="fdr")
    }
    sigmask=which(pvec<0.05)
    len=length(sigmask)
    if(len>0){
        compdnames=paste0(omicsvec[sigmask],featvec[sigmask])
        dirc=ifelse(detltavec[sigmask]>0,"spikh","spikl")
        listnames=paste0(dirc,foodtypes[foodsigi])
        for(listname_i in seq(length(listnames))){
            listname=listnames[listname_i]
            list_sig[[listname]]=c(list_sig[[listname]],compdnames[listname_i])
        }
    }
    stattab=rbind(stattab,data.frame(food=rep(foodtypes[foodsigi],times=len),feature=featvec[sigmask],omics=omicsvec[sigmask],delta=detltavec[sigmask],pval=pvec[sigmask],padj=padj[sigmask],statistics=statvec[sigmask],ci=civec[sigmask],Df=Dfvec[sigmask],effectsz=effectvec[sigmask]))
}
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
for(grouphere in names(list_sig)){
    locvec=list_sig[[grouphere]]
    ind_prot=str_detect(string=locvec,pattern="proteomics")
    ind_prot_select=locvec %in% prot_sele_list
    list_sig[[grouphere]]=locvec[(!ind_prot)|(ind_prot_select)]
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
save(list_sig,annolist,stattab,file="correlated_features_gp.RData")
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
                    # pathgenes=pathgenes[pathgenes!=""]
                    # idlist=str_split(string=id_df[[idtype]],pattern="\\|")
                    # locind=sapply(idlist,function(x) any(any(x %in% pathgenes)))
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
save(clust_res_anno,file="enrichment_gptest.RData")
# save for output 
enrich_protein=clust_res_anno[["enrich"]][["nonmeta"]]
enrich_metabolomics=clust_res_anno[["enrich"]][["meta"]]
enrich_lipids=clust_res_anno[["enrich"]][["lipid"]]
save(enrich_protein,enrich_metabolomics,enrich_lipids,file="enrichment_summary_gptest.RData")
write.table(stattab,file="omics_gp_comp.tsv",row.names=FALSE)
# stat test comparison bread
center_gr="Bread"
compfeat="N1-Methyladenosine"
foodsvis=groupvec[seq(5)]
plotdf=testtab[,c(compfeat,foodsvis)]
worstfoodvec=rep(NA,times=nrow(plotdf))
for(foodhere in foodsvis){
    worstfoodvec[which(plotdf[,foodhere]==1)]=foodhere
}
plotdf$worstfood=worstfoodvec
colnames(plotdf)[1]="value"
comparisons=list()
for(food in setdiff(foodsvis,center_gr)){
    comparisons[[length(comparisons)+1]]=c(food,center_gr)
}
p<-ggboxplot(plotdf,x="worstfood",y="value",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
ggsave(plot=p,paste0("boxplot_spiketype_meta",compfeat,"_compare_",center_gr,".pdf"))
# lipids bar 0.2 fdr cut off
lipidlist=c("PALMITOLEIC ACID","MYRISTIC ACID","TAG56:3-FA16:0","TAG52:0-FA16:0","C12:1,DC FA(4)")
foodpair=c("Potatoes","Grapes")
plotdf=testtab[,c(lipidlist,foodpair)]
worstfoodvec=rep(NA,times=nrow(plotdf))
for(foodhere in foodpair){
    worstfoodvec[which(plotdf[,foodhere]==1)]=foodhere
}
worstfoodvec[is.na(worstfoodvec)]="Others"
plotdf$worstfood=worstfoodvec
comparisons=list("1"=c("Potatoes","Grapes"),"2"=c("Potatoes","Others"),"3"=c("Grapes","Others"))
for(colfeat in lipidlist){
    plotdftemp=plotdf[,c(colfeat,"worstfood")]
    colnames(plotdftemp)[1]="value"
    p<-ggboxplot(plotdftemp,x="worstfood",y="value",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)+ggtitle(colfeat)
    ggsave(plot=p,paste0("boxplot_spiketype_meta",colfeat,"_compare_.pdf"))
}
# plot showing each groups
for(thecent in c("Potatoes","Grapes")){
    for(colfeat in lipidlist){
        plotdf=testtab[,c(colfeat,foodsvis)]
        worstfoodvec=rep(NA,times=nrow(plotdf))
        for(foodhere in foodsvis){
            worstfoodvec[which(plotdf[,foodhere]==1)]=foodhere
        }
        plotdf$worstfood=worstfoodvec
        colnames(plotdf)[1]="value"
        comparisons=list()
        for(food in setdiff(foodsvis,thecent)){
            comparisons[[length(comparisons)+1]]=c(food,thecent)
        }
        p<-ggboxplot(plotdf,x="worstfood",y="value",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
        ggsave(plot=p,paste0("boxplot_spiketype_meta",colfeat,"_compare_",thecent,".pdf"))
    }
}
compdlist=c("N1-Methyladenosine","PALMITOLEIC ACID","MYRISTIC ACID","TAG56:3-FA16:0","TAG52:0-FA16:0")
savdf=testtab[,c(foodsvis,compdlist)]
worstfoodvec=rep(NA,times=nrow(savdf))
for(foodhere in foodsvis){
    worstfoodvec[which(savdf[,foodhere]==1)]=foodhere
}
savdf$carb_response_type=worstfoodvec
savdf=savdf[,c("carb_response_type",compdlist)]
savdf=savdf[!is.na(savdf[,2]),]
