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
omicsvec=c(rep("metabolomics",times=dim(annolist[[1]])[[1]]),rep("lipidomics",times=dim(annolist[[2]])[[1]]),rep("proteomics",times=dim(annolist[[3]])[[1]]))
featvec=c(annolist[[1]][,"Compound.name"],annolist[[2]][,"Compound.name"],annolist[[3]][,"UniProt"])
save(omicsmat,omicsvec,featvec,file="omics_data_clean.RData")
# select correlated features
food_sig_mat=cbind(mat_quantile,mat_mitigator_effect)
# potato vs grapes
# food_sig_mat$"potato_vs_grape"=unlist(food_sig_mat[,"Potatoes"]/food_sig_mat[,"Grapes"])
# 
nomics=dim(omicsmat)[2]
nfoodsig=dim(food_sig_mat)[2]
foodtypes=colnames(food_sig_mat)
# 
cor_mat=matrix(NA,nrow=nomics,ncol=nfoodsig)
colnames(cor_mat)=foodtypes
p_mat=cor_mat
list_sig=vector(mode="list",length=nfoodsig*2)
names(list_sig)=paste0(rep(c("posicorr","negacorr"),each=nfoodsig),rep(foodtypes,times=2))
stattab=c()
for(omicsi in seq(nomics)){
    for(foodsigi in seq(nfoodsig)){
        xvec=food_sig_mat[,foodsigi]
        yvec=omicsmat[,omicsi]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        cor_mat[omicsi,foodsigi]=corres$estimate
        p_mat[omicsi,foodsigi]=corres$p.value
        if(corres$p.value<0.05){
            compdname=paste0(omicsvec[omicsi],featvec[omicsi])
            dirc=ifelse(corres$estimate>0,"posicorr","negacorr")
            listname=paste0(dirc,foodtypes[foodsigi])
            list_sig[[listname]]=c(list_sig[[listname]],compdname)
            stattab=rbind(stattab,data.frame(food=foodtypes[foodsigi],feature=featvec[omicsi],omics=omicsvec[omicsi],corr=corres$estimate,pval=corres$p.value))
        }
    }
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
# update the significant feature list accoding to used protein annotation table
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
save(list_sig,annolist,stattab,file="correlated_features.RData")
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
save(clust_res_anno,file="enrichment.RData")
# save for output 
enrich_protein=clust_res_anno[["enrich"]][["nonmeta"]]
enrich_metabolomics=clust_res_anno[["enrich"]][["meta"]]
enrich_lipids=clust_res_anno[["enrich"]][["lipid"]]
save(enrich_protein,enrich_metabolomics,enrich_lipids,file="enrichment_summary.RData")
# Metabolomcis: heatmap, rank by fiber effect on rice spike and plot non lipid metabite and non-medicine (mostly amino acids)
# plotlist=c("Threonine","O-Succinyl-L-homoserine","N-Acetyl-D-galactosaminitol","Phenylalanyl-Tryptophan","D(-)-Fructose","Glycine","Homoarginine","L-Proline","Ethopabate","Glucaric acid","Cortisol","Thiabendazole","Indole-3-carboxaldehyde","Hippuric acid","3-carboxy-4-methyl-5-propyl-2-furanpropanoate (CMPF)","Hydroxyprogesterone","4-formyl Indole(1)","Citramalic acid","Cys-Gly or Gly-Cys","D-TRYPTOPHAN")
plotlist=unlist(list_sig[c("posicorrRice+Fiber","negacorrRice+Fiber")])
plotlist=plotlist[str_detect(string=plotlist,pattern="metabolomics")]
plotlist=str_remove(string=plotlist,pattern="metabolomics")
subj_order=names(sort(food_sig_mat[,"Rice+Fiber"]))#negative to positive
omicsind=sapply(plotlist,function(x) which(featvec==x))
showmat=omicsmat[subj_order,omicsind]
colnames(showmat)=plotlist
h=Heatmap(t(showmat),name="omicsfiber",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE)#
pdf(paste0("pathway_metabolomics_fiber_heatmap.pdf"))
draw(h,main_heatmap="omicsfiber")
dev.off()
write.table(showmat,file="fiber_heatmap.txt")
# Beans positive and negative assocaiton, realted to Histidine metaboism 
meta_enrich_tab=clust_res_anno[["enrich"]][["meta"]]
inds=str_which(string=meta_enrich_tab[,"cluster_id"],pattern="Beans")
plotvec=c()
namesvec=c()
for(ind in inds){
    loctab=meta_enrich_tab[ind,]
    clustereles<-loctab[,"mapped_id"]%>%str_split(string=.,pattern=";")%>%unlist()%>%unique()
    clusterty=ifelse(loctab[,"cluster_id"]=="posicorrBeans",1,-1)
    plotvec=c(plotvec,rep(clusterty,times=length(clustereles)))
    namesvec=c(namesvec,clustereles)
}
names(plotvec)=namesvec
pv.out<-pathview(cpd.data=plotvec,pathway.id="00340",species="hsa",out.suffix="00340",kegg.native=TRUE,cpd.idtype="kegg",same.layer=TRUE)
genelist=pv.out$plot.data.gene[,"kegg.names"]
genevec=rep(1,times=length(genelist))
names(genevec)=genelist
pv.out<-pathview(gene.data=genevec,cpd.data=plotvec,pathway.id="00340",species="hsa",out.suffix="00340",kegg.native=TRUE,cpd.idtype="kegg",same.layer=TRUE)
# heatmap beans spike related
# plotlist=c("Ergothioneine","3-Methylhistidine","L-Glutamic acid","Methylimidazoleacetic acid","L-HISTIDINE","Aspartate","1-Acetylimidazole","Beta-Citryl-L-glutamic acid","(S)-(+)-2-Hydroxy-3-methylbutyric acid","Hydroxybutyric acid(1)","Alpha-ketoisovaleric acid","gamma-glutamylthreonine(1)")
plotlist=unlist(list_sig[c("posicorrBeans","negacorrBeans")])
plotlist=plotlist[str_detect(string=plotlist,pattern="metabolomics")]
plotlist=str_remove(string=plotlist,pattern="metabolomics")
plotlist=plotlist[!plotlist%in%c("HISTIDINE","Histidine")]
subj_order=names(sort(food_sig_mat[,"Beans"]))#negative to positive
omicsind=sapply(plotlist,function(x) which(featvec==x))
showmat=omicsmat[subj_order,omicsind]
# 
# load(paste0(pardir,"result/meal_prop/meal_prop.RData"))
# mealprpfeatlist=c("proteins")
# tab_prop=tab_food_prop
# tab_prop=tab_prop[,c("subject_id",mealprpfeatlist)]
# tab_prop[,"subject_id"]<-tab_prop[,"subject_id"]%>%str_remove("STUDYID-")%>%as.numeric()%>%as.character()
# tab_prop=tab_prop[tab_prop[,"subject_id"]!="SOMESUBJECTID",]
# rownames(tab_prop)=tab_prop[,"subject_id"]
# showmat=cbind(showmat,scale(tab_prop[subj_order,mealprpfeatlist]))
# colnames(showmat)=c(plotlist,"protein prop daily")
colnames(showmat)=c(plotlist)
h=Heatmap(t(showmat),name="omicsbean",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE)#
pdf(paste0("pathway_metabolomics_bean_heatmap.pdf"))
draw(h,main_heatmap="omicsbean")
dev.off()
write.table(showmat,file="bean_heatmap.txt")
# heatmap rice spike related
vec_carnitine=featvec[str_detect(string=featvec,pattern="[Cc]arnitin")]
plotlist=vec_carnitine
names(plotlist)=vec_carnitine
subj_order=names(sort(food_sig_mat[,"Rice"]))#negative to positive
omicsind=sapply(plotlist,function(x) which(featvec==x))
showmat=omicsmat[subj_order,omicsind]
colnames(showmat)=names(plotlist)
# ethnicity
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
metafeatlist=c("ethnicity")
metadata=metadata[,c("study_id",metafeatlist)]
metadata[,"study_id"]<-metadata[,"study_id"]%>%str_remove("STUDYID-")%>%as.numeric()%>%as.character()
metadata=metadata[metadata[,"study_id"]!="SOMESUBJECTID",]
rownames(metadata)=metadata[,"study_id"]
metatab=metadata[subj_order,metafeatlist,drop=FALSE]
# 
h=Heatmap(t(showmat),name="omicsrice",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE)#
colors=structure(1:4,names=c("White","Asian","Asian, White","Hispanic or Latino")) 
h_metadata_ch=Heatmap(t(metatab),name="metadata_ch",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=FALSE,show_column_names=TRUE,col=colors)#
pdf(paste0("hormone_metabolomics_rice_heatmap.pdf"))
draw(h %v% h_metadata_ch,main_heatmap="omicsrice")
dev.off()
write.table(showmat,file="rice_heatmap.txt")
# heatmap potatoes spiking 
plotlist=c("Taurocholic acid","Tauroursodeoxycholic acid","P41159","P09619")
names(plotlist)=c("Taurocholic acid","Tauroursodeoxycholic acid","LEP","PDGFRB")
subj_order=rownames(food_sig_mat)
foods=colnames(food_sig_mat)
omicsind=sapply(plotlist,function(x) which(featvec==x))
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
metafeatlist=c("bmi_avg_all")
metadata=metadata[,c("study_id",metafeatlist)]
metadata[,"study_id"]<-metadata[,"study_id"]%>%str_remove("STUDYID-")%>%as.numeric()%>%as.character()
metadata=metadata[metadata[,"study_id"]!="SOMESUBJECTID",]
rownames(metadata)=metadata[,"study_id"]
showmat=cbind(food_sig_mat,omicsmat[,omicsind])
showmat=cbind(showmat,scale(metadata[subj_order,metafeatlist]))
colnames(showmat)=c(foods,names(plotlist),"BMI")
omicsvec=c(names(plotlist),"BMI")
cor_mat=matrix(NA,nrow=length(omicsvec),ncol=length(foods))
colnames(cor_mat)=foods
rownames(cor_mat)=omicsvec
p_mat=cor_mat
for(omics in omicsvec){
    for(food in foods){
        xvec=showmat[,omics]
        yvec=showmat[,food]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        cor_mat[omics,food]=corres$estimate
        p_mat[omics,food]=corres$p.value
    }
}
h=Heatmap(cor_mat,name="corrmat",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,
        cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat[i,j]<0.05){
                grid.text("*",x,y)
            }
        },
        col=circlize::colorRamp2(c(-1,0,1), c("blue","white","red")))#
pdf(paste0("omics_corr_potato_heatmap.pdf"))
draw(h,main_heatmap="corrmat")
dev.off()
write.table(cor_mat,file="potato_cor_heatmap.txt")
# count the number of unsaturated lipids for each group
food_dirc=names(list_sig)
plotdf=c()
for(fd_name in food_dirc){
    parts=str_split(string=fd_name,pattern="corr")[[1]]
    featlist=list_sig[[fd_name]]
    lipidlist=featlist[str_which(string=featlist,pattern="lipidomics")]
    c_n_tab=lipidlist%>%str_extract(string=.,pattern="TAG\\d+:\\d+")%>%str_remove(string=.,pattern="TAG")%>%str_split(string=.,pattern=":",simplify=TRUE)
    if(dim(c_n_tab)[1]<1){
        next
    }
    c_n_tab=c_n_tab[!is.na(c_n_tab[,1]),,drop=FALSE]
    class(c_n_tab)="numeric"
    if(dim(c_n_tab)[1]<1){
        next
    }
    unsat_ratio_vec=c_n_tab[,2]/c_n_tab[,1]
    nlip=length(unsat_ratio_vec)
    plotdf=rbind(plotdf,data.frame(food=rep(parts[2],times=nlip),direction=rep(parts[1],times=nlip),prop=unsat_ratio_vec))
}
p<-ggplot(plotdf,aes(x=food,y=prop,color=direction))+geom_boxplot()+ggtitle("proportion of unsaturation in TAG")+theme(strip.background=element_blank(),strip.text.x=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))#
ggsave("TAG_unsaturate_prop.pdf",plot=p)
# count numbers of correlated lipids in each class
food_dirc=names(list_sig)
lipid_interest_list=c("TAG","PE","PC","FFA","DAG")
plotdf=c()
for(fd_name in food_dirc){
    parts=str_split(string=fd_name,pattern="corr")[[1]]
    featlist=list_sig[[fd_name]]
    lipidlist=featlist[str_which(string=featlist,pattern="lipidomics")]
    lipidlist=str_remove(string=lipidlist,pattern="lipidomics")
    for(lipidele in lipid_interest_list){
        n=sum(str_detect(string=lipidlist,pattern=paste0("^",lipidele)))
        plotdf=rbind(plotdf,data.frame(food=parts[2],direction=parts[1],class=lipidele,n=n))
    }
}
p<-ggplot(plotdf,aes(x=food,y=n,fill=class))+geom_bar(stat="identity")+facet_wrap(~direction,ncol=1)+theme(strip.background=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))#
ggsave(paste0("lipid_n.pdf"),plot=p)
# correlation heatmap of omics X food+mitigators
for(dirc in c("positive","negative")){
    if(dirc=="positive"){
        loctab=stattab[stattab[,"corr"]>0,]
    }else{
        loctab=stattab[stattab[,"corr"]<0,]
    }
    for(omics in c("metabolomics","proteomics","lipidomics")){
        loctab2=loctab[loctab[,"omics"]==omics,]
        feat_all=unique(loctab2[,"feature"])
        showmat=matrix(0.0,nrow=length(feat_all),ncol=length(foodtypes))
        rownames(showmat)=feat_all
        colnames(showmat)=foodtypes
        for(rowi in seq(dim(loctab2)[1])){
            showmat[loctab2[rowi,"feature"],loctab2[rowi,"food"]]=loctab2[rowi,"corr"]
        }
        h=Heatmap(showmat,name="omicscorr",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,col=circlize::colorRamp2(c(-1,0,1), c("blue","white","red")))#
        pdf(paste0("pathway_",dirc,omics,"_heatmap.pdf"),height=20)
        draw(h,main_heatmap="omicscorr")
        dev.off()
    }
}
save(list=ls(all.names=TRUE),file="pathway_heatmap_plot.RData")
# selected lists for heatmap, metabolomics
lists_compd=list(centr_meta=c("Glucose","Glycerol 3-phosphate","FUMARATE","L-Malic acid","DL-Lactic acid"),aa=c("Threonine","HISTIDINE","Serine","DL-Glutamate","Aspartate","D-TRYPTOPHAN","Glycine","L-Proline","L-Tyrosine","L-Leucine","L-ISOLEUCINE","L-Valine"),keto=c("Hydroxybutyric acid(1)","(R)-3-Hydroxybutyric acid","Alpha-ketoisovaleric acid","(S)-(+)-2-Hydroxy-3-methylbutyric acid"),bileacids=c("Tauroursodeoxycholic acid"),carnitine=c(),others=c("Uridine"))
allfeat=unique(stattab[,"feature"])
lists_compd[["carnitine"]]=allfeat[str_detect(string=allfeat,pattern=fixed("carnitine",ignore_case=TRUE))]
featplot=unlist(lists_compd)
omicsind=sapply(featplot,function(x) which(featvec==x))
# meta data
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
metafeatlist=c("a1c_avg_all","sspg_avg_all")
metadata=metadata[,c("study_id",metafeatlist)]
metadata[,"study_id"]<-metadata[,"study_id"]%>%str_remove("STUDYID-")%>%as.numeric()%>%as.character()
metadata=metadata[metadata[,"study_id"]!="SOMESUBJECTID",]
rownames(metadata)=metadata[,"study_id"]
# 
showmat=cbind(food_sig_mat,omicsmat[,omicsind])
showmat=cbind(showmat,scale(metadata[subj_order,metafeatlist]))
colnames(showmat)=c(foodtypes,featplot,"a1c","sspg")
# 
plottypes=c(foodtypes,"a1c","sspg")
cor_mat=matrix(NA,nrow=length(featplot),ncol=length(plottypes))
colnames(cor_mat)=plottypes
rownames(cor_mat)=featplot
p_mat=cor_mat
for(omics in featplot){
    for(food in plottypes){
        xvec=showmat[,omics]
        yvec=showmat[,food]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        cor_mat[omics,food]=corres$estimate
        p_mat[omics,food]=corres$p.value
    }
}
splitn=sapply(lists_compd,length)
split=rep(names(lists_compd),times=splitn)
# 
h=Heatmap(cor_mat,name="corrmat",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,row_split=split,
        cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat[i,j]<0.05){
                grid.text("*",x,y)
            }
        },
        col=circlize::colorRamp2(c(-1,0,1),c("blue","white","red")))#
pdf(paste0("omics_corr_heatmap_main.pdf"))
draw(h,main_heatmap="corrmat")
dev.off()
write.table(cor_mat,file="fullselect_cor_heatmap.txt")
outputatab=c()
for(rowele in rownames(cor_mat)){
    for(colele in colnames(cor_mat)){
        outputatab=rbind(outputatab,data.frame(metabolites=rowele,features=colele,correlation=cor_mat[rowele,colele],pvalue=p_mat[rowele,colele]))
    }
}
outputatab$adjust=p.adjust(outputatab$pvalue,method="fdr")
write.table(outputatab,file="metabolomics_select_cor_heatmap_tab.txt",row.names=FALSE)