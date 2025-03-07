# the enrichment function for lipidomics
lipid_enrich<-function(query,universe){
    # code adapated from https://github.com/PNNL-Comp-Mass-Spec/Rodin/blob/master/Example%20script/Example_of_use_target_list_mode.R
    cleanedquery=clean.lipid.list(query)
    cleaneduniverse=clean.lipid.list(universe)
    lipid.miner(cleanedquery,name="Query",TGcollapse.rm=TRUE)
    lipid.miner(cleaneduniverse,name="Universe",TGcollapse.rm=TRUE)
    # different enrichment runs
    intact_fisher_cat_result=intact.fisher(Query.intact$Category,Universe.intact$Category)
    intact_fisher_main_result=intact.fisher(Query.intact$"Main class",Universe.intact$"Main class")
    intact_fisher_sub_result=intact.fisher(Query.intact$"Sub class",Universe.intact$"Sub class")
    # chain_fisher_result=chain.fisher(Query.chain,Universe.chain)
    allchains_fisher_result=allchains.fisher(Query.allchains,Universe.allchains)
    total_carbon_cat_result=total.carbon.cat(Query.intact,Universe.intact,enrich=FALSE)
    total_carbon_main_result=total.carbon.main(Query.intact,Universe.intact,enrich=FALSE)
    # total_carbon_sub_result=total.carbon.sub(Query.intact,Universe.intact,enrich=FALSE)
    total_DB_cat_result=total.DB.cat(Query.intact,Universe.intact,enrich=FALSE)
    total_DB_main_result=total.DB.main(Query.intact,Universe.intact,enrich=FALSE)
    # total_DB_sub_result=total.DB.sub(Query.intact,Universe.intact,enrich=FALSE)
    allchains_cat_result=allchains.cat(Query.intact,Universe.intact,enrich=FALSE)
    allchains_main_result=allchains.main(Query.intact,Universe.intact,enrich=FALSE)
    # allchains_sub_result=allchains.sub(Query.intact,Universe.intact,enrich=FALSE)
    subclass_mainclass_result=subclass.mainclass(Query.intact,Universe.intact,test="Fisher",enrich=FALSE)
    #reformat the tables to add the type of test
    intact_fisher_cat_result<-cbind(Type=(rep("Intact")),intact_fisher_cat_result)
    intact_fisher_main_result<-cbind(Type=(rep("Main Class")),intact_fisher_main_result)
    intact_fisher_sub_result<-cbind(Type=(rep("Sub Class")),intact_fisher_sub_result)
    # chain_fisher_result<-cbind(Type=(rep("Chain characteristics")),chain_fisher_result)
    allchains_fisher_result<-cbind(Type=(rep("Specific chain")),allchains_fisher_result)
    total_carbon_cat_result<-cbind(Type=(rep("TotalCarbon by Cat")),total_carbon_cat_result)
    total_carbon_main_result<-cbind(Type=(rep("TotalCarbon by Main")),total_carbon_main_result)
    # total_carbon_sub_result<-cbind(Type=(rep("TotalCarbon by Sub")),total_carbon_sub_result)
    total_DB_cat_result<-cbind(Type=(rep("total DB by Cat")),total_DB_cat_result)
    total_DB_main_result<-cbind(Type=(rep("total DB by Main")),total_DB_main_result)
    # total_DB_sub_result<-cbind(Type=(rep("total DB by Sub")),total_DB_sub_result)
    allchains_cat_result<-cbind(Type=(rep("Specific Chain by Cat")),allchains_cat_result)
    allchains_main_result<-cbind(Type=(rep("Specific Chain by Main")),allchains_main_result)
    # allchains_sub_result<-cbind(Type=(rep("Specific Chain by Sub")),allchains_sub_result)
    subclass_mainclass_result<-cbind(Type=(rep("Subclass by Mainclass")),subclass_mainclass_result)
    # place all the tables in a globaloutput table
    globaloutput<-rbind(intact_fisher_cat_result,intact_fisher_main_result,intact_fisher_sub_result,total_carbon_cat_result,total_carbon_main_result,total_DB_cat_result,total_DB_main_result,allchains_cat_result,allchains_main_result,subclass_mainclass_result, allchains_fisher_result)#chain_fisher_result,total_carbon_sub_result,total_DB_sub_result,allchains_sub_result
    return(globaloutput)
}
