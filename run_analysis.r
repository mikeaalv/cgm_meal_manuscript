rm(list=ls())
pathdir="PATH"
ee<-new.env()
# preprocessing
# step 2 of metadata preprocessing (no manual steps afterwards)
sys.source(paste0(pathdir,"metadata/preprocess_step2.r"),ee)
# omics preprocessing 
sys.source(paste0(pathdir,"omics/preprocess.r"),ee)
# cgm preprocesing preprocess.py (run in jupyternotebook) -> preprocess.r -> qc_plot.m -> manual_feature.py (run in terminal)
# qc_plot.m: 1C, S1, 3A, 4C, S4, S5, 
sys.source(paste0(pathdir,"cgm/preprocessing/preprocess.r"),ee)
# qc heatmap replot
sys.source(paste0(pathdir,"cgm/preprocessing/qc_heatmap.r"),ee) 
# replicates patterns (Fig S2)
# cgm/replicates/compare_replicates.py
# cgm/replicates/compare_replicates.py

# meta data
# statistics of sample numbers (Fig 1B)
sys.source(paste0(pathdir,"metadata/statistics_sample.r"),ee)
# daily meal information from cronometer
sys.source(paste0(pathdir,"metadata/daily_meal_prop.r"),ee)
# cohort statistic table (table 1)
sys.source(paste0(pathdir,"metadata/cohort_stat_tab.r"),ee)

# cgm feature plot (Fig 1D, 1E, S3)
# ./cgm/visualization/figures_cgm_features.py
# mean curves of classes (Fig 2D, 4B, S8)
# ./cgm/visualization/mean_curv_classes.m
# patterns between spike and mean nutrient in food (Fig 1F)
sys.source(paste0(pathdir,"cgm/visualization/cgm_mean_nutrient.r"),ee)
# cgm and meta data realted plot (Fig 2A, 2C, 3B, 3D, 4A, S7, S10, S19)(Table S2, Table S8,)
sys.source(paste0(pathdir,"cgm/visualization/individual_foodtypes.r"),ee)
# cgm and meta data realted plot for AUC, time to peak, and glucose at 2h
sys.source(paste0(pathdir,"cgm/visualization/individual_foodtypes_expd.r"),ee)
# cgm data correlation heatmap (Fig 2B, S6, )
sys.source(paste0(pathdir,"cgm/visualization/plot_corr.r"),ee)
# cgm data correlation heatmap for AUC and time to peak (Fig S9)
sys.source(paste0(pathdir,"cgm/visualization/plot_corr_expd.r"),ee)
# cgm data compare between groups (Fig 3C, 4D, 4E)
sys.source(paste0(pathdir,"cgm/visualization/group_compare_plot.r"),ee)
# a different (quantile 1d) definination of spikers
sys.source(paste0(pathdir,"cgm/visualization/group_compare_plot_thres.r"),ee)
# test between cgm curves to different food (table S1, table S9)
sys.source(paste0(pathdir,"cgm/visualization/cgm_feature_test.r"),ee)
# linear model check for a few targetted statistic test (table S3)
sys.source(paste0(pathdir,"cgm/visualization/covar_check.r"),ee)

# metabolic lipidomics data assocation with metadata (Fig 5A, 5B, S11)
sys.source(paste0(pathdir,"omics/association_metadata.r"),ee)
# pathway enrichment (Fig S14, S15, S16) (Data S1)
sys.source(paste0(pathdir,"omics/pathway_enrich.r"),ee)
# annova test between subtype for each omics features (Fig S12, S13)
sys.source(paste0(pathdir,"omics/metab_subtype_association.r"),ee)
# omics features that associated with a few different spikes or mitigation effects
sys.source(paste0(pathdir,"omics/omics_repeat_heatmap.r"),ee)
# additional figures for lipidomics
sys.source(paste0(pathdir,"omics/lipids_target_vis.r"),ee)
# associaiton between omics and food spikes
sys.source(paste0(pathdir,"omics/asso_oimcs_spike.r"),ee)
# group wise test  (Fig 5C)(Table S4)
sys.source(paste0(pathdir,"omics/pathway_enrich_gp_test.r"),ee)

# microbiome sample data matching
sys.source(paste0(pathdir,"microbiome/idmatch.r"),ee)
# microbiome data formating and initial visualization 
sys.source(paste0(pathdir,"microbiome/quan_eda.r"),ee)
# diversity calculation
# ./microbiome/diversity_calculation.py
# beta diversity pcoa scatter plot
sys.source(paste0(pathdir,"microbiome/beta_diversity_pcoa.r"),ee)
# microbiome association with metabolic subtypes (Fig 5D, S13, )
sys.source(paste0(pathdir,"microbiome/micro_assoc_subtype.R"),ee)
# microbiome association with food spikes (Fig S18)
sys.source(paste0(pathdir,"microbiome/micro_assoc_spike.R"),ee)
# microbiome association with metabolomics
sys.source(paste0(pathdir,"microbiome/micro_assoc_otheromics.R"),ee)
# lm between food spiker types
sys.source(paste0(pathdir,"microbiome/micro_ap_comp_spike.R"),ee)

# mitigator group comparison in clinical and omics data
sys.source(paste0(pathdir,"cgm/visualization/mitigation_g_comp.r"),ee)
# mean spikes and mitigators associated with microbiome (Table S6)
sys.source(paste0(pathdir,"microbiome/meanspike_omics.r"),ee)

# mediation analysis (Fig 5E) (Table S5)
sys.source(paste0(pathdir,"mediation/mediat_micr_metlip_clispi.r"),ee)

# correlation network and clustering
# sys.source(paste0(pathdir,"corr_network/corrnet.r"),ee)
# partial correaltion for subtype and spikes 
sys.source(paste0(pathdir,"omics/subtype_spike_pcorr.r"),ee)

# prediction model (Fig S20)
# predictive_model

# food analysis
sys.source(paste0(pathdir,"foodana/starch_comp.r"),ee)