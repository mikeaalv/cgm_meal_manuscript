#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time 200:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH -D ANOTHERPATH
#SBATCH -J huamnn_join
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=FAIL
#SBATCH -p batch
#SBATCH --account=USERNAME
#SBATCH -o PATH/THEFILE_huamnn_join_%A.log
#SBATCH -e PATH/THEFILE_huamnn_join_%A.err

module load anaconda/3_5.0.1_20180125
module load humann3

merge_metaphlan_tables.py -o ../humann/species_full.tsv ../humann/*bugs_list* 
humann_join_tables --input ../humann/ --output ../humann/genefamilies_full.tsv --file_name _genefamilies_cpm.tsv
humann_join_tables --input ../humann/ --output ../humann/kegg_full.tsv --file_name _kegg-cpm.tsv
humann_join_tables --input ../humann/ --output ../humann/pfam_full.tsv --file_name pfam-cpm.tsv
humann_join_tables --input ../humann/ --output ../humann/go_full.tsv --file_name go-cpm.tsv
humann_join_tables --input ../humann/ --output ../humann/rxn_full.tsv --file_name rxn-cpm.tsv
humann_join_tables --input ../humann/ --output ../humann/pathway_full.tsv --file_name pathabundance.tsv
humann_join_tables --input ../humann/ --output ../humann/pathway_full_cpm.tsv --file_name pathabundance_addon_cpm.tsv