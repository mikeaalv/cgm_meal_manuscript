#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time 200:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH -D ANOTHERPATH
#SBATCH -J huamnn_idconv
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=FAIL
#SBATCH -p batch
#SBATCH --account=USERNAME
#SBATCH -o PATH/THEFILE_humann_idconv_%A.log
#SBATCH -e PATH/THEFILE_humann_idconv_%A.err

module load anaconda/3_5.0.1_20180125
module load humann3
# convert ids for go, kegg, pfam, rxn
humann_rename_table --input go_full.tsv --output go_full_idconv.tsv --custom ../datamap/map_go_name.txt.gz
humann_rename_table --input kegg_full.tsv --output kegg_full_idconv.tsv --custom ../datamap/map_ko_name.txt.gz
humann_rename_table --input pfam_full.tsv --output pfam_full_idconv.tsv --custom ../datamap/map_pfam_name.txt.gz
humann_rename_table --input rxn_full.tsv --output rxn_full_idconv.tsv --names metacyc-rxn