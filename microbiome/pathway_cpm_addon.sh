#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time 200:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH -D ANOTHERPATH
#SBATCH -J pathway_addon
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=FAIL
#SBATCH -p batch
#SBATCH --account=USERNAME
#SBATCH -o PATH/THEFILE_pathwayaddon_%A.log
#SBATCH -e PATH/THEFILE_pathwayaddon_%A.err

module load anaconda/3_5.0.1_20180125
module load humann3

for thefile in *_pathabundance.tsv; do
    [[ $thefile =~ (.*)_pathabundance.tsv ]];
    echo ${BASH_REMATCH[1]}
    humann_renorm_table --input ../humann/${BASH_REMATCH[1]}_pathabundance.tsv --output ../humann/${BASH_REMATCH[1]}_pathabundance_addon_cpm.tsv --units cpm --update-snames
done
