#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time 200:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH -D ANOTHERPATH
#SBATCH -J huamnn
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=FAIL
#SBATCH -p batch
#SBATCH --account=USERNAME
#SBATCH -o PATH/THEFILE_humann_%A.log
#SBATCH -e PATH/THEFILE_humann_%A.err

module load anaconda/3_5.0.1_20180125
module load humann3

# work on temp folder for the huge temp files (even though deleted later)
humann -i ../trim/THEFILE_clean.fastq.gz -o $TMPDIR --threads 12 -r -v --metaphlan-options "-t rel_ab --bowtie2db ../database/ --index mpa_v30_CHOCOPhlAn_201901" --nucleotide-database /labs/mpsnyder/brooksa/db_humann/chocophlan/chocophlan/ --protein-database /labs/mpsnyder/brooksa/db_humann/uniref90/uniref/  --nucleotide-identity-threshold 0.0 --nucleotide-query-coverage-threshold 90.0 --nucleotide-subject-coverage-threshold 50.0 --o-log ../humann/log/THEFILE.txt --output-basename THEFILE

cp $TMPDIR/THEFILE* ../humann
tempfoldname=$(ls -d $TMPDIR/*/ |head -n 1)
cp ${tempfoldname}/*bugs_list.tsv ../humann

# normalizaiton
humann_renorm_table --input ../humann/THEFILE_genefamilies.tsv --output ../humann/THEFILE_genefamilies_cpm.tsv --units cpm --update-snames
# each quantificaiton types
humann_regroup_table --input ../humann/THEFILE_genefamilies_cpm.tsv --output ../humann/THEFILE_rxn-cpm.tsv --groups uniref90_rxn
humann_regroup_table --input ../humann/THEFILE_genefamilies_cpm.tsv --output ../humann/THEFILE_go-cpm.tsv -c ../datamap/map_go_uniref90.txt.gz
humann_regroup_table --input ../humann/THEFILE_genefamilies_cpm.tsv --output ../humann/THEFILE_pfam-cpm.tsv -c ../datamap/map_pfam_uniref90.txt.gz
humann_regroup_table --input ../humann/THEFILE_genefamilies_cpm.tsv --output ../humann/THEFILE_kegg-cpm.tsv -c ../datamap/map_ko_uniref90.txt.gz

# format conversion
module load qiime2/2019.7
grep -v "|" < ../humann/THEFILE_genefamilies.tsv > ../humann/THEFILE_genefamilies_unstratified.tsv
grep -v "|" < ../humann/THEFILE_pathabundance.tsv > ../humann/THEFILE_pathabundance_unstratified.tsv
# 
biom convert -i ../humann/THEFILE_genefamilies_unstratified.tsv -o ../humann/THEFILE_genefamilies.biom --table-type="Gene table" --to-hdf5
biom convert -i ../humann/THEFILE_pathabundance_unstratified.tsv -o ../humann/THEFILE_pathabundance.biom --table-type="Pathway table" --to-hdf5

