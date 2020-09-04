#! /bin/bash
#$ -q iblm.q
#$ -V
#$ -o /iblm/netapp/home/jezhou/sge/out
#$ -e /iblm/netapp/home/jezhou/sge/err
#$ -cwd

#Rscript $HOME/power_analysis/scripts/estimate_nb_params.R \
#--count_mtx /iblm/netapp/home/jezhou/power_analysis/data/snRNA/Rat_Opioid_HS_1_premrna/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
#--sample Rat_Opioid_HS_1_premrna 


Rscript $HOME/power_analysis/scripts/estimate_nb_params.R 
