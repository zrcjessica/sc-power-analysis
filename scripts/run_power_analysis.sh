#! /bin/bash
#$ -q iblm.q
#$ -V
#$ -o /iblm/netapp/home/jezhou/sge/out
#$ -e /iblm/netapp/home/jezhou/sge/err
#$ -cwd

proj_dir=$HOME/power_analysis
Rscript $proj_dir/scripts/power_analysis.R --replicates 2 --effectSize 3