#!/bin/bash
#$ -V
#$ -o /iblm/netapp/home/jezhou/sge/out
#$ -e /iblm/netapp/home/jezhou/sge/err
#$ -cwd

cellranger aggr --id=rat_snrnaseq_agg --csv=$HOME/power_analysis/rat_snrnaseq_agg.csv

