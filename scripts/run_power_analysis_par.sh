#! /bin/bash
#$ -q iblm.q
#$ -V
#$ -o /iblm/netapp/home/jezhou/sge/out
#$ -e /iblm/netapp/home/jezhou/sge/err
#$ -cwd
#$ -t 1-20
#$ -tc 5

proj_dir=$HOME/power_analysis

sample="Rat_Amygdala_787A_all_seq"

params=`awk -v line=$SGE_TASK_ID 'NR == line' $proj_dir/outs/simulation_params.csv`
repls=`echo $params | cut -d',' -f1`
effect=`echo $params | cut -d',' -f2`

outfh=${proj_dir}/outs/${sample}_${repls}reps_${effect}effect.txt

echo "running power analysis with replicates = $repls and effect size = $effect"
echo "writing output to $outfh"
Rscript $proj_dir/scripts/power_analysis.R --replicates $repls --effectSize $effect --out $outfh