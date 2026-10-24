#!/bin/bash


#SBATCH --job-name=Heatmap_para_16
#SBATCH --mail-user=q.guilloit@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-16
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%A_%a.out 
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00-72:00:00

## Load software modules
module load julia/1.10.5

# Execute
arg2=4
arg3=0.0001
arg4="fixeSimple"
arg5=50
arg6=0
arg7="Heatmap_para_16"


if (( $SLURM_ARRAY_TASK_ID % arg2 == 1 )); then
    julia --project=mon_environnement Main_Modele0_MeSU.jl $SLURM_ARRAY_TASK_ID $arg2 $arg3 $arg4 $arg5 $arg6 $arg7
elif (( $SLURM_ARRAY_TASK_ID % arg2 == 2 )); then
    julia --project=mon_environnement Main_Modele1_DD_MeSU.jl $SLURM_ARRAY_TASK_ID $arg2 $arg3 $arg4 $arg5 $arg6 $arg7
elif (( $SLURM_ARRAY_TASK_ID % arg2 == 3 )); then
    julia --project=mon_environnement Main_Modele2_DD_MeSU.jl $SLURM_ARRAY_TASK_ID $arg2 $arg3 $arg4 $arg5 $arg6 $arg7
elif (( $SLURM_ARRAY_TASK_ID % arg2 == 0 )); then
    julia --project=mon_environnement Main_Modele3_DD_MeSU.jl $SLURM_ARRAY_TASK_ID $arg2 $arg3 $arg4 $arg5 $arg6 $arg7
else
    echo "Index de tâche inattendu: $SLURM_ARRAY_TASK_ID"
    exit 1
fi



## Finish gracefully
exit 0




