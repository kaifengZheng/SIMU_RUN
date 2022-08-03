#!/usr/bin/env bash

#SBATCH --job-name=array_test
#SBATCH --output=array_test.%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH -N 1
#SBATCH -p short-96core
#SBATCH -t 04:00:00
#SBATCH --array=1-3


mkdir temp
module load python
module load gnu-parallel/6.0

START=$SLURM_ARRAY_TASK_ID
NUMLINES=2
STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

scontrol show hostnames >> hostname_$SLURM_ARRAY_TASK_ID

echo "START=$START"
echo "STOP=$STOP"

for (( N = $START; N <= $STOP; N++ ))
do
    LINE=$(sed -n "$N"p indexes.txt)
    echo $LINE >> temp/tasks_${START}_${STOP}
done

#cat temp/tasks_${START}_${STOP}

cat temp/tasks_${START}_${STOP} | parallel --verbose --jobs 2 python main_nocif.py -np 24 -in {} --hostfile hostname_$SLURM_ARRAY_TASK_ID --host 0

# clean up the temp files
rm temp/tasks_${START}_${STOP}
