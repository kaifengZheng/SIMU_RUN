#!/usr/bin/env bash
#SBATCH --nodes=4
#SBATCH --time=10:59:59
#SBATCH --partition=long-24core
#SBATCH --output=array_test.%N_%a.log
#SBATCH --job-name=
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

module load slurm
module load  python/3.9.7
##module load  mpi/2021.4.0
##module load intel/oneAPI/2021.4

##module load gnu-parallel

NUM_NODES=4
CPU_PER_TASK=24
CPU_PER_PARTITION=24



nodeindex=0
scontrol show hostnames >> hostname

filenum=$(ls -1 "input" | wc -l)

let numcpus=$((NUM_NODES*CPU_PER_PARTITION))

let dividcpus=$((numcpus/CPU_PER_TASK))

for ((i=0;i<$filenum;i++))
  do
    if (( $nodeindex>=$NUM_NODES )); then
        let nodeindex=0
    fi

    if [ $((i%dividcpus)) -eq 0 -a $i -ne 0 ];then
        wait
    fi
    echo $i 
    echo "using the $nodeindex th node" 
    python main_cif.py -np $CPU_PER_TASK -in $i --hostfile hostname --host $nodeindex &  
    let "modu = ((((i+1)*CPU_PER_TASK)%CPU_PER_PARTITION))"
    if [ $modu -eq 0 -a $i ];then
       (( nodeindex += 1 ))
    fi
  done
wait

