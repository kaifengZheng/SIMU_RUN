#!/usr/bin/env bash
rm -f *.log
rm -rf /gpfs/scratch/kaifzheng/fdmnes
rm -rf js
rm -f hostname
mkdir /gpfs/scratch/kaifzheng/fdmnes
mkdir /gpfs/scratch/kaifzheng/fdmnes/calculation
mkdir js

sbatch run_bash.sh

