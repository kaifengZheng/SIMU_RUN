#!/usr/bin/env bash
rm -f *.log
rm -rf /gpfs/scratch/user/fdmnes
rm -rf js
rm -f hostname
mkdir /gpfs/scratch/user/fdmnes
mkdir /gpfs/scratch/user/fdmnes/calculation
mkdir js

sbatch run_bash.sh

