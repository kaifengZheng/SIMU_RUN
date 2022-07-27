#!/usr/bin/env bash
rm -f *.log
rm -f *.log
rm -rf ./fdmnes
rm -rf js
rm -f hostname
mkdir ./fdmnes
mkdir ./fdmnes/calculation
mkdir js

sbatch run_bash.sh

