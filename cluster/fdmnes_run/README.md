# Running XAS simulations on cluster using full speed
## Chapter I FDMNES simulation(test successfully on Seawulf cluster(Stony Brook University))
FDMNES, for Finite Difference Method Near Edge Strucuture uses DFT to calculate K edge and L<sub>23</sub> edge XANES of all elements. The Project link: https://fdmnes.neel.cnrs.fr/<br><br>
This folder is for calculating mulptiple FDMNES calculations. In this package, there are two mode. One mode is for calculation of crystalized materials, which need cif input files. The second mode is for clusters or molecules, which don't have periodic boundary condition. The input files are multiple fdmnes.inps.<br><br>
When use this package, please do the following check before submit your job
**TODO LIST**
[] Check the fdmnes.inp file (edge, absorber...)
[] Check the main_cif.py or main_nocif.py
* **main_cif.py(main_nocif.py)**
  * line 61-64(line 65-68) change the dir to the workspace to preform all the simulations in. It needs enough big for all the calculation.
```python
    try:
        fdmnes_scratch_path = js['fdmnes_scratch_path']
    except:
        fdmnes_scratch_path = '/gpfs/scratch/user/fdmnes/calculation'
```
  * line 71-74(line 75-78) change the dir to the place stored fdmnes executable files
```python
    try:
        exe_path = js['exe_path']
    except:
        exe_path = '/gpfs/home/user/software/parallel_fdmnes'
```
  * line 132(line 136) Running FDMNES using 4 cores for MUMPS calculation, user can change this value, but the total number of cpus should be the multiplier of this number. 
```python
 _ = subprocess.run(['HOST_NUM_FOR_MUMPS=4  bash %s/%s -np %d --hostfile %s --host %s >> fdmnes.out ' %(exe_path,e,ncores,'hostfile',host)],shell=True)
```
  * line 283(line 285) The dir user stores stuff in
```python
path='/gpfs/home/user/FDMNES_cal/Fe_proj/'
```
* **run_bash.sh**
  * line 1-8 configuration
```bash
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --time=167:59:59
#SBATCH --partition=extended-96core
#SBATCH --output=array_test.%N_%a.log
#SBATCH --job-name=Athena
#SBATCH --mail-user=<user email>
#SBATCH --mail-type=ALL
```
  * line 10-11 module
```bash
module load slurm
module load  python/3.9.7
```
  * line 17-19 how to specify cpus for several jobs after cpu allocation. The total number of cpus should smaller than the total cpus allocated.
```bash
NUM_NODES=1
CPU_PER_TASK=12
CPU_PER_PARTITION=96
```
* **start_bash.sh**
   * check the dirs
