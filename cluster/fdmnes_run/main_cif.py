import numpy as np
import glob
import os
import shutil
import subprocess
import random
import string
import numpy as np
import time,datetime
import json
import argparse

import sys




        
#inputs = sys.argv     
#ncores  = int(inputs[1])
#ind   = int(inputs[2])
#print(ind)


parser = argparse.ArgumentParser(description='calculation configuration')
parser.add_argument('-np','--ncores',type=int,help="number of cores")
parser.add_argument('-in','--inputfile',type=int,help="index of specific input file")
parser.add_argument('--hostfile',type=str,help="The location of hostfile")
parser.add_argument('--host',type=int,help="assign specific host to run the program")
args = parser.parse_args()
ncores = args.ncores
ind = args.inputfile
hostname=args.hostfile
host_ind=args.host








def fdmnes_calculator_mpi(js,ncores,host_ind):

    '''
    # USAGE:
    with open('/data/software/FDMNES/fdmnes_tests/4_nospin/fdmnes.inp') as fi:
        fdmnesinp = fi.readlines()
    js =  {'task':'fdmnes_run', 'fdmnes_inp': fdmnesinp, 'ncores':12} 
    js_out = fdmnes_calculator_mpi(js)
    js_out
    '''
    
    #runtime parameters
    #try:
    #    ncores = js['ncores']
    #except:
    #    ncores=4

        
    try:
        fdmnes_scratch_path = js['fdmnes_scratch_path']
    except:
        fdmnes_scratch_path = '/gpfs/scratch/kaifzheng/fdmnes/calculation'
        
    # try:
    #     mpirun_cmd = js['mpirun_cmd']
    # except:
    #     mpirun_cmd = '/opt/intel/oneapi/mpi/2021.5.1/bin/mpirun'
    
    try:
        exe_path = js['exe_path']
    except:
        exe_path = '/gpfs/home/kaifzheng/software/parallel_fdmnes'

    try:
        cleanup = js['cleanup']
    except:
        cleanup = 'true'
        
    
    
    
    try:

        os.chdir(fdmnes_scratch_path)
        uid = 'fdmnes_'+''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=10))
        os.makedirs(uid,exist_ok=True)
        os.chdir(uid)
        
        print('Running FDMNES calculation at \n %s/%s \n using %d cores\n'%(fdmnes_scratch_path,uid,ncores))
        
        fi = open("fdmnes.inp", "w")
        for ii in js['fdmnes_inp']:
            fi.write(ii)
        fi.close()
        

        try:
            js['cif']
            cif = open("structure.cif", "w")
            for cc in js['cif']:
                cif.write(cc)
            cif.close()
        except:
            pass
        




        #fdmnes wants this
        fi = open("fdmfile.txt", "w")
        fi.write('1\n')
        fi.write('fdmnes.inp\n')
        fi.close()          
        fi = open("hostfile","w")
        for hh in js['hostfile']:
            fi.write(hh)
        fi.close()
        exe_list = [
            'mpirun_fdmnes'                               #'fdmnes_mpi_linux64',
        ]

        start_time = time.time()  

        begin_time = time.time()
        
        host = js["hostfile"][host_ind][:-1]
        print(host)
        for e in exe_list:
            _ = subprocess.run(['HOST_NUM_FOR_MUMPS=4  bash %s/%s -np %d --hostfile %s --host %s >> fdmnes.out ' %(exe_path,e,ncores,'hostfile',host)],shell=True)
            #_ = subprocess.run(['srun -n %d %s/%s >> fdmnes.out' %(ncores,exe_path,e)],shell=True)
           # _ = subprocess.run(['srun --nodelist=%s -N 1 -n %d %s/%s >> fdmnes.out' %(node,ncores,exe_path,e)],shell=True)
            if _.returncode > 0:
                print('error at %s'%e)
                print(_)
                break
        finish_time = time.time()
        print(f"running time:{finish_time-begin_time}")
        print("===================FINISHED======================\n\n\n")

        
      

        with open('fdmnes.inp') as f:
            fdmnesinp = f.readlines()
        with open('fdmnes.out') as f:
            fdmnesout = f.readlines()
        with open('fdmnes_bav.txt') as f:
            fdmnesbav = f.readlines()

        js = {
            'start_time': datetime.datetime.fromtimestamp(start_time).strftime('%Y-%m-%d__%H:%M:%S'),
            'finish_time': datetime.datetime.fromtimestamp(finish_time).strftime('%Y-%m-%d__%H:%M:%S'),
            'time_elapsed': finish_time - start_time,
            'fdmnesinp': fdmnesinp,
            'fdmnesout': fdmnesout,
            'fdmnesbav': fdmnesbav,
            'ncores': ncores,
            'exe_path': exe_path,
            }            
            


        # single site
        if os.path.isfile('fdmnes.txt'):
            with open('fdmnes.txt') as f:
                fdmnestxt = f.readlines()
            # non-magnetic calculation
            if len(fdmnestxt[-1].split()) == 2:
                js_new = {
                    'fname_fdmnestxt': 'fdmnes.txt',
                    'header_fdmnestxt': fdmnestxt[0:2],
                    'e':     [float(i.split()[0]) for i in fdmnestxt[2:]],
                    'mu':    [float(i.split()[1]) for i in fdmnestxt[2:]],
                    }
                js.update(js_new)
            # magnetic calculation
            elif len(fdmnestxt[-1].split()) == 3:
                js_new = {
                    'fname_fdmnestxt': 'fdmnes.txt',
                    'header_fdmnestxt': fdmnestxt[0:2],
                    'e':     [float(i.split()[0]) for i in fdmnestxt[2:]],
                    'mu_up':  [float(i.split()[1]) for i in fdmnestxt[2:]],
                    'mu_dw':  [float(i.split()[2]) for i in fdmnestxt[2:]],
                    'mu':  [(float(i.split()[1])+float(i.split()[2])) for i in fdmnestxt[2:]],
                    }
                js.update(js_new)


        # multiple sites
        else:

            import glob 
            outputs = glob.glob('fdmnes_*.txt')

            for i in outputs:
                try:
                    t = int(i.split('_')[1].split('.txt')[0])
                    
                    with open('fdmnes_%d.txt'%t) as f:
                        fdmnestxt = f.readlines()
                            
                    # non-magnetic calculation
                    if len(fdmnestxt[-1].split()) == 2:
                        js_new = {
                            'fname_fdmnestxt_%d'%t: i,
                            'header_fdmnestxt_%d'%t: fdmnestxt[0:2],
                            'e_%d'%t:     [float(i.split()[0]) for i in fdmnestxt[2:]],
                            'mu_%d'%t:    [float(i.split()[1]) for i in fdmnestxt[2:]],
                            }
                        js.update(js_new)

                    # magnetic calculation
                    elif len(fdmnestxt[-1].split()) == 3:
                        js_new = {
                            'fname_fdmnestxt_%d'%t: i,
                            'header_fdmnestxt_%d'%t: fdmnestxt[0:2],
                            'e_%d'%t:     [float(i.split()[0]) for i in fdmnestxt[2:]],
                            'mu_up_%d'%t:  [float(i.split()[1]) for i in fdmnestxt[2:]],
                            'mu_dw_%d'%t:  [float(i.split()[2]) for i in fdmnestxt[2:]],
                            'mu_%d'%t:  [(float(i.split()[1])+float(i.split()[2])) for i in fdmnestxt[2:]],
                            }
                        js.update(js_new)  

                except:
                    pass            





            
        if os.path.isfile('fdmnes_conv.txt'):
            with open('fdmnes_conv.txt') as f:
                fdmnestxt_conv = f.readlines()
            js_new = {
                'fname_fdmnesconvtxt': 'fdmnes_conv.txt',
                'header_fdmnesconvtxt': fdmnestxt[0:1],
                'e_conv':     [float(i.split()[0]) for i in fdmnestxt_conv[1:]],
                'mu_conv':    [float(i.split()[1]) for i in fdmnestxt_conv[1:]],
                }
            js.update(js_new)


        os.chdir('..')

        if cleanup == 'true':
            shutil.rmtree(uid)
            
            
        return js


    except Exception as exc:
        print('something is wrong\n')
        with open("fdmnes_error.txt") as f:
             error=f.realines()
        print(error)
        print(exc)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
import glob
        
path='/gpfs/home/kaifzheng/FDMNES_cal/Fe_proj/'
#cif_files=glob.glob(path+"scaled_structures/*.cif")

os.chdir(path)
cif_files=glob.glob("input/*.cif")


with open(cif_files[ind]) as fi:
    cif = fi.readlines()

with open(hostname) as fi:
    hn=fi.readlines() 

with open(path+'fdmnes.inp') as fi:
    fdmnesinp = fi.readlines()

        
js_in  = {'fdmnes_inp': fdmnesinp,'cif': cif,"hostfile":hn,'cleanup':'false'}        


print(f"Running calculation using {cif_files[ind]}")
js_out =fdmnes_calculator_mpi(js_in,ncores,host_ind)


os.chdir(path)
pathlen = 6

with open(path+'js/js_'+cif_files[ind][pathlen:-4]+'.json', 'w') as f:
    json.dump(js_out, f)


        
        
        
        
        
        
        
        
        
        
        
        
        



