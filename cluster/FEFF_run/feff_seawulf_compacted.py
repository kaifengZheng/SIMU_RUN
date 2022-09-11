#!/data/software/miniconda3/bin/python

import matplotlib
matplotlib.use('Agg') # see http://matplotlib.org/faq/usage_faq.html#what-is-a-backend


import numpy as np
import os
import subprocess
from scipy.spatial import distance_matrix
from pymatgen.core.periodic_table import Element
import shutil
from scipy.interpolate import InterpolatedUnivariateSpline
import argparse
from collections import Counter
from collections import Counter
from pymatgen.core.periodic_table import Element
from scipy.spatial import distance_matrix
import numpy as np
import pandas as pd
import time

import pyarrow as pa
import pandas as pd
import numpy as np
import json
import pyarrow.parquet as pq



from sys import argv

import subprocess

from glob import glob
filename=glob("input/*.xyz")

parser = argparse.ArgumentParser(description='calculation configuration')
parser.add_argument('-np','--ncores',type=int,help="number of cores")
parser.add_argument('-in','--inputfile',type=int,help="index of specific input file")
parser.add_argument('-hf','--hostfile',type=str,default="hostfile",help="The name of hostfile(for cluster user), by default, \"hostfile\" ")
parser.add_argument('-host',type=int,help="assign specific host to run the program(for cluster user)")
parser.add_argument('-ca','--center_atom',type=str,help="the element of absorber atom")
parser.add_argument('-edge',type=str,default='K',help="edge of measuring atom, by default, edge=K")
parser.add_argument('-rFMS','-rfms',default=7.2,type=float,help="rfms, by default, rfms=7.2")
parser.add_argument('-rSCF','-rscf',default=5.6,type=float,help="rscf, by default escf=5.6")
parser.add_argument('--corehole','-ch',default='RPA',type=str,help="corhole, by default, corhole=RPA")
parser.add_argument('--nmax','-nm',default=20,type=int,help="nmax, by default nmax=20")
parser.add_argument('--dmax','-dm',default=12,type=float,help="dmax, by default dmax=12")
parser.add_argument('-S02', default=0,type=float,help="amplitude reduction factor, by default S02=0(calculate automatically by FEFF)")
parser.add_argument('-show','--ShowInputFiles',action='store_true',help="get the input file names and the indexes of files")
parser.add_argument('-wm','--WorkMachine',type=int,default=0,help="work place:(1): workstation(a single linux system) or (0):cluster(Seawulf)")
#parser.add_argument('--site_file','-sites',type=int,default=0,help="index of the site file in the particle folder")

args = parser.parse_args()
ncores = args.ncores
ind = args.inputfile
ca=args.center_atom
edge=args.edge

hostname=args.hostfile
host_ind=args.host
S02=args.S02
rFMS=args.rFMS
rSCF=args.rSCF
nmax=args.nmax
dmax=args.dmax
corehole=args.corehole
#site_file=args.site_file


if args.ShowInputFiles==True:
    for i in range(len(filename)):
        print(f"{i}\t{filename[i]}")
    exit()
if args.WorkMachine!=0 and args.WorkMachine!=1:
    print("argument value error, the value should be eith 0(on cluster) or 1(on workstation)")
    exit()
# inputs = sys.argv
# for n,i in enumerate(inputs):
#     s=i.split('=')
#     if s[0] == 'ca' : ca = s[1]
#     if s[0] == 'rFMS' or s[0] == 'rfms' : rFMS = float(s[1])
#     if s[0] == 'rSCF' or s[0] == 'rscf' : rSCF = float(s[1])
#     if s[0] == 'corehole' or s[0] == 'ch' : corehole = s[1]
#     if s[0] == 'nmax' or s[0] == 'nm' : nmax = int(s[1])
#     if s[0] == 'dmax' or s[0] == 'dm' : dmax = float(s[1])    



def feff_runner(ncores,hostfile=None,host_ind=None):

    exe_list = [
        'rdinp',
        'dmdw',
        'atomic',
        'pot',
        'ldos',
        'screen',
        'crpa',
        'opconsat',
        'xsph',
        'fms',
        'mkgtr',
        'path',
        'genfmt',
        'ff2x',
        'sfconv',
        'compton',
        'eels',
        'rhorrp'
    ]
    if args.WorkMachine==0:
        mpirun_cmd='mpirun'
        exe_path='/gpfs/projects/FrenkelGroup/softwares/FEFF10/bin/MPI'
    if args.WorkMachine==1:
        mpirun_cmd = '/opt/intel/oneapi/mpi/2021.5.1/bin/mpirun'
        exe_path =  '/data/software/FEFF/bin/mpi/'
    open('feff.out', 'a' ).close()
    if args.WorkMachine==0:
        #work on cluster
        if args.host==None or args.hostfile==None:
          for e in exe_list:
             _ = subprocess.run(['%s -np %d %s/%s'%(mpirun_cmd,ncores,exe_path,e)],shell=True)
        if args.host !=None and args.hostfile!=None:
          
          #fi = open("../../hostname","r")
          #with open(hostfile) as fi:
          #   hn=fi.readlines()
          #fi.clost()
          #start_time = time.time()
          #host = hn[host_ind][:-1]

          for e in exe_list:
             _ = subprocess.run(['%s -np %d --hostfile %s --host %s %s/%s ' %(mpirun_cmd,ncores,"../../hostname",host_ind,exe_path,e)],shell=True)
             subprocess.run(["wait"])

    if args.WorkMachine==1:
        #work on the workstation
        for e in exe_list:
             _ = subprocess.run(['%s -np %d %s/%s >> feff.out' %(mpirun_cmd,ncores,exe_path,e)],shell=True)# >> feff.out
        
   
    
def del_emptyline(data):
    return [line for line in data if line.strip() != '' and len(line.split())==4]
    
def read_xyz(filename):
    with open(filename,'r') as f1:
        xyzdata=f1.readlines()
    xyzdata=del_emptyline(xyzdata[1:])
    labels=[]
    positions=[]
    for i in range(len(xyzdata)):
        lines=xyzdata[i].split()
        #print(lines)
        labels.append(lines[0])
        positions.append([float(lines[1]),float(lines[2]),float(lines[3])])
    natoms=dict(Counter(labels))
    return labels,natoms,positions
   
def equ_sites(ca_atom:str,labels,natoms,positions,cutoff,randomness=4):
    """
    :param positions: coordinates
    :param cutoff:    cutoff distance defined by mutiple scattering radius
    :return:          non-equ position indexes
    """
    # cutoff method
    def duplicates(lst, item):
        """
        :param lst: the whole list
        :param item: item which you want to find the duplication in the list
        :return: the indexes of duplicate items
        """
        return [i for i, x in enumerate(lst) if x == item]
    
    caatoms=[i for i, e in enumerate(labels) if e == ca_atom]

    dis_all =  np.around(distance_matrix(np.array(positions),np.array(positions)[caatoms]),decimals=randomness)
    dis_all.sort(axis=1)
    dis_cut = [list(dis_all[i][dis_all[i] < cutoff]) for i in range(len(dis_all))]
    dup = []
    for i in range(len(dis_cut)):
        dup.append(duplicates(dis_cut, dis_cut[i])[0])
    #unique_index = list(set(dup))  # set can delete all duplicated items in a list
    unique_index = dict()
    for i in range(len(dup)):
        if dup[i] in unique_index:
            unique_index[dup[i]].append(i)
        else:
            unique_index.update({dup[i]:[i]})
    # formula_text
    formula = '$'
    for i,(k,v) in enumerate(natoms.items()):
        formula = formula + k + '_{' + str(int(v)) + '}'
    formula = formula + '$'
    # sort it using sorted method. Do not use list.sort() method, because it returns a nonetype.
    #unique_index = np.array(sorted(unique_index))
    #print("number of atoms: {}".format(len(positions)))
    #print("number of unique atoms: {}".format(len(atom_index))) #
    
    return unique_index,formula  #keys are those unique sites, values are the cooresponding equ-sites for those unique_sites

# this writes feff.inp


    
    
def sorted_coord(ca,pos,labels,natoms):
    dis=distance_matrix(pos,[pos[ca]])
    dis_dict=dict()
    #create a dataframe, which is easy to sort and reorder
    for i in range(len(pos)):
        dis_dict[i]={'dis':float(dis[i]),'labels':labels[i],'x':pos[i][0],'y':pos[i][1],'z':pos[i][2]}
    dis_dict=pd.DataFrame(dis_dict).T
    center_atom=pos[ca]
    for i in range(len(dis_dict)):
        dis_dict['x'][i]=dis_dict['x'][i]-center_atom[0]
        dis_dict['y'][i]=dis_dict['y'][i]-center_atom[1]
        dis_dict['z'][i]=dis_dict['z'][i]-center_atom[2]
    dis_dict=dis_dict.sort_values(by='dis')
    dis_dict=dis_dict.reset_index(drop=True)
    

    dis_array=[]
    positions=[]
    labels=[]
    for i in range(len(dis_dict)):
        positions.append([dis_dict['x'][i],dis_dict['y'][i],dis_dict['z'][i]])
        dis_array.append(dis_dict['dis'][i])
        labels.append(dis_dict['labels'][i])
    index_array=[]
    for i in range(len(labels)):
        if i==0:
            index_array.append(0)
        else:
            index_array.append(list(natoms.keys()).index(labels[i])+1)
            
    positions=np.array(positions)
    dis_array=np.array(dis_array)
    index_array=np.array(index_array)
    return positions,dis_array,labels,index_array
        
    
   

def write_feffinp(ca,dmax,rFMS,rSCF,corehole,edge,S02,filename): 
    

    [labels, natoms, positions] = read_xyz(filename)
    #print("=======================================\n")
    #print(natoms)
    #print("\n========================================\n")
    #natoms=list(natoms.values())
    labels_short = []
    for i in labels:
        if i not in labels_short:
            labels_short.append(i)
    positions_modified,dis_array,labels,index=sorted_coord(ca,positions,labels,natoms)
    
 
    f=open('feff.inp',"w")    
    f.write("""TITLE xx             
                                                   
EDGE      %(edge)s
S02       %(S02)2.1f
COREHOLE  %(corehole)s                                  
CONTROL   1 1 1 1 1 1                               
                                                   
XANES 6 0.05 0.1
                             
                                                   
FMS       %(rFMS)2.1f                      
          
EXCHANGE  0                                 
SCF       %(rSCF)2.1f  0 100 0.1 15                                   
RPATH     -1 
    
""" % vars())
    
    Sca =  labels[ca]      
    el = Element(Sca); d = el.data; Zca = d['Atomic no']    
    
    f.write("""POTENTIALS
*   ipot   Z      element   l_scmt   l_fms   stoichiometry
    0      %(Zca)i     %(Sca)s        -1       -1      0.001 """ % vars())     
    
    
    for i in range(len(labels_short)):
        n = (i+1)
        s = labels_short[i]
        el = Element(s); d = el.data; z = d['Atomic no'] 
        st = list(natoms.values())[i]
        f.write("""
    %(n)i      %(z)i     %(s)s        -1       -1      %(st)i """ % vars())         
        
    f.write("\n \n")  
    
    f.write("ATOMS\n")
    #f.write("      -0.000000     0.000000     0.000000     0    0.0\n")    
    for i in range(len(positions_modified)):      
        f.write('  %13.6f%13.6f%13.6f   %3d   %6.3f\n' % (positions_modified[i][0], positions_modified[i][1], positions_modified[i][2], index[i], dis_array[i]))
    f.write("END\n")    
    f.close()       
    
    for_dist_plot = []
    for i in range(len(positions_modified)):
        for_dist_plot.append([dis_array[i],index[i]])
    return for_dist_plot, labels_short



#from pymatgen.io.xyz import XYZ
if len(filename)==1 or type(filename)==str:
    ind=0
#scratchdir="scratch/.."
labels,natoms,positions=read_xyz(filename[ind])
unique_index,formula = equ_sites(ca,labels,natoms,positions,10) #10 is cutoff distance for equivalent site calculation.
#print("===================================================\n")
#print(formula)
#print("\n===================================================\n")


all_list=list(unique_index.keys())


# do feff calculation for inequivalent ions   
occupancy=[]
xmus_calculated = []
atoms_calculated = []
print("=========================Start calculation================================\n")
print(f"calculating {filename[ind][6:-4]}           Number of inequivalent sites: {len(all_list)}\n")
start_time=time.time()
for i,(k,v) in enumerate(unique_index.items()):       
    print(f"=========================site.{k+1}===============================\n")
    f = f'{filename[ind][6:-4]}/at-{k+1:03d}'
    os.makedirs(f,exist_ok=True)
    os.chdir(f)
    #print(os.getcwd())
    for_dist_plot, labels_short = write_feffinp(ca=k,dmax=dmax,rFMS=rFMS,rSCF=rSCF,corehole=corehole,edge=edge,S02=S02,filename="../../"+filename[ind])
    print('running feff calculation at '+os.getcwd())
    feff_runner(ncores,hostname,host_ind)
    xmu = np.loadtxt('xmu.dat', unpack=True, comments='#', usecols=(0,3), skiprows=0)
    xmus_calculated.append(xmu)    
    atoms_calculated.append([for_dist_plot, labels_short, k+1])        
    occupancy.append(len(v))
    os.chdir('../..') 


#remapping spectra into the same grids   
# find boundaries
e_lowers =[]
e_uppers =[]
for i in xmus_calculated:
    e_lowers.append(i[0][0])
    e_uppers.append(i[0][-1]) 
e_lower = max(e_lowers)
e_upper = min(e_uppers)
e_int = np.arange(e_lower, e_upper, 0.1)  


# interpolate xmus
#to write in json the value in list need to be converted to list
#to write in parquet the value in list need to be converted to dict
xmus_calculated_int = dict()
xmus_calculated_int["E"]=e_int.tolist()
j=0
for d in xmus_calculated:
    f = InterpolatedUnivariateSpline(d[0],d[1])
    d_int = f(e_int)
    xmus_calculated_int[f"site_{atoms_calculated[j][2]}_mu"]=d_int.tolist()
    j=j+1
    
   
# weight-averaged total XAS
xmu_total = 0
for i in range(len(occupancy)):
    xmu_total = xmu_total + np.array(xmus_calculated_int[f"site_{atoms_calculated[i][2]}_mu"])*occupancy[i] 
xmu_total = xmu_total/sum(occupancy)

xmus_calculated_int["total_mu"]=xmu_total.tolist()    






#============================================================================================================

custom_meta_key=f"{filename[ind][6:-4]}"
#print(xmus_calculated_int)
site_spec=json.dumps(xmus_calculated_int)
if not os.path.exists("spectrum.parquet"):
    total_spectrum=dict()
    total_spectrum["E"]=e_int
    total_spectrum[f"{filename[ind][6:-4]}_total_mu"]=xmu_total
    total_spectrum=pd.DataFrame(total_spectrum)
    table = pa.Table.from_pandas(total_spectrum)
    
    existing_meta = table.schema.metadata
    combined_meta = {
        custom_meta_key.encode() : site_spec.encode(),
        **existing_meta
    }
    table = table.replace_schema_metadata(combined_meta)
    pq.write_table(table, 'spectrum.parquet', compression='GZIP')
else:
    restored_table = pq.read_table('spectrum.parquet')
    restored_df = restored_table.to_pandas()
    print(xmu_total[:len(restored_df)])
    restored_table=restored_table.append_column(f"{filename[ind][6:-4]}_total_mu",pa.array(xmu_total[:len(restored_df)]))
    
    #table = pa.Table.from_pandas(restored_df)
    existing_meta =  restored_table.schema.metadata
    combined_meta = {
        custom_meta_key.encode() :site_spec.encode(),
        **existing_meta
    }
    table = restored_table.replace_schema_metadata(combined_meta)
    pq.write_table(table, 'spectrum.parquet', compression='GZIP')
    
    
# shutil.rmtree(f'feff_{ind}', ignore_errors=True)
shutil.rmtree(f'{filename[ind][6:-4]}',ignore_errors=True)



# out = np.column_stack((e_int,xmu_total))


# np.savetxt(f'{filename[ind][6:-4]}_out/feff.dat',out)
finish_time=time.time()
subprocess.run(["wait"])
print(f"\n running time:{finish_time-start_time}")
print("=========================FINISHED==============\n\n\n")
