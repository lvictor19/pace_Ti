import numpy as np
import ase
import os
from pyace import PyACECalculator
from ase import Atoms,Atom
import pymatgen
from tqdm import tqdm
import pandas as pd
import json
import glob
from benchmark_utils import *
import argparse

dir_name = os.path.basename(os.path.abspath('.'))
outputname=dir_name.strip('regularization_').replace('order_','')

dataset=read_json('.','always_include_0410.json')
df=convert_to_pickle(dataset)
potfiles=glob.glob('*ladder_step*')
for potfile in potfiles:
    num=potfile.split('.')[0].split('_')[-1]
    p = PyACECalculator(potfile)
    energies=[]
    forces=[]
    for i in tqdm(np.arange(df.shape[0])):
        atoms=df.iloc[i]['ase_atoms'].copy()
        atoms.calc = p
        energies.append(atoms.get_potential_energy())
        forces.append(atoms.get_forces())
    df['pace_energy_'+outputname+'_'+num]=energies
    df['pace_forces_'+outputname+'_'+num]=forces
    
df.to_pickle('always_include_'+outputname+'.pckl.gzip', compression='gzip', protocol=4)

dataset=read_json('.','not_always_include_0410.json')
df=convert_to_pickle(dataset)
potfiles=glob.glob('*ladder_step*')
for potfile in potfiles:
    num=potfile.split('.')[0].split('_')[-1]
    p = PyACECalculator(potfile)
    energies=[]
    forces=[]
    for i in tqdm(np.arange(df.shape[0])):
        atoms=df.iloc[i]['ase_atoms'].copy()
        atoms.calc = p
        energies.append(atoms.get_potential_energy())
        forces.append(atoms.get_forces())
    df['pace_energy_'+outputname+'_'+num]=energies
    df['pace_forces_'+outputname+'_'+num]=forces
    
df.to_pickle('not_always_include_'+outputname+'.pckl.gzip', compression='gzip', protocol=4)

dataset=read_json('.','selected_0410.json')
df=convert_to_pickle(dataset)
potfiles=glob.glob('*ladder_step*')
for potfile in potfiles:
    num=potfile.split('.')[0].split('_')[-1]
    p = PyACECalculator(potfile)
    energies=[]
    forces=[]
    for i in tqdm(np.arange(df.shape[0])):
        atoms=df.iloc[i]['ase_atoms'].copy()
        atoms.calc = p
        energies.append(atoms.get_potential_energy())
        forces.append(atoms.get_forces())
    df['pace_energy_'+outputname+'_'+num]=energies
    df['pace_forces_'+outputname+'_'+num]=forces
    
df.to_pickle('selected_'+outputname+'.pckl.gzip', compression='gzip', protocol=4)

