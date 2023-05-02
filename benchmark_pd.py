import numpy as np
import os
from pyace import PyACECalculator
from tqdm import tqdm
import pandas as pd
import glob

dir_name = os.path.basename(os.path.abspath('.'))
outputname=dir_name.strip('regularization_').replace('order_','')

df = pd.read_pickle("Ti_data_0410+volumrtric.pckl.gzip", compression="gzip")
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
    
df.to_pickle('always_include_+volumetric_'+outputname+'.pckl.gzip', compression='gzip', protocol=4)
