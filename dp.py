import numpy as np
import ase
import os
from ase.io import read
from deepmd.calculator import DP
from ase import Atoms
import pymatgen
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm
import json
import gc
from glob import glob
from monty.json import MontyDecoder,MontyEncoder


p=DP(model=glob(r'*.pb')[0])
filename=glob(r'*.json')[0]
with open(filename) as f:
    data=json.loads(f.read(),cls=MontyDecoder)
for i in tqdm(np.arange(len(data))):
    try:
        atoms=AseAtomsAdaptor.get_atoms(data[i]['structure'])
    except:
        atoms=AseAtomsAdaptor.get_atoms(data[i]['calc_properties'][0]['structure'])
    atoms.calc = p
    data[i]["2021"]={}
    data[i]["2021"]["energy"]=atoms.get_potential_energy()
    data[i]["2021"]["forces"]=atoms.get_forces()
    del atoms
    if i%100==0:
        gc.collect()

with open(filename.split('.')[0]+"_2021"+".json","w") as f:
    json.dump(data,f,cls=MontysEncoder)
    
