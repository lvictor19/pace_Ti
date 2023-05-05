import pandas as pd
from pymatgen.io.ase import AseAtomsAdaptor


def convert_to_pickle(dataset,pace=False):
    '''
    create pickle file for a json file
    Parameters
    ----------
    dataset : json format
    ref_energy : reference energy
    output : output path
    '''
    
    energies=[]
    forces=[]
    ase_atoms=[]
    calc=[]
    metadata=[]
    perturbation=[]
    if pace:
        pace_energy=[]
        pace_forces=[]
    for data in dataset:
        energies.append(data['energy'])
        forces.append(data['forces'])
        ase_atoms.append(AseAtomsAdaptor.get_atoms(data['structure']))
        calc.append(data['calc'])
        perturbation.append(data['metadata']['perturbation'])
        try:
            metadata.append(data['metadata'])
        except:
            datametadata_new=data['metadata'].copy()
            datametadata_new['reference_cell']=datametadata_new['reference_cell'].as_dict()
            metadata.append(datametadata_new)
        if pace:
            pace_energy.append(data['pace']['energy'])
            pace_forces.append(data['pace']['forces'])
    datadict={
    'energy':energies,
    'forces':forces,
    'ase_atoms':ase_atoms,
    'calc':calc,
    'perturbation':perturbation,
    'metadata':metadata,
    }
    if pace:
        datadict['pace_energy']=pace_energy
        datadict['pace_forces']=pace_forces
    df=pd.DataFrame(data=datadict)
    return df
