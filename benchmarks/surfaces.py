import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas.plotting import table
from .general import get_pert_data, get_icsd_ref_energy_multiple, gentable

def surfaces(surface_data: pd.DataFrame):
    surface_data=get_pert_data(surface_data,'surfaces')
    return surface_data

def get_E_sf(surface_data: pd.DataFrame,key:str):
    '''
    get surface energy
    '''
    return np.array(surface_data[key])

def get_N_sf(surface_data: pd.DataFrame):
    '''
    get number of atoms in cell
    '''
    return np.array(surface_data["ase_atoms"].map(lambda x: len(x)))


def get_area(surface_data:pd.DataFrame):
    '''
    get area
    '''
    lattice_a=list(surface_data["ase_atoms"].map(lambda x: x.cell[0]))
    lattice_b=list(surface_data["ase_atoms"].map(lambda x: x.cell[1]))
    area=np.linalg.norm(np.cross(np.array(lattice_a),np.array(lattice_b),axisa=1,axisb=1),axis=1)
    return area

def analysis(surface_data: pd.DataFrame,dataset,key):
    """
    surfaces data analysis
    """
    protos=surface_data['metadata'].map(lambda x: x['proto'])
    planes=surface_data["metadata"].map(lambda x: x['miller'])
    N_sf = get_N_sf(surface_data)
    E_sf=get_E_sf(surface_data,'energy')
    E_icsd = np.array(get_icsd_ref_energy_multiple(dataset,protos,'energy'))
    E_sf_pred=get_E_sf(surface_data,key)
    E_icsd_pred=np.array(get_icsd_ref_energy_multiple(dataset,protos,key))
    area=get_area(surface_data)
    deltaE_sf = (E_sf - E_icsd * N_sf ) / 2 / area
    deltaEpred_sf = (E_sf_pred - E_icsd * N_sf ) / 2 / area
    discrepancy = abs((deltaEpred_sf - deltaE_sf) / deltaE_sf) > 0.15
    charlist=['+' if value else '' for value in discrepancy]
    d = {
        "proto": protos,
        "plane": planes,
        r"$E^{DFT}(eV/Å^{2})$": np.round(deltaE_sf, 5),
        r"$E^{pred}(eV/Å^{2})$": np.round(deltaEpred_sf, 5),
        "discrepancy": charlist,
    }
    df = pd.DataFrame(data=d)
    return df
