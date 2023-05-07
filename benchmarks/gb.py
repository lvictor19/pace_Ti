import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas.plotting import table
from .general import get_pert_data,get_icsd_ref_energy_multiple,gentable


def grain_boundary(dataset: pd.DataFrame):
    gb_data=get_pert_data(dataset,'gb')
    return gb_data

def get_E_gb(gb_data: pd.DataFrame,key:str):
    '''
    get grain boudary energy
    '''
    return np.array(gb_data[key])

def get_N_gb(gb_data: pd.DataFrame):
    '''
    get number of atoms in cell
    '''
    return np.array(gb_data["ase_atoms"].map(lambda x: len(x)))


def get_area(gb_data:pd.DataFrame):
    '''
    get area
    '''
    lattice_a=list(gb_data["ase_atoms"].map(lambda x: x.cell[0]))
    lattice_b=list(gb_data["ase_atoms"].map(lambda x: x.cell[1]))
    area=np.linalg.norm(np.cross(np.array(lattice_a),np.array(lattice_b),axisa=1,axisb=1),axis=1)
    return area

def analysis(gb_data:pd.DataFrame,dataset:pd.DataFrame,key):
    '''
    grain boudary data analysis
    '''
    protos=gb_data["metadata"].map(lambda x: x['proto'])
    planes=gb_data["metadata"].map(lambda x: x['plane'])
    sigmas=gb_data["metadata"].map(lambda x: x['sigma'])

    E_gb_dft=get_E_gb(gb_data,'energy')
    E_gb_pred=get_E_gb(gb_data,key)
    N_gb = get_N_gb(gb_data)
    E_icsd = np.array(get_icsd_ref_energy_multiple(dataset,protos,'energy'))
    E_icsd_pred = np.array(get_icsd_ref_energy_multiple(dataset,protos,key))
    area=get_area(gb_data)
    deltaE_gb = (E_gb_dft - np.array(E_icsd) * N_gb ) / 2 / area
    deltaEpred_gb = (E_gb_pred - np.array(E_icsd_pred) * N_gb ) / 2 / area
    discrepancy = abs((deltaEpred_gb - deltaE_gb) / deltaE_gb) > 0.15
    charlist=['+' if value else '' for value in discrepancy]

    d = {
        "proto": protos,
        "plane": planes,
        "sigma": sigmas,
        r"$E^{DFT}(eV/Å^{2})$": np.round(deltaE_gb, 5),
        r"$E^{pred}(eV/Å^{2})$": np.round(deltaEpred_gb, 5),
        "discrepancy": charlist,
    }
    df = pd.DataFrame(data=d)

    return df