import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas.plotting import table
from .general import get_pert_data,get_icsd_ref_energy,gentable


def grain_boundary(data_collection: pd.DataFrame):
    """
    Grain boundary energy analysis
    input data_collection::list
    output a picture of the table of the grain boundary energies
    """
    gb_data=get_pert_data(data_collection,'gb')
    return gb_data

def get_deltaE(row,data_collection,key):
    '''
    get grain boundary energy
    '''
    a=row["ase_atoms"].cell[0]
    b=row["ase_atoms"].cell[1]
    N_gb=len(row["ase_atoms"])
    A = np.linalg.norm(np.cross(a,b))
    E_gb=row[key]
    ref_elem=row['metadata']['ref_elem']
    ref_proto_dict={"Mo":"bcc","Ti":"hcp"}
    proto=ref_proto_dict[ref_elem]['ref_elem']
    E_icsd=get_icsd_ref_energy(data_collection,proto,'energy')
    deltaE=(E_gb - E_icsd * N_gb ) / 2 / A
    return deltaE

def analysis(gb_data:pd.DataFrame,data_collection:pd.DataFrame,key):
    '''
    grain boudary data analysis
    '''
    ref_elem=gb_data['metadata'].map(lambda x: x['ref_elem'])
    ref_proto_dict={"Mo":"bcc","Ti":"hcp"}
    protos=[ref_proto_dict[elem] for elem in ref_elem]
    planes=gb_data["metadata"].map(lambda x: x['plane'])
    sigmas=gb_data["metadata"].map(lambda x: x['sigma'])

    gb_data['deltaE_gb']=gb_data.apply(lambda x: get_deltaE(x,data_collection,'energy'),axis=1)
    gb_data['deltaEpred_gb']=gb_data.apply(lambda x: get_deltaE(x,data_collection,key),axis=1)

    Edft=np.array(gb_data['deltaE_gb'])
    Epred=np.array(gb_data['deltaEpred_gb'])

    discrepancy = abs((Epred - Edft) / Edft) > 0.15
    charlist=['+' if value else '' for value in discrepancy]
    d = {
        "proto": protos,
        "plane": planes,
        "sigma": sigmas,
        r"$E^{DFT}(eV/Å^{2})$": np.round(Edft, 3),
        r"$E^{pred}(eV/Å^{2})$": np.round(Epred, 3),
        "discrepancy": charlist,
    }
    df = pd.DataFrame(data=d)

    return df