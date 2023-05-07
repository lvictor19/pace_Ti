import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas.plotting import table
from .general import get_pert_data, get_icsd_ref_energy_multiple, gentable,get_number_of_atoms,get_area

def surfaces(surface_data: pd.DataFrame):
    surface_data=get_pert_data(surface_data,'surfaces')
    return surface_data

def analysis(surface_data: pd.DataFrame,dataset,key):
    """
    surfaces data analysis
    """
    protos=surface_data['metadata'].map(lambda x: x['proto'])
    planes=surface_data["metadata"].map(lambda x: x['miller'])

    N_sf = get_number_of_atoms(surface_data)
    E_sf=surface_data['energy']
    E_icsd = np.array(get_icsd_ref_energy_multiple(dataset,protos,'energy'))
    E_sf_pred=surface_data[key]
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
