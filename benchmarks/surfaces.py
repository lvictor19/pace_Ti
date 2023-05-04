import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas.plotting import table
from .general import get_pert_data, get_icsd_ref_energy, gentable

def savetable(data: dict, name: str, dpi: int):
    """
    input a dictionary of a table; the name and dpi of output picture
    output a picture of the table
    """
    df = pd.DataFrame(data=data)
    fig = plt.figure(figsize=(5, 6))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    table(ax, df, loc="center")
    plt.savefig(name, dpi=dpi)

def surfaces(data_collection: pd.DataFrame):
    """
    Surface boundary energy analysis
    input data_collection::list
    output a picture of the table of the surface energies
    """
    surface_data=get_pert_data(data_collection,'surfaces')
    protos=surface_data['metadata'].map(lambda x: x['proto'])
    planes=surface_data["metadata"].map(lambda x: x['miller'])
    N_sf = np.array(surface_data["ase_atoms"].map(lambda x: len(x)))
    E_sf = np.array(surface_data["energy"])
    keys=list(data_collection.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    E_icsd = [get_icsd_ref_energy(data_collection,proto,'energy') for proto in protos]
    for key in keys:
        E_sf_pred = surface_data[key]
        E_icsd_pred = [get_icsd_ref_energy(data_collection,proto,key) for proto in protos]
        a = list(surface_data["ase_atoms"].map(lambda x: x.cell[0]))
        b = list(surface_data["ase_atoms"].map(lambda x: x.cell[1]))
        A = np.linalg.norm(np.cross(np.array(a),np.array(b),axisa=1,axisb=1),axis=1)
        deltaE_sf = (E_sf - np.array(E_icsd) * N_sf ) / 2 / A
        deltaEpred_sf = (E_sf_pred - np.array(E_icsd) * N_sf ) / 2 / A
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
        fig = plt.figure(figsize=(5, 6))
        ax = fig.add_subplot(111, frame_on=False)
        ax=gentable(ax,df)
        plt.savefig("Surfaces/Surfaces_s{0}.png".format(key.split('_')[-1]), dpi=200)
        plt.close()