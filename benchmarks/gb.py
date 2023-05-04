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
    ref_elem=gb_data['metadata'].map(lambda x: x['ref_elem'])
    ref_proto_dict={"Mo":"bcc","Ti":"hcp"}
    protos=[ref_proto_dict[elem] for elem in ref_elem]
    planes=gb_data["metadata"].map(lambda x: x['plane'])
    sigmas=gb_data["metadata"].map(lambda x: x['sigma'])
    N_gb = np.array(gb_data["ase_atoms"].map(lambda x: len(x)))
    E_gb = np.array(gb_data["energy"])
    keys=list(data_collection.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    E_icsd = [get_icsd_ref_energy(data_collection,proto,'energy') for proto in protos]
    for key in keys:
        E_gb_pred = gb_data[key]
        E_icsd_pred = [get_icsd_ref_energy(data_collection,proto,key) for proto in protos]
        a = list(gb_data["ase_atoms"].map(lambda x: x.cell[0]))
        b = list(gb_data["ase_atoms"].map(lambda x: x.cell[1]))
        A = np.linalg.norm(np.cross(np.array(a),np.array(b),axisa=1,axisb=1),axis=1)
        deltaE_gb = (E_gb - np.array(E_icsd) * N_gb ) / 2 / A
        deltaEpred_gb = (E_gb_pred - np.array(E_icsd) * N_gb ) / 2 / A
        discrepancy = abs((deltaEpred_gb - deltaE_gb) / deltaE_gb) > 0.15
        charlist=['+' if value else '' for value in discrepancy]
        d = {
            "proto": protos,
            "plane": planes,
            "sigma": sigmas,
            r"$E^{DFT}(eV/Å^{2})$": np.round(np.array(deltaE_gb), 3),
            r"$E^{pred}(eV/Å^{2})$": np.round(np.array(deltaEpred_gb), 3),
            "discrepancy": charlist,
        }
        df = pd.DataFrame(data=d)
        fig = plt.figure(figsize=(5, 6))
        ax = fig.add_subplot(111, frame_on=False)
        ax=gentable(ax,df)
        plt.savefig("gb/Grain_boundary_s{0}.png".format(key.split('_')[-1]), dpi=200)