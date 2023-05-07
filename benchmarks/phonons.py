import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from .plotting_settings import colors,markers
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.io.phonopy import get_pmg_structure, get_phonopy_structure
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from .general import get_pert_data


def calculate_phonons_pd(calc_df: pd.DataFrame,key) -> Phonopy:

    phonopy_supercell = get_phonopy_structure(calc_df.iloc[0]["metadata"]["reference_cell"])
    phonon = Phonopy(
        unitcell=phonopy_supercell,
        primitive_matrix="auto",
    )

    phonon_dataset = {"natom": len(calc_df.iloc[0]["metadata"]["reference_cell"]), "first_atoms": []}
    for index in range(calc_df.shape[0]):
       calc_dataset=calc_df.iloc[index]['metadata']['phonopy']
       calc_dataset["forces"]=np.array(calc_df.iloc[index][key])
       phonon_dataset["first_atoms"].append(calc_dataset)
    phonon.dataset = phonon_dataset
    phonon.produce_force_constants()
    return phonon

def phonons(data_collection: pd.DataFrame):
    """
    analysis for the phonon dispersion curves
    """
    phonon_data=get_pert_data(data_collection,'phonons')
    return phonon_data

def analysis(phonon_data:pd.DataFrame,proto:str,key:str):
    proto_df=phonon_data[phonon_data['metadata'].map(lambda x: x['proto'])==proto]
    phonon=calculate_phonons_pd(proto_df,key=key)
    band_fig = phonon.auto_band_structure(plot=True)
    for ax in band_fig.gcf().axes:
        ax.tick_params(axis="both", labelsize=12.0)
        for line_idx in range(len(ax.get_lines())):
            ax.get_lines()[line_idx].set_color("black")
    band_fig.gcf().axes[0].set_ylabel(r"Frequency (THz)", fontsize=14.0)
    return band_fig
