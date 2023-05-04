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

def phonons(data_collection: list):
    """
    analysis for the phonon dispersion curves
    """
    phonon_data=get_pert_data(data_collection,'phonons')
    protos=set(phonon_data['metadata'].map(lambda x: x['proto']))
    keys=list(phonon_data.keys())
    keys=[x for x in keys if x.startswith('pace_forces')]
    for proto in protos:
        proto_df=phonon_data[phonon_data['metadata'].map(lambda x: x['proto'])==proto]
        phonon=calculate_phonons_pd(proto_df,key='forces')
        band_fig = phonon.auto_band_structure(plot=True)
        for ax in band_fig.gcf().axes:
            ax.tick_params(axis="both", labelsize=12.0)
            for line_idx in range(len(ax.get_lines())):
                ax.get_lines()[line_idx].set_color("black")
        band_fig.gcf().axes[0].set_ylabel(r"Frequency (THz)", fontsize=14.0)
        band_fig.savefig(os.path.join("Phonons","phonon_{0}.png".format(proto)), transparent=True)
        plt.close()

        for key in keys:
            proto_df=phonon_data[phonon_data['metadata'].map(lambda x: x['proto'])==proto]
            phonon=calculate_phonons_pd(proto_df,key=key)
            band_fig = phonon.auto_band_structure(plot=True)
            for ax in band_fig.gcf().axes:
                ax.tick_params(axis="both", labelsize=12.0)
                for line_idx in range(len(ax.get_lines())):
                    ax.get_lines()[line_idx].set_color("black")
            band_fig.gcf().axes[0].set_ylabel(r"Frequency (THz)", fontsize=14.0)
            band_fig.savefig(os.path.join("Phonons","phonon_{0}_s{1}.png".format(proto,key.split('_')[-1])), transparent=True)
            plt.close()