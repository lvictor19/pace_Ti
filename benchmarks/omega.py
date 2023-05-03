import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
import scienceplots
from general import get_icsd_ref_energy, get_pathway_data, reference_energies,scatter_data,plot_data


def benchmark_bcc_omega_pathway(data_collection: pd.DataFrame):
    """
    analysis for the pathway connecting bcc and omega
    ref_omega_dft::float   reference energy (per atom) for omega structure from dft calculation
    ref_omega_pace::float  reference energy (per atom) for omega structure from pace potential
    """
    bcc_omega_data=get_pathway_data(data_collection,'bcc_omega')
    values=np.array(bcc_omega_data['metadata'].apply(lambda x: x['value']))
    Natom=np.array(bcc_omega_data['ase_atoms'].map(lambda x: len(x)))
    keys=list(bcc_omega_data.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    Energies_dft=np.array(bcc_omega_data['energy'])/Natom
    index = np.argsort(values)
    values = values[index]
    Energies_dft = Energies_dft[index]
    ref_omega_dft=get_icsd_ref_energy(data_collection,'omega','energy')
    print(keys)
    for key in keys:
        Energies_pace=np.array(bcc_omega_data[key])/Natom
        Energies_pace = Energies_pace[index]
        ref_omega_pace=get_icsd_ref_energy(data_collection,'omega',key)
        with plt.style.context("science"):
            fig=plt.figure(figsize=(5, 4))
            ax=fig.add_subplot(111)
            ax=scatter_data(ax,Energies_dft-ref_omega_dft,c="black",s=4,label="dft")
            ax=plot_data(ax,Energies_pace-ref_omega_pace,c="black",label="pace")
            ax.set_xlabel("Omega Transformation")
            ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
            ax.tick_params(axis='both', which='both', direction='in')
            ax.set_xticks([])
            plt.legend()
            plt.savefig("Omega/Omegas{0}.png".format(key.split('_')[-1]), dpi=150)
            plt.close()