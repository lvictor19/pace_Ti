import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
import scienceplots
from .general import get_pathway_data,get_icsd_ref_energy,scatter_data,plot_data

def split_hcp_omega_path(data_collection:pd.DataFrame):
    '''
    split hcp_omega path to SSNEB relaxed and not SSNEB relaxed
    '''
    ssNEB=data_collection['metadata'].map(lambda x: x['SSNEB']==True)
    nssNEB=data_collection['metadata'].map(lambda x: x['SSNEB']==False)

    return data_collection[ssNEB],data_collection[nssNEB]


def benchmark_hcp_omega_pathway(data_collection: pd.DataFrame):
    """
    analysis for the pathways connecting hcp and omega
    """
    hcp_omega_data=get_pathway_data(data_collection,'hcp_omega')
    ssNEB,static=split_hcp_omega_path(hcp_omega_data)
    hcp_energy_dft = get_icsd_ref_energy(data_collection,'hcp','energy')
    omega_energy_dft = get_icsd_ref_energy(data_collection,'omega','energy')
    pathways=static["metadata"].apply(lambda x: x["pathway_number"])
    keys=list(hcp_omega_data.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    for key in keys:
        hcp_energy_pace=get_icsd_ref_energy(data_collection,'hcp',key)
        omega_energy_pace=get_icsd_ref_energy(data_collection,'omega',key)
        for pathway in pathways:
            ispathway=static['metadata'].map(lambda x: x["pathway_number"]==pathway)
            ispathway_ssNEB=ssNEB['metadata'].map(lambda x: x["pathway_number"]==pathway)
            pathway_static=static[ispathway]
            pathway_ssNEB=ssNEB[ispathway_ssNEB]
            images_static=pathway_static['metadata'].map(lambda x: int(x["image"]))
            images_ssNEB=pathway_ssNEB['metadata'].map(lambda x: int(x["image"]))
            Natom=len(pathway_static.iloc[0]['ase_atoms'])
            energies_dft=np.array(pathway_static['energy'])/Natom
            energies_dft_NEB=np.array(pathway_ssNEB['energy'])/Natom
            energies_pace=np.array(pathway_static[key])/Natom
            energies_pace_NEB=np.array(pathway_ssNEB[key])/Natom
            index_static=np.argsort(images_static)
            index_ssNEB=np.argsort(images_ssNEB)
            
            # reference the energies and add the energies of hcp and omega
            energies_dft=np.insert(np.append(energies_dft[index_static]-omega_energy_dft,hcp_energy_dft-omega_energy_dft),0,0)
            energies_dft_NEB=np.insert(np.append(energies_dft_NEB[index_ssNEB]-omega_energy_dft,hcp_energy_dft-omega_energy_dft),0,0)
            energies_pace=np.insert(np.append(energies_pace[index_static]-omega_energy_pace,hcp_energy_pace-omega_energy_pace),0,0)
            energies_pace_NEB=np.insert(np.append(energies_pace_NEB[index_ssNEB]-omega_energy_pace,hcp_energy_pace-omega_energy_pace),0,0)
            with plt.style.context("science"):
                fig=plt.figure(figsize=(5, 4))
                ax=fig.add_subplot(111)
                ax=scatter_data(ax,energies_dft, c="k", label="dft")
                ax=plot_data(ax,energies_pace, c="k", label="pace")
                ax=scatter_data(ax,energies_dft_NEB, c="r", label="dft_SSNEB")
                ax=plot_data(ax, energies_pace_NEB, c="r", label="pace_SSNEB")
                ax.set_xlabel("Transformation Coordinates")
                ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
                ax.tick_params(axis='both', which='both', direction='in')
                ax.set_xticks([0, 6], ["omega", "hcp"])
                plt.legend()
                plt.savefig("hcp_omega/{0}_s{1}.png".format(pathway,key.split('_')[-1]))
                plt.close()