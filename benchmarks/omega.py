import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
import scienceplots
from .general import get_icsd_ref_energy, get_pathway_data, reference_energies,unroll_pathway_metadata,plot_pathway_data


def plot_omega(ax,bcc_omega_data:pd.DataFrame,key):
    '''
    plot burger_bain pathway energy profile
    key : key in the dataframe for the pace energy
    '''
    bcc_omega_data.sort_values(by='value', inplace=True)
    plotlists=[]
    plotlists.append(('scatter',np.array(bcc_omega_data['value']),np.array(bcc_omega_data['energy'])/np.array(bcc_omega_data['ase_atoms'].map(lambda x:len(x))), {"c":"black","s":4,"label":"dft"}))
    plotlists.append(('plot',np.array(bcc_omega_data['value']), np.array(bcc_omega_data[key])/np.array(bcc_omega_data['ase_atoms'].map(lambda x:len(x))), {"c":"black","label":"pace"}))
    ax=plot_pathway_data(ax,plotlists)
    return ax

def benchmark_bcc_omega_pathway(data_collection: pd.DataFrame):
    """
    analysis for the pathway connecting bcc and omega
    ref_omega_dft::float   reference energy (per atom) for omega structure from dft calculation
    ref_omega_pace::float  reference energy (per atom) for omega structure from pace potential
    """
    bcc_omega_data=unroll_pathway_metadata(get_pathway_data(data_collection,'bcc_omega'))
    ref_omega_dft=get_icsd_ref_energy(data_collection,'omega','energy')
    keys=list(bcc_omega_data.keys())
    ref_dict={'energy':ref_omega_dft}
    for key in [x for x in keys if x.startswith('pace_energy')]:
        omegaref_pace=get_icsd_ref_energy(data_collection,'omega',key)
        ref_dict[key]=omegaref_pace
    bcc_omega_data=reference_energies(bcc_omega_data,ref_dict)
    
    return bcc_omega_data
