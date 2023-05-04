import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
from .general import get_icsd_ref_energy,reference_energies,get_pathway_data,unroll_pathway_metadata,plot_pathway_data
import scienceplots


def split_Burgers_Bain_path(data_collection:pd.DataFrame):
    '''
    split Burgers_Bain path to shuffled and unshuffled
    '''
    BBshuffled=data_collection[data_collection['pathway_number']==1]
    BBunshuffled=data_collection[data_collection['pathway_number']==0]

    return BBshuffled,BBunshuffled

def plot_burgers_bain(ax,shuffled_referenced:pd.DataFrame,unshuffled_referenced:pd.DataFrame,key):
    '''
    plot burger_bain pathway energy profile
    key : key in the dataframe for the pace energy
    '''
    shuffled_referenced.sort_values(by='value', inplace=True)
    unshuffled_referenced.sort_values(by='value', inplace=True)
    plotlists=[]
    plotlists.append(('scatter',np.array(shuffled_referenced['value']),np.array(shuffled_referenced['energy'])/np.array(shuffled_referenced['ase_atoms'].map(lambda x:len(x))), {"s":6, "c":"r", "label":"shuffled_dft"}))
    plotlists.append(('scatter',np.array(unshuffled_referenced['value']),np.array(unshuffled_referenced['energy'])/np.array(unshuffled_referenced['ase_atoms'].map(lambda x:len(x))), {"s":6, "c":"b", "label":"unshuffled_dft"}))
    plotlists.append(('plot',np.array(shuffled_referenced['value']), np.array(shuffled_referenced[key])/np.array(shuffled_referenced['ase_atoms'].map(lambda x:len(x))), {"c":"r", "label":"shuffled_pace"}))
    plotlists.append(('plot',np.array(unshuffled_referenced['value']),np.array(unshuffled_referenced[key])/np.array(unshuffled_referenced['ase_atoms'].map(lambda x:len(x))), {"c":"b", "label":"unshuffled_pace"}))
    ax=plot_pathway_data(ax,plotlists)
    return ax

def benchmark_bain_burgers_pathway(data_collection:pd.DataFrame):
    """
    analysis for Burger_Bains pathways connecting fcc-bcc-hcp
    """
    Burgers_Bain_data=unroll_pathway_metadata(get_pathway_data(data_collection,'Burgers_Bain'))
    bccref_dft=get_icsd_ref_energy(data_collection,'bcc','energy')
    keys=list(data_collection.keys())
    ref_dict={'energy':bccref_dft}
    for key in [x for x in keys if x.startswith('pace_energy')]:
        bccref_pace=get_icsd_ref_energy(data_collection,'bcc',key)
        ref_dict[key]=bccref_pace

    Burgers_Bain_data=reference_energies(Burgers_Bain_data,ref_dict)
    shuffled,unshuffled=split_Burgers_Bain_path(Burgers_Bain_data)
    return shuffled,unshuffled

