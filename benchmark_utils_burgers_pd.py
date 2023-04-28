import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
import scienceplots

def get_icsd_ref_energy(dataset:pd.DataFrame,protoname:str,which)->float:
    '''
    get the reference energy per atom for a prototype
    '''
    isfinal=dataset['calc'].map(lambda x: x=='final')
    isicsd=dataset['perturbation'].map(lambda x: x=='icsd')
    isproto=dataset['metadata'].map(lambda x: 'proto' in x.keys() and x['proto']==protoname)
    data=dataset[isfinal&isicsd&isproto]
    index=data['ase_atoms'].keys()[0]
    return data.loc[index][which]/len(data.loc[index]['ase_atoms'])

def get_pathway_data(data_collection:pd.DataFrame,pathwaylabel):
    '''
    get pathway data in a pandas dataframe
    '''
    ispathway=data_collection['perturbation'].map(lambda x: x=='pathways')
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    isBurgers_Bain=data_collection['metadata'].map(lambda x: 'pathway_label' in x.keys() and x['pathway_label']==pathwaylabel)
    Burgers_Bain_data=data_collection[ispathway&isfinal&isBurgers_Bain].copy()
    return Burgers_Bain_data

def split_Burgers_Bain_path(data_collection:pd.DataFrame):
    '''
    split Burgers_Bain path to shuffled and unshuffled
    '''
    isshuffled=data_collection['metadata'].map(lambda x: 'shuffled' in x.keys() and x['shuffled']==True)
    isnotshuffled=isshuffled.map(lambda x: not x)
    BBshuffled=data_collection[isshuffled]
    BBunshuffled=data_collection[isnotshuffled]

    return BBshuffled,BBunshuffled

def sort_list(sorting_list,*lists):
    '''
    sort *lists with sorting index identical to sorting the sorting_list
    '''
    for list in lists:
        assert len(sorting_list)==len(list)

    sorting_list=np.array(sorting_list)
    sorting_index=np.argsort(sorting_list)

    sorted_list=sorting_list[sorting_index].tolist()
    yield sorted_list
    for list in lists:
        yield np.array(list)[sorting_index].tolist()

def reference_energies(data_collection:pd.DataFrame,reference_energies:dict):
    '''
    reference the energies in a dataframe according to a list of reference energies
    '''
    for ref_e in reference_energies.keys():
        data_collection[ref_e]=data_collection[ref_e]/data_collection['ase_atoms'].map(lambda x:len(x))-reference_energies[ref_e]
    return data_collection

def plot_burgers_bain(ax,shuffled_referenced:pd.DataFrame,unshuffled_referenced:pd.DataFrame,key):
    '''
    plot burger_bain pathway energy profile
    key : key in the dataframe for the pace energy
    '''
    shuffled_referenced['value']=shuffled_referenced['metadata'].map(lambda x: x['value'])
    unshuffled_referenced['value']=unshuffled_referenced['metadata'].map(lambda x: x['value'])
    shuffled_referenced.sort_values(by='value', inplace=True)
    unshuffled_referenced.sort_values(by='value', inplace=True)

    ax.scatter(
        shuffled_referenced['value'], shuffled_referenced['energy'], s=6, c="r", label="shuffled_dft"
    )
    ax.scatter(
        unshuffled_referenced['value'],
        unshuffled_referenced['energy'],
        s=6,
        c="b",
        label="unshuffled_dft",
    )
    ax.plot(shuffled_referenced['value'], shuffled_referenced[key], c="r", label="shuffled_pace")
    ax.plot(
        unshuffled_referenced['value'],unshuffled_referenced[key], c="b", label="unshuffled_pace"
    )
    return ax

def benchmark_bain_burgers_pathway(data_collection:pd.DataFrame):
    """
    analysis for Burger_Bains pathways connecting fcc-bcc-hcp
    """
    Burgers_Bain_data=get_pathway_data(data_collection,'Burgers_Bain')
    bccref_dft=get_icsd_ref_energy(data_collection,'bcc','energy')
    keys=list(data_collection.keys())
    ref_dict={'energy':bccref_dft}
    for key in [x for x in keys if x.startswith('pace_energy')]:
        bccref_pace=get_icsd_ref_energy(data_collection,'bcc',key)
        ref_dict[key]=bccref_pace

    Burgers_Bain_data=reference_energies(Burgers_Bain_data,ref_dict)
    shuffled,unshuffled=split_Burgers_Bain_path(Burgers_Bain_data)
    return shuffled,unshuffled

