import pandas as pd
import numpy as np
from pandas.plotting import table

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

def reference_energies(data_collection:pd.DataFrame,reference_energies:dict):
    '''
    reference the energies in a dataframe according to a list of reference energies
    '''
    for ref_e in reference_energies.keys():
        data_collection[ref_e]=data_collection[ref_e]/data_collection['ase_atoms'].map(lambda x:len(x))-reference_energies[ref_e]
    return data_collection

def get_pert_data(data_collection:pd.DataFrame,pertname):
    '''
    get the strain data
    '''
    ispert=data_collection['perturbation'].map(lambda x: x==pertname)
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    pert_data=data_collection[ispert&isfinal].copy()
    return pert_data

def get_proto_data(data_collection:pd.DataFrame,pertname,protoname):
    '''
    get the strain data for a prototype
    '''
    if not data_collection['perturbation'].eq(pertname).all():
        raise ValueError('Perturbation of input data not all strain for getting strain data belonging to a prototype')
    isproto=data_collection["metadata"].map(lambda x: x["proto"]==protoname)
    strain_data=data_collection[isproto].copy()
    return strain_data

def get_pathway_data(data_collection:pd.DataFrame,pathwaylabel):
    '''
    get pathway data in a pandas dataframe
    '''
    ispathway=data_collection['perturbation'].map(lambda x: x=='pathways')
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    isBurgers_Bain=data_collection['metadata'].map(lambda x: 'pathway_label' in x.keys() and x['pathway_label']==pathwaylabel)
    Burgers_Bain_data=data_collection[ispathway&isfinal&isBurgers_Bain].copy()
    return Burgers_Bain_data

def gentable(ax,data: dict):
    """
    input a dictionary of a table; the name and dpi of output picture
    output a picture of the table
    """
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax=table(ax, data, loc="center")
    return ax

def scatter_data(ax,Energies_dft:np.ndarray,x=None,**kwargs):
    if x==None:
        x=np.arange(len(Energies_dft))
    ax.scatter(np.arange(len(Energies_dft)),Energies_dft,**kwargs)
    return ax

def plot_data(ax,Energies_pace:np.ndarray,x=None,**kwargs):
    if x==None:
        x=np.arange(len(Energies_pace))
    ax.scatter(np.arange(len(Energies_pace)),Energies_pace,**kwargs)
    return ax
