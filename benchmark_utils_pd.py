import numpy as np
import json
from uberplot import *
import pandas as pd
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.io.phonopy import get_pmg_structure, get_phonopy_structure
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
import scienceplots


def get_reference_energy(dataset:pd.DataFrame)->float:
    '''
    get the reference energy per atom from omega structure
    '''
    isfinal=dataset['calc'].map(lambda x: x=='final')
    isicsd=dataset['perturbation'].map(lambda x: x=='icsd')
    isomega=dataset['metadata'].map(lambda x: 'proto' in x.keys() and x['proto']=='omega')
    data=dataset[isfinal&isicsd&isomega]
    index=data['ase_atoms'].keys()[0]
    return data.loc[index]['energy']/len(data.loc[index]['ase_atoms'])

def get_reference_energy_pace(dataset:pd.DataFrame)->float:
    '''
    get the reference pace energy per atom from omega structure
    '''
    isfinal=dataset['calc'].map(lambda x: x=='final')
    isicsd=dataset['perturbation'].map(lambda x: x=='icsd')
    isomega=dataset['metadata'].map(lambda x: 'proto' in x.keys() and x['proto']=='omega')
    data=dataset[isfinal&isicsd&isomega]
    index=data['ase_atoms'].keys()[0]
    return data.loc[index]['pace_energy']/len(data.loc[index]['ase_atoms'])

def get_bcc_ref(data_collection):
    '''
    get the reference energy for bcc structure for dft calculation and pace prediction
    '''
    isicsd=data_collection['perturbation'].map(lambda x: x=='icsd')
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    isbcc=data_collection['metadata'].map(lambda x: 'proto' in x.keys() and x['proto']=='bcc')
    where=data_collection[isicsd&isfinal&isbcc]['ase_atoms'].keys()[0]
    bccref_dft=data_collection.loc[where]['energy']/len(data_collection.loc[where]['ase_atoms'])
    bccref_pace=data_collection.loc[where]['pace_energy']/len(data_collection.loc[where]['ase_atoms'])

    return bccref_dft,bccref_pace

def sort_list(sorting_list,*lists):
    '''
    sort lists with regards to the sorting_list accordingly
    '''
    for list in lists:
        assert len(sorting_list)==len(list)

    sorting_list=np.array(sorting_list)
    sorting_index=np.argsort(sorting_list)

    sorted_list=sorting_list[sorting_index].tolist()
    yield sorted_list
    for list in lists:
        yield np.array(list)[sorting_index].tolist()

    

def burgers_bain(data_collection:pd.DataFrame,step):
    """
    analysis for Burger_Bains pathways connecting fcc-bcc-hcp
    """
    ispathway=data_collection['perturbation'].map(lambda x: x=='pathways')
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    isBurgers_Bain=data_collection['metadata'].map(lambda x: 'pathway_label' in x.keys() and x['pathway_label']=='Burgers_Bain')

    Burgers_Bain_data=data_collection[ispathway&isfinal&isBurgers_Bain]

    bccref_dft,bccref_pace = get_bcc_ref(data_collection)
    isshuffled=data_collection['metadata'].map(lambda x: 'shuffled' in x.keys() and x['shuffled']==True)
    isnotshuffled=isshuffled.map(lambda x: not x)
    BBs=Burgers_Bain_data[isshuffled]
    BBns=Burgers_Bain_data[isnotshuffled]
    
    shuffled_values=list(BBs['metadata'].map(lambda x: x['value']))
    shuffled_energies_dft=list(np.array(BBs['energy'])/np.array(BBs['ase_atoms'].apply(lambda x: len(x)))-bccref_dft)
    shuffled_energies_pace=list(np.array(BBs['pace_energy'])/np.array(BBs['ase_atoms'].apply(lambda x: len(x)))-bccref_pace)

    unshuffled_values=list(BBns['metadata'].map(lambda x: x['value']))
    unshuffled_energies_dft=list(np.array(BBns['energy'])/np.array(BBns['ase_atoms'].apply(lambda x: len(x)))-bccref_dft)
    unshuffled_energies_pace=list(np.array(BBns['pace_energy'])/np.array(BBns['ase_atoms'].apply(lambda x: len(x)))-bccref_pace)

    unshuffled_values,unshuffled_energies_dft,unshuffled_energies_pace=sort_list(unshuffled_values,unshuffled_energies_dft,unshuffled_energies_pace)
    shuffled_values,shuffled_energies_dft,shuffled_energies_pace=sort_list(shuffled_values,shuffled_energies_dft,shuffled_energies_pace)

    with plt.style.context("science"):
        plt.figure(figsize=(5, 4))
        plt.scatter(
            shuffled_values, shuffled_energies_dft, s=6, c="r", label="shuffled_dft"
        )
        plt.scatter(
            unshuffled_values,
            unshuffled_energies_dft,
            s=6,
            c="b",
            label="unshuffled_dft",
        )
        plt.plot(shuffled_values, shuffled_energies_pace, c="r", label="shuffled_pace")
        plt.plot(
            unshuffled_values, unshuffled_energies_pace, c="b", label="unshuffled_pace"
        )
        plt.xlabel("Strain Magnitude")
        plt.ylabel(r"$\Delta E \enspace (eV/atom)$")
        plt.tick_params(axis="y", direction="in")
        plt.tick_params(axis="x", direction="in")
        plt.legend()
        plt.savefig("Ti_Burgers_strain"+"s"+step+".png")
        plt.close()
