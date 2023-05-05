import numpy as np
import json
import pandas as pd
from matplotlib import pyplot as plt
import scienceplots
from .general import get_pathway_data,get_icsd_ref_energy,unroll_pathway_metadata,plot_pathway_data

def split_hcp_omega_path(data_collection:pd.DataFrame):
    '''
    split hcp_omega path to SSNEB relaxed and not SSNEB relaxed
    '''
    ssNEB=data_collection[(data_collection['SSNEB']==True)&(data_collection['calc']=='final')]
    nssNEB=data_collection[(data_collection['SSNEB']==False)&(data_collection['calc']=='final')]

    return ssNEB,nssNEB

def get_hcp_omega_reference_energy(data_collection: pd.DataFrame,protoname:str,key):
    '''
    get the reference energy for hcp or omega for the hcp-omega pathway
    '''
    return data_collection[data_collection['proto']==protoname].iloc[0][key]

def plot_hcp_omega(ax,ssNEB:pd.DataFrame,static:pd.DataFrame,pathway:str,key:str):
    pathway_ssNEB=ssNEB[(ssNEB['pathway_number']==pathway)&(ssNEB['value']!=0)]
    pathway_static=static[(static['pathway_number']==pathway)&(static['value']!=0)]
    images_static=np.array(pathway_static['value']).astype(int)
    images_ssNEB=np.array(pathway_ssNEB['value']).astype(int)    
    Natom=len(pathway_static.iloc[0]['ase_atoms'])
    energies_dft=np.array(pathway_static['energy'])/Natom
    energies_dft_NEB=np.array(pathway_ssNEB['energy'])/Natom
    energies_pace=np.array(pathway_static[key])/Natom
    energies_pace_NEB=np.array(pathway_ssNEB[key])/Natom

    hcp_energy_dft=get_hcp_omega_reference_energy(static,'hcp','energy')/Natom
    omega_energy_dft=get_hcp_omega_reference_energy(static,'omega','energy')/Natom
    hcp_energy_pace=get_hcp_omega_reference_energy(static,'hcp',key)/Natom
    omega_energy_pace=get_hcp_omega_reference_energy(static,'omega',key)/Natom

    index_static=np.argsort(images_static)
    index_ssNEB=np.argsort(images_ssNEB)
    # reference the energies and add the energies of hcp and omega
    energies_dft=np.insert(np.append(energies_dft[index_static]-omega_energy_dft,hcp_energy_dft-omega_energy_dft),0,0)
    energies_dft_NEB=np.insert(np.append(energies_dft_NEB[index_ssNEB]-omega_energy_dft,hcp_energy_dft-omega_energy_dft),0,0)
    energies_pace=np.insert(np.append(energies_pace[index_static]-omega_energy_pace,hcp_energy_pace-omega_energy_pace),0,0)
    energies_pace_NEB=np.insert(np.append(energies_pace_NEB[index_ssNEB]-omega_energy_pace,hcp_energy_pace-omega_energy_pace),0,0)

    index_static=np.insert(np.append(images_static[index_static],6),0,0)
    index_ssNEB=np.insert(np.append(images_ssNEB[index_ssNEB],6),0,0)

    plotlists=[]
    plotlists.append(('scatter',index_static,energies_dft,{"c":"k", "label":"dft"}))
    plotlists.append(('plot',index_static,energies_pace,{"c":"k", "label":"pace"}))
    plotlists.append(('scatter',index_ssNEB,energies_dft_NEB,{"c":"r", "label":"dft_SSNEB"}))
    plotlists.append(('plot',index_ssNEB,energies_pace_NEB,{"c":"r", "label":"pace_SSNEB"}))
    ax=plot_pathway_data(ax,plotlists)
    
    return ax

def benchmark_hcp_omega_pathway(data_collection: pd.DataFrame):
    """
    analysis for the pathways connecting hcp and omega
    """
    hcp_omega_data=unroll_pathway_metadata(get_pathway_data(data_collection,'hcp_omega'))
    ssNEB,static=split_hcp_omega_path(hcp_omega_data)
    return ssNEB,static