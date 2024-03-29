import pandas as pd
import numpy as np
from pandas.plotting import table

def get_icsd_ref_energy(dataset:pd.DataFrame,protoname:str,key:str)->float:
    '''
    get the reference energy per atom for a prototype
    '''
    isfinal=dataset['calc'].map(lambda x: x=='final')
    isicsd=dataset['perturbation'].map(lambda x: x=='icsd')
    isproto=dataset['metadata'].map(lambda x: 'proto' in x.keys() and x['proto']==protoname)
    data=dataset[isfinal&isicsd&isproto]
    index=data['ase_atoms'].keys()[0]
    return data.loc[index][key]/len(data.loc[index]['ase_atoms'])

def get_number_of_atoms(dataset: pd.DataFrame):
    '''
    get number of atoms in cell
    '''
    return dataset["ase_atoms"].map(lambda x: len(x))

def get_icsd_ref_energy_multiple(dataset:pd.DataFrame,protonames:pd.Series,key:str)->pd.Series:
    '''
    get the reference energy per atom for a prototype
    '''
    isfinal=dataset['calc'].map(lambda x: x=='final')
    isicsd=dataset['perturbation'].map(lambda x: x=='icsd')
    data=dataset[isfinal&isicsd]
    data['proto']=data['metadata'].map(lambda x: x['proto'])
    filtered_df = data[data['proto'].isin(protonames)]
    print(protonames)
    ordered_df = protonames.to_frame(name='proto').merge(data, on='proto', how='left')
    print(ordered_df)
    return ordered_df[key]/ordered_df['ase_atoms'].map(lambda x: len(x))

def get_area(dataset:pd.DataFrame):
    '''
    get area
    '''
    lattice_a=list(dataset["ase_atoms"].map(lambda x: x.cell[0]))
    lattice_b=list(dataset["ase_atoms"].map(lambda x: x.cell[1]))
    area=np.linalg.norm(np.cross(np.array(lattice_a),np.array(lattice_b),axisa=1,axisb=1),axis=1)
    return area

def reference_energies(data_collection:pd.DataFrame,reference_energies:dict):
    '''
    reference the energies in a dataframe according to a list of reference energies
    '''
    for ref_e in reference_energies.keys():
        data_collection[ref_e]=data_collection[ref_e]-data_collection['ase_atoms'].map(lambda x:len(x))*reference_energies[ref_e]
    return data_collection

def get_pert_data(data_collection:pd.DataFrame,pertname:str,is_final=True):
    '''
    get the data for a perturbation
    '''
    ispert=data_collection['perturbation'].map(lambda x: x==pertname)
    if is_final:
        isfinal=data_collection['calc'].map(lambda x: x=='final')
    else:
        isfinal=data_collection['calc'].map(lambda x: True)
    pert_data=data_collection[ispert&isfinal].copy()
    return pert_data

def unroll_pathway_metadata(data_collection:pd.DataFrame):
    '''
    unroll metadata and write in the dataframe in seperate columns
    '''
    if not data_collection['perturbation'].eq('pathways').all():
        raise ValueError('Perturbations of input data are not all pathways')
    for key in data_collection.iloc[0]['metadata']:
        data_collection[key]=data_collection['metadata'].map(lambda x: x[key])
    return data_collection

def get_proto_data(data_collection:pd.DataFrame,protoname):
    '''
    get the data for a prototype
    '''
    isproto=data_collection["metadata"].map(lambda x: x["proto"]==protoname)
    strain_data=data_collection[isproto].copy()
    return strain_data

def get_pathway_data(data_collection:pd.DataFrame,pathwaylabel):
    '''
    get pathway data in a pandas dataframe
    '''
    ispathway=data_collection['perturbation'].map(lambda x: x=='pathways')
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    ispathway=data_collection['metadata'].map(lambda x: 'pathway_label' in x.keys() and x['pathway_label']==pathwaylabel)
    pathway_data=data_collection[ispathway&isfinal&ispathway].copy()
    return pathway_data

def gentable(ax,data: dict):
    """
    input a dictionary of a table; the name and dpi of output picture
    output a picture of the table
    """
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax=table(ax, data, loc="center")
    return ax

def plot_pathway_data(ax,plotlist):
    '''
    plot pathway data,plotlist contains tuples in which the entry is either 'scatter' or 'plot',
    the second and third entries are the x-y values, the third entry is a dictionary of plotting args
    '''
    for plot in  plotlist:
        if plot[0]=='scatter':
            ax.scatter(plot[1],plot[2],**plot[3])
        if plot[0]=='plot':
            ax.plot(plot[1],plot[2],**plot[3])
    return ax