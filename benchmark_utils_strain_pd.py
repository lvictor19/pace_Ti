import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from plotting_settings import colors,markers

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

def strain(data_collection: pd.DataFrame):
    """
    Strain-Energy relations
    input
        data_collection::list
        ref_omega_dft::float   reference energy (per atom) for omega structure from dft calculation
        ref_omega_pace::float  reference energy (per atom) for omega structure from pace potential
    output pictures of the Energy-Strain curves (a picture contains different strain directions for the same prototype)
    """
    ref_omega_dft=get_icsd_ref_energy(data_collection,'omega','energy')
    strain_data=get_pert_data(data_collection,'strain')
    protos=set(strain_data['metadata'].map(lambda x: x['proto']))
    keys=list(data_collection.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    for proto in protos:
        strain_proto_data=get_proto_data(strain_data,'strain',proto)
        Es=np.array(list(strain_proto_data['metadata'].map(lambda x : x['special_direction'])))
        groups = []
        already = np.zeros(len(Es))
        isidentical = abs(np.dot(Es, Es.T) - 1) < 1e-12
        for i in np.arange(len(Es)):
            if not already[i] == 1:
                thisgroup = np.where(isidentical[i] == True)
                groups.append(thisgroup[0].tolist())
                already[thisgroup] = 1
        for key in keys:
            i=0
            for group in groups:
                ref_omega_pace=get_icsd_ref_energy(data_collection,'omega',key)
                Natom=np.array(strain_proto_data.iloc[group]['ase_atoms'].map(lambda x: len(x)))
                magnitudes_=np.array(strain_proto_data.iloc[group]['metadata'].map(lambda x: x['magnitude']))
                Energies_=np.array(strain_proto_data.iloc[group]['energy'])/Natom
                Energies_pred_=np.array(strain_proto_data.iloc[group][key])/Natom
                magnitudes = []
                Energies = []
                Energies_pred = []
                for checkindex in np.arange(len(magnitudes_)):
                    samemag = np.where(abs(magnitudes_ - magnitudes_[checkindex]) < 1e-2)[0]
                    if len(samemag) == 1:
                        magnitudes.append(magnitudes_[checkindex])
                        Energies.append(Energies_[checkindex])
                        Energies_pred.append(Energies_pred_[checkindex])
                    else:
                        flag = 1
                        for mag in samemag:
                            if Energies_[checkindex] > Energies_[mag]:
                                flag = 0
                        if flag == 1:
                            magnitudes.append(magnitudes_[checkindex])
                            Energies.append(Energies_[checkindex])
                            Energies_pred.append(Energies_pred_[checkindex])
                magnitudes = np.array(magnitudes)
                Energies = np.array(Energies)
                Energies_pred = np.array(Energies_pred)
                order = magnitudes.argsort()
                plt.scatter(
                    magnitudes[order],
                    Energies[order] - ref_omega_dft,
                    c=colors[i],
                    marker=markers[i],
                    s=2
                )
                plt.plot(
                    magnitudes[order],
                    Energies_pred[order] - ref_omega_pace,
                    c=colors[i],
                )
                plt.xlabel("Strain Magnitude")
                plt.ylabel(r"$\Delta E \enspace (eV/atom)$")
                i += 1
                
            plt.savefig(os.path.join("Strains","Strains_" + proto +"s"+key.split('_')[-1]+".png"), dpi=200)
            plt.close()
