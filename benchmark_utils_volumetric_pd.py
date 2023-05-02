import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from plotting_settings import colors,markers,renaming

def get_pert_data(data_collection:pd.DataFrame,pertname):
    '''
    get the data for a perturbation
    '''
    ispert=data_collection['perturbation'].map(lambda x: x==pertname)
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    pert_data=data_collection[ispert&isfinal].copy()
    return pert_data

def get_proto_data(data_collection:pd.DataFrame,pertname,protoname):
    '''
    get the data for a prototype
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

def volumetric_deformation_df(data_collection: pd.DataFrame):
    """
    Volumetric deformation energy analysis
    input data_collection::pd.DataFrame
    """
    ref_omega_dft=get_icsd_ref_energy(data_collection,'omega','energy')
    volumetric_data=get_pert_data(data_collection,'volumetric')
    protos=set(volumetric_data['metadata'].map(lambda x: x['proto']))
    keys=list(volumetric_data.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    print(protos)
    for key in keys:
        ref_omega_pace=get_icsd_ref_energy(data_collection,'omega',key)
        i=0
        for proto in protos:
            proto_data=get_proto_data(volumetric_data,'volumetric',proto)
            Natom=np.array(proto_data['ase_atoms'].map(lambda x: len(x)))
            Vs_=np.array(proto_data['metadata'].map(lambda x: x['volume']))/Natom
            Es_=np.array(proto_data['energy'])/Natom
            Epreds_=np.array(proto_data[key])/Natom
            Vs = []
            Es = []
            Epreds = []
            for checkindex in np.arange(len(Vs_)):
                samemag = np.where(
                    abs(Vs_ - Vs_[checkindex]) < 5e-2
                )[0]
                if len(samemag) == 1:
                    Vs.append(Vs_[checkindex])
                    Es.append(Es_[checkindex])
                    Epreds.append(Epreds_[checkindex])
                else:
                    flag = 1
                    for mag in samemag:
                        if Es_[checkindex] > Es_[mag]:
                            flag = 0
                    if flag == 1:
                        Vs.append(Vs_[checkindex])
                        Es.append(Es_[checkindex])
                        Epreds.append(Epreds_[checkindex])
            Vs = np.array(Vs)
            Es = np.array(Es) - ref_omega_dft
            Epreds = np.array(Epreds) - ref_omega_pace
            order = np.argsort(Vs)
            label = proto
            if label in renaming:
                label = renaming[label]
            plt.scatter(
                Vs[order],
                Es[order],
                c=colors[i],
                marker=markers[i],
                s=0.2,
                label=label,
            )
            plt.plot(Vs[order], Epreds[order], c=colors[i])
            i+=1
        plt.xlabel(r"$\overline{V} \enspace (Å^{3}/atom)$")
        plt.ylabel(r"$\Delta E (eV/atom)$")
        plt.legend()
        plt.savefig(os.path.join("Volumetric","Volumetric_deformation_s{0}.png".format(key.split('_')[-1])), dpi=200)
        plt.close()
       