import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pandas.plotting import table
from plotting_settings import colors,markers,renaming
from general import get_pert_data,get_proto_data,gentable

def loc_proto_ideal(vacancy_data:pd.DataFrame,vacancy_ideal:pd.DataFrame):
    '''
    return the location of a prototype from vacancy_data in vacancy_ideal DataFrame
    '''
    locs=[]
    for index,row in vacancy_data.iterrows():
        proto=row['metadata']['proto']
        print(proto)
        for i,r in vacancy_ideal.iterrows():
            #  print(r['metadata']['proto'])
             if r['metadata']['proto']==proto:
                  locs.append(i)
    print(locs)
    return locs

def vacancy_formation(data_collection: pd.DataFrame):
    """
    Vacancy formation energy analysis
    input data_collection::list
    output a picture of the table of the vacancy formation energies
    """
    vacancy_data=get_pert_data(data_collection,'vacancies')
    vacancy_ideal=get_pert_data(data_collection,'vacancies_ideal')
    protos=list(vacancy_data['metadata'].map(lambda x: x['proto']))
    keys=list(vacancy_data.keys())
    keys=[x for x in keys if x.startswith('pace_energy')]
    for key in keys:
        deltaEs_collect=[]
        deltaEs_pred_collect=[]
        discrepancies_collect=[]
        protos_=[]
        for proto in protos:
            if proto in []:#['casi-alth_3','mp-865373POSCAR']:
                continue
            else:
                protos_.append(proto)
            vacancy_proto_data=get_proto_data(vacancy_data,'vacancies',proto)
            locs=loc_proto_ideal(vacancy_proto_data,vacancy_ideal)
            Evac=np.array(vacancy_proto_data['energy'])
            Eideal=np.array([vacancy_ideal.loc[i]['energy'] for i in locs])
            N_proto = np.array([len(vacancy_ideal.loc[i]['ase_atoms']) for i in locs])
            Evac_pred=np.array(vacancy_proto_data[key])
            Eideal_pred=[vacancy_ideal.loc[i][key] for i in locs]
            deltaE = np.round(Evac - (N_proto - 1) * Eideal / N_proto,3)
            deltaE_pred = np.round(Evac_pred - (N_proto - 1) * Eideal_pred / N_proto,3)
            deltaEs_collect.append(", ".join([str(f) for f in deltaE]))
            deltaEs_pred_collect.append(", ".join([str(f) for f in deltaE_pred]))
            dis = (deltaE_pred-deltaE)/deltaE > 0.15
            mapping = {False: "", True: "+"}
            dis = [mapping[item] for item in dis]
            dis = ",".join(dis)
            if dis == ",":
                dis = ""
            discrepancies_collect.append(dis)
        d = {
        "proto": protos_,
        r"$E^{DFT} (eV)$": deltaEs_collect,
        r"$E^{pred} (eV)$": deltaEs_pred_collect,
        "discrepancy": discrepancies_collect,
        }
        df = pd.DataFrame(data=d)
        fig = plt.figure(figsize=(5, 6))
        ax = fig.add_subplot(111, frame_on=False)
        ax=gentable(ax,df)
        plt.savefig("Vacancies/Vacancy_formation_s{0}.png".format(key.split('_')[-1]), dpi=200)