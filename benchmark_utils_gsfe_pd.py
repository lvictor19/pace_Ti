import numpy as np
from matplotlib import pyplot as plt
import json
import pandas as pd
import os


def in_plane_vectors(record):
    a,b=record["grid"]
    au,bu=record["shift_units"]
    avec,bvec=np.array(au)*a,np.array(bu)*b
    return avec,bvec

def transform_energy_to_surface_energy(E,avec,bvec):
    area=np.linalg.norm(np.cross(avec,bvec))
    ev_to_joules=1.602176565e-19
    angstrom_to_meter=1e-10
    return ev_to_joules*E/(area*angstrom_to_meter*angstrom_to_meter)

def get_whichlines() -> dict:
    '''
    get the gsfe linescans
    0_ for a horizotal linescan
    __ for diagonal
    _-_ for back diagonal
    '''
    whichlines = {}
    whichlines["bcc100"] = [("0_", r"$[100]$"), ("__", r"$[110]$")]
    whichlines["bcc110"] = [
        ("0_", r"$1/2 \enspace [1 \overline{1} 1]$"),
        ("__", r"$[100]$"),
        ("_-_", r"$[0 \overline{1} 1]$"),
    ]
    whichlines["fcc100"] = [("0_", r"$1/2 \enspace [110]$"), ("__", r"$[010]$")]
    whichlines["fcc111"] = [
        ("0_", r"$[101]$"),
        ("__", r"$[112]$"),
        ("_-_", r"$[1 \overline{1} 0]$"),
    ]
    whichlines["basal"] = [
        ("0_", r"$1/3 \enspace [1 1 \overline{2} 0]$"),
        ("__", r"$1/3 \enspace [\overline{2} 4 \overline{2} 0]$"),
        ("_-_", r"$1/3 \enspace [4 2 \overline{2} 0]$"),
    ]
    whichlines["prismatic"] = [
        ("0_", r"$1/3 \enspace [1 1 \overline{2} 0]$"),
        ("_0", r"$[0001]$"),
        ("__", r"$1/3 \enspace [1 1 \overline{2} 3]$"),
    ]
    whichlines["pyramidal"] = [
        ("_0", r"$1/3 \enspace [\overline{2} 1 1 0]$"),
        ("special", r"$1/3 \enspace [0 \overline{1} 1 2]$"),
    ]
    whichlines["pyramidal2nd"] = [
        ("0_", r"$[0 0 1 1]$"),
    ]
    return whichlines


def get_pert_data(data_collection:pd.DataFrame,pertname):
    '''
    get the strain data
    '''
    ispert=data_collection['perturbation'].map(lambda x: x==pertname)
    isfinal=data_collection['calc'].map(lambda x: x=='final')
    pert_data=data_collection[ispert&isfinal].copy()
    return pert_data
    
def get_gsfe_plane_data(data_collection:pd.DataFrame,protoname,plane):
    '''
    get the gsfe data for a plane
    '''
    if not data_collection['perturbation'].eq('gsfe').all():
        raise ValueError('Perturbation of input data not all gsfe for getting gsfe plane data')
    isproto=data_collection["metadata"].map(lambda x: x["proto"]==protoname)
    isplane=data_collection["metadata"].map(lambda x: x["plane"]==plane)
    gsfe_plane_data=data_collection[isproto&isplane].copy()
    return gsfe_plane_data

def plotlinescan(datapace: np.ndarray, datadft: np.ndarray, name: str, axisname: str):
    """
    plot stacking fault energy linescan
    input datapace
    """
    with plt.style.context("science"):
        plt.figure(figsize=(5, 4))
        plt.scatter(
            np.arange(0, 1.01, 1 / 12), datadft, marker=".", color="k", label="dft"
        )
        plt.plot(np.arange(0, 1.01, 1 / 12), datapace, color="k", label="pace")
        plt.xlabel(axisname)
        plt.ylabel(r"$Excess \enspace Energy \enspace (mJ/m^2)$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join('linescan',name))
        plt.close()

def gsfe(data_collection: pd.DataFrame):
    faces = ["bcc100", "bcc110", "fcc100", "fcc111", "basal", "prismatic", "pyramidal","pyramidal2nd"]
    protos = ["bcc", "bcc", "fcc", "fcc", "hcp", "hcp", "hcp","hcp"]
    planes = ["100", "110", "100", "111", "basal", "prismatic", "pyramidal","pyramidal2nd"]
    whichlines = get_whichlines()
    for whichface in range(len(faces)):
        face = faces[whichface]
        planedata=get_gsfe_plane_data(data_collection,protos[whichface],planes[whichface])
        with open("record/"+face + "record.json") as f:
            record = json.loads(f.read())
        equivalents = []
        for equi1 in record["equivalents"]:
            equi_orbit = []
            for equi2 in equi1:
                if equi2.endswith("0.000000"):
                    a, b, c = equi2.split(":")
                    equi_orbit.append(a + "." + b)
            if not equi_orbit == []:
                equivalents.append(equi_orbit)
        keys=list(data_collection.keys())
        for key in [x for x in keys if x.startswith('pace_energy')]:
            values_pace = np.zeros((13, 13), dtype=float)
            values_dft = np.zeros((13, 13), dtype=float)
            for i in np.arange(len(equivalents)):
                first = equivalents[i][0]
                isshift0=planedata["metadata"].map(lambda x: x["shifts"][0]==int(first.split(".")[0]))
                isshift1=planedata["metadata"].map(lambda x: x["shifts"][1]==int(first.split(".")[1]))
                value_pace=0
                value_dft=0
                try:
                    value_dft=planedata[isshift0&isshift1].iloc[0]['energy']
                    value_pace=planedata[isshift0&isshift1].iloc[0][key]
                except:
                    pass
                avec, bvec = in_plane_vectors(record)
                value_dft=transform_energy_to_surface_energy(value_dft, avec, bvec)
                value_pace=transform_energy_to_surface_energy(value_pace, avec, bvec)
                for j in np.arange(len(equivalents[i])):
                    which = equivalents[i][j].split(".")
                    values_pace[int(which[0]), int(which[1])] = value_pace
                    values_dft[int(which[0]), int(which[1])] = value_dft
            values_dft[:,12]=values_dft[:,0]
            values_dft[12,:]=values_dft[0,:]
            values_pace[:,12]=values_pace[:,0]
            values_pace[12,:]=values_pace[0,:]
            values_dft-=values_dft[0,0]
            values_pace-=values_pace[0,0]
            for line in whichlines[faces[whichface]]:
                if line[0] == "0_":
                    plotlinescan(
                        values_pace[0, :],
                        values_dft[0, :],
                        faces[whichface] + "-" + "linescan0_"+"s"+key.split('_')[-1]+".png",
                        line[1],
                    )
                elif line[0] == "_0":
                    plotlinescan(
                        values_pace[:, 0],
                        values_dft[:, 0],
                        faces[whichface] + "-" + "linescan_0"+"s"+key.split('_')[-1]+".png",
                        line[1],
                    )
                elif line[0] == "__":
                    plotlinescan(
                        np.diag(values_pace),
                        np.diag(values_dft),
                        faces[whichface] + "-" + "linescan__"+"s"+key.split('_')[-1]+".png",
                        line[1],
                    )
                elif line[0] == "_-_":
                    plotlinescan(
                        np.diag(np.flipud(values_pace)),
                        np.diag(np.flipud(values_dft)),
                        faces[whichface] + "-" + "linescan_-_"+"s"+key.split('_')[-1]+".png",
                        line[1],
                    )
                elif line[0] == "special":
                    plotlinescan(
                        values_pace[
                            np.arange(13),
                            np.array([0, 2, 4, 6, 8, 10, 12, 2, 4, 6, 8, 10, 12]),
                        ],
                        values_dft[
                            np.arange(13),
                            np.array([0, 2, 4, 6, 8, 10, 12, 2, 4, 6, 8, 10, 12]),
                        ],
                        faces[whichface] + "-" + "linescanspecial"+"s"+key.split('_')[-1]+".png",
                        line[1],
                    )