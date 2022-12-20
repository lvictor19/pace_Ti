# -*- coding: utf-8 -*-

import os
import json
import glob
import pandas as pd
from pandas.plotting import table
from monty.json import MontyDecoder, MontyEncoder
from sklearn.metrics import mean_absolute_error, mean_squared_error
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import scienceplots
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from uberplot import (
    load_record,
    tabulize_record,
    in_plane_vectors,
    transform_energy_to_surface_energy,
    uber_fit,
    uber,
)
from pptx import Presentation
from pptx.util import Inches, Pt

import json
from monty.json import MontyDecoder, MontyEncoder
from sklearn.metrics import mean_absolute_error, mean_squared_error

######setup plotting#####
colors = [
    "black",
    "dimgray",
    "lightcoral",
    "brown",
    "red",
    "orangered",
    "chocolate",
    "burlywood",
    "orange",
    "darkgoldenrod",
    "gold",
    "olive",
    "darkolivegreen",
    "chartreuse",
    "green",
    "turquoise",
    "darkcyan",
    "deepskyblue",
    "dodgerblue",
    "midnightblue",
    "blue",
    "blueviolet",
    "violet",
    "deeppink",
]
index = np.array(
    [
        0,
        4,
        8,
        12,
        16,
        20,
        1,
        5,
        9,
        13,
        17,
        21,
        2,
        6,
        10,
        14,
        18,
        22,
        3,
        7,
        11,
        15,
        19,
        23,
    ]
)
colors_ = []
for i in index:
    colors_.append(colors[i])
colors = colors_
markers = [
    ".",
    "o",
    "v",
    "^",
    "<",
    ">",
    "s",
    "*",
    "D",
    "1",
    "+",
    "x",
    ".",
    "o",
    "v",
    "^",
    "<",
    ">",
    "s",
    "*",
    "D",
    "1",
    "+",
    "x",
]


#####save a picture of a table#####
def savetable(data: dict, name: str, dpi: int):
    """
    input a dictionary of a table; the name and dpi of output picture
    output a picture of the table
    """
    df = pd.DataFrame(data=data)
    fig = plt.figure(figsize=(5, 6))
    ax = fig.add_subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    table(ax, df, loc="center")
    plt.savefig(name, dpi=dpi)


#####get id for prototypes in GB,surface calculations#####
def get_ref_id(data_collection: list):
    """
    get the ids for prototypes bcc, fcc, hcp, omega
    """
    ref_id = {}
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "icsd"
            and data_collection[i]["calc"] == "final"
            and "proto" in data_collection[i]["metadata"]
        ):
            if data_collection[i]["metadata"]["proto"] in [
                "bcc.vasp",
                "fcc.vasp",
                "hcp.vasp",
                "omega.vasp",
            ]:
                ref_id[data_collection[i]["metadata"]["proto"]] = i
    return ref_id


#####plot surface heatmap#####
def plot_surface_heatmap(ax, X, Y, Z):
    pc = ax.pcolormesh(X, Y, Z)
    ax.set_aspect("equal")
    im_ratio = (max(Y.ravel()) - min(Y.ravel())) / (max(X.ravel()) - min(X.ravel()))
    cbar = plt.colorbar(pc, ax=ax, fraction=0.047 * im_ratio, pad=0.04)
    return ax, cbar


####### Vacancy formation energy #######
########################################
def Vacancy_formation(data_collection: list):
    """
    Vacancy formation energy analysis
    input data_collection::list
    output a picture of the table of the vacancy formation energies
    """
    vacancy_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "vacancies"
            and data_collection[i]["calc"] == "final"
            and not data_collection[i]["metadata"]["proto"] == "casi-alth"
            and not data_collection[i]["metadata"]["proto"] == "mp-865373POSCAR"
        ):
            vacancy_id.append(i)
    vacancy_ideal_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "vacancies_ideal"
            and data_collection[i]["calc"] == "final"
        ):
            vacancy_ideal_id.append(i)
    protos = []
    protos_ = []
    for id in vacancy_id:
        protos.append(
            data_collection[id]["metadata"]["proto"]
            + "_"
            + str(data_collection[id]["metadata"]["deleted_site"])
            + "|"
            + str(id)
        )
        protos_.append(
            data_collection[id]["metadata"]["proto"]
            + "_"
            + str(data_collection[id]["metadata"]["deleted_site"])
        )
    protos_ideal = {}
    for id in vacancy_ideal_id:
        protos_ideal[data_collection[id]["metadata"]["proto"]] = id
    deltaEs = []
    deltaEs_pred = []
    for proto in list(protos):
        delid = int(proto.split("|")[1])
        idealid = protos_ideal[proto.split("_")[0]]
        N_proto = len(data_collection[idealid]["structure"].species)
        Eideal = data_collection[idealid]["energy"]
        Evac = data_collection[delid]["energy"]
        Eideal_pred = data_collection[idealid]["pace"]["energy"]
        Evac_pred = data_collection[delid]["pace"]["energy"]
        deltaE = Evac - (N_proto - 1) * Eideal / N_proto
        deltaE_pred = Evac_pred - (N_proto - 1) * Eideal_pred / N_proto
        deltaEs.append(deltaE)
        deltaEs_pred.append(deltaE_pred)
    d = {
        "Structure_site": protos_,
        r"$E^{DFT} (eV)$": np.round(np.array(deltaEs), 5),
        r"$E^{pred} (eV)$": np.round(np.array(deltaEs_pred), 5),
    }
    savetable(d, "Vacancy_formation_pace.png", dpi=200)
    return


####### Volumetric deformation #######
#############################
def Volumetric_deformation(
    data_collection: list, ref_omega_dft: float, ref_omega_pace: float
):
    """
    Volumetric deformation energy analysis
    input data_collection::list
    ref_omega_dft::float   reference energy (per atom) for omega structure from dft calculation
    ref_omega_pace::float  reference energy (per atom) for omega structure from pace potential
    output pictures of Energy curves with deformation
    """
    volumetric_id = []
    for i in np.arange(len(data_collection)):
        if data_collection[i]["metadata"]["perturbation"] == "volumetric":
            volumetric_id.append(i)
    protos = set()
    for i in volumetric_id:
        protos.add(data_collection[i]["metadata"]["proto"])
    with plt.style.context("science"):
        plt.figure()
        for i in np.arange(len(protos)):
            Vs = []
            Es = []
            Epreds = []
            for j in volumetric_id:
                if data_collection[j]["metadata"]["proto"] == list(protos)[i]:
                    Vs.append(
                        data_collection[j]["metadata"]["volume"]
                        / len(data_collection[j]["structure"].species)
                    )
                    Es.append(
                        data_collection[j]["energy"]
                        / len(data_collection[j]["structure"].species)
                    )
                    Epreds.append(
                        data_collection[j]["pace"]["energy"]
                        / len(data_collection[j]["structure"].species)
                    )
            Vs = np.array(Vs)
            Es = np.array(Es) - ref_omega_dft
            Epreds = np.array(Epreds) - ref_omega_pace
            order = np.argsort(Vs)
            plt.scatter(
                Vs[order],
                Es[order],
                c=colors[i],
                marker=markers[i],
                s=0.2,
                label=list(protos)[i],
            )
            plt.plot(Vs[order], Epreds[order], c=colors[i])
        plt.xlabel(r"$\overline{V} \enspace (Å^{3}/atom)$")
        plt.ylabel(r"$\Delta E (eV/atom)$")
        plt.legend(bbox_to_anchor=(0.5, -0.05))
        plt.savefig("Volumetric_deformation_pace.png", dpi=200)
        plt.close()
    for iter in np.arange(4):
        with plt.style.context("science"):
            plt.figure()
            for i in np.arange(6):
                Vs = []
                Es = []
                Epreds = []
                for j in volumetric_id:
                    if (
                        data_collection[j]["metadata"]["proto"]
                        == list(protos)[i + iter * 6]
                    ):
                        Vs.append(
                            data_collection[j]["metadata"]["volume"]
                            / len(data_collection[j]["structure"].species)
                        )
                        Es.append(
                            data_collection[j]["energy"]
                            / len(data_collection[j]["structure"].species)
                        )
                        Epreds.append(
                            data_collection[j]["pace"]["energy"]
                            / len(data_collection[j]["structure"].species)
                        )
                Vs = np.array(Vs)
                Es = np.array(Es) - ref_omega_dft
                Epreds = np.array(Epreds) - ref_omega_pace
                order = np.argsort(Vs)
                plt.scatter(
                    Vs[order],
                    Es[order],
                    c=colors[i],
                    marker=markers[i],
                    s=3,
                    label=list(protos)[i + iter * 6],
                )
                plt.plot(Vs[order], Epreds[order], c=colors[i])
            plt.xlabel(r"$\overline{V} \enspace (Å^{3}/atom)$")
            plt.ylabel(r"$\Delta E \enspace (eV/atom)$")
            plt.legend(bbox_to_anchor=(0.5, -0.15))
            plt.savefig("Volumetric_deformation_" + str(iter) + "_pace.png", dpi=250)
            plt.close()


####### Grain boudary #######
#############################
def Grain_boundary(data_collection: list):
    """
    Grain boundary energy analysis
    input data_collection::list
    output a picture of the table of the grain boundary energies
    """
    protoids = get_ref_id(data_collection)
    print(protoids)
    gb_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "gb"
            and data_collection[i]["calc"] == "final"
        ):
            gb_id.append(i)
    icsd_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "icsd"
            and data_collection[i]["calc"] == "final"
        ):
            icsd_id.append(i)
    deltaEs_gb = []
    deltaEspred_gb = []
    protos = []
    for id in gb_id:
        protoid = 0
        if data_collection[id]["metadata"]["ref_elem"] == "Mo":
            protoid = protoids["bcc.vasp"]
            protos.append("bcc")
        if data_collection[id]["metadata"]["ref_elem"] == "Ti":
            protoid = protoids["hcp.vasp"]
            protos.append("hcp")
        N_gb = len(data_collection[id]["structure"].species)
        M_icsd = len(data_collection[protoid]["structure"].species)
        E_gb = data_collection[id]["energy"]
        E_gb_pred = data_collection[id]["pace"]["energy"]
        E_icsd = data_collection[protoid]["energy"]
        E_icsd_pred = data_collection[protoid]["pace"]["energy"]
        lattice = data_collection[id]["structure"].lattice.matrix
        A = np.linalg.norm(np.cross(lattice[0, :], lattice[1, :]))
        deltaE_gb = (E_gb - E_icsd * N_gb / M_icsd) / 2 / A
        deltaEpred_gb = (E_gb_pred - E_icsd_pred * N_gb / M_icsd) / 2 / A
        deltaEs_gb.append(deltaE_gb)
        deltaEspred_gb.append(deltaEpred_gb)
    planes = [data_collection[gb_id[i]]["metadata"]["plane"] for i in np.arange(8)]
    sigmas = [data_collection[gb_id[i]]["metadata"]["sigma"] for i in np.arange(8)]
    d = {
        "proto": protos,
        "plane": planes,
        "sigma": sigmas,
        r"$E^{DFT}(eV/Å^{2})$": np.round(np.array(deltaEs_gb), 5),
        r"$E^{pred}(eV/Å^{2})$": np.round(np.array(deltaEspred_gb), 5),
    }
    savetable(d, "Grain_boundary_pace.png", 200)


####### Surface #######
#######################
def Surface(data_collection: list):
    """
    Surface boundary energy analysis
    input data_collection::list
    output a picture of the table of the surface energies
    """
    surface_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "surfaces"
            and data_collection[i]["calc"] == "final"
        ):
            surface_id.append(i)
    icsd_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "icsd"
            and data_collection[i]["calc"] == "final"
        ):  # and data_collection[i]['metadata']['proto']=='bcc.vasp' or data_collection[i]['metadata']['proto']=='hcp.vasp':
            icsd_id.append(i)
    protoids = get_ref_id(data_collection)
    deltaEs_sf = []
    deltaEspred_sf = []
    for id in surface_id:
        protoid = protoids[data_collection[id]["metadata"]["proto"]]
        N_sf = len(data_collection[id]["structure"].species)
        M_icsd = len(data_collection[protoid]["structure"].species)
        E_sf = data_collection[id]["energy"]
        E_sf_pred = data_collection[id]["pace"]["energy"]
        E_icsd = data_collection[protoid]["energy"]
        E_icsd_pred = data_collection[protoid]["pace"]["energy"]
        lattice = data_collection[id]["structure"].lattice.matrix
        A = np.linalg.norm(np.cross(lattice[0, :], lattice[1, :]))
        deltaE_sf = (E_sf - E_icsd * N_sf / M_icsd) / 2 / A
        deltaEpred_sf = (E_sf_pred - E_icsd_pred * N_sf / M_icsd) / 2 / A
        deltaEs_sf.append(deltaE_sf)
        deltaEspred_sf.append(deltaEpred_sf)
    protos = [
        data_collection[id]["metadata"]["proto"].split(".")[0] for id in surface_id
    ]
    planes = [data_collection[id]["metadata"]["miller"] for id in surface_id]
    d = {
        "proto": protos,
        "plane": planes,
        r"$E^{DFT}(eV/Å^{2})$": np.round(np.array(deltaEs_sf), 5),
        r"$E^{pred}(eV/Å^{2})$": np.round(np.array(deltaEspred_sf), 5),
    }
    savetable(d, "Surfaces_pace.png", 200)


####### Strain #######
######################
def Strain(data_collection: list, ref_omega_dft: float, ref_omega_pace):
    """
    Strain-Energy relations
    input
        data_collection::list
        ref_omega_dft::float   reference energy (per atom) for omega structure from dft calculation
        ref_omega_pace::float  reference energy (per atom) for omega structure from pace potential
    output pictures of the Energy-Strain curves (a picture contains different strain directions for the same prototype)
    """
    strain_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "strain"
            and data_collection[i]["calc"] == "final"
        ):
            strain_id.append(i)

    protos = set()
    for id in strain_id:
        protos.add(data_collection[id]["metadata"]["proto"])
    protos = list(protos)
    for proto in protos:
        strain_id = []
        for i in np.arange(len(data_collection)):
            if (
                data_collection[i]["metadata"]["perturbation"] == "strain"
                and data_collection[i]["calc"] == "final"
                and data_collection[i]["metadata"]["proto"] == proto
            ):
                strain_id.append(i)
        Es = np.zeros((len(strain_id), 6))
        i = 0
        for id in strain_id:
            Es[i] = np.array(data_collection[id]["metadata"]["special_direction"])
            i += 1
        groups = []
        already = np.zeros(len(Es))
        isidentical = abs(np.dot(Es, Es.T) - 1) < 1e-12
        for i in np.arange(len(Es)):
            if not already[i] == 1:
                thisgroup = np.where(isidentical[i] == True)
                groups.append(thisgroup[0].tolist())
                already[thisgroup] = 1
        i = 0
        with plt.style.context("science"):
            plt.figure()
            for group in groups:
                magnitudes_ = []
                Energies_ = []
                Energies_pred_ = []
                for groupid in group:
                    Natom = len(
                        data_collection[strain_id[groupid]]["structure"].species
                    )
                    magnitudes_.append(
                        data_collection[strain_id[groupid]]["metadata"]["magnitude"]
                    )
                    Energies_.append(
                        data_collection[strain_id[groupid]]["energy"] / Natom
                    )
                    Energies_pred_.append(
                        data_collection[strain_id[groupid]]["pace"]["energy"] / Natom
                    )
                magnitudes = []
                Energies = []
                Energies_pred = []
                magnitudes_ = np.array(magnitudes_)
                for checkindex in np.arange(len(magnitudes_)):
                    samemag = np.where(
                        abs(magnitudes_ - magnitudes_[checkindex]) < 1e-2
                    )[0]
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
                    s=2,
                    label=np.round(
                        data_collection[strain_id[group[0]]]["metadata"][
                            "special_direction"
                        ],
                        4,
                    ),
                )
                plt.plot(
                    magnitudes[order],
                    Energies_pred[order] - ref_omega_pace,
                    c=colors[i],
                )
                # plt.legend(bbox_to_anchor=(0.5, -0.05))
                plt.xlabel("Strain Magnitude")
                plt.ylabel(r"$\Delta E \enspace (eV/atom)$")
                i += 1
            plt.savefig("Strains_" + proto + "_pace.png", dpi=200)
            plt.close()


def plotlinescan(datapace: np.ndarray, datadft: np.ndarray, name: str, axisname: str):
    """
    plot stacking fault energy linescan
    input datapace
    """
    with plt.style.context("science"):
        plt.figure(figsize=(5, 4))
        plt.scatter(np.arange(0, 1.01, 1 / 12), datadft, marker=".", color="k")
        plt.plot(np.arange(0, 1.01, 1 / 12), datapace, color="orangered")
        plt.xlabel(axisname)
        plt.ylabel(r"$Excess \enspace Energy \enspace (mJ/m^2)$")
        plt.tight_layout()
        plt.savefig(name)
        plt.close()


####### Generalized stacking fault energy #######
################ without uberfit ################
def Gsfe_withoutUber(data_collection: list):
    """
    Generalized stacking fault energy analysis without UBER fittting
    input data_collection::list
    output pictures of the gamma-surfaces
    """

    faces = ["bcc100", "bcc110", "fcc100", "fcc111", "basal", "prismatic", "pyramidal"]
    protos = ["bcc", "bcc", "fcc", "fcc", "hcp", "hcp", "hcp"]
    planes = ["100", "110", "100", "111", "basal", "prismatic", "pyramidal"]

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

    for whichface in np.arange(7):
        face = faces[whichface]
        with open(face + "record.json") as f:
            record = json.loads(f.read())
        equivalents = []
        for equi1 in record["equivalents"]:
            equi_orbit = []
            for equi2 in equi1:
                if equi2.endswith("0.000000"):
                    a, b, c = equi2.split(":")
                    equi_orbit.append(a + "." + b)
                    # print(equi_orbit)
            if not equi_orbit == []:
                equivalents.append(equi_orbit)

        values_pace = np.zeros((13, 13), dtype=float)
        values_dft = np.zeros((13, 13), dtype=float)
        for i in np.arange(len(equivalents)):
            first = equivalents[i][0]
            eqe = 0
            for data in data_collection:
                if (
                    data["calc"] == "final"
                    and data["metadata"]["perturbation"] == "gsfe"
                    and data["metadata"]["proto"] == protos[whichface]
                    and data["metadata"]["plane"] == planes[whichface]
                    and data["metadata"]["shift0"] == 0
                    and data["metadata"]["shift1"] == 0
                    and data["metadata"]["cleavage"] == 0
                ):
                    eqe = data["pace"]["energy"]
            value_pace = 0
            value_dft = 0
            for data in data_collection:
                if (
                    data["calc"] == "final"
                    and data["metadata"]["perturbation"] == "gsfe"
                    and data["metadata"]["proto"] == protos[whichface]
                    and data["metadata"]["plane"] == planes[whichface]
                    and data["metadata"]["shift0"] == int(first.split(".")[0])
                    and data["metadata"]["shift1"] == int(first.split(".")[1])
                    and data["metadata"]["cleavage"] == 0
                ):
                    value_pace = data["pace"]["energy"]
                    value_dft = data["energy"]

            for j in np.arange(len(equivalents[i])):
                which = equivalents[i][j].split(".")
                values_pace[int(which[0]), int(which[1])] = value_pace
                values_dft[int(which[0]), int(which[1])] = value_dft

        values_pace[12, :] = values_pace[0, :]
        values_pace[:, 12] = values_pace[:, 0]
        values_pace = (values_pace - values_pace[0, 0]) * 1000

        values_dft[12, :] = values_dft[0, :]
        values_dft[:, 12] = values_dft[:, 0]
        values_dft = (values_dft - values_dft[0, 0]) * 1000

        def plot_surface_heatmap(ax, X, Y, Z):
            pc = ax.pcolormesh(X, Y, Z)
            ax.set_aspect("equal")
            im_ratio = (max(Y.ravel()) - min(Y.ravel())) / (
                max(X.ravel()) - min(X.ravel())
            )
            cbar = plt.colorbar(pc, ax=ax, fraction=0.047 * im_ratio, pad=0.04)

            return ax, cbar

        values_pace = (values_pace - values_pace[0, 0]) * 1000
        values_dft = (values_dft - values_dft[0, 0]) * 1000

        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = np.arange(0, 1 + 1e-7, 1 / 12)
        y = np.arange(0, 1 + 1e-7, 1 / 12)
        X, Y = np.meshgrid(x, y)
        ax, cbar = plot_surface_heatmap(ax, X, Y, values_pace)
        plt.tight_layout()
        plt.savefig("2d_" + faces[whichface] + "_pace.png")
        plt.close()
        for line in whichlines[faces[whichface]]:
            if line[0] == "0_":
                plotlinescan(
                    values_pace[0, :],
                    values_dft[0, :],
                    faces[whichface] + "-" + "linescan0_.png",
                    line[1],
                )
            elif line[0] == "_0":
                plotlinescan(
                    values_pace[:, 0],
                    values_dft[:, 0],
                    faces[whichface] + "-" + "linescan_0.png",
                    line[1],
                )
            elif line[0] == "__":
                plotlinescan(
                    np.diag(values_pace),
                    np.diag(values_dft),
                    faces[whichface] + "-" + "linescan__.png",
                    line[1],
                )
            elif line[0] == "_-_":
                plotlinescan(
                    np.diag(np.flipud(values_pace)),
                    np.diag(np.flipud(values_dft)),
                    faces[whichface] + "-" + "linescan_-_.png",
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
                    faces[whichface] + "-" + "linescan_-_.png",
                    line[1],
                )


#####modified functions from uber.py#####
def shift_energies_to_surface_energy(unwinded, eqe):
    gamma = surface_energy(unwinded, eqe)
    unwinded["energy"] -= 2 * gamma + eqe
    return unwinded


def surface_energy(unwinded, eqe):
    maxcleave = max(unwinded["cleavage"])

    surface_energies = (
        unwinded.loc[unwinded["cleavage"] == maxcleave]["energy"] - eqe
    ) * 0.5
    if max(surface_energies) - min(surface_energies) > 0.001:
        raise ValueError(
            "Surface energies don't match at different shifts. Did you separate enough?"
        )

    return np.mean(surface_energies)


####### Generalized stacking fault energy #######
################ with uberfit ###################
def gsfe(data_collection):
    """
    Generalized stacking fault energy analysis with UBER fittting
    input data_collection::list
    output pictures of the gamma-surfaces
    """
    faces = ["bcc100", "bcc110", "fcc100", "fcc111", "basal", "prismatic", "pyramidal"]
    protos = ["bcc", "bcc", "fcc", "fcc", "hcp", "hcp", "hcp"]
    planes = ["100", "110", "100", "111", "basal", "prismatic", "pyramidal"]
    for whichface in np.arange(7):
        face = faces[whichface]
        with open(face + "record.json") as f:
            record = json.loads(f.read())
        equivalents = []
        for equi1 in record["equivalents"]:
            equi_orbit = []
            for equi2 in equi1:
                if equi2.endswith("0.000000"):
                    a, b, c = equi2.split(":")
                    equi_orbit.append(a + "." + b)
                    # print(equi_orbit)
            if not equi_orbit == []:
                equivalents.append(equi_orbit)

        values = np.zeros((13, 13), dtype=float)
        for i in np.arange(len(equivalents)):
            first = equivalents[i][0]
            eqe = 0
            for data in data_collection:
                if (
                    data["calc"] == "final"
                    and data["metadata"]["perturbation"] == "gsfe"
                    and data["metadata"]["proto"] == protos[whichface]
                    and data["metadata"]["plane"] == planes[whichface]
                    and data["metadata"]["shift0"] == 0
                    and data["metadata"]["shift1"] == 0
                    and data["metadata"]["cleavage"] == 0
                ):
                    eqe = data["pace"]["energy"]
            cleaves = []
            cleave_energies = []
            for data in data_collection:
                if (
                    data["calc"] == "final"
                    and data["metadata"]["perturbation"] == "gsfe"
                    and data["metadata"]["proto"] == protos[whichface]
                    and data["metadata"]["plane"] == planes[whichface]
                    and data["metadata"]["shift0"] == int(first.split(".")[0])
                    and data["metadata"]["shift1"] == int(first.split(".")[1])
                ):
                    cleaves.append(data["metadata"]["cleavage"])
                    cleave_energies.append(data["pace"]["energy"])
            value = 0
            cleaves = np.array(cleaves)
            cleave_energies = np.array(cleave_energies)
            index = np.argsort(cleaves)
            energies = cleave_energies[index]
            record = load_record("./" + faces[whichface] + "_cleave_record.json")
            unwinded = tabulize_record(record)
            avec, bvec = in_plane_vectors(record)
            unwinded["energy"] = energies
            unwinded_slice = shift_energies_to_surface_energy(unwinded, eqe)
            popt, pcov = uber_fit(unwinded_slice)
            value = -2 * popt[1] + popt[3]
            for j in np.arange(len(equivalents[i])):
                which = equivalents[i][j].split(".")
                values[int(which[0]), int(which[1])] = value
        values[12, :] = values[0, :]
        values[:, 12] = values[:, 0]
        values = (values - values[0, 0]) * 1000

        def plot_surface_heatmap(ax, X, Y, Z):
            pc = ax.pcolormesh(X, Y, Z)
            ax.set_aspect("equal")
            im_ratio = (max(Y.ravel()) - min(Y.ravel())) / (
                max(X.ravel()) - min(X.ravel())
            )
            cbar = plt.colorbar(pc, ax=ax, fraction=0.047 * im_ratio, pad=0.04)

            return ax, cbar

        values = (values - values[0, 0]) * 1000

        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = np.arange(0, 1 + 1e-7, 1 / 12)
        y = np.arange(0, 1 + 1e-7, 1 / 12)
        X, Y = np.meshgrid(x, y)
        ax, cbar = plot_surface_heatmap(ax, X, Y, values)
        plt.tight_layout()
        plt.savefig("2d_" + faces[whichface] + "_pace.png")
        plt.close()


def Burgers_Bains(data_collection: list):
    """
    analysis for Burger_Bains pathways connecting fcc-bcc-hcp
    """
    B_B_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "Burgers_Bains"
            and data_collection[i]["metadata"]["value"][:3] in ["sh_", "sp_", "sn_"]
            # and data_collection[i]["calc"] == "final"
        ):
            B_B_id.append(i)
    bccref_pace = 0
    bccref_dft = 0
    for data in data_collection:
        if (
            data["metadata"]["perturbation"] == "icsd"
            and data["calc"] == "final"
            and data["metadata"]["proto"] == "bcc.vasp"
        ):
            bccref_pace = data["pace"]["energy"] / len(data["structure"])
            bccref_dft = data["energy"] / len(data["structure"])

    bccref_pace
    unshuffled_values = []
    unshuffled_energies_dft = []
    unshuffled_energies_pace = []

    shuffled_values = []
    shuffled_energies_dft = []
    shuffled_energies_pace = []

    for id in B_B_id:
        if "sh" in data_collection[id]["metadata"]["value"]:
            symbol = 1 if "p" in data_collection[id]["metadata"]["value"] else -1
            shuffled_values.append(
                symbol * float(data_collection[id]["metadata"]["value"].split("_")[-1])
            )
            shuffled_energies_dft.append(
                float(data_collection[id]["energy"])
                / len(data_collection[id]["structure"])
                - bccref_dft
            )
            shuffled_energies_pace.append(
                float(data_collection[id]["pace"]["energy"])
                / len(data_collection[id]["structure"])
                - bccref_pace
            )
        else:
            symbol = 1 if "p" in data_collection[id]["metadata"]["value"] else -1
            unshuffled_values.append(
                symbol * float(data_collection[id]["metadata"]["value"].split("_")[-1])
            )
            unshuffled_energies_dft.append(
                float(data_collection[id]["energy"])
                / len(data_collection[id]["structure"])
                - bccref_dft
            )
            unshuffled_energies_pace.append(
                float(data_collection[id]["pace"]["energy"])
                / len(data_collection[id]["structure"])
                - bccref_pace
            )
    unshuffled_values = np.array(unshuffled_values)
    unshuffled_energies_dft = np.array(unshuffled_energies_dft)
    unshuffled_energies_pace = np.array(unshuffled_energies_pace)

    shuffled_values = np.array(shuffled_values)
    shuffled_energies_dft = np.array(shuffled_energies_dft)
    shuffled_energies_pace = np.array(shuffled_energies_pace)

    unshuffled_index = np.argsort(unshuffled_values)
    shuffled_index = np.argsort(shuffled_values)

    unshuffled_values = unshuffled_values[unshuffled_index]
    unshuffled_energies_dft = unshuffled_energies_dft[unshuffled_index]
    unshuffled_energies_pace = unshuffled_energies_pace[unshuffled_index]

    shuffled_values = shuffled_values[shuffled_index]
    shuffled_energies_dft = shuffled_energies_dft[shuffled_index]
    shuffled_energies_pace = shuffled_energies_pace[shuffled_index]

    with plt.style.context("science"):
        plt.figure()
        plt.scatter(shuffled_values, shuffled_energies_dft, s=6, c="r")
        plt.scatter(unshuffled_values, unshuffled_energies_dft, s=6, c="b")
        plt.plot(shuffled_values, shuffled_energies_pace, c="r")
        plt.plot(unshuffled_values, unshuffled_energies_pace, c="b")
        plt.ylabel("Energy (eV/atom)", fontsize=12)
        plt.tick_params(axis="y", direction="in")
        plt.tick_params(axis="x", direction="in")
        plt.savefig("Ti_Burgers_strain.png")
        plt.close()


def bcc_omega(data_collection: list, ref_omega_dft: float, ref_omega_pace: float):
    """
    analysis for the pathway connecting bcc and omega
    ref_omega_dft::float   reference energy (per atom) for omega structure from dft calculation
    ref_omega_pace::float  reference energy (per atom) for omega structure from pace potential
    """
    b_o_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "bcc_omega"
            and data_collection[i]["calc"] == "final"
        ):
            b_o_id.append(i)
    values = []
    Energies_dft = []
    Energies_pace = []
    for id in b_o_id:
        values.append(data_collection[id]["metadata"]["value"])
        Energies_pace.append(
            data_collection[id]["pace"]["energy"]
            / len(data_collection[id]["structure"])
        )
        Energies_dft.append(
            data_collection[id]["energy"] / len(data_collection[id]["structure"])
        )
    values = np.array(values)
    Energies_dft = np.array(Energies_dft)
    Energies_pace = np.array(Energies_pace)
    index = np.argsort(values)
    values = values[index]
    Energies_dft = Energies_dft[index]
    Energies_pace = Energies_pace[index]
    with plt.style.context("science"):
        plt.figure()
        plt.plot(
            np.arange(len(Energies_pace)), Energies_pace - ref_omega_pace, c="black"
        )
        plt.scatter(
            np.arange(len(Energies_dft)), Energies_dft - ref_omega_dft, c="black", s=4
        )

        plt.xlabel("Omega Transformation")
        plt.xticks([])
        plt.ylabel("Energy/Atom (eV)")
        plt.legend(fontsize=8)
        plt.savefig("Omega.png", dpi=150)
        plt.close()


def hcp_omega(data_collection: list):
    """
    analysis for the pathways connecting hcp and omega
    """
    h_o_static_id = []
    for i in np.arange(len(data_collection)):
        if data_collection[i]["metadata"]["perturbation"] == "hcp_omega":
            try:
                print(data_collection[i]["metadata"])
                if (
                    "NEB" in data_collection[i]["metadata"].keys()
                    and data_collection[i]["metadata"]["NEB"] == "False"
                ):
                    h_o_static_id.append(i)
            except:
                pass
    h_o_NEB_id = []
    for i in np.arange(len(data_collection)):
        if data_collection[i]["metadata"]["perturbation"] == "hcp_omega":
            if (
                "NEB" in data_collection[i]["metadata"].keys()
                and data_collection[i]["metadata"]["NEB"] == "True"
            ):
                h_o_NEB_id.append(i)
    hcp_energy_pace = 0
    omega_energy_pace = 0
    hcp_energy_dft = 0
    omega_energy_dft = 0
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
            and data_collection[i]["metadata"]["proto"] == "hcp"
        ):
            hcp_energy_dft = data_collection[i]["energy"] / len(
                data_collection[i]["structure"]
            )
            hcp_energy_pace = data_collection[i]["pace"]["energy"] / len(
                data_collection[i]["structure"]
            )
        if (
            data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
            and data_collection[i]["metadata"]["proto"] == "omega"
        ):
            omega_energy_dft = data_collection[i]["energy"] / len(
                data_collection[i]["structure"]
            )
            omega_energy_pace = data_collection[i]["pace"]["energy"] / len(
                data_collection[i]["structure"]
            )
    pathways = set()
    for id in h_o_static_id:
        pathways.add(data_collection[id]["metadata"]["#pathway"])
    pathways = list(pathways)

    for pathway in pathways:
        energies_dft = np.zeros(7)
        energies_dft[0] = omega_energy_dft
        energies_dft[6] = hcp_energy_dft
        energies_pace = np.zeros(7)
        energies_pace[0] = omega_energy_pace
        energies_pace[6] = hcp_energy_pace
        for image in np.arange(1, 6, 1):
            for i in np.arange(len(data_collection)):
                if (
                    "NEB" in data_collection[i]["metadata"]
                    and data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
                    and data_collection[i]["metadata"]["NEB"] == "True"
                    and data_collection[i]["metadata"]["#pathway"] == pathway
                    and int(data_collection[i]["metadata"]["image"]) == image
                ):
                    energies_pace[image] = data_collection[i]["pace"]["energy"] / len(
                        data_collection[i]["structure"]
                    )
                    energies_dft[image] = data_collection[i]["energy"] / len(
                        data_collection[i]["structure"]
                    )
        energies_dft = energies_dft - energies_dft[0]
        energies_pace = energies_pace - energies_pace[0]
        with plt.style.context("science"):
            plt.figure(figsize=(4, 3.2))
            plt.scatter(np.arange(7), energies_dft)
            plt.plot(np.arange(7), energies_pace)
            plt.xlabel("Transformation Coordinates")
            plt.ylabel(r"$Excess \enspace Energy \enspace (eV/Å)$")
            plt.xticks([0, 6], ["omega", "hcp"])
            plt.savefig(str(pathway) + "_NEB_pace.png")
            plt.close()

        energies_dft = np.zeros(7)
        energies_dft[0] = omega_energy_dft
        energies_dft[6] = hcp_energy_dft
        energies_pace = np.zeros(7)
        energies_pace[0] = omega_energy_pace
        energies_pace[6] = hcp_energy_pace
        for image in np.arange(1, 6, 1):
            for i in np.arange(len(data_collection)):
                if (
                    "NEB" in data_collection[i]["metadata"]
                    and data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
                    and data_collection[i]["metadata"]["NEB"] == "False"
                    and data_collection[i]["metadata"]["#pathway"] == pathway
                    and int(data_collection[i]["metadata"]["image"]) == image
                ):
                    energies_pace[image] = data_collection[i]["pace"]["energy"] / len(
                        data_collection[i]["structure"]
                    )
                    energies_dft[image] = data_collection[i]["energy"] / len(
                        data_collection[i]["structure"]
                    )
        energies_dft = energies_dft - energies_dft[0]
        energies_pace = energies_pace - energies_pace[0]
        with plt.style.context("science"):
            plt.figure(figsize=(4, 3.2))
            plt.scatter(np.arange(7), energies_dft)
            plt.plot(np.arange(7), energies_pace)
            plt.xlabel("Transformation Coordinates")
            plt.ylabel(r"$Excess \enspace Energy \enspace (eV/Å)$")
            plt.xticks([0, 6], ["omega", "hcp"])
            plt.savefig(str(pathway) + "_static_pace.png")
            plt.close()


def collect_in_powerpoint(
    data_collection, RMSE_E: float, MAE_E: float, RMSE_F: float, MAE_F: float
):
    """
    collect things in a powerpoint file
    input   RMSE_E::float Energy RMSE error
            MAE_E::float   Energy MAE error
            RMSE_F::float  Force RMSE error
            MAE_F::float   Force MAE error
    """
    prs = Presentation()
    blank_slide_layout = prs.slide_layouts[6]
    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = ""
    p = tf.add_paragraph()
    p.text = "Energy"
    p.font.size = Pt(30)
    p0 = tf.add_paragraph()
    p0.text = ""
    p1 = tf.add_paragraph()
    p1.text = "RMSE=%.5f" % RMSE_E
    p2 = tf.add_paragraph()
    p2.text = "MAE=%.5f" % MAE_E
    img_path = "pace_energy.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)
    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = ""
    p = tf.add_paragraph()
    p.text = "Force"
    p.font.size = Pt(30)
    p0 = tf.add_paragraph()
    p0.text = ""
    p1 = tf.add_paragraph()
    p1.text = "RMSE=%.5f" % RMSE_F
    p2 = tf.add_paragraph()
    p2.text = "MAE=%.5f" % MAE_F
    img_path = "pace_force.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)
    for i in np.arange(4):
        slide = prs.slides.add_slide(blank_slide_layout)
        left = top = width = height = Inches(1)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        tf.text = "Volumetric_" + str(i)
        img_path = "Volumetric_deformation_" + str(i) + "_pace.png"
        picleft = Inches(3)
        picheight = Inches(5.5)
        pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)
    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = "Grain Boundary"
    img_path = "Grain_boundary_pace.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = "Surfaces"
    img_path = "Surfaces_pace.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = "Vacancy formation"
    img_path = "Vacancy_formation_pace.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = "Burgers_Bains"
    img_path = "Ti_Burgers_strain.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    slide = prs.slides.add_slide(blank_slide_layout)
    left = top = width = height = Inches(1)
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.text = "Omega transaformation"
    img_path = "Omega.png"
    picleft = Inches(3)
    picheight = Inches(5.5)
    pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    strain_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "strain"
            and data_collection[i]["calc"] == "final"
        ):
            strain_id.append(i)

    protos = set()
    for id in strain_id:
        protos.add(data_collection[id]["metadata"]["proto"])
    protos = list(protos)
    for proto in protos:
        slide = prs.slides.add_slide(blank_slide_layout)
        left = top = width = height = Inches(1)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        tf.text = proto
        img_path = "Strains_" + proto + "_pace.png"
        picleft = Inches(3)
        picheight = Inches(5.5)
        pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    pathnum = [
        "0_static",
        "0_NEB",
        "1_static",
        "1_NEB",
        "67_static",
        "67_NEB",
        "297_static",
        "297_NEB",
        "357_static",
        "357_NEB",
        "499_static",
        "499_NEB",
    ]
    for ptn in pathnum:
        slide = prs.slides.add_slide(blank_slide_layout)
        left = top = width = height = Inches(1)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        tf.text = ptn
        img_path = ptn + "_pace.png"
        picleft = Inches(3)
        picheight = Inches(5.5)
        pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    faces = ["bcc100", "bcc110", "fcc100", "fcc111", "basal", "prismatic", "pyramidal"]
    for face in faces:
        slide = prs.slides.add_slide(blank_slide_layout)
        left = top = width = height = Inches(1)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        tf.text = face
        img_path = "2d_" + face + "_pace.png"
        picleft = Inches(3)
        picheight = Inches(5.5)
        pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    linescans_ = glob.glob("*linescan*")
    linescans1 = [l for l in linescans_ if "bcc100" in l]
    linescans2 = [l for l in linescans_ if "bcc110" in l]
    linescans3 = [l for l in linescans_ if "fcc100" in l]
    linescans4 = [l for l in linescans_ if "fcc111" in l]
    linescans5 = [l for l in linescans_ if "basal" in l]
    linescans6 = [l for l in linescans_ if "prismatic" in l]
    linescans7 = [l for l in linescans_ if "pyramidal" in l]
    linescans = (
        linescans1
        + linescans2
        + linescans3
        + linescans4
        + linescans5
        + linescans6
        + linescans7
    )

    for linescan in linescans:
        slide = prs.slides.add_slide(blank_slide_layout)
        left = top = width = height = Inches(1)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        tf.text = linescan.split(".")[0]
        img_path = linescan
        picleft = Inches(3)
        picheight = Inches(5.5)
        pic = slide.shapes.add_picture(img_path, picleft, top, height=picheight)

    prs.save("0000.pptx")


def main():
    # filename=glob(r'*.json')[0]
    filename = "Ti_data_pace.json"
    with open(filename) as f:
        data_collection_ = json.loads(f.read(), cls=MontyDecoder)
    data_collection = []
    for data in data_collection_:
        if data["metadata"]["perturbation"] != "pairs":
            data_collection.append(data)
    del data_collection_

    ######get reference energies#####
    ref_omega_dft = 0
    for data in data_collection:
        if (
            data["metadata"]["perturbation"] == "icsd"
            and data["calc"] == "final"
            and data["metadata"]["proto"] == "omega.vasp"
        ):
            ref_omega_dft = data["energy"]
    ref_omega_dft = ref_omega_dft / 3

    ref_omega_pace = 0
    for data in data_collection:
        if (
            data["metadata"]["perturbation"] == "icsd"
            and data["calc"] == "final"
            and data["metadata"]["proto"] == "omega.vasp"
        ):
            ref_omega_pace = data["pace"]["energy"]
    ref_omega_pace = ref_omega_pace / 3

    #####evaluate energy and force prediction#####
    energies_dft = []
    forces_dft = []
    energies_pace = []
    forces_pace = []

    energies_dft_stdt = []
    forces_dft_stdt = []
    energies_pace_stdt = []
    forces_pace_stdt = []
    for i in tqdm(np.arange(len(data_collection))):
        if not "standout" in data_collection[i].keys():
            atom_num = len(data_collection[i]["structure"].species)
            energies_dft.append(data_collection[i]["energy"] / atom_num - ref_omega_dft)
            forces_dft += data_collection[i]["forces"]
            energies_pace.append(
                data_collection[i]["pace"]["energy"] / atom_num - ref_omega_pace
            )
            forces_pace += data_collection[i]["pace"]["forces"].tolist()
        else:
            atom_num = len(data_collection[i]["structure"].species)
            energies_dft_stdt.append(
                data_collection[i]["energy"] / atom_num - ref_omega_dft
            )
            forces_dft_stdt += data_collection[i]["forces"]
            energies_pace_stdt.append(
                data_collection[i]["pace"]["energy"] / atom_num - ref_omega_pace
            )
            forces_pace_stdt += data_collection[i]["pace"]["forces"].tolist()
    with plt.style.context("science"):
        plt.figure()
        plt.scatter(energies_dft, energies_pace, c="k", s=0.3)
        plt.scatter(energies_dft_stdt, energies_pace_stdt, c="r", s=0.3)
        plt.xlabel(r"$\Delta E^{DFT} \enspace (eV/atom)$")
        plt.ylabel(r"$\Delta E^{pred} \enspace (eV/atom)$")
        plt.savefig("pace_energy.png", dpi=200)
        plt.close()

    energies_dft += energies_dft_stdt
    energies_pace += energies_pace_stdt
    RMSE_E = np.sqrt(mean_squared_error(energies_dft, energies_pace))
    MAE_E = mean_absolute_error(energies_dft, energies_pace)

    forces_sca_dft = np.sqrt(np.sum(np.array(forces_dft) ** 2, axis=1))
    forces_sca_pace = np.sqrt(np.sum(np.array(forces_pace) ** 2, axis=1))
    try:
        forces_sca_dft_stdt = np.sqrt(np.sum(np.array(forces_dft_stdt) ** 2, axis=1))
        forces_sca_pace_stdt = np.sqrt(np.sum(np.array(forces_pace_stdt) ** 2, axis=1))
    except:
        pass

    with plt.style.context("science"):
        plt.figure()
        plt.scatter(forces_sca_dft, forces_sca_pace, c="k", s=0.3, label="pace")
        try:
            plt.scatter(
                forces_sca_dft_stdt, forces_sca_pace_stdt, c="r", s=0.3, label="pace"
            )
        except:
            pass
        plt.xlabel(r"$|F^{i}|^{DFT} \enspace (eV/Å)$")
        plt.ylabel(r"$|F^{i}|^{pred} \enspace (eV/Å)$")
        plt.savefig("pace_force.png", dpi=200)
        plt.close()

    try:
        forces_sca_dft += forces_sca_dft_stdt
        forces_sca_pace += forces_sca_pace_stdt
    except:
        pass
    
    RMSE_F = np.sqrt(mean_squared_error(forces_sca_dft, forces_sca_pace))
    MAE_F = mean_absolute_error(forces_sca_dft, forces_sca_pace)
    Vacancy_formation(data_collection)
    Volumetric_deformation(data_collection, ref_omega_dft, ref_omega_pace)
    Grain_boundary(data_collection)
    Surface(data_collection)
    Strain(data_collection, ref_omega_dft, ref_omega_pace)
    Gsfe_withoutUber(data_collection)
    Burgers_Bains(data_collection)
    bcc_omega(data_collection, ref_omega_dft, ref_omega_pace)
    hcp_omega(data_collection)
    Strain(data_collection, ref_omega_dft, ref_omega_pace)
    collect_in_powerpoint(data_collection, RMSE_E, MAE_E, RMSE_F, RMSE_F)


if __name__ == "__main__":
    main()
