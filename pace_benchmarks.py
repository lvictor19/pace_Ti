# -*- coding: utf-8 -*-

import os
import json
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
        lattice = data_collection[id]["structure"].lattice.matrix
        A = np.linalg.norm(np.cross(lattice[0, :], lattice[1, :]))
        deltaE_gb = (E_gb - E_icsd * N_gb / M_icsd) / 2 / A
        deltaEpred_gb = (E_gb_pred - E_icsd * N_gb / M_icsd) / 2 / A
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
        lattice = data_collection[id]["structure"].lattice.matrix
        A = np.linalg.norm(np.cross(lattice[0, :], lattice[1, :]))
        deltaE_sf = (E_sf - E_icsd * N_sf / M_icsd) / 2 / A
        deltaEpred_sf = (E_sf_pred - E_icsd * N_sf / M_icsd) / 2 / A
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
                magnitudes = []
                Energies = []
                Energies_pred = []
                for groupid in group:
                    Natom = len(
                        data_collection[strain_id[groupid]]["structure"].species
                    )
                    magnitudes.append(
                        data_collection[strain_id[groupid]]["metadata"]["magnitude"]
                    )
                    Energies.append(
                        data_collection[strain_id[groupid]]["energy"] / Natom
                    )
                    Energies_pred.append(
                        data_collection[strain_id[groupid]]["pace"]["energy"] / Natom
                    )
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
                plt.legend(bbox_to_anchor=(0.5, -0.05))
                i += 1
            plt.savefig("Strains_" + proto + "_pace.png", dpi=200)
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
            value = 0
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
                    value = data["pace"]["energy"]

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
        plt.savefig("2d_" + faces[whichface] + ".png")
        plt.close()


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
def gmsf(data_collection):
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
        plt.savefig("2d_" + faces[whichface] + ".png")
        plt.close()


def Burgers_Bains(data_collection: list):
    """
    analysis for Burger_Bains pathways connecting fcc-bcc-hcp
    """
    B_B_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "Burgers_Bains"
            and data_collection[i]["calc"] == "final"
        ):
            B_B_id.append(i)
    bccref_pace = 0
    for data in data_collection:
        if (
            data["metadata"]["perturbation"] == "icsd"
            and data["calc"] == "final"
            and data["metadata"]["proto"] == "bcc.vasp"
        ):
            bccref_pace = data["pace"]["energy"] / len(data["structure"])
    bccref_pace
    unshuffled_values = []
    unshuffled_energies = []
    shuffled_values = []
    shuffled_energies = []
    for id in B_B_id:
        if "sh" in data_collection[id]["metadata"]["value"]:
            symbol = 1 if "p" in data_collection[id]["metadata"]["value"] else -1
            shuffled_values.append(
                symbol * float(data_collection[id]["metadata"]["value"].split("_")[-1])
            )
            shuffled_energies.append(
                float(data_collection[id]["pace"]["energy"])
                / len(data_collection[id]["structure"])
                - bccref_pace
            )
        else:
            symbol = 1 if "p" in data_collection[id]["metadata"]["value"] else -1
            unshuffled_values.append(
                symbol * float(data_collection[id]["metadata"]["value"].split("_")[-1])
            )
            unshuffled_energies.append(
                float(data_collection[id]["pace"]["energy"])
                / len(data_collection[id]["structure"])
                - bccref_pace
            )
    unshuffled_values = np.array(unshuffled_values)
    unshuffled_energies = np.array(unshuffled_energies)
    shuffled_values = np.array(shuffled_values)
    shuffled_energies = np.array(shuffled_energies)
    unshuffled_index = np.argsort(unshuffled_values)
    shuffled_index = np.argsort(shuffled_values)
    unshuffled_values = unshuffled_values[unshuffled_index]
    unshuffled_energies = unshuffled_energies[unshuffled_index]
    shuffled_values = shuffled_values[shuffled_index]
    shuffled_energies = shuffled_energies[shuffled_index]
    with plt.style.context("science"):
        plt.figure()
        plt.scatter(shuffled_values, shuffled_energies, s=6, c="r")
        plt.scatter(unshuffled_values, unshuffled_energies, s=6, c="b")
        plt.ylabel("Energy (eV/atom)", fontsize=12)
        plt.tick_params(axis="y", direction="in")
        plt.tick_params(axis="x", direction="in")
        plt.savefig("Ti_Burgers_strain.png")
        plt.close()


def bcc_omega(data_collection: list):
    """
    analysis for the pathway connecting bcc and omega
    """
    b_o_id = []
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "bcc_omega"
            and data_collection[i]["calc"] == "final"
        ):
            b_o_id.append(i)
    values = []
    Energies = []
    for id in b_o_id:
        values.append(data_collection[id]["metadata"]["value"])
        Energies.append(
            data_collection[id]["pace"]["energy"]
            / len(data_collection[id]["structure"])
        )
    values = np.array(values)
    Energies = np.array(Energies)
    index = np.argsort(values)
    values = values[index]
    Energies = Energies[index]
    with plt.style.context("science"):
        plt.figure()
        plt.plot(np.arange(len(Energies)), Energies - bp, c="black")
        plt.scatter(np.arange(len(Energies)), Energies - bp, c="black", s=4)

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
    hcp_energy = 0
    omega_energy = 0
    for i in np.arange(len(data_collection)):
        if (
            data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
            and data_collection[i]["metadata"]["proto"] == "hcp"
        ):
            hcp_energy = data_collection[i]["pace"]["energy"] / len(
                data_collection[i]["structure"]
            )
            omega_energy = data_collection[i]["pace"]["energy"] / len(
                data_collection[i]["structure"]
            )
    pathways = set()
    for id in h_o_static_id:
        pathways.add(data_collection[id]["metadata"]["#pathway"])
    pathways = list(pathways)

    for pathway in pathways:
        energies = np.zeros(7)
        energies[0] = omega_energy
        energies[6] = hcp_energy
        for image in np.arange(1, 6, 1):
            for i in np.arange(len(data_collection)):
                if (
                    "NEB" in data_collection[i]["metadata"]
                    and data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
                    and data_collection[i]["metadata"]["NEB"] == "True"
                    and data_collection[i]["metadata"]["#pathway"] == pathway
                    and int(data_collection[i]["metadata"]["image"]) == image
                ):
                    energies[image] = data_collection[i]["pace"]["energy"]
        energies = energies - energies[0]
        with plt.style.context("science"):
            plt.figure()
            plt.plot(np.arange(7), energies)
            plt.savefig(str(pathway) + "_NEB.png")
            plt.close()

        energies = np.zeros(7)
        energies[0] = omega_energy
        energies[6] = hcp_energy
        for image in np.arange(1, 6, 1):
            for i in np.arange(len(data_collection)):
                if (
                    "NEB" in data_collection[i]["metadata"]
                    and data_collection[i]["metadata"]["perturbation"] == "hcp_omega"
                    and data_collection[i]["metadata"]["NEB"] == "False"
                    and data_collection[i]["metadata"]["#pathway"] == pathway
                    and int(data_collection[i]["metadata"]["image"]) == image
                ):
                    energies[image - 1] = data_collection[i]["pace"]["energy"]
        energies = energies - energies[0]
        with plt.style.context("science"):
            plt.figure()
            plt.plot(np.arange(7), energies)
            plt.savefig(str(pathway) + "_static.png")
            plt.close()


def collect_in_powerpoint(RMSE_E: float, MAE_E: float, RMSE_F: float, MAE_F: float):
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

    ############################### hcp 2 omega ################################

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
            ref_omega_pace = data["energy"]
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
    energies_dft = np.zeros(len(data_collection))
    forces_dft = []
    energies_pace = np.zeros(len(data_collection))
    forces_pace = []
    for i in tqdm(np.arange(len(data_collection))):
        atom_num = len(data_collection[i]["structure"].species)
        energies_dft[i] = data_collection[i]["energy"] / atom_num - ref_omega_dft
        forces_dft += data_collection[i]["forces"]
        energies_pace[i] = (
            data_collection[i]["pace"]["energy"] / atom_num - ref_omega_pace
        )
        forces_pace += data_collection[i]["pace"]["forces"].tolist()

    with plt.style.context("science"):
        plt.figure()
        plt.scatter(energies_dft, energies_pace, c="r", s=0.3)
        plt.legend()
        plt.xlabel(r"$\Delta E^{DFT} \enspace (eV/atom)$")
        plt.ylabel(r"$\Delta E^{pred} \enspace (eV/atom)$")
        plt.xlim([0, 6])
        plt.ylim([0, 6])
        plt.savefig("pace_energy.png", dpi=200)
        plt.close()
    RMSE_E = np.sqrt(mean_squared_error(energies_dft, energies_pace))
    MAE_E = mean_absolute_error(energies_dft, energies_pace)

    forces_sca_dft = np.sqrt(np.sum(np.array(forces_dft) ** 2, axis=1))
    forces_sca_pace = np.sqrt(np.sum(np.array(forces_pace) ** 2, axis=1))

    with plt.style.context("science"):
        plt.figure()
        plt.scatter(forces_sca_dft, forces_sca_pace, c="r", s=0.3, label="pace")
        plt.xlabel(r"$|F^{i}|^{DFT} \enspace (eV/Å)$")
        plt.ylabel(r"$|F^{i}|^{pred} \enspace (eV/Å)$")
        plt.savefig("pace_force.png", dpi=200)
        plt.close()

    RMSE_F = np.sqrt(mean_squared_error(forces_sca_dft, forces_sca_pace))
    MAE_F = mean_absolute_error(forces_sca_dft, forces_sca_pace)
    Vacancy_formation(data_collection)
    Grain_boundary(data_collection)
    Surface(data_collection)
    Strain(data_collection, ref_omega_dft, ref_omega_pace)
    Gsfe_withoutUber(data_collection)
    Burgers_Bains(data_collection)
    bcc_omega(data_collection)
    hcp_omega(data_collection)
    Strain(data_collection, ref_omega_dft, ref_omega_pace)
    collect_in_powerpoint(RMSE_E, MAE_E, RMSE_F, RMSE_F)


if __name__ == "__main__":
    main()
