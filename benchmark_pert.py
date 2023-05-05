import benchmarks
import pandas as pd
import glob
import scienceplots
from matplotlib import pyplot as plt
import pandas as pd
import os

statements=['os.mkdir("Burgers")',
            'os.mkdir("hcp_omega")',
            'os.mkdir("Omega")',
            'os.mkdir("Phonons")',
            'os.mkdir("linescan")',
            'os.mkdir("Strains")',
            'os.mkdir("Surfaces")',
            'os.mkdir("Vacancies")',
            'os.mkdir("Volumetric")',
            'os.mkdir("gb")']

for statement in statements:
    try:
        exec(statement)
    except Exception:
        pass

dataset=pd.read_pickle(glob.glob('all_perts*gzip')[0],compression='gzip')
keys=list(dataset.keys())
keys=[x for x in keys if x.startswith('pace_energy')]
##############Burgers#################
print("Working on Burgers_Bain")
shuffled_referenced,unshuffled_referenced=benchmarks.burgers.benchmark_bain_burgers_pathway(dataset)
for key in keys:
    with plt.style.context("science"):
        fig=plt.figure(figsize=(5, 4))
        ax=fig.add_subplot(111)
        ax=benchmarks.burgers.plot_burgers_bain(ax,shuffled_referenced,unshuffled_referenced,key)
        ax.set_xlabel("Strain Magnitude")
        ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
        ax.tick_params(axis='both', which='both', direction='in')
        plt.legend()
        plt.savefig("Burgers/Ti_Burgers_strains{0}.png".format(key.split('_')[-1]))
        plt.close()

##############omega#################
print("Working on bcc_omega")
bcc_omega_data=benchmarks.omega.benchmark_bcc_omega_pathway(dataset)
for key in keys:
    with plt.style.context("science"):
        fig=plt.figure(figsize=(5, 4))
        ax=fig.add_subplot(111)
        ax=benchmarks.omega.plot_omega(ax,bcc_omega_data,key)
        ax.set_xlabel("Omega Transformation")
        ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
        ax.tick_params(axis='both', which='both', direction='in')
        ax.set_xticks([0,0.167],['bcc','omega'])
        plt.legend()
        plt.savefig("Omega/Omegas{0}.png".format(key.split('_')[-1]))
        plt.close()

##############hcp_omega################
print("Working on hcp_omega")
ssNEB,static=benchmarks.hcpomega.benchmark_hcp_omega_pathway(dataset)
pathways=list(set(static["pathway_number"]))
keys=list(static.keys())
keys=[x for x in keys if x.startswith('pace_energy')]
for key in keys:
    for pathway in pathways:
        fig=plt.figure(figsize=(5, 4))
        ax=fig.add_subplot(111)
        ax=benchmarks.hcpomega.plot_hcp_omega(ax,ssNEB,static,pathway,key)
        ax.set_xlabel("Transformation Coordinates")
        ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
        ax.tick_params(axis='both', which='both', direction='in')
        ax.set_xticks([0, 6], ["omega", "hcp"])
        plt.legend()
        plt.savefig("hcp_omega/{0}_s{1}.png".format(pathway,key.split('_')[-1]))
        plt.close()

##############phonons#################
print("Working on phonons")
benchmarks.phonons.phonons(dataset)

##############gsfe#################
print("Working on gsfe")
gsfedata=benchmarks.get_pert_data(dataset,'gsfe')
benchmarks.gsfe.gsfe(gsfedata)

##############strain#################
print("Working on strains")
benchmarks.strain.strain(dataset)

##############gb#################
print("Working on grain boundaries")
benchmarks.gb.grain_boundary(dataset)

##############surfaces#################
print("Working on surfaces")
benchmarks.surfaces.surfaces(dataset)

#############Vacancies#################
print("Working on vacancies")
benchmarks.vacancies.vacancy_formation(dataset)

##############Volumetric#################
print("Working on volumetric deformations")
benchmarks.volumetric.volumetric_deformation(dataset)

