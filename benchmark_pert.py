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

dataset=pd.read_pickle(glob.glob('always_include_*_*gzip')[0],compression='gzip')


##############Burgers#################
print("Working on Burgers_Bain")
shuffled_referenced,unshuffled_referenced=benchmarks.burgers.benchmark_bain_burgers_pathway(dataset)
keys=list(shuffled_referenced.keys())
for key in [x for x in keys if x.startswith('pace_energy')]:
    with plt.style.context("science"):
        fig=plt.figure(figsize=(5, 4))
        ax=fig.add_subplot(111)
        ax=benchmarks.burgers.plot_burgers_bain(ax,shuffled_referenced,unshuffled_referenced,key)
        ax.set_xlabel("Strain Magnitude")
        ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
        ax.tick_params(axis='both', which='both', direction='in')
        plt.legend()
        plt.savefig("Burgers/Ti_Burgers_strain"+"s"+key.split('_')[-1]+".png")
        plt.close()

##############hcp_omega################
print("Working on hcp_omega")
benchmarks.hcpomega.benchmark_hcp_omega_pathway(dataset)

##############omega#################
print("Working on bcc_omega")
benchmarks.omega.benchmark_bcc_omega_pathway(dataset)

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
dataset=pd.read_pickle(glob.glob('all_perts*gzip')[0],compression='gzip')
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

