import benchmarks
import pandas as pd
import glob
import scienceplots
from matplotlib import pyplot as plt
import pandas as pd
import os

##############Burgers#################
dataset=pd.read_pickle(glob.glob('always_include_*_*gzip')[0],compression='gzip')
shuffled_referenced,unshuffled_referenced=benchmarks.benchmark_bain_burgers_pathway(dataset)

keys=list(shuffled_referenced.keys())
for key in [x for x in keys if x.startswith('pace_energy')]:
    with plt.style.context("science"):
        fig=plt.figure(figsize=(5, 4))
        ax=fig.add_subplot(111)
        ax=benchmarks.plot_burgers_bain(ax,shuffled_referenced,unshuffled_referenced,key)
        ax.set_xlabel("Strain Magnitude")
        ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
        ax.tick_params(axis='both', which='both', direction='in')
        plt.legend()
        plt.savefig("Burgers/Ti_Burgers_strain"+"s"+key.split('_')[-1]+".png")
        plt.close()

##############gsfe#################
try:
    os.mkdir('linescan')
except:
    pass

gsfedata=benchmarks.get_pert_data(dataset,'gsfe')
benchmarks.gsfe(gsfedata)

##############hcp_omega################
try:
    os.mkdir('hcp_omega')
except:
    pass

benchmarks.benchmark_hcp_omega_pathway(dataset)

##############omega#################
try:
    os.mkdir('Omega')
except:
    pass

benchmarks.benchmark_bcc_omega_pathway(dataset)

##############omega#################
try:
    os.mkdir('Phonons')
except:
    pass

benchmarks.phonons(dataset)

##############strain#################
try:
    os.mkdir('Strains')
except:
    pass

benchmarks.strain(dataset)

##############gb#################
try:
    os.mkdir('gb')
except:
    pass

dataset=pd.read_pickle(glob.glob('all_perts*gzip')[0],compression='gzip')
benchmarks.grain_boundary(dataset)

##############surfaces#################
try:
    os.mkdir('Surfaces')
except:
    pass

dataset=pd.read_pickle(glob.glob('all_perts*gzip')[0],compression='gzip')
benchmarks.surface(dataset)

##############Vacancies#################
try:
    os.mkdir('Vacancies')
except:
    pass

benchmarks.vacancy_formation(dataset)

##############Volumetric#################
try:
    os.mkdir('Volumetric')
except:
    pass

benchmarks.volumetric_deformation_df(dataset)

