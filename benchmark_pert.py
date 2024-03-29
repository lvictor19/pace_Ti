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
#############Burgers#################
def burgers_bain():
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
def bcc_omega():
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

# ##############hcp_omega################
def hcp_omega():
    print("Working on hcp_omega")
    ssNEB,static=benchmarks.hcpomega.benchmark_hcp_omega_pathway(dataset)
    pathways=list(set(static["pathway_number"]))
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
def phonons():
    print("Working on phonons")
    phonon_data=benchmarks.phonons.phonons(dataset)
    protos=set(phonon_data['metadata'].map(lambda x: x['proto']))
    keys=list(phonon_data.keys())
    keys=[x for x in keys if x.startswith('pace_forces')]
    for proto in protos:
        band_fig=benchmarks.phonons.analysis(phonon_data,proto,'forces')
        band_fig.savefig(os.path.join("Phonons","phonon_{0}.png".format(proto)), transparent=True)
        plt.close()
        for key in keys:
            band_fig=benchmarks.phonons.analysis(phonon_data,proto,key)
            band_fig.savefig(os.path.join("Phonons","phonon_{0}.png".format(proto)), transparent=True)
            plt.close()


##############gsfe#################
def gsfe():
    print("Working on gsfe")
    benchmarks.gsfe.gsfe(dataset)

##############strain#################
def strains():
    print("Working on strains")
    benchmarks.strain.strain(dataset)


##############gb#################
def gb():
    print("Working on grain boundaries")
    gb_data=benchmarks.gb.grain_boundary(dataset)
    for key in keys:
        with plt.style.context("science"):
            fig = plt.figure(figsize=(5, 6))
            ax = fig.add_subplot(111, frame_on=False)
            df=benchmarks.gb.analysis(gb_data,dataset,key)
            ax=benchmarks.gentable(ax,df)
            plt.savefig("gb/Grain_boundary_s{0}.png".format(key.split('_')[-1]),dpi=200)
            plt.close()

#############surfaces#################
def surfaces():
    print("Working on surfaces")
    surface_data=benchmarks.surfaces.surfaces(dataset)
    for key in keys:
        with plt.style.context("science"):
            fig = plt.figure(figsize=(5, 6))
            ax = fig.add_subplot(111, frame_on=False)
            df=benchmarks.surfaces.analysis(surface_data,dataset,key)
            ax=benchmarks.gentable(ax,df)
            plt.savefig("Sufaces/Surfaces_s{0}.png".format(key.split('_')[-1]),dpi=200)
            plt.close()

#############Vacancies#################
def vacancies():
    print("Working on vacancies")
    benchmarks.vacancies.vacancy_formation(dataset)

##############Volumetric#################
def volumetric():
    print("Working on volumetric deformations")
    benchmarks.volumetric.volumetric_deformation(dataset)

if __name__=="__main__":
    burgers_bain()
    bcc_omega()
    hcp_omega()
    phonons()
    gsfe()
    strains()
    gb()
    surfaces()
    vacancies()
    volumetric()


