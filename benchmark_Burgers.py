import pandas as pd
from benchmark_utils_pd import *
import glob
import scienceplots
from matplotlib import pyplot as plt
import pandas as pd


nums=len(glob.glob('interim*ladder_step*yaml'))
dataset=pd.read_pickle(glob.glob('always_include_*_*gzip')[0],compression='gzip')
shuffled_referenced,unshuffled_referenced=burgers_bain(dataset)

keys=list(shuffled_referenced.keys())
for key in [x for x in keys if x.startswith('pace_energy')]:
    with plt.style.context("science"):
        fig=plt.figure(figsize=(5, 4))
        ax=fig.add_subplot(111)
        ax=plot_burgers_bain(ax,shuffled_referenced,unshuffled_referenced,key)
        ax.set_xlabel("Strain Magnitude")
        ax.set_ylabel(r"$\Delta E \enspace (eV/atom)$")
        ax.tick_params(axis='both', which='both', direction='in')
        plt.legend()
        plt.savefig("Ti_Burgers_strain"+"s"+key.split('_')[-1]+".png")
        plt.close()
