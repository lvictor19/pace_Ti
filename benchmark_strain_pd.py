import pandas as pd
from benchmark_utils_strain_pd import *
import glob
import scienceplots
from matplotlib import pyplot as plt
import pandas as pd
import os

try:
    os.mkdir('Strains')
except:
    pass

dataset=pd.read_pickle(glob.glob('always_include_*_*gzip')[0],compression='gzip')
strain(dataset)