import pandas as pd
from benchmark_utils_volumetric_pd import *
import glob
import scienceplots
from matplotlib import pyplot as plt
import pandas as pd
import os

try:
    os.mkdir('Volumetric')
except:
    pass

dataset=pd.read_pickle(glob.glob('*+volumetric*gzip')[0],compression='gzip')
volumetric_deformation_df(dataset)