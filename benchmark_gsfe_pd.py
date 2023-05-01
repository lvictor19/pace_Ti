import pandas as pd
from benchmark_utils_gsfe_pd import *
import glob
import scienceplots
from matplotlib import pyplot as plt
import pandas as pd
import os

try:
    os.mkdir('linescan')
except:
    pass

dataset=pd.read_pickle(glob.glob('always_include_*_*gzip')[0],compression='gzip')
gsfedata=get_pert_data(dataset,'gsfe')
gsfe(gsfedata)