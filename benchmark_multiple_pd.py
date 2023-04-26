import pandas as pd
from benchmark_utils_pd import *

potfiles=glob.glob('*interim*ladder_step*yaml')
nums=[filename.split('.')[0].split('_')[-1] for filename in potfiles]
for i in nums:
    dataset=pd.read_pickle('Ti_data_pace'+i+'.pckl.gzip',compression='gzip')
    burgers_bain(dataset,i)



    







    

