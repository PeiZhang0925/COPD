#!/usr/bin/env python
# coding: utf-8

# In[7]:


# Imports
import math, argparse
import numpy as np
import pandas as pd
from numpy.random.mtrand import sample
from tigramite import data_processing as pp
from tigramite.toymodels import structural_causal_processes as toys
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr,GPDC, CMIknn, CMIsymb
from causallearn.search.ConstraintBased.PC import pc
from causallearn.utils.cit import fisherz
from causallearn.search.ScoreBased.GES import ges
from causallearn.search.Granger.Granger import Granger
from matplotlib import pyplot as plt
import sys
sys.path.append("")
from multiprocessing import Pool
import igraph as ig
import random
from causalnex.structure.dynotears import from_pandas_dynamic
from causalnex.plots import plot_structure, NODE_STYLE, EDGE_STYLE
from IPython.display import Image



# In[8]:


data_all = pd.read_csv("20230222_COPD_TTTF.csv")

data_all = data_all.loc[:, ["copd","temp","rh","o3h8max","co","fsp","no2"]]
data_all.rename(columns={'copd': "COPD",
                          'o3h8max': 'O3',
                          'fsp': 'PM2.5',
                          'no2': 'NO2',
                         'co': 'CO',
                         'temp': 'Temp',
                         'rh': 'Humid'}, inplace=True)
df = data_all.to_numpy()

data_all


# In[9]:


g_learnt,w,a = from_pandas_dynamic(data_all,p=6,lambda_w=0.1,lambda_a=0.1,w_threshold=0.05)
res=pd.DataFrame(g_learnt.edges(data="weight"))
res=res[res[1].str.contains('COPD')]
res


# In[10]:


g_learnt,w,a = from_pandas_dynamic(data_all,p=6,lambda_w=0.1,lambda_a=0.1,w_threshold=0.01)
res=pd.DataFrame(g_learnt.edges(data="weight"))
res=res[res[1].str.contains('COPD')]
res


# In[ ]:




