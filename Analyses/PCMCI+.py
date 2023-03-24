#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import pyarrow.feather as feather
from matplotlib import pyplot as plt
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr
import numpy as np

data_all = pd.read_csv("20230222_COPD_TTTF.csv")
data_all.rename(columns={'copd': r'$COPD$',
                           'co': r'$CO$',
                          'o3h8max': r'$O_{3}$',
                          'fsp': r'$PM_{2.5}$',
                          'no2': r'$NO_{2}$',
                         'temp': r'$Temp.$',
                         'rh': r'$Humid.$'}, inplace=True)
data_all



# Function for preparing the values of the variables specified by 'var_names'
def get_data(data_all,  var_names):

    data_out = data_all

    # Select the columns
    data_out = data_out[var_names]

    # Turn into numpy array
    data_out = data_out.values

    # Return
    return data_out
    # Return
    return data_out


def get_data_and_mask(data_all, season, var_names):
    if season == "Year":
        data_out = get_data(data_all,  var_names)
        return data_out, None


    data_with_season = data_all.copy()

    # Select the state and columns, including 'season'
    data_with_season = data_with_season[var_names + [r'$Season$']]

    # Data without season
    data_out = data_with_season[var_names]
    data_out = data_out.values

    # Mask
    if season == "Summer":
        # Mask the winter months ('True' means masked, 'False' means not masked)
        data_with_season.loc[data_with_season[r'$Season$'] == 0, var_names] = True
        data_with_season.loc[data_with_season[r'$Season$'] == 1, var_names] = False
    elif season == "Winter":
        # Mask the sommer months ('True' means masked, 'False' means not masked)
        data_with_season.loc[data_with_season[r'$Season$'] == 0, var_names] = False
        data_with_season.loc[data_with_season[r'$Season$'] == 1, var_names] = True
    else:
        raise ValueError("Season must be in ['Year', 'Summer', 'Winter'].")
    mask_out = data_with_season[var_names]
    mask_out = mask_out.values

    # Return
    return data_out, mask_out

def get_selected_links(var_names, tau_min, tau_max):
    
    # Get index of the season variable, if it exists
    if r'$Season$' in var_names:
        season_idx = np.argwhere(np.array(var_names) == r'$Season$')[0, 0]
    else:
        season_idx = None

    # Get index of the temperature variable, if it exists
    if r'$Temp.$' in var_names:
        temp_idx = np.argwhere(np.array(var_names) == r'$Temp.$')[0, 0]
    else:
        temp_idx = None

    # Get index of the humidity variable, if it exists
    if r'$Humid.$' in var_names:
        humid_idx = np.argwhere(np.array(var_names) == r'$Humid.$')[0, 0]
    else:
        humid_idx = None

    # Build dictionary
    selected_links = {}
    
    for idx, var in enumerate(var_names):

        if var == r'$COPD$':
            # COPD may be influenced by all variables other than season at lags
            selected_links[idx] = [(other_idx, -tau) for other_idx, other_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1) if other_var != r'$Season$']
            
        elif var == r'$PM_{2.5}$':
            # pm2.5 may be influenced by all variables other than copd
            selected_links[idx] = [(other_idx, -tau) for other_idx, other_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1) if other_var != r'$COPD$']
            
            
        elif var == r'$Humid.$':
            # Humiditiy may be influenced by itself at all lags
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]

            # Humidity may also be influenced by temperature
            selected_links[idx] = [(temp_idx, -tau) for tau in range(max(1, tau_min), tau_max + 1)]

        elif var == r'$Temp.$':
            # Temperature may be influenced by itself at all lags
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]

            # Temperature may also be influenced by humidity 
            selected_links[idx] = [(humid_idx, -tau) for tau in range(max(1, tau_min), tau_max + 1)]

        elif var == r'$Season$':
            # Season may be influenced by itself at all lags
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]

        else:
            # All other variables,may be influenced by all variables other than
            # COPD and pm2.5
            selected_links[idx] = [(other_idx, -tau) for other_idx, other_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1) if
                                   ((other_var != r'$COPD$') & (other_var != r'$PM_{2.5}$'))]
            
        # Season may influence all variables at lag tau_min
        if var != r'$Season$' and season_idx is not None:
            selected_links[idx].append((season_idx, -tau_min))

    # Return
    return selected_links

def apply_pcmci(data_all,
                var_names,
                season,
                tau_min,
                tau_max,
                pc_alpha,
                verbosity,
               name):
    # Get the data and mask
    data, mask = get_data_and_mask(data_all=data_all,
                                   season=season,
                                   var_names=var_names)

    # Prepare the DataFrame object
    dataframe = pp.DataFrame(data,
                             mask=mask,
                             var_names=var_names,
                             missing_flag=999.)

    # Prepare the independence test and PCMCI object
    if season == "Year":
        parcorr = ParCorr()
    else:
        parcorr = ParCorr(mask_type='y')
    pcmci = PCMCI(dataframe=dataframe,
                  cond_ind_test=parcorr,
                  verbosity=verbosity)

    # Get the selected_links arguement
    selected_links = get_selected_links(var_names,
                                        tau_min,
                                        tau_max)

    # Run PCMCI^+ with these parameters
    results = pcmci.run_pcmciplus(tau_min=tau_min,
                                  tau_max=tau_max,
                                  pc_alpha=pc_alpha,
                                  selected_links=selected_links)

    tp.plot_graph(
        arrow_linewidth=7.0,
        figsize=(10, 5),
        vmin_edges=-0.5,
        vmax_edges=0.5,
        node_label_size=13,
        link_label_fontsize=7,
        val_matrix=results['val_matrix'],
        graph=results['graph'],
        var_names=var_names,
        link_colorbar_label='cross-MCI (edges)',
        node_colorbar_label='auto-MCI (nodes)',
        node_aspect=1,
        node_size=0.35,
        label_fontsize=11,
        network_lower_bound=0.2,
        show_colorbar=1
    );

    plt.suptitle(name,
                 size=14,
                 weight="semibold",
                 color="gray",
                 verticalalignment="top",
                 horizontalalignment="center")
    plt.show()
    
    #Plot time series graph    
    tp.plot_time_series_graph(
    figsize=(12, 6),
    node_size=0.05,
    vmin_edges=-0.5,
    vmax_edges=0.5,
    val_matrix=results['val_matrix'],
    graph=results['graph'],
    var_names=var_names,
    link_colorbar_label='MCI',
    )
    plt.show()

    return results


# In[2]:


pc_alpha = 0.001
results = apply_pcmci(data_all=data_all,
                      var_names=[r'$COPD$', r'$O_{3}$',r'$PM_{2.5}$',r'$NO_{2}$',r'$CO$',
                                 r'$Temp.$',r'$Humid.$'
                                    ],
                      season="Year",
                      tau_min=0,
                      tau_max=6,
                      pc_alpha=pc_alpha,
                      verbosity=0,
                      name="Causal network"
                      )
print(pc_alpha)


# In[3]:


pc_alpha = 0.01
results = apply_pcmci(data_all=data_all,
                      var_names=[r'$COPD$', r'$O_{3}$',r'$PM_{2.5}$',r'$NO_{2}$',r'$CO$',
                                 r'$Temp.$',r'$Humid.$'
                                    ],
                      season="Year",
                      tau_min=0,
                      tau_max=6,
                      pc_alpha=pc_alpha,
                      verbosity=0,
                      name="Causal network"
                      )
print(pc_alpha)

