#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 16:47:11 2020

@author: attari.v
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import re

#font = {'family': 'serif',
#        'color':  'darkred',
#        'weight': 'normal',
#        'size': 16,
#        }

params = {"text.color" : "white"}
plt.rcParams.update(params)

def read_vahid_file(filename):
    file_obj = open(filename,'r')
    zone_dict = {}
    array_list = []
    for line in file_obj:
        if 'variables' in line:
            re_finds_obj = re.findall(r'["](.*?)["]',line)
            headers = [find for find in re_finds_obj]
            
        elif 'ZONE' in line:
            # pass
            re_finds_obj = re.findall(r'[I][=][\s]*[0-9]+',line)
            I_val = int(re_finds_obj[0].split('=')[1].strip())
            #print(re_finds_obj)
            re_finds_obj = re.findall(r'[J][=][\s]*[0-9]+',line)
            J_val = int(re_finds_obj[0].split('=')[1].strip())
            zone_dict['I'] = I_val
            zone_dict['J'] = J_val
        elif line:
            line = line.replace('D','E')
            line = re.sub('^[\s]+','',line)
            line = re.sub('[\s]+',' ', line)
            array_list.append(np.fromstring(line,sep = ' '))
        else:
            pass
    #variables =  "IG"  "JG"  "PHIMAX"  "C1"  "PHI1_2"
    #file_obj.readlines()
    file_obj.close()
    df_out = pd.DataFrame(np.array(array_list),columns=headers)

    return df_out, zone_dict

def read_mts_pars(filename):
    file_obj = open(filename,'r')
    zone_dict = {}
    array_list = []

    with open(filename) as f_in:
        lines = (line.rstrip() for line in f_in) 
        lines = list(line for line in lines if line) # Non-blank lines in a list
    
    for line in lines:
        if 'variables' in line:
            re_finds_obj = re.findall(r'["](.*?)["]',line)
            headers = [find for find in re_finds_obj]
        elif line:
            l = [x.strip() for x in line.split('=',2)]
            array_list.append(l)
        else:
            pass
    file_obj.close
    df_out = pd.DataFrame(array_list)
    df_out = pd.DataFrame(df_out).T
    df_out.columns = df_out.iloc[0]
    df_out = df_out.iloc[1: , :]
    df_out = df_out.astype(float)

    return df_out

## plot 2 --- phase-field order parameter
filename = sys.argv[-1]
filename = filename;
df2, zd = read_vahid_file(filename)

#df2, _ = read_vahid_file('../images/phi_con000000000.plt')
#df2, _ = read_vahid_file('../results/flux_pars000000000.plt')

a = df2.PHI.values
phi = np.reshape(a,[zd['J'],zd['I']])
phi = np.transpose(phi)
phi = np.flip(phi)
#########

fig3, ax3 = plt.subplots()

axins1 = inset_axes(ax3,
                    width="40%",  # width = 50% of parent_bbox width
                    height="2%",  # height : 5%
                    loc='upper right')

im1 = ax3.imshow(phi, cmap='binary_r', interpolation='nearest')
ax3.set(xticks=[], yticks=[])
fig3.colorbar(im1, cax=axins1, orientation="horizontal", ticks=[0, 0.5, 1])
axins1.tick_params(labelsize=7) 
axins1.tick_params(colors='white') 
axins1.tick_params(grid_color='white') 


#ax.imshow(weights, **imshow_kwargs)
plt.savefig(filename+'_order_p.jpg',dpi=200, bbox_inches='tight',pad_inches = 0)




