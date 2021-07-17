#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 12:07:39 2021

@author: seblohier
"""




import matplotlib.pyplot as plt
from pycalphad import Database, binplot, variables as v

# Load database and choose the phases that will be plotted
db_nbre = Database('Mg-Bi.tdb')
my_phases_nbre = ['LIQUID','RHOMBO_A7','HCP_A3','BCT_A5','BI2MG3_L','BI2MG3_H','MG2SN']

# Create a matplotlib Figure object and get the active Axes
fig = plt.figure(figsize=(9,6))
axes = fig.gca()

# Plot the phase diagram on the existing axes using the `plot_kwargs={'ax': axes}` keyword argument
binplot(db_nbre, ['MG', 'BI','VA'] , my_phases_nbre, {v.X('MG'): (0,1,0.01), v.T: (300, 1200, 20), v.P:101325}, plot_kwargs={'ax': axes})

plt.show()