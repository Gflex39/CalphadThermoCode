#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 12:07:39 2021

@author: seblohier
"""




import matplotlib.pyplot as plt
from pycalphad import Database, binplot, variables as v

# Load database and choose the phases that will be plotted
db_nbre = Database('Fe-Co.tdb')
my_phases_nbre = ['LIQUID','B2_BCC','FCC_A1','HCP_A3']

# Create a matplotlib Figure object and get the active Axes
fig = plt.figure(figsize=(9,6))
axes = fig.gca()

# Plot the phase diagram on the existing axes using the `plot_kwargs={'ax': axes}` keyword argument
binplot(db_nbre, ['Co', 'Fe','VA'] , my_phases_nbre, {v.X('FE'): (0,1,0.01), v.T: (0, 2000, 20), v.P:101325}, plot_kwargs={'ax': axes})
axes.set_title("Phase diagram of Cu-Ni System", fontsize=16)
axes.set_xlabel("X(Ni)")


plt.show()