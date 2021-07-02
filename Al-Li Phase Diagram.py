#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 12:07:39 2021

@author: seblohier
"""




import matplotlib.pyplot as plt
from pycalphad import Database, binplot, variables as v

# Load database and choose the phases that will be plotted
db_nbre = Database('Al-Li.tdb')
my_phases_nbre = ['LIQUID','BCC_A2','FCC_A1','BCC_DIS','BCC_B2','AL2LI3','AL4LI9']

# Create a matplotlib Figure object and get the active Axes
fig = plt.figure(figsize=(9,6))
axes = fig.gca()

# Plot the phase diagram on the existing axes using the `plot_kwargs={'ax': axes}` keyword argument
binplot(db_nbre, ['AL', 'LI'] , my_phases_nbre, {v.X('LI'): (0,1,0.01), v.T: (300, 1000, 10), v.P:101325}, plot_kwargs={'ax': axes})


plt.show()