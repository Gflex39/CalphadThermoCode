#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:30:26 2021

@author: seblohier
"""
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import background
import csv
import os
import pandas as pd
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np
import alphashape
from scipy.interpolate import interp1d
import itertools
def potfinder(temper,plot):
    # Load database and choose the phases that will be plotted
    db_nbre = Database('Al-Li.tdb')
    my_phases_nbre = ['LIQUID','BCC_A2','FCC_A1','BCC_DIS','BCC_B2','AL2LI3','AL4LI9']
    testx= [i for i in np.around(np.arange(0.01,.99,.0001),2)]
    # Get the colors that map phase names to colors in the legend
    legend_handles, color_dict = phase_legend(my_phases_nbre)

    fig = plt.figure(figsize=(9,6))
    ax = fig.gca()
    DataValues=[]
    df = pd.DataFrame(list())
    phase_potentials={}
    hullval=[]
    # Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
    for phase_name in my_phases_nbre:
        result = calculate(db_nbre, ['AL', 'LI','VA'], phase_name, P=101325, T=temper, output='GM')
        print(phase_name)

        temp=result.GM.values
        tempx=result.X.sel(component='LI').values

    
        x = tempx.squeeze() 
        y = temp.squeeze()
        # ax.scatter(x,y)
        #Here the (x,y) pairs of points are but into a more simple and callable form
        #We use a try except case here because one of the phases has only one data point. Meaning when
        #try to ge the length it gives a type eror
        try:
            #This portion of the code needs to be changed but all it does is get the hull of the gibbs surface
            if phase_name=='BCC_B2':
                minihull=[]
                hullx=[]
                hully=[]
                count=0
                
            for i in range(len(x)):
                app=[x[i],y[i]]
                if phase_name == 'BCC_B2':
                    if i%400==0:
                        if x[i] not in hullx:
                            count+=1
                            minihull.append(app)

                            hullx.append(x[i])
                            hully.append(y[i])
                    if count>=100:
                        break
                hullval.append(app)
            if phase_name=='BCC_B2':
                minihull=np.array(minihull)
                alpha = 0.95 * alphashape.optimizealpha(minihull)
                print("optimized")
                hull = alphashape.alphashape(minihull, alpha)
                hull_pts = hull.exterior.coords.xy
                
                x=hull_pts[0]
                y=hull_pts[1]
            #From this point all the code makes a polynomial estimate of the specific phase curve 
            #then find the derivative
            if phase_name=='BCC_B2':
                x, y = zip(*sorted(zip(x, y)))
                f = interp1d(x, y)
                estimate=np.polyfit(testx,f(testx),8)
                testy=np.polyval(estimate,testx)
                print(len(testx))
                print(len(testy))
                plt.plot(testx,testy, color=color_dict[phase_name],linewidth=4)
            else:
                estimate=np.polyfit(x,y,3)
            # testy=np.polyval(estimate,testx)
            # plt.plot(testx,testy)
            derivative=background.find_derivative(estimate)
            alphapot=[]
            betapot=[]
            for i in range(len(x)):
                #This derivative is used to find the potential of the phase at each composition
                slope=derivative(x[i])
                potential=background.find_potentials(slope, x[i], y[i], True)
                alphapot.append(potential[0])
                betapot.append(potential[1])
            print(len(x))
            phase_potentials[phase_name]={'alpha':alphapot,'beta':betapot,'x':x}
        except TypeError:
            print("We good")
        #The points are the plotted
        if plot:
            if phase_name!='BCC_B2':
                ax.scatter(x,y, marker='.', s=5, color=color_dict[phase_name])

        #This commented portion takes the values attained and turns them into a csv file labeling what it is
        # background.csvform("Al", "Li", phase_name, temper)
            

    # Format the plot
    if plot:
        ax.set_xlabel('X(Li)')
        ax.set_ylabel('GM')
        ax.set_xlim((0, 1))
        ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))

        pplot=plt.figure(2)
        ax2=pplot.gca()
        phase='BCC_A2'

        ax2.set_xlabel('Mole Percent of Li',fontsize=16)
        ax2.set_ylabel('Chemical Potential')
        ax2.set_title('Chemical Potential of the '+phase+' phase at '+str(temper)+" K")
        ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['alpha'],color='red',label='Pure Al potential')
        ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['beta'],color='blue',label='Pure Li potential')
        ax2.legend(loc='best', scatterpoints=100)

        plt.show()
    return phase_potentials

temps=[]
pots=[]

for i in range(1216,1226):
    temps.append(i)
    curr=potfinder(i, False)
    minimum=min(curr['LIQUID']['x'])
    for j in range(len(curr['LIQUID']['x'])):
        if curr['LIQUID']['x'][j]==minimum:
            if i==1221:
                Gibb=curr['LIQUID']['alpha'][j]
            pots.append(curr['LIQUID']['alpha'][j])
            break
potestimate=np.polyfit(temps,pots,4)
potderivative=background.find_derivative(potestimate)
entropy=-potderivative(1221)
enthalpy=Gibb+entropy*1221
print("Enthalpy="+str(enthalpy))
print("Entropy="+str(entropy))
