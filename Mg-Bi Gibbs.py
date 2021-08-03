#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:30:26 2021

@author: seblohier
"""
import background
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import os
import csv
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np
import alphashape

def potfinder(temper,plot):
    # Load database and choose the phases that will be plotted
    db_nbre = Database('Mg-Bi.tdb')
    my_phases_nbre = ['LIQUID', 'BI2MG3_H', 'BI2MG3_L']

    # Get the colors that map phase names to colors in the legend
    legend_handles, color_dict = phase_legend(my_phases_nbre)
    hullval=[]
    fig = plt.figure(num=1,figsize=(9,6))
    ax = fig.gca()

    # Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
    phase_potentials={}
    for phase_name in my_phases_nbre:
        result = calculate(db_nbre, ['BI','MG','VA'], phase_name, P=101325, T=temper, output='GM')
        
        
            #ax.scatter(result.X.sel(component='MG'), result.GM, marker='.', s=5, color=color_dict[phase_name])
        
        print(phase_name)
        temp=result.GM.values
        tempx=result.X.sel(component='MG').values
        
        x = tempx.squeeze()     
        y = temp.squeeze()

        if phase_name=='LIQUID':
            newx=[]
            newy=[]
            for i in np.around(np.arange(0,1,.05),3):
                lowx=[]
                lowy=[]
                
                for j in range(len(x)):
                    if abs(x[j]-i)<0.0005:
                        lowx.append(x[j])
                        lowy.append(y[j])
                bot=min(lowy)
                newx.append(i)
                newy.append(bot)
            x=newx
            y=newy

        #Here the (x,y) pairs of points are but into a more simple and callable form
        #This portion of the code needs to be changed but all it does is get the hull of the gibbs surface
        for i in range(len(x)):
            app=[x[i],y[i]]
            hullval.append(app)
        #From this point all the code makes a polynomial estimate of the specific phase curve 
        #then find the derivative
        estimate=np.polyfit(x,y,8)
        derivative=background.find_derivative(estimate)
        alphapot=[]
        betapot=[]
        for i in range(len(x)):
            slope=derivative(x[i])
            potential=background.find_potentials(slope, x[i], y[i], True)
            alphapot.append(potential[0])
            betapot.append(potential[1])
        print(len(x))
        phase_potentials[phase_name]={'alpha':alphapot,'beta':betapot,'x':x}

        
    
        if plot:
            if phase_name!='LIQUID':
                ax.scatter(x,y, marker='.', s=5, color=color_dict[phase_name])
            else:
                plt.plot(x,y,color=color_dict[phase_name])


    if plot:    
            #This commented portion takes the values attained and turns them into a csv file labeling what it is
            # background.csvform("Mg", "Bi", phase_name, temper)
        #This finds the convex hull of the plot
        hullval=np.array(hullval)
        hull = ConvexHull(hullval)

        # plt.plot(hullval[:][0], hullval[:][1], 'o')
        #for simplex in hull.simplices:

            #plt.plot(hullval[simplex,0], hullval[simplex,1], 'k-')
        #plt.plot(hullval[hull.vertices,0], hullval[hull.vertices,1], 'r--', lw=2)
        #plt.plot(hullval[hull.vertices[0],0], hullval[hull.vertices[0],1], 'ro')

        # Format the plot   
        ax.set_xlabel('Mole Percent of Mg',fontsize=16)
        ax.set_ylabel('Molar Gibbs Energy',fontsize=16)
        ax.set_title(' Molar Gibbs Energy of Mg-Bi System at '+str(temper)+"\xb0"+" K")
        ax.set_xlim((0, 1))
        ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))

        pplot=plt.figure(2)
        ax2=pplot.gca()
        phase='LIQUID'
        ax2.set_xlabel('Mole Percent of Mg',fontsize=16)
        ax2.set_ylabel('Chemical Potential')
        ax2.set_title('Chemical Potential of the '+phase+' phase at '+str(temper)+" K")
        ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['alpha'],color='red',label='Potential of pure Mg ')
        ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['beta'],color='blue',label='Potential of pure Bi')
        ax2.legend(loc='best', scatterpoints=100)

        
    return phase_potentials
temps=[]
pots1=[]
pots2=[]
pots3=[]
pots4=[]

for i in range(539,549):
    check1=True
    check2=True
    check3=True
    check4=True
    temps.append(i)
    plot=False
    if i == 544:
        plot=True
    curr=potfinder(i, plot)
    minimum=min(curr['LIQUID']['x'])
    for j in range(len(curr['LIQUID']['x'])):
        if curr['LIQUID']['x'][j]==minimum and check1:
            if i==544:
                Gibb1=curr['LIQUID']['alpha'][j]
            pots1.append(curr['LIQUID']['alpha'][j])
            check1=False
        if abs(curr['LIQUID']['x'][j]-.25)<.01 and check2:
            if i==544:
                Gibb2=curr['LIQUID']['alpha'][j]
            pots2.append(curr['LIQUID']['alpha'][j])
            check2=False
        if abs(curr['LIQUID']['x'][j]-.5)<.01 and check3:
            if i==544:
                Gibb3=curr['LIQUID']['alpha'][j]
            pots3.append(curr['LIQUID']['alpha'][j])
            check3=False
        if abs(curr['LIQUID']['x'][j]-.75)<.01 and check4:
            if i==544:
                Gibb4=curr['LIQUID']['alpha'][j]
            pots4.append(curr['LIQUID']['alpha'][j])
            check4=False
            
potestimate1=np.polyfit(temps,pots1,4)
potderivative1=background.find_derivative(potestimate1)
entropy1=-potderivative1(544)
enthalpy1=Gibb1+entropy1*544
print("Enthalpy at 0="+str(enthalpy1))
print("Entropy at 0="+str(entropy1))

potestimate2=np.polyfit(temps,pots2,4)
potderivative2=background.find_derivative(potestimate2)
entropy2=-potderivative2(544)
enthalpy2=Gibb1+entropy2*544
print("Enthalpy at .25="+str(enthalpy2))
print("Entropy at .25="+str(entropy2))

potestimate3=np.polyfit(temps,pots3,4)
potderivative3=background.find_derivative(potestimate3)
entropy3=-potderivative3(544)
enthalpy3=Gibb1+entropy3*544
print("Enthalpy at .5="+str(enthalpy3))
print("Entropy at .5="+str(entropy3))

potestimate4=np.polyfit(temps,pots4,4)
potderivative4=background.find_derivative(potestimate4)
entropy4=-potderivative4(544)
enthalpy4=Gibb1+entropy4*544
print("Enthalpy at .75="+str(enthalpy4))
print("Entropy at .75="+str(entropy4))

plt.show()


