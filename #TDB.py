#TDB
import background
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import os
import csv
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np
import alphashape
import mendeleev as me

def tdb(filename,comps,temp,surfaces,phases):
    easy=np.around(np.arange(0,1.05,.075),1)
    def check(comp):
        if comp != 'VA' and '%'  not in comp: 
            return me.element(comp).name 
        else: 
            return 'VA'

    graphcomp=[check(i) for i in comps]
    comps=[i.upper() for i in comps]
    # easy=np.append(easy,[.995])
    # print(easy)

    def potfinder(temper,plot):
        # Load database and choose the phases that will be plotted
        # db_nbre = Database('Mg-Bi.tdb')
        db_nbre = Database(filename)
        my_phases_nbre = phases

        # Get the colors that map phase names to colors in the legend
        legend_handles, color_dict = phase_legend(my_phases_nbre)
        hullval=[]
        fig = plt.figure(num=1,figsize=(9,6))
        ax = fig.gca()

        # Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
        phase_potentials={}
        for phase_name in my_phases_nbre:
            result = calculate(db_nbre, comps, phase_name, P=101325, T=temper, output='GM')
            
            
                #ax.scatter(result.X.sel(component='MG'), result.GM, marker='.', s=5, color=color_dict[phase_name])
            
            print(phase_name)
            temp=result.GM.values
            tempx=result.X.sel(component=comps[1]).values
            
            x = tempx.squeeze()

            y = temp.squeeze()

            if phase_name in surfaces:
                newx=[]
                newy=[]
                for i in easy:
                    lowx=[]
                    lowy=[]
                    
                    for j in range(len(x)):
                        if abs(x[j]-i)<0.005:
                            lowx.append(x[j])
                            lowy.append(y[j])
                    bot=min(lowy)
                    newx.append(i)
                    newy.append(bot)
                x=newx
                y=newy

            #Here the (x,y) pairs of points are but into a more simple and callable form
            #This portion of the code needs to be changed but all it does is get the hull of the gibbs surface
            
            try:
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
                phase_potentials[phase_name]={'alpha':alphapot,'beta':betapot,'x':x,'y':y}
            except TypeError:
                print("We good")
                phase_potentials[phase_name]={'x':x,'y':y}

            
        
            if plot:
                if phase_name not in surfaces:
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
            ax.set_xlabel('Mole Percent of '+graphcomp[1],fontsize=16)
            ax.set_ylabel('Molar Gibbs Energy (J/mol)',fontsize=16)
            ax.set_title(' Molar Gibbs Energy of'+graphcomp[1]+'-'+graphcomp[0]+' System at '+str(temper)+"\xb0"+" K")
            ax.set_xlim((0, 1))
            ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))

            pplot=plt.figure(2)
            ax2=pplot.gca()
            phase='LIQUID'
            ax2.set_xlabel('Mole Percent of '+graphcomp[1],fontsize=20)
            ax2.set_ylabel('Chemical Potential (J/mol)',fontsize=20 )
            ax2.set_title('Chemical Potential of the '+phase+' phase at '+str(temper)+" K",fontsize=25)
            ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['alpha'],color='red',label='Potential of '+graphcomp[1])
            ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['beta'],color='blue',label='Potential of '+graphcomp[0])
            ax2.legend(loc='best', scatterpoints=100)

            
        return phase_potentials
    temps=[]
    vals=[[] for i in range(len(easy))]
    Gibbs=[0 for i in range(len(easy))]
    # potsb=[[] for i in range(40)]
    # Gibbsb=[0 for i in range(40)]


    for i in np.around(np.arange(temp-2,temp+2,.25),3):
        print('-----------------------------------------')
        print(str(i)+' K')
        print('-----------------------------------------')

        check=[True for i in range(len(easy))]
        # checkb=[True for i in range(40)]
        temps.append(i)
        plot=False
        if i == temp:
            plot=True
        curr=potfinder(i, plot)
        # minimum=min(curr['LIQUID']['x'])
        for j in range(len(curr['LIQUID']['x'])):
            for k in range(len(easy)):
                if abs(curr['LIQUID']['x'][j]-easy[k])<.005 and check[k]:
                    if i==temp:
                        print("Gibbbs added: "+str(curr['LIQUID']['y'][j]))
                        Gibbs[k]=curr['LIQUID']['y'][j]
                    vals[k].append(curr['LIQUID']['y'][j])
                    check[k]=False

        # for j in range(len(curr['LIQUID']['x'])):
        #     for k in range(len(easy)):
        #         if abs(curr['LIQUID']['x'][j]-easy[k])<.005 and checkb[k]:
        #             if i==544:
        #                 Gibbsb[k]=curr['LIQUID']['y'][j]
        #             potsb[k].append(curr['LIQUID']['y'][j])
        #             checkb[k]=False
            # if curr['LIQUID']['x'][j]==minimum and check1:
            #     if i==544:
            #         Gibb1=curr['LIQUID']['alpha'][j]
            #     pots1.append(curr['LIQUID']['alpha'][j])
            #     check1=False
            # if abs(curr['LIQUID']['x'][j]-.25)<.01 and check2:
            #     if i==544:
            #         Gibb2=curr['LIQUID']['alpha'][j]
            #     pots2.append(curr['LIQUID']['alpha'][j])
            #     check2=False
            # if abs(curr['LIQUID']['x'][j]-.5)<.01 and check3:
            #     if i==544:
            #         Gibb3=curr['LIQUID']['alpha'][j]
            #     pots3.append(curr['LIQUID']['alpha'][j])
            #     check3=False
            # if abs(curr['LIQUID']['x'][j]-.75)<.01 and check4:
            #     if i==544:
            #         Gibb4=curr['LIQUID']['alpha'][j]
            #     pots4.append(curr['LIQUID']['alpha'][j])
            #     check4=False
        
    potestimate=[]
    potderivative=[]
    entropy=[]
    enthalpy=[]
    for i in range(len(easy)):
        potestimate.append(np.polyfit(temps,vals[i],4))
        potderivative.append(background.find_derivative(potestimate[-1]))
        entropy.append(-potderivative[-1](temp))
        enthalpy.append(Gibbs[i]+entropy[-1]*temp)
        print(str(enthalpy[-1]/entropy[-1]))
    
    mixe=[]
    for i in range(len(easy)):
        mixe.append(entropy[i]-(easy[i]*entropy[-1]+(1-easy[i])*entropy[0]))
        # print(str(entropy[i])+"-"+(str(easy[i])+"*"+str(entropy[-1])+"+"+(str(1-easy[i]))+"*"+str(entropy[0]))+" = "+str(entropy[i]-(easy[i]*entropy[-1]+(1-easy[i])*entropy[0])))

        # mixe.append(enthalpy[i]-(easy[i]*enthalpy[-1]+(1-easy[i])*enthalpy[0]))
        # print(str(enthalpy[i])+"-"+(str(easy[i])+"*"+str(enthalpy[-1])+"+"+(str(1-easy[i]))+"*"+str(enthalpy[0]))+" = "+str(enthalpy[i]-(easy[i]*enthalpy[-1]+(1-easy[i])*enthalpy[0])))
        # mixe.append(Gibbs[i]-(easy[i]*Gibbs[-1]+(1-easy[i])*Gibbs[0]))
    mplot=plt.figure(3)
    ax3=mplot.gca()
    ax3.scatter(easy,mixe)
    ax3.set_xlabel('Mole Percent of '+graphcomp[1],fontsize=16)
    ax3.set_ylabel('Molar Enthalpy of Mixing (J)',fontsize=16)
    ax3.set_title(' Molar Enthalpy of Mixing of '+graphcomp[1]+'-'+graphcomp[0]+' System at '+str(temp)+"\xb0"+" K")
    ax3.set_xlim((0, 1))
    plt.show()


# tdb('Al-Li.tdb',['Al', 'Li','VA'],1250,['BCC_B2'],['LIQUID','BCC_A2','FCC_A1','BCC_DIS','BCC_B2','AL2LI3','AL4LI9'])
# tdb('Cu-Ni.tdb',['Cu', 'Ni','CU%','VA%','NI%'],1800,[],['LIQUID','FCC_A1'])
# tdb('Ca-Na.tdb',['Ca','Na','VA'],1500,[],['LIQUID','FCC_A1','BCC_A2','HCP_A3'])
tdb('Mg-Bi.tdb',['Bi','Mg','VA'],1400,['LIQUID'],['LIQUID', 'BI2MG3_H', 'BI2MG3_L'])
# tdb('Ba-Ga.tdb',['Ba','Ga', 'VA'],1200,[],["LIQUID","BCC","ORTHO","BA10GA","BA8GA7","BAGA2","BAGA4","BA5GA6"])