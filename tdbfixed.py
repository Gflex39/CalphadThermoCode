#TDB
import background
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
from pycalphad import Database, calculate, variables as v
from pycalphad.plot.utils import phase_legend
import numpy as np
import os
import matplotlib as mpl
import time
# cur=os.path.dirname(__file__)
# print(cur)
# print("---------------------------------")
import mendeleev as me
class tdb:
    def __init__(self,filename,comps,temp,surfaces,phases,full=True,catalog=False,interest='LIQUID'):
        # if cur[0]:
            # print(cur[1])
            # os.chdir(cur[1])
        self.interest=interest
        self.inverse={}
        self.catalog=catalog
        self.hullval=[]
        self.filename=filename
        self.comps=comps
        self.comps1=comps
        self.temp=temp
        self.surfaces=surfaces
        self.phases=phases
        self.easy=np.around(np.arange(0,1.05,.05),2)
        def checker(comp):
            if comp != 'VA' and '%'  not in comp: 
                return me.element(comp).name 
            else: 
                return 'VA'

        self.graphcomp=[checker(i) for i in comps]
        self.comps=[i.upper() for i in comps]



        self.temps=[]
        self.vals=[[] for i in range(len(self.easy))]
        self.Gibbs=[0 for i in range(len(self.easy))]
      


        for i in np.around(np.arange(temp-2,temp+2,1),3):
            print('-----------------------------------------')
            print(str(i)+' K')
            print('-----------------------------------------')

            self.check=[True for i in range(len(self.easy))]
            
            self.temps.append(i)
            plot=False
            if i == temp and full:
                plot=True

            self.curr=self.potfinder(i, plot)
            
            if i ==temp:
                self.act=self.curr
            bottom=min(self.curr[self.interest]['x'])
            top=max(self.curr[self.interest]['x'])
            for j in range(len(self.curr[self.interest]['x'])):
                if self.curr[self.interest]['x'][j]==bottom:
                    bottomy=self.curr[self.interest]['y'][j]
                elif self.curr[self.interest]['x'][j]==top:
                    topy=self.curr[self.interest]['y'][j]

            
            for j in range(len(self.curr[self.interest]['x'])):
                for k in range(len(self.easy)):
                    # print(abs(self.curr['LIQUID']['x'][j]-self.easy[k]))
                    # print(self.easy[k])
                    # print(self.check[k])

                    if abs(self.curr[self.interest]['x'][j]-self.easy[k])<.005 and self.check[k]:
                        if i==temp:
                            # print("Gibbbs added: "+str(self.curr['LIQUID']['y'][j]))
                            self.Gibbs[k]=self.curr[self.interest]['y'][j]
                        self.vals[k].append(self.curr[self.interest]['y'][j])
                        self.check[k]=False
            which=True
            for j in range(len(self.Gibbs)):
                if self.easy[j]<.5:
                    which=False
                if self.Gibbs[j]==0:
                    if i==temp:
                        # print("Gibbbs added: "+str(self.curr['LIQUID']['y'][j]))
                        self.Gibbs[j]=bottomy
                        if which:  
                            self.vals[j]=[bottomy for i in np.around(np.arange(temp-2,temp+2,1),3)]
                        else:
                            self.vals[j]=[topy for i in np.around(np.arange(temp-2,temp+2,1),3)]

                    
                    


                        

                
            


        if full:   
            self.potestimate=[]
            self.potderivative=[]
            self.entropy=[]
            self.enthalpy=[]
            self.Igibb=[background.ideal_gibbs(i, self.temp) for i in self.easy]
            # for i in range(len(self.Igibb)):
            #     print([self.easy[i],self.Igibb[i]])
            self.Ientropy=[background.ideal_entropy(i) for i in self.easy]
            self.Egibb=[]
            self.Eentropy=[]
            for i in range(len(self.easy)):
                self.potestimate.append(np.polyfit(self.temps,self.vals[i],4))
                self.potderivative.append(background.find_derivative(self.potestimate[-1]))
                self.entropy.append(-self.potderivative[-1](temp))
                self.enthalpy.append(self.Gibbs[i]+self.entropy[-1]*temp)

                    

        
        #This commented portion takes the values attained and turns them into a csv file labeling what it is
        # background.csvform("Mg", "Bi", phase_name, temper)
        #This finds the convex hull of the plot
        self.hullval=np.array(self.hullval)
        self.hull = ConvexHull(self.hullval)
        

        # plt.plot(hullval[:][0], hullval[:][1], 'o')
        #for simplex in hull.simplices:

        #plt.plot(hullval[simplex,0], hullval[simplex,1], 'k-')
        #plt.plot(hullval[hull.vertices,0], hullval[hull.vertices,1], 'r--', lw=2)
        #plt.plot(hullval[hull.vertices[0],0], hullval[hull.vertices[0],1], 'ro')
        
    def potfinder(self,temper,plot):

        
        # Load database and choose the phases that will be plotted
        # db_nbre = Database('Mg-Bi.tdb')
        # print(os.getcwd())
        db_nbre = Database(self.filename)
        my_phases_nbre = self.phases

        # Get the colors that map phase names to colors in the legend
        self.legend_handles, self.color_dict = phase_legend(my_phases_nbre)
        fig = plt.figure(num=1,figsize=(9,6))
        ax = fig.gca()

        # Loop over phases, calculate the Gibbs energy, and scatter plot GM vs. X(RE)
        phase_potentials={}
        for phase_name in my_phases_nbre:
            # if True:
            try:
                # assert False
                # print(phase_name)
                
                phase_potentials[phase_name]=background.checkdata(self.graphcomp[1], self.graphcomp[0], phase_name, temper)

                x=phase_potentials[phase_name]['x']
                y=phase_potentials[phase_name]['y']                
                for i in range(len(phase_potentials[phase_name]['x'])):

                    app=[x[i],y[i]]
                    if temper==self.temp:
                        self.hullval.append(app)
                    

                
            except:
                result = calculate(db_nbre, self.comps, phase_name, P=101325, T=temper, output='GM')
                print(phase_name)
                temp=result.GM.values
                tempx=result.X.sel(component=self.comps[1]).values
                
                x = tempx.squeeze()

                y = temp.squeeze()

                if phase_name in self.surfaces:
                    newx=[]
                    newy=[]
                    for i in self.easy:
                        lowx=[]
                        lowy=[]
                        
                        for j in range(len(x)):
                            if abs(x[j]-i)<0.01:
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
                        if temper==self.temp:
                            self.hullval.append(app)
                    #From this point all the code makes a polynomial estimate of the specific phase curve 
                    #then find the derivative
                    estimate=np.polyfit(x,y,15)
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
                    empty=[None for i in range(len(x))]
                    print("We good")
                    phase_potentials[phase_name]={'alpha':empty,'beta':empty,'x':x,'y':y}
                if self.catalog:
                    background.csvform(self.graphcomp[1], self.graphcomp[0], phase_name, temper,phase_potentials[phase_name],['alpha','beta','x','y'])

                
            
            if plot:
                if phase_name not in self.surfaces:
                    ax.scatter(x,y, marker='.', s=5, color=self.color_dict[phase_name])
                else:
                    plt.figure(1)
                    plt.plot(x,y,color=self.color_dict[phase_name])
                ax.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
                ax.set_ylabel('Molar Gibbs Energy (J/mol)',fontsize=16)
                ax.set_title(' Molar Gibbs Energy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(temper)+"\xb0"+" K")
                ax.set_xlim((0, 1))
                ax.legend(handles=self.legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))
   


                # Format the plot  
                if phase_name=='LIQUID':

                    fig2=plt.figure(2)
                    ax2=fig2.gca()
                    phase='LIQUID'
                    ax2.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=20)
                    ax2.set_ylabel('Chemical Potential (J/mol)',fontsize=20 )
                    ax2.set_title('Chemical Potential of the '+phase+' phase at '+str(temper)+" K",fontsize=25)
                    ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['alpha'],color='red',label='Potential of '+self.graphcomp[1])
                    ax2.scatter(phase_potentials[phase]['x'],phase_potentials[phase]['beta'],color='blue',label='Potential of '+self.graphcomp[0])

                    ax2.legend(loc='best', scatterpoints=100)

            
        return phase_potentials
    def graphdata(self,name,catalog):
        """Can graph data by using keys
        Entropy,GibbsM,Enthalpy,EntropyM,EnthalpyM,EntropyE,GibbE,EntropyI,GibbI """
        entropymix=[]
        enthalpymix=[]
        gibbmix=[]
        for i in range(len(self.easy)):
            entropymix.append(self.entropy[i]-(self.easy[i]*self.entropy[-1]+(1-self.easy[i])*self.entropy[0]))
            # print(str(entropy[i])+"-"+(str(easy[i])+"*"+str(entropy[-1])+"+"+(str(1-easy[i]))+"*"+str(entropy[0]))+" = "+str(entropy[i]-(easy[i]*entropy[-1]+(1-easy[i])*entropy[0])))

            enthalpymix.append(self.enthalpy[i]-(self.easy[i]*self.enthalpy[-1]+(1-self.easy[i])*self.enthalpy[0]))
            # print(str(enthalpy[i])+"-"+(str(easy[i])+"*"+str(enthalpy[-1])+"+"+(str(1-easy[i]))+"*"+str(enthalpy[0]))+" = "+str(enthalpy[i]-(easy[i]*enthalpy[-1]+(1-easy[i])*enthalpy[0])))
            gibbmix.append(self.Gibbs[i]-(self.easy[i]*self.Gibbs[-1]+(1-self.easy[i])*self.Gibbs[0]))
            self.Egibb.append(gibbmix[i]-self.Igibb[i])
            self.Eentropy.append(entropymix[i]-self.Ientropy[i])
        
        fig3=plt.figure(3)
        ax3=fig3.gca()
        if name=='Entropy':
            plt.plot(self.easy,self.entropy)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Molar Entropy (J/K)',fontsize=16)
            ax3.set_title(' Molar Entropy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='GibbsM':
            plt.plot(self.easy,gibbmix)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Molar Gibbs Energy of Mixing (J)',fontsize=16)
            ax3.set_title(' Molar Gibbs Energy of Mixing of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='Enthalpy':
            plt.plot(self.easy,self.enthalpy)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Molar Enthalpy (J)',fontsize=16)
            ax3.set_title(' Molar Enthalpy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='EntropyM':
            plt.plot(self.easy,entropymix)
            if catalog:
                background.csvform(self.graphcomp[1], self.graphcomp[0]+" Entropy", self.interest, self.temp,{'x':self.easy,'y':entropymix,'temp':[self.temp for i in self.easy]},['x','y','temp'])
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Molar Entropy of Mixing (J/K)',fontsize=16)
            ax3.set_title(' Molar Entropy of Mixing of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='EnthalpyM':
            plt.plot(self.easy,enthalpymix)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Molar Enthalpy of Mixing (J)',fontsize=16)
            ax3.set_title(' Molar Enthalpy of Mixing of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='EntropyE':
            plt.plot(self.easy,self.Eentropy)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Excess Molar Entropy (J/K)',fontsize=16)
            ax3.set_title('Excess Molar Entropy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='GibbE':
            plt.plot(self.easy,self.Egibb)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Excess Molar Gibbs Energy (J)',fontsize=16)
            ax3.set_title('Excess Molar Gibbs Energy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='EntropyI':
            plt.plot(self.easy,self.Ientropy)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Ideal Molar Entropy (J/K)',fontsize=16)
            ax3.set_title('Ideal Molar Entropy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))
        elif name=='GibbI':
            plt.plot(self.easy,self.Igibb)
            ax3.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
            ax3.set_ylabel('Ideal Molar Gibbs Energy (J)',fontsize=16)
            ax3.set_title('Ideal Molar Gibbs Energy of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
            ax3.set_xlim((0, 1))

        else:
            print("Dont have that data")
    def composition(self, comp):
        # print(comp)
        mode=[]
        modey=[]
        phase=[]
        for simplex in self.hull.simplices:
            sides=(self.hullval[simplex,0])
            
            if comp > min(sides) and comp < max(sides):
                if min(sides)> .0000000000001 and max(sides)!=1:
                    mode=sides
                    modey=self.hullval[simplex,1]

                    
            
            # print((self.hullval[simplex,0], self.hullval[simplex,1]))
            
        # plt.plot(self.hullval[simplex,0], self.hullval[simplex,1], 'k-')
        
        for i in range(2):
            breaker=False
            for j in self.act.keys():
                if breaker:
                    break
                for k in range(len(self.act[j]['x'])):
                    # print(mode[i]-self.act[j]['x'][k])
                    # print(modey[i]-self.act[j]['y'][k])
                    if mode[i]==self.act[j]['x'][k] and modey[i]==self.act[j]['y'][k]:
                        phase.append(j)
                        breaker==True
                        
        if mode[1]<mode[0]:
            mode=[mode[1],mode[0]]
            modey=[modey[1],modey[0]]
            # print(phase)
            phase=[phase[1],phase[0]]
        if phase[0]==phase[1]:
            # print(phase[0])
            return phase[0]
        else:
            range_=mode[1]-mode[0]
            diff=comp-mode[0]
            # print(str(np.round(diff/range_*100,2))+"% "+phase[0]+" and "+str(np.round((1-diff/range_)*100,2))+"% "+phase[1])
            return None
    def phase_diagram(self):
        points=[]
        layers=[]
        for i in np.arange(self.temp-300,self.temp+300,1):
            layers.append(tdb(self.filename,self.comps1,i,self.surfaces,self.phases,full=False,catalog=True))
        for i in range(len(layers)):
            if i%100==0:
                print(i)
            # print("#####################################################################################")
            for j in self.easy[1:-1]:
                # print(self.temp-30+i)
                # print(j)
                # print("*******************************************************************************")
                points.append([j,self.temp-150+i,layers[i].composition(j)])
        fig4=plt.figure(4)
        ax4=fig4.gca()
        for i in points:
            if i[2] in self.color_dict.keys():
                ax4.scatter(i[0],i[1],color=self.color_dict[i[2]],s=30)
                # ax4.add_patch(mpl.patches.Circle((i[0],i[1]), radius=0.5,color=self.color_dict[i[2]]))
            else:
                # ax4.add_patch(mpl.patches.Circle((i[0],i[1]), radius=0.5,color='hotpink'))
                ax4.scatter(i[0],i[1],color='hotpink',s=30)
        ax4.legend(handles=self.legend_handles, loc='center left', bbox_to_anchor=(1, 0.6))
        ax4.set_xlabel('Mole Percent of '+self.graphcomp[1],fontsize=16)
        ax4.set_ylabel('Temperature (K)',fontsize=16)
        ax4.set_title(' Phase Diagram of '+self.graphcomp[1]+'-'+self.graphcomp[0]+' System at '+str(self.temp)+"\xb0"+" K")
        ax4.set_xlim((0, 1))

        

            





        
        




        

            


# AlLi=tdb('Al-Li.tdb',['Al', 'Li','VA'],1250,['BCC_B2'],['LIQUID','BCC_A2','FCC_A1','BCC_DIS','BCC_B2','AL2LI3','AL4LI9'])

#CaNa=tdb('Ca-Na.tdb',['Ca','Na','VA'],1500,[],['LIQUID','FCC_A1','BCC_A2','HCP_A3'])
# MgBi=tdb('Mg-Bi.tdb',['Bi','Mg','VA'],673,['LIQUID'],['LIQUID', 'BI2MG3_H', 'BI2MG3_L'])

# BaGa=tdb('Ba-Ga.tdb',['Ba','Ga', 'VA'],1200,[],["LIQUID","BCC","ORTHO","BA10GA","BA8GA7","BAGA2","BAGA4","BA5GA6"])
# CuNi.composition(.4)
# CuNi.phase_diagram()
# if True:
for i in range(1):
    FeNi=tdb('Fe-Ni.tdb',['Ni', 'Fe','VA'],500+10*i,['BCC_A2','BCC4','FCC4'],['LIQUID','BCC_A2','FCC_A1','FCC4','BCC4'],catalog=True,interest="FCC4")
    FeNi.graphdata("EntropyM",True)

# CuNi=tdb('Cu-Ni.tdb',['Cu', 'Ni','CU%','VA%','NI%'],1150+10*i,[],['LIQUID','FCC_A1'],catalog=True,interest='FCC_A1')
# CuNi.graphdata('EntropyM',True)
# FeCo=tdb('Fe-Ni.tdb',['Ni', 'Fe','VA'],950,[],['LIQUID','BCC_A2','FCC_A1'])
# FeCr=tdb('Fe-Cr.tdb',['Cr', 'Fe','VA'],950,[],['LIQUID','BCC_A2','FCC_A1','SIGMA'])
clock=time.time()
# CuNi.phase_diagram()
print(time.time()-clock)
plt.show()