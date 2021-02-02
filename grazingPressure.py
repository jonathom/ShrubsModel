from pcraster import *
from pcraster.framework import *
import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import transpose

class ShrubManage(DynamicModel):
    def __init__(self):
        DynamicModel.__init__(self)
        # need to fix single shrub patches on bottom margin in initial map
        setclone('initialStateA.map')
    def initial(self):
        self.b1 = 6.8
        self.b2 = 0.387
        self.b3 = 38.8
        self.s1 = 0.109
        self.s2 = 0.061
        self.bg = 0.5
        self.c = 0.0015
        self.ds = 0.028
        self.dg = 0.125
        self.theta = 0.8

        #initializing parameters for management practices
        self.year = 0
        print(cell)
        self.h = cell[0]    # grazing pressure (0 to 1)
        self.n = cell[1]    # removal event period (in years)
        self.f = cell[2]    # fraction of shrub area removed

        #2=shrubs, 1=grass, 0=empty
        self.biotop = self.readmap('initialStateA')
    def dynamic(self):

        #adding smallest positive float number to shrubDensity to avoid division by zero error in log() later on
        #and minimize deviance in result
        prevShrubTotal = cellvalue(maptotal(scalar(self.biotop==2)),0)[0]
        shrubDensity.append([prevShrubTotal/(200*200),self.currentTimeStep()])
        self.report(self.biotop, "biotop")
        random = uniform(1)
        self.grass = self.biotop==1
        weg = empty2grass(self)
        wge = grass2empty(self)
        #grassGrowth contains all the newly grown grass
        self.grassGrowth = ifthenelse(weg>random,self.biotop==0,self.grass)
        self.grassGrowth = ifthenelse(self.grassGrowth==True, self.grassGrowth ^ self.grass, False)
        #grassDeath true when a grass cell has died
        self.grassDeath = ifthenelse(wge>random,self.grass,False)
        #this biotop calculation is only for grass
        self.biotop = ifthenelse(self.grassGrowth,1,self.biotop)
        self.biotop = ifthenelse(self.grassDeath,0,self.biotop)
        self.shrub = self.biotop == 2
        wes = empty2shrub(self)
        wgs = grass2shrub(self)
        wse = shrub2empty(self)
        #shrubGrowth contains newly grown shrubs (from empty and from grass)
        shrubGrowthFromEmpty = ifthenelse(wes > random, self.biotop == 0, self.shrub)
        shrubGrowthFromGrass = ifthenelse(wgs > random, self.biotop == 1, self.shrub)
        self.shrubGrowth = pcror(shrubGrowthFromEmpty, shrubGrowthFromGrass)
        self.shrubGrowth = ifthenelse(self.shrubGrowth == True, self.shrubGrowth ^ self.shrub, False)
        #shrubDeath true when shrub cell dies
        self.shrubDeath = ifthenelse(wse > random, self.shrub, False)
        #update biotop for shrub
        self.biotop = ifthenelse(self.shrubGrowth, 2, self.biotop)
        self.biotop = ifthenelse(self.shrubDeath, 0, self.biotop)

        # mechanical removal event every nth year with fraction f being removed
        # TODO: Removed shrub cells selected randomly within patch, not one clump at the edge. How to do this?
        self.year = self.year + 1
        realization = uniform(1) < self.f
        if (self.year == self.n):
            self.shrubRemoved = pcrand(realization, self.shrub)
            self.biotop = ifthenelse(self.shrubRemoved, 0, self.biotop)
            self.year = 0


def empty2grass(self):
    #used von Neumann neighbourhood for all functions (as defined on p.162, 2.2)
    qge = window4total(scalar(self.grass))/4
    #200x200 is the size of the test biotop
    pg = maptotal(scalar(self.grass))/(200*200)
    weg = self.bg * qge + self.theta * pg
    return weg
def grass2empty(self):
    return self.dg
def empty2shrub(self):
    qse = window4total(scalar(self.shrub))/4
    wes = ((1 - qse) * self.b1 + self.s1) * qse
    return wes
def grass2shrub(self):
    # including grazing pressure self.h (last sentence in 2.7)
    qsg = window4total(scalar(self.shrub))/4
    wgs = ((1 - qsg) * self.b2 * (1 - self.h) + self.s2) * qsg
    return wgs
def shrub2empty(self):
    qss = window4total(scalar(self.shrub))/4
    wse = self.ds + self.c * qss
    return wse

# in here the shrub density is recorded for every timestep in the model run
shrubDensity = []

# this array sets up the run only for different grazing pressure, without any mechanical removal.
# [<grazing pressure (0 to 1)>, <removal event period (in years)>, <fraction of shrub area removed>]
grazingArray=np.array([
    [[0.0, 0, 0.0],[0.1, 0, 0.0],[0.2, 0, 0.0],[0.3, 0, 0.0],[0.4, 0, 0.0],[0.5, 0, 0.0],
     [0.6, 0, 0.0],[0.7, 0, 0.0],[0.8, 0, 0.0],[0.9, 0, 0.0],[1.0, 0, 0.0]]])
# initializing grazing array and fraction
grazing = []
grazingPressure = 0.0


# number of timesteps for every run
nrOfTimeSteps=20
# model setup
myModel = ShrubManage()

dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)


# looping through grazingArray
for idx, row in enumerate(grazingArray, start=0):
    for idy, cell in enumerate(row, start=0):
        print(cell)
        dynamicModel.run()
        slope, intercept = np.polyfit(np.log([item[1] for item in shrubDensity]),
                                      np.log([item[0] for item in shrubDensity]), 1)
        grazingArray[idx][idy] = slope
        grazing.append([grazingArray[idx][idy], grazingPressure])
        grazingPressure = grazingPressure + 0.1


#plot grazing array without mechanical removal, grazing pressure fraction (x) and slope of growth rate (y)
plt.plot([item[1] for item in grazing], [item[0] for item in grazing])
plt.xlabel("grazing pressure (fraction)", size=10)
plt.ylabel("slope of shrub growth", size=10)
plt.show()