from pcraster import *
from pcraster.framework import *
import math
import matplotlib.pyplot as plt
import numpy as np

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
        self.h = cell[0]          # grazing pressure (0 to 1)
        self.n = cell[1]    # removal event period (in years)
        self.f = cell[2]    # fraction of shrub area removed

        #2=shrubs, 1=grass, 0=empty
        self.biotop = self.readmap('initialStateA')


    def dynamic(self):
        prevShrubTotal = cellvalue(maptotal(scalar(self.biotop==2)),0)[0]
        shrubDensity.append([prevShrubTotal/(200*200),self.currentTimeStep()])

        #mechanical removal event every nth year with fraction f being removed
        #TODO: Removed shrub cells selected randomly within patch, not one clump at the edge. How to do this?
        self.year = self.year + 1
        realization = uniform(1) < self.f
        if (self.year == self.n):
            self.shrubRemoved = pcrand(realization, self.shrub)
            self.biotop = ifthenelse(self.shrubRemoved,0,self.biotop)
            self.year = 0


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
# variable gp (grazing pressure) can be changed manually, e.g. for 3 total model runs (gp = 0, gp = 0.5, gp = 1)
gp = 0
# this is a 2D array that stores a parameter configuration for every run
# [<grazing pressure (0 to 1)>, <removal event period (in years)>, <fraction of shrub area removed>]
parameterArray=np.array([
    [[gp, 2, 0.1],[gp, 6, 0.1]],
    [[gp, 2, 0.8],[gp, 6, 0.8]]])
# this will store the resulting slope (initially filled with infinite values)
resultArray=np.full((parameterArray.shape[0],parameterArray.shape[1]), np.inf)
# storing 2d array with 0 (negative) and 1 (positive values) for visualization, derived from resultArray
visualizationArray=np.full((parameterArray.shape[0],parameterArray.shape[1]), np.inf)
# number of timesteps for every run
nrOfTimeSteps=100
# model setup
myModel = ShrubManage()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
# looping through the parameterArray, to run the model for every parameter configuration
for idx, row in enumerate(parameterArray, start=0):
    for idy, cell in enumerate(row, start=0):
        print(cell)
        dynamicModel.run()
        #if shrubDensity is above 0.0, slope calculation with log(), if it becomes 0.0 at some point,
        #resultArray is assigned -1 to avoid division by zero error
        if all(item[0] for item in shrubDensity):
            slope, intercept = np.polyfit(np.log([item[1] for item in shrubDensity]),
                               np.log([item[0] for item in shrubDensity]), 1)
            resultArray[idx][idy] = slope
        else:
            resultArray[idx][idy] = -1
        ''' to visualize the shrub density slope for this run
        plt.loglog([item[1] for item in shrubDensity],
           [item[0] for item in shrubDensity], '--')
        plt.show()
        '''
        shrubDensity = []
        #filling visualizationArray with 0 or 1
        if resultArray[idx][idy] > 0:
            visualizationArray[idx][idy] = 1
        else:
            visualizationArray[idx][idy] = 0

print()
print(resultArray)
print()
print(visualizationArray)

#to be modified...
plt.imshow(visualizationArray, interpolation='none', cmap=plt.get_cmap('gray'))
plt.xlabel("Fraction removed (in percent)", size=10)
plt.ylabel("removal period (in years)", size=10)
plt.xlim([1, 2])
plt.ylim([1, 2])
plt.show()
