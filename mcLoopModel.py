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

# variable gp (grazing pressure) can be changed manually, e.g. for 3 total model runs (gp = 0, gp = 0.5, gp = 1)
gp = 1
# this is a 2D array that stores a parameter configuration for every run
# [<grazing pressure (0 to 1)>, <removal event period (in years)>, <fraction of shrub area removed>]
parameterArray=np.array([
    [[gp, 1, 0.0],[gp, 1, 0.1],[gp, 1, 0.2],[gp, 1, 0.3],[gp, 1, 0.4],[gp, 1, 0.5],[gp, 1, 0.6],[gp, 1, 0.7],[gp, 1, 0.8],[gp, 1, 0.9],[gp, 1, 1.0]],
    [[gp, 2, 0.0],[gp, 2, 0.1],[gp, 2, 0.2],[gp, 2, 0.3],[gp, 2, 0.4],[gp, 2, 0.5],[gp, 2, 0.6],[gp, 2, 0.7],[gp, 2, 0.8],[gp, 2, 0.9],[gp, 2, 1.0]],
    [[gp, 3, 0.0],[gp, 3, 0.1],[gp, 3, 0.2],[gp, 3, 0.3],[gp, 3, 0.4],[gp, 3, 0.5],[gp, 3, 0.6],[gp, 3, 0.7],[gp, 3, 0.8],[gp, 3, 0.9],[gp, 3, 1.0]],
    [[gp, 4, 0.0],[gp, 4, 0.1],[gp, 4, 0.2],[gp, 4, 0.3],[gp, 4, 0.4],[gp, 4, 0.5],[gp, 4, 0.6],[gp, 4, 0.7],[gp, 4, 0.8],[gp, 4, 0.9],[gp, 4, 1.0]],
    [[gp, 5, 0.0],[gp, 5, 0.1],[gp, 5, 0.2],[gp, 5, 0.3],[gp, 5, 0.4],[gp, 5, 0.5],[gp, 5, 0.6],[gp, 5, 0.7],[gp, 5, 0.8],[gp, 5, 0.9],[gp, 5, 1.0]],
    [[gp, 6, 0.0],[gp, 6, 0.1],[gp, 6, 0.2],[gp, 6, 0.3],[gp, 6, 0.4],[gp, 6, 0.5],[gp, 6, 0.6],[gp, 6, 0.7],[gp, 6, 0.8],[gp, 6, 0.9],[gp, 6, 1.0]],
    [[gp, 7, 0.0],[gp, 7, 0.1],[gp, 7, 0.2],[gp, 7, 0.3],[gp, 7, 0.4],[gp, 7, 0.5],[gp, 7, 0.6],[gp, 7, 0.7],[gp, 7, 0.8],[gp, 7, 0.9],[gp, 7, 1.0]],
    [[gp, 8, 0.0],[gp, 8, 0.1],[gp, 8, 0.2],[gp, 8, 0.3],[gp, 8, 0.4],[gp, 8, 0.5],[gp, 8, 0.6],[gp, 8, 0.7],[gp, 8, 0.8],[gp, 8, 0.9],[gp, 8, 1.0]],
    [[gp, 9, 0.0],[gp, 9, 0.1],[gp, 9, 0.2],[gp, 9, 0.3],[gp, 9, 0.4],[gp, 9, 0.5],[gp, 9, 0.6],[gp, 9, 0.7],[gp, 9, 0.8],[gp, 9, 0.9],[gp, 9, 1.0]],
    [[gp, 10,0.0],[gp,10, 0.1],[gp,10, 0.2],[gp,10, 0.3],[gp,10, 0.4],[gp,10, 0.5],[gp,10, 0.6],[gp,10, 0.7],[gp,10, 0.8],[gp,10, 0.9],[gp,10, 1.0]]])

# this array sets up the run only for different grazing pressure, without any mechanical removal.
# [<grazing pressure (0 to 1)>, <removal event period (in years)>, <fraction of shrub area removed>]
grazingArray=np.array([
    [[0.0, 0, 0.0],[0.1, 0, 0.0],[0.2, 0, 0.0],[0.3, 0, 0.0],[0.4, 0, 0.0],[0.5, 0, 0.0],
     [0.6, 0, 0.0],[0.7, 0, 0.0],[0.8, 0, 0.0],[0.9, 0, 0.0],[1.0, 0, 0.0]]])
# initializing grazing array and fraction
grazing = []
grazingPressure = 0.0

# this will store the resulting slope (initially filled with infinite values)
resultArray=np.full((parameterArray.shape[0],parameterArray.shape[1]), np.inf)
# storing 2d array with 0 for negative and 1 for positive values (for visualization) derived from resultArray
visualizationArray=np.full((parameterArray.shape[0],parameterArray.shape[1]), np.inf)


# number of timesteps for every run
nrOfTimeSteps=50
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

        shrubDensity = []
        #filling visualizationArray with zero or one
        if resultArray[idx][idy] > 0:
            visualizationArray[idx][idy] = 1
        else:
            visualizationArray[idx][idy] = 0

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

print()
print('Slope of shrub growth')
print(resultArray)
print()
print('Shrub expansion: 0 = controlled, 1 = not controlled')
print(visualizationArray)

#plot grazing array without mechanical removal, grazing pressure fraction (x) and slope of growth rate (y)
plt.plot([item[1] for item in grazing], [item[0] for item in grazing])
plt.xlabel("grazing pressure (fraction)", size=10)
plt.ylabel("slope of shrub growth", size=10)
plt.show()

#plot visualization array, mechanical removal fraction (x) and period (y) for given grazing pressure
plt.imshow(transpose(visualizationArray), interpolation='none', cmap=plt.get_cmap('gray_r'),
            origin='lower', extent=[0.0, 10.0, 0.0, 10.0])
plt.xlabel("Fraction removed (in percent)", size=10)
plt.ylabel("removal period (in years)", size=10)
plt.xticks(np.arange(11), ('0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'))
plt.xlim([0, 10])
plt.ylim([0, 10])
plt.show()
