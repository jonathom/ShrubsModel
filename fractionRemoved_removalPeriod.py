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

# variable gp (grazing pressure) can be changed manually, e.g. for 3 total model runs (gp = 0, gp = 0.5, gp = 1)
gp = 1

nrOfTimeSteps=50
# model setup
myModel = ShrubManage()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)

print("[INFO] running model loop for variable mechanical removal and grazing")
grazing = [gp, gp+0.01, 1]
interval = [1, 10.01, 1]
mechanical = [0.1, 1.01, 0.1]

parameters = np.mgrid[grazing[0]:(grazing[1]):grazing[2], interval[0]:(interval[1]):interval[2], mechanical[0]:(mechanical[1]):mechanical[2]].reshape(3,-1).T
parameterArray=np.reshape(parameters, (10, 10, 3))
# this will store the resulting slope (initially filled with infinite values)
resultArray=np.full((parameterArray.shape[0],parameterArray.shape[1]), np.inf)
# storing 2d array with 0 (negative) and 1 (positive values) for visualization, derived from resultArray
visualizationArray=np.full((parameterArray.shape[0],parameterArray.shape[1]), np.inf)
# number of timesteps for every run
# looping through the parameterArray, to run the model for every parameter configuration
for idx, row in enumerate(parameterArray, start=0):
    for idy, cell in enumerate(row, start=0):
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

print()
print('Slope of shrub growth')
print(resultArray)
print()
print('Shrub expansion: 0 = controlled, 1 = not controlled')
print(visualizationArray)


#plot visualization array, mechanical removal fraction (x) and period (y) for given grazing pressure
plt.imshow(transpose(visualizationArray), interpolation='none', cmap=plt.get_cmap('gray_r'),
            origin='lower', extent=[0.0, 10.0, 0.0, 10.0])
plt.xlabel("Fraction removed (in percent)", size=10)
plt.ylabel("removal period (in years)", size=10)
plt.xticks(np.arange(11), ('0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'))
plt.xlim([0, 10])
plt.ylim([0, 10])
plt.show()
