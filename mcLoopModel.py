from pcraster import *
from pcraster.framework import *
import math
import matplotlib.pyplot as plt
import numpy as np

class ShrubManage(DynamicModel, MonteCarloModel):
    def __init__(self):
        DynamicModel.__init__(self)
        # need to fix single shrub patches on bottom margin in initial map
        setclone('initialStateA.map')
        MonteCarloModel.__init__(self)
        DynamicModel.__init__(self)

    def premcloop(self):
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
        self.h = 0    # grazing pressure (0 to 1)
        self.n = 5    # removal event period (in years)
        self.f = 0.99  # fraction of shrub area removed

        self.initialMap = self.readmap('initialStateA')

        #2=shrubs, 1=grass, 0=empty
        self.biotop = self.initialMap

        self.shrubExpansion = []
        self.shrubDensity = []

    def initial(self):
        self.biotop = self.initialMap

    def dynamic(self):
        prevShrubTotal = cellvalue(maptotal(scalar(self.biotop==2)),0)[0]
        self.shrubDensity.append([prevShrubTotal/(200*200),self.currentTimeStep(),self.currentSampleNumber()])

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

        postShrubTotal = cellvalue(maptotal(scalar(self.biotop==2)),0)[0]
        # this corresponds to chapter 3.0 "logarithm of shrub cover in year t+1 divided by shrub cover in year t"
        intrinsicGrowth = math.log(postShrubTotal/prevShrubTotal)
        self.shrubExpansion.append([intrinsicGrowth,self.currentTimeStep(),self.currentSampleNumber()])

    def postmcloop(self):
        ''' this is to visualize the shrub growth per timestep (like figure 2b)
        for i in range(nrOfSamples):
            plt.plot([item[1] for item in self.shrubExpansion if item[2] == i],[item[0] for item in self.shrubExpansion if item[2] == i])
        '''
        slope, intercept = np.polyfit(np.log([item[1] for item in self.shrubDensity if item[2] == 1]),
                                      np.log([item[0] for item in self.shrubDensity if item[2] == 1]), 1)
        # this is the slope of the  shrub Expansion
        # if this is negative the shrubs are not expanding
        print(slope)
        plt.loglog([item[1] for item in self.shrubDensity if item[2] == 1],
                   [item[0] for item in self.shrubDensity if item[2] == 1], '--')
        plt.show()

        names = ['biotop']
        sampleNumbers = self.sampleNumbers()
        timesteps = self.timeSteps()
        percentiles = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        mcpercentiles(names, percentiles, sampleNumbers, timesteps)
        mcaveragevariance(names, sampleNumbers, timesteps)


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


nrOfSamples = 2
nrOfTimeSteps = 100
myModel = ShrubManage()
dynamicModel = DynamicFramework(myModel, nrOfTimeSteps)
mcModel = MonteCarloFramework(dynamicModel, nrOfSamples)
mcModel.run()
