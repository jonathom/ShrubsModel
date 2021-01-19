from pcraster import *
from pcraster.framework import *

class ShrubManage(DynamicModel, MonteCarloModel):
    def __init__(self):
        # this is just a test grid (the dimensions are not right)
        setclone('initialStateA.map')
        MonteCarloModel.__init__(self)
        DynamicModel.__init__(self)
    
    def premcloop(self):
        #2=shrubs, 1=grass, 0=empty
        self.biotop = self.readmap('initialStateA')
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
        
    def initial(self):
        pass
        
    def dynamic(self):
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
        
    def postmcloop(self):
        names = ['biotop']
        sampleNumbers = self.sampleNumbers()
        timesteps = self.timeSteps()
        percentiles = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        mcpercentiles(names, percentiles, sampleNumbers, timesteps)
        mcaveragevariance(names, sampleNumbers, timesteps)
        
        
def empty2grass(self):
    #TODO: discuss what kind of neighborhood we will use
    qge = window4total(scalar(self.grass))/4
    #30x40 is the size of the test biotop
    pg = maptotal(scalar(self.grass))/(30*40)
    weg = self.bg * qge + self.theta * pg
    return weg

def grass2empty(self):
    return self.dg

def empty2shrub(self):
    #used von Neumann neighbourhood type (as defined on p.162, 2.2)
    qse = window4total(scalar(self.shrub))/4
    wes = ((1 - qse) * self.b1 + self.s1) * qse
    return wes

def grass2shrub(self):
    qsg = window4total(scalar(self.shrub))/4
    wgs = ((1 - qsg) * self.b2 + self.s2) * qsg
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
