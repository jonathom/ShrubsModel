from pcraster import *
from pcraster.framework import *

class ShrubManage(DynamicModel):
    def __init__(self):
        DynamicModel.__init__(self)
        # this is just a test grid (the dimensions are not right)
        setclone('test.map')
        
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
        
        #this is just for testing until we have the actual initial maps
        #2=shrubs, 1=grass, 0=empty
        random = uniform(1)
        self.biotop = ifthenelse(0.1 > random,scalar(2),ifthenelse(0.2 > random, scalar(1),scalar(0)))
        
        
    def dynamic(self):
        self.report(self.biotop, "biotop")
        random = uniform(1)
        self.grass = self.biotop==1
        weg = empty2grass(self)
        wge = grass2empty(self)
        
        #grassGrowth contains all the newly grown grass
        self.grassGrowth = ifthenelse(weg>random,self.biotop==0,self.grass)
        self.grassGrowth = ifthenelse(self.grassGrowth==True, self.grassGrowth ^ self.grass, False)
        #grass can only grow in cells that have grass as a neighbor
        hasGrassNeighbor = window4total(scalar(self.grass))
        self.grassGrowth = ifthenelse(hasGrassNeighbor>0,self.grassGrowth,False)
        #grassDeath true when a grass cell has died
        self.grassDeath = ifthenelse(wge>random,self.grass,False)
        #this biotop calculation is only for grass
        self.biotop = ifthenelse(self.grassGrowth,1,self.biotop)
        self.biotop = ifthenelse(self.grassDeath,0,self.biotop)
        
def empty2grass(self):
    #TODO: discuss what kind of neighborhood we will use
    qge = window4total(scalar(self.grass))/4
    #30x40 is the size of the test biotop
    pg = maptotal(scalar(self.grass))/(30*40)
    weg = self.bg * qge + self.theta * pg
    return weg

def grass2empty(self):
    return self.dg

nrOfTimeSteps=100
myModel = ShrubManage()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()