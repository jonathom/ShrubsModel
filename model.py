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
        
        #this biotop only contains true and false values (so that I can test grass)
        #true=grass, false=empty
        self.grass = 0.1 > uniform(1)
        
        
    def dynamic(self):
        self.report(self.grass, "grass")
        random = uniform(1)
        weg = empty2grass(self)
        wge = grass2empty(self)
        
        ###### this calculation is only for testing the grass #######
        #grassGrowth contains all the newly grown grass
        self.grassGrowth = ifthenelse(weg>random,True,self.grass)
        self.grassGrowth = self.grassGrowth ^ self.grass
        #grassDeath contains all the grass from the previous time that hasnt died this year
        self.grassDeath = ifthenelse(wge>random,False,self.grass)
        #grass at the next timestep will be everywhere where grass has newly grown or hasnt died
        self.grass = self.grassGrowth | self.grassDeath;
        
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