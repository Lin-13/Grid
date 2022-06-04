import NonliearPack as nlp
import scipy.optimize as opt
from main import *
def calcV1(R):
    [x,y,theta]=R
    return nlp.calcV_2d(grid1,x,y,theta)
print(calcV1([1,1,-0.5]))
a=opt.fmin(calcV1,[0,4,5])
print(a)