import grid_body_2d_nonlinear as gb2d_nl
import scipy.optimize as opt
from main import *
def calcV1(R):
    [x,y,theta]=R
    return gb2d_nl.calcV(grid1,x,y,theta)
print(calcV1([1,1,-0.5]))
a=opt.fmin(calcV1,[0,4,5])
print(a)