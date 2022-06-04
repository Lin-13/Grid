import NonliearPack as nlp
import scipy.optimize as opt
from my_grid import *
def calcV1(R):
    [x,y,theta]=R
    return nlp.calcV_2d(grid1,x,y,theta)
print(calcV1([1,1,-0.5]))
#非线性方法
a=opt.fmin(calcV1,[0,4,5])

print(a)
#线性方法
b=np.dot(np.linalg.inv(grid1.coef),grid1.calc_loadMat())
print(b)