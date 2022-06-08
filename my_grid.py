
import numpy as np
from scipy import stats as st # 导入统计函数计算概率

import grid_body_2d as gb2d
import pole2d as pl2
grid1=gb2d.Grid()
pole1=pl2.Pole(theta=0)#默认 E=1,A=1,L=1,enable=True
pole2=pl2.Pole(E=1,theta=np.pi/2)
pole3=pl2.Pole(theta=np.pi)
pole4=pl2.Pole(theta=-np.pi/2,E=1)
grid1.add_anchor(gb2d.Anchor(0,0,[pole1,pole2,pole3,pole4]))
grid1.add_anchor(gb2d.Anchor(1,0,[pole1,pole2,pole3,pole4]))
grid1.print_anchors()
grid1.print_coef()
F1=gb2d.Load(0.01,0,0,0,loadtype='F')
M1=gb2d.Load(0.01,loadtype='M')
grid1.add_load(F1)
grid1.add_load(M1)
mat=grid1.calc_loadMat()#计算加载矩阵，[Fx,Fy,M]
'''
print(mat)
print("anchor0:{},{}\nanchor1:{},{}"\
    .format(grid1.anchors[0].x,grid1.anchors[0].y\
        ,grid1.anchors[1].x,grid1.anchors[1].y))
print("center:{}".format(grid1.center))

print(grid1.anchors[0].matrix.mat_K)
print(grid1.anchors[0].matrix.mat_G)
print(grid1.anchors[1].matrix.mat_K)
print(grid1.anchors[1].matrix.mat_G)
'''
#K=grid1.anchors[0].matrix.mat_K
#G=grid1.anchors[0].matrix.mat_G
#print(np.dot(G,np.dot(K,G.T)))

#print(grid1.coef)
import matplotlib.pyplot as plt
st.norm.cdf(0.5)
import scipy.optimize as opt
#opt.minimize()