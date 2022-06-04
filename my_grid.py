
import numpy as np
from scipy import stats as st # 导入统计函数计算概率

import grid_body_2d as gb2d
import pole2d as pl2
grid1=gb2d.Grid()
pole1=pl2.pole(theta=0)#默认 E=1,A=1,L=1,enable=True
pole2=pl2.pole(E=1,theta=np.pi/2)
pole3=pl2.pole(theta=np.pi)
pole4=pl2.pole(theta=-np.pi/2,E=1)
grid1.add_anchor(gb2d.anchor(0,0,[pole1,pole2,pole3,pole4]))
grid1.add_anchor(gb2d.anchor(1,0,[pole1,pole2,pole3,pole4]))
grid1.print_anchors()
grid1.print_coef()
F1=gb2d.load(1,0,0,0,loadtype='F')
M1=gb2d.load(1,loadtype='M')
grid1.add_load(F1)
#grid1.add_load(M1)
mat=grid1.calc_loadMat()#计算加载矩阵，[Fx,Fy,M]
np.dot(np.linalg.inv(grid1.coef),mat)
'''
print(grid1.anchors[0].matrix.mat_K)
print(grid1.anchors[0].matrix.mat_G)
print(grid1.anchors[1].matrix.mat_K)
print(grid1.anchors[1].matrix.mat_G)
K=grid1.anchors[0].matrix.mat_K
G=grid1.anchors[0].matrix.mat_G
print(np.dot(G,np.dot(K,G.T)))
'''
#print(grid1.coef_matrix)
import matplotlib.pyplot as plt
st.norm.cdf(0.5)
import scipy.optimize as opt
#opt.minimize()