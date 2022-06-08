from cmath import pi
import numpy as np
from scipy import stats as st # 导入统计函数计算概率

import grid_body_3d as gb3d
import pole3d as pl3
grid1=gb3d.Grid()
pole1=pl3.Pole(E=1,A=1,L=1,theta=np.pi/2,phi=pi/3,enable=True,)#默认 E=1,A=1,L=1,enable=True
pole2=pl3.Pole(E=1,A=1,L=1,theta=np.pi/2,phi=0)
pole3=pl3.Pole(E=1,A=1,L=1,theta=np.pi,phi=pi/3)
pole4=pl3.Pole(E=1,A=1,L=1,theta=-np.pi/2,phi=0)
pole5=pl3.Pole(E=1,A=1,L=1,theta=0,phi=pi/2)
pole6=pl3.Pole(E=1,A=1,L=1,theta=0,phi=0)
grid1.add_anchor(gb3d.Anchor(0,0,0,[pole1,pole2,pole3,pole4,pole5,pole6]))
grid1.add_anchor(gb3d.Anchor(1,1,1,[pole1,pole2,pole3,pole4,pole5,pole6]))
grid1.add_anchor(gb3d.Anchor(1,0,1,[pole1,pole2,pole3,pole4,pole5,pole6]))
grid1.add_anchor(gb3d.Anchor(0,0,1,[pole1,pole2,pole3,pole4,pole5,pole6]))

#grid1.print_anchors()
#grid1.print_coef()
F1=gb3d.Load(1,0,0,0,0,pi/2,loadtype='F')
M1=gb3d.Load(0.000001,x=1,phi=0,loadtype='M')
grid1.add_load(F1)
grid1.add_load(M1)
mat=grid1.calc_loadMat()  #计算加载矩阵，[Fx,Fy,M]
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