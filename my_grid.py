import numpy as np
import grid_body_2d as gb2d
import pole2d as pl2
grid1=gb2d.Grid()
pole1=pl2.Pole(theta=0)#默认 E=1,A=1,L=1,enable=True
pole2=pl2.Pole(theta=np.pi/2)
pole3=pl2.Pole(theta=np.pi)
pole4=pl2.Pole(theta=-np.pi/2)
grid1.add_anchor(gb2d.Anchor(0,0,[pole1,pole2,pole3,pole4]))
grid1.add_anchor(gb2d.Anchor(1,0,[pole1,pole2,pole3,pole4]))
grid1.print_anchors()
grid1.print_coef()
F1=gb2d.Load(0.01,0,0,np.pi/2,loadtype='F')
M1=gb2d.Load(-0.01,loadtype='M')
grid1.add_load(F1)
grid1.add_load(M1)
mat=grid1.calc_loadMat()#计算加载矩阵，[Fx,Fy,M]
print(mat)