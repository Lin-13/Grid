import NonliearPack as nlp
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from my_grid import *
def calcV1(R):
    [x,y,theta]=R
    return nlp.calcV_2d(grid1,x,y,theta)
print(calcV1([1,1,-0.5]))
#非线性方法
a=opt.fmin(calcV1,[0,0,0])

print(a)
#线性方法
b=np.dot(np.linalg.inv(grid1.coef),grid1.calc_loadMat())
print(b)
#动态仿真：
grid1.set(1,1,[5,0])
next_fram=np.array([0,0,0,0,0,0])
it=2000
delta_t=0.01
frames=np.zeros((it+1,6))
frames[0]=next_fram
for i in range(it):
    next_fram=nlp.next_frame_2d(grid1,next_fram,t=delta_t,b=1,grad_step=0.001)
    frames[i+1]=next_fram
print("通过动力学仿真(衰减因子b=0.05,m=1,J=1),100s后位置如下:\n[x,y,theta,vx,vy,omega]:")
print(frames[-1])
#print(frames)
print(np.shape(frames))
title=['x','y','theta','Vx','Vy','omega']
y_axis=['mm','mm','rad','mm/s','mm/s','rad/s']
plt.figure(1)
for i in range(6):
    plt.subplot(2,3,i+1)
    plt.title(title[i])
    y=np.zeros((it+1))
    for j in range(it+1):
        y[j]=frames[j][i]
    t=np.linspace(0,delta_t*it,it+1,endpoint=True)
    plt.plot(t,y)
    plt.xlabel('time(s)')
    plt.ylabel(y_axis[i])
plt.show()