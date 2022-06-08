import NonliearPack as nlp
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from my_grid import *

def calcV1(R):
    [x,y,theta]=R
    return nlp.calcV_2d(grid1,x,y,theta)
print("-----------------------------------------------------")
#非线性方法
print("能量法")
a=opt.fmin(calcV1,[0,0,0])
print(a)                            #结果
#线性方法
print("线性方法")
b=np.dot(np.linalg.inv(grid1.coef),grid1.calc_loadMat())
print(b)#结果
#动态仿真：
#参数设置
grid1.set(1,1,[0.5,0])                #设置质量，矩
next_fram=np.array([0,0,0,0,0,0])   #初始状态
it=2000                             #迭代次数
delta_t=0.01                        #时间步长
frames=np.zeros((it+1,6))
b=1                                 #衰减因子
frames[0]=next_fram
#开始迭代仿真
for i in range(it):
    next_fram=nlp.next_frame_2d(grid1,next_fram,t=delta_t,b=b,grad_step=0.001)
    frames[i+1]=next_fram
#打印仿真结果
print("通过动力学仿真\n(衰减因子b={},m={},J={}),100s后位置如下:\n[x,y,theta,vx,vy,omega]:".format(b,grid1.m,grid1.J))
print(frames[-1])
#print(frames)
#plot
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