import NonliearPack as nlp
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from my_grid_3d import *
import tqdm
def calcV1(R):
    [x,y,z,appha,beta,gamma]=R
    return nlp.calcV_3d(grid1,x,y,z,appha,beta,gamma)
print("-----------------------------------------------------")
#非线性方法
print("能量法")
a=opt.fmin(calcV1,[0,0,0,0,0,0])
print(a)                            #结果
#线性方法
print("线性方法")
b=np.dot(np.linalg.inv(grid1.coef),grid1.calc_loadMat())
print(b)#结果
#动态仿真：
#参数设置
print("仿真")
grid1.set(M=10,J=[10,10,10],center=[0,0,0])                #设置质量，矩
#next_fram=np.array([0,0,0,0,0,0,0,0,0,0,0,0])   #初始状态
next_fram=np.zeros((12))
it=10000                             #迭代次数
delta_t=0.1                         #时间步长
frames=np.zeros((it+1,12))
b=4                                 #衰减因子
frames[0]=next_fram
#开始迭代仿真
bar=tqdm.tqdm(total=it)#进度条
for i in range(it):
    next_fram=nlp.next_frame_3d(grid1,next_fram,t=delta_t,b=b,grad_step=0.001)
    frames[i+1]=next_fram
    bar.update(1)
bar.close()
#打印仿真结果
print("通过动力学仿真\n(衰减因子b={},m={},J={}),100s后位置如下:\n[x,y,theta,vx,vy,omega]:".format(b,grid1.m,grid1.J))
print(frames[-1])
#print(frames)
#plot
title=['x','y','z','alpha','beta','gamma','Vx','Vy','Vz','Ox','Oy','Oz']
y_axis=['mm','mm','mm','rad','rad','rad','mm/s','mm/s','mm/s','rad/s','rad/s','rad/s']
y=np.zeros((it+1,12))
y=frames.T
plt.figure("迭代次数{},阻尼{}".format(it,b))
for i in range(12):
    plt.subplot(4,3,i+1)
    plt.title(title[i])
    t=np.linspace(0,delta_t*it,it+1,endpoint=True)
    plt.plot(t,y[i])
    plt.xlabel('time(s)')
    plt.ylabel(y_axis[i])
plt.show()