import pole2d as pl2
import grid_body_2d as gb2d
import numpy as np
import scipy.optimize as opt
def calcV(Grid:gb2d.Grid,x,y,theta):#刚体位移广义坐标
    V=0.0
    for anchor in Grid.anchors:
        K=anchor.matrix.mat_K
        rotate_mat=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        old=np.array([[anchor.x],[anchor.y]],dtype=np.float64)
        new=np.dot(rotate_mat,old)+np.array([[x],[y]],dtype=np.float64)
        delta=new-old
        V=V+np.dot(delta.T,np.dot(K,delta))
        return V