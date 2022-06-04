import pole2d as pl2
import grid_body_2d as gb2d
import numpy as np
import scipy.optimize as opt
def calcV(Grid:gb2d.Grid,x,y,theta)->float:#刚体位移广义坐标
    V=0.0
    for anchor in Grid.anchors:
        K=anchor.matrix.mat_K
        rotate_mat=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        old=np.array([[anchor.x],[anchor.y]],dtype=np.float64)
        new=np.dot(rotate_mat,old)+np.array([[x],[y]],dtype=np.float64)
        delta=new-old
        grid_loadMat=Grid.calc_loadMat()
        V_load=-np.dot(grid_loadMat.T,np.array([[x],[y],[theta]]))
        V=V+np.dot(delta.T,np.dot(K,delta))+V_load
        
    return V