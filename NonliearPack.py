import grid_body_2d as gb2d
import grid_body_3d as gb3d
import numpy as np
from numpy import cos,sin
def calcV_2d(Grid:gb2d.Grid,x,y,theta)->float:#刚体位移广义坐标
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
def calcV_3d(Grid:gb3d.Grid,x,y,z,alpha,beta,gamma)->float:
    V=0
    for anchor in Grid.anchors:
        K=anchor.matrix.mat_K
        rotate_mat=np.array(\
            [[cos(gamma)*cos(alpha)-sin(alpha)*sin(beta)*sin(gamma),\
                -cos(beta)*sin(alpha)-sin(gamma)*sin(beta),\
                    -sin(gamma)*cos(beta)],\
            [cos(beta)*sin(alpha),\
                cos(beta)*cos(alpha),\
                    -sin(beta)],\
            [sin(gamma)*cos(alpha)+sin(alpha)*sin(beta)*cos(gamma),\
                -sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)*sin(beta),\
                    cos(gamma)*cos(beta)]],dtype=np.float64)
        old=np.array([[anchor.x],[anchor.y],[anchor.z]],dtype=np.float64)
        new=np.dot(rotate_mat,old)+np.array([[x],[y],[z]],dtype=np.float64)
        delta=old-new
        grid_loadMat=Grid.calc_loadMat()
        V_load=-np.dot(grid_loadMat.T,\
            np.array([[x],[y],[z],[alpha],[beta],[gamma]]))
        V=V+np.dot(delta.T,np.dot(K,delta))+V_load
    return V