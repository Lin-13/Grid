import scipy.stats as st
import numpy as np
import grid_body_2d as gb2d
import grid_body_3d as gb3d
#广义坐标下anchor中pole的受力
def calc_pole_force_2d(grid:gb2d.Grid,x:float,y:float,theta:float):
    r=np.array([[x],[y],[theta]])#广义坐标向量
    for anchor in grid.anchors:
        if type(anchor) is not gb2d.Anchor:
            raise Exception('grid.anchors must be a list of gb2d.Anchor')
        index=0
        F_k=np.zeros((anchor.count,2))#2d，两个方向
        sigma_N=np.zeros((anchor.count))#正应力
        for pole in anchor.poles:
            if type(pole)is not gb2d.pl2.Pole:
                raise Exception('pole: type error')
            F_k[index]=np.dot(anchor.mat_pole_K[index],\
                np.dot(anchor.matrix.mat_G.T,r))
            sigma_N[index]=np.sqrt(np.dot(F_k[index],F_k[index]))/pole.A
            index=index+1
        anchor.info['F_k']=F_k
        anchor.info['sigma_N']=sigma_N
def calc_pole_force_3d(grid:gb3d.Grid,x:float,y:float,\
    z:float,alpha:float,beta:float,gamma:float):
    r=np.array([[x],[y],[z],[alpha],[beta],[gamma]])#广义坐标向量
    for anchor in grid.anchors:
        if type(anchor) is not gb2d.Anchor:
            raise Exception('grid.anchors must be a list of gb2d.Anchor')
        index=0
        F_k=np.zeros((anchor.count,3))#3d，3个方向
        sigma_N=np.zeros((anchor.count))#正应力
        for pole in anchor.poles:
            if type(pole)is not gb3d.pl3.Pole:
                raise Exception('pole 3d: type error')
            F_k[index]=np.dot(anchor.mat_pole_K[index],\
                np.dot(anchor.matrix.mat_G.T,r))
            sigma_N[index]=np.sqrt(np.dot(F_k[index],F_k[index]))/pole.A
            index=index+1
        anchor.info['F_k']=F_k
        anchor.info['sigma_N']=sigma_N