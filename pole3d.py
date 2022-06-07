## pole3d.py
import numpy as np
import scipy
from math import *
#类和函数
class Pole:
    E=0
    A=0
    L=0
    theta=0
    phi=0
    k=0
    enable=True
    def __init__(self,E,A,L,theta,phi,enable=True):
        self.theta=theta
        self.E=E
        self.A=A
        self.L=L
        self.k=E*A/L
        self.enable=enable
        self.phi=phi
def GenMatrix(poles):
    ret=np.zeros((3,3),dtype=np.float64)
    for p in poles:
        if type(p) is not Pole:
            raise Exception("p is not a pole")
        if p.enable==False:
            continue
        [xi,yi,zi]=[cos(p.phi)*cos(p.theta),cos(p.phi)*sin(p.theta),sin(p.phi)]
        ret=ret+np.array([p.k],dtype=np.float64)*\
            np.array([[xi**2,xi*yi,xi*zi],[xi*yi,yi**2,yi*zi],[xi*zi,yi*zi,zi**2]],dtype=np.float64)
    return ret
def calc_forces_mat(poles:list):
    forces_mat=list()
    for p in poles:
        if type(p) is not Pole:
            raise Exception("p is not a pole")
        if p.enable==False:
            continue
        [xi,yi,zi]=[cos(p.phi)*cos(p.theta),cos(p.phi)*sin(p.theta),sin(p.phi)]
        mat=-np.array([p.k],dtype=np.float64)*\
            np.array([[-xi**2,xi*yi,xi*zi],[xi*yi,yi**2,yi*zi],[xi*zi,yi*zi,zi**2]],dtype=np.float64)
        forces_mat.append(mat)
    return forces_mat    
def CalR(F,mat,error=0.01):
    mat_inv=np.linalg.inv(mat)
    cond=np.linalg.cond(mat,np.inf)
    R_error=cond*error/(1-cond*error)
    #print("cond(mat)={}".format(cond))
    ret=np.dot(mat_inv,F)
    return [ret,R_error]