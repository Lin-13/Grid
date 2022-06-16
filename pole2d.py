from typing import overload
import numpy as np
class Pole:
    def __init__(self,theta,E=1,A=1,L=1,enable=True):
        self.theta=theta
        self.E=E
        self.A=A
        self.L=L
        self.k=E*A/L
        self.enable=enable
    def __repr__(self) -> str:
        return 'pole(theta=%f,E=%f,A=%f,L=%f,enable=%s)'%(self.theta,self.E,self.A,self.L,self.enable)
    
def GenMatrix(poles:list):
    sum1=0.0
    sum2=0.0
    sum3=0.0
    for p in poles:
        if p.enable==False:
            continue
        sum1=sum1+p.k*np.cos(p.theta)**2
        sum2=sum2+p.k*np.cos(p.theta)*np.sin(p.theta)
        sum3=sum3+p.k*np.sin(p.theta)**2
    ret=np.array([[sum1,sum2],[sum2,sum3]],dtype=np.float64)
    return ret
#计算各个杆力的矩阵，并返回一个列表,第k个杆的力Fk=force_mat[k]*[x;y]
def calc_forces_mat(poles:list):
    forces_mat=list()
    for p in poles:
        if type(p) is not Pole:
            raise Exception("p is not a pole")        
        if p.enable==False:
            continue
        a11=p.k*np.cos(p.theta)**2
        a12=p.k*np.sin(p.theta)*np.sin(p.theta)
        a22=p.k*np.cos(p.theta)**2
        mat=-np.array([[a11,a12],[a12,a22]],dtype=np.float64)
        forces_mat.append(mat)
    return forces_mat
#计算二力杆组在受到外力时的移动
def CalR(F,mat,error=0.01):
    mat_inv=np.linalg.inv(mat)
    cond=np.linalg.cond(mat,np.inf)
    R_error=cond*error/(1-cond*error)
    #print("cond(mat)={}".format(cond))
    ret=np.dot(mat_inv,F)
    return [ret,R_error]