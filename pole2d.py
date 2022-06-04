from typing import overload
import numpy as np
class pole:
    E=0
    A=0
    L=0
    theta=0
    k=0
    enable=True
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