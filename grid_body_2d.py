import numpy as np
from pole2d import *
# 定义二力杆与刚体的约束
class load:
    def __init__(self,m=0,x=0,y=0,theta=0,loadtype='F'):
        self.x=x
        self.y=y
        self.m=m
        self.theta=theta
        self.loadtype=loadtype
        self.Fx,self.Fy=0,0
        if loadtype=='F':
            self.Fx=self.m*np.cos(self.theta)
            self.Fy=self.m*np.sin(self.theta)
        elif loadtype=='M':
            self.Fx=0
            self.Fy=0
        else:
            raise Exception('type accept F or M')
        return
class anchor:
    
    # Brief: 反映二力杆的物理参数以及二力杆对于刚体的影响的矩阵列表
    # K : 二维弹力矩阵 F=K*x
    # G :反映二力杆和刚体的几何联系的矩阵
    class Matrix:
        mat_K=np.array([[0,0],[0,0]],dtype=np.float64)
        mat_G=np.zeros((3,2),dtype=np.float64)
        def __init__(self) -> None:
            pass
    matrix=Matrix()
    a=0
    def __init__(self,x:float,y:float,poles:list,mode='xy'):
        self.x=x
        self.y=y
        self.poles=poles
        self.mode=mode
        if(self.mode=='xy'):
            self.z=0
        self.matrix.mat_K=GenMatrix(self.poles)

class Grid:
    anchors=list()
    loadlist=list()#负载
    center=[0,0]
    m=0
    J=0
    coef_matrix=list()
    def __init__(self,anchors:list,mode='auto'):
        self.anchors=anchors
        self.mode=mode
        self.__calcCenter()
        self.__clacCoefMatrix()
        if mode=='simple':
            self.m=0
            self.J=0
        else:
            raise Exception('mode accept simple not supported till now')
    def __init__(self,mode='auto') -> None:
        self.anchors=list()
        self.center=[0,0]
        self.mode=mode
        if mode=='auto' or mode=='manual':
            self.m=0
            self.J=0
        else:
            raise Exception('mode manual not implemented')
        self.coef_matrix=list()
        self.cnt=0
    def __calcCenter(self):
        if self.mode=='auto':
            sum_x=0
            sum_y=0
            cnt=0
            for i in self.anchors:
                sum_x=sum_x+i.x
                sum_y=sum_y+i.y
                cnt=cnt+1
            self.center[0]=sum_x/cnt
            self.center[1]=sum_y/cnt
        else:
            self.center[0],self.center[1]=0,0
    # 计算系数矩阵C,[Fx,Fy,Mc]=C*[x,y,theta]
    def __clacCoefMatrix(self):
        for anchor in self.anchors:
            anchor.matrix.mat_G=np.array([[1,0],[0,1],[-anchor.y+self.center[1],anchor.x-self.center[0]]],dtype=np.float64)
        self.coef_matrix.append(np.dot(anchor.matrix.mat_G,np.dot(anchor.matrix.mat_K,anchor.matrix.mat_G.T)))
        self.coef=np.sum(self.coef_matrix,axis=0)
    def add_anchor(self,anchor:anchor) -> None:
        self.anchors.append(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
    def add_load(self,load:load) -> None:
        self.loadlist.append(load)   
    def remove_load(self,load:load) -> None:
        self.loadlist.remove(load)   
    def remove_anchor(self,anchor:anchor) -> None:
        self.anchors.remove(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
    def setCenter(self,x:float,y:float) -> None:
        self.center[0]=x
        self.center[1]=y
        self.mode='manual'
        self.__clacCoefMatrix()
    def calc(self,arr:np.array)->np.array:
        return np.dot(self.coef,arr)
    def calc_loadMat(self) -> np.array: #计算负载向量,[Fx,Fy,Mz]
        loadMat=np.zeros((3,1),dtype=np.float64)
        for load in self.loadlist:
            if load.loadtype=='F':
                loadMat[0]=loadMat[0]+load.Fx
                loadMat[1]=loadMat[1]+load.Fy
                loadMat[2]=(load.x-self.center[0])*load.Fy-(load.y-self.center[1])*load.Fx
            elif load.loadtype=='M':
                loadMat[0]=loadMat[0]+load.Fx
                loadMat[1]=loadMat[1]+load.Fy
                loadMat[2]=loadMat[2]+load.m
        return loadMat
                
    def print_anchors(self) -> None:
        for i in self.anchors:
            print("x={0},y={1},poles:{2}".format(i.x,i.y,i.poles))
    def print_coef(self) -> None:
        print(self.coef)