import numpy as np
import pole2d as pl2
# 定义二力杆与刚体的约束
class Load:
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
class Anchor:
    
    # Brief: 反映二力杆的物理参数以及二力杆对于刚体的影响的矩阵列表
    # K : 二维弹力矩阵 F=K*x
    # G :反映二力杆和刚体的几何联系的矩阵
    class Matrix:
        mat_K=np.array([[0,0],[0,0]],dtype=np.float64)
        mat_G=np.zeros((3,2),dtype=np.float64)
        mat_pole_K=None
        def __init__(self):
            pass
    count=0
    matrix=Matrix()
    info=dict()#为其他算法提供空间
    poles=list()
    def __init__(self,x:float,y:float,poles:list):
        self.x=x
        self.y=y
        for pole in poles:
            if type(pole) is not pl2.Pole:
                raise Exception("pole:type error")
            if pole.enable is True:
                self.poles.append(pole)
        self.matrix.mat_K=pl2.GenMatrix(self.poles)
        self.mat_pole_K=pl2.calc_forces_mat(self.poles)
        self.count=len(self.poles)

class Grid:
    anchors=list()
    loadlist=list()#负载
    center=[0,0]
    m=0
    J=0
    mode='auto'
    coef_matrix=list()
    def __init__(self,anchors:list):
        self.anchors=anchors
        self.__calcCenter()
        self.__clacCoefMatrix()
    def __init__(self):
        pass
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
    # 计算系数矩阵C,[Fx,Fy,Mc]=C*[x,y,theta]
    def __clacCoefMatrix(self):
        self.coef_matrix=list()
        for anchor in self.anchors:
            if type(anchor) is not Anchor:
                raise Exception('anchor : type error')
            anchor.matrix.mat_G=np.array([[1,0],[0,1],[-anchor.y+self.center[1],anchor.x-self.center[0]]],dtype=np.float64)
            self.coef_matrix.append(np.dot(anchor.matrix.mat_G,np.dot(anchor.matrix.mat_K,anchor.matrix.mat_G.T)))
        self.coef=np.sum(self.coef_matrix,axis=0)
    def add_anchor(self,anchor:Anchor):
        self.anchors.append(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
    def add_load(self,load:Load):
        self.loadlist.append(load)   
    def remove_load(self,load:Load):
        self.loadlist.remove(load)   
    def remove_anchor(self,anchor:Anchor):
        self.anchors.remove(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
######
    def __setCenter(self,x:float,y:float):
        self.center[0]=x
        self.center[1]=y
        self.__clacCoefMatrix()
    def __setJ(self,J:float):
        self.J=J
    def __setM(self,M:float):
        self.m=M
    def set(self,M:float,J:float,center:list):
        self.__setM(M)
        self.__setJ(J)
        self.__setCenter(center[0],center[1])
        self.mode='manual'
    def calc(self,arr:np.array):
        return np.dot(self.coef,arr)
    def calc_loadMat(self): #计算负载向量,[Fx,Fy,Mz]
        loadMat=np.zeros((3,1),dtype=np.float64)
        for load in self.loadlist:
            if type(load) is not Load:
                raise Exception('load : type error')
            if load.loadtype=='F':
                loadMat[0]=loadMat[0]+load.Fx
                loadMat[1]=loadMat[1]+load.Fy
                loadMat[2]=(load.x-self.center[0])*load.Fy-(load.y-self.center[1])*load.Fx
            elif load.loadtype=='M':
                loadMat[0]=loadMat[0]+load.Fx
                loadMat[1]=loadMat[1]+load.Fy
                loadMat[2]=loadMat[2]+load.m
        return loadMat
                
    def print_anchors(self):
        for i in self.anchors:
            if type(i) is not Anchor:
                raise Exception('anchor : type error')
            print("x={0},y={1},poles:{2}".format(i.x,i.y,i.poles))
    def print_coef(self):
        print(self.coef)