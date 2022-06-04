## grid_body_3d.py
import numpy as np
import pole3d as pl3
# 定义二力杆与刚体的约束
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
    def __init__(self,x:float,y:float,z:float,poles:list,mode='xy'):
        self.x=x
        self.y=y
        self.z=z
        self.poles=poles
        self.mode=mode
        if(self.mode=='xy'):
            self.z=0
        self.matrix.mat_K=pl3.GenMatrix(self.poles)

class Grid:
    anchors=list()
    center=[0,0]
    m=0
    J=0
    mode='auto'
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
    def __init__(self):
        self.anchors=list()
        self.center=[0,0,0]
        self.m=0
        self.J=0
        self.coef_matrix=list()
        self.cnt=0
    def __calcCenter(self):
        sum_x=0
        sum_y=0
        sum_z=0
        cnt=0
        for i in self.anchors:
            sum_x=sum_x+i.x
            sum_y=sum_y+i.y
            sum_z=sum_z+i.z
            cnt=cnt+1
        self.center[0]=sum_x/cnt
        self.center[1]=sum_y/cnt
        self.center[2]=sum_z/cnt
    # 计算系数矩阵C:[Fx,Fy,Fz,Mcx,Mcy,Mcz]=C*[x,y,z,alpha,beta,gamma]
    def __clacCoefMatrix(self):
        for anchor in self.anchors:
            anchor.matrix.mat_G=\
                np.array([[1,0,0],                                                  #Fx
                        [0,1,0],                                                  #Fy
                        [0,0,1],                                                  #Fz
                        [0,-anchor.z+self.center[2],anchor.y-self.center[1]],     #Mcx
                        [-anchor.z+self.center[2],0,anchor.x-self.center[0]],     #Mcy
                        [-anchor.y+self.center[1],anchor.x-self.center[0],0]],    #Mcz
                        dtype=np.float64)
        self.coef_matrix.append(np.dot(anchor.matrix.mat_G,np.dot(anchor.matrix.mat_K,anchor.matrix.mat_G.T)))
        self.coef=np.sum(self.coef_matrix,axis=0)
    def add_anchor(self,anchor:anchor) -> None:
        self.anchors.append(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
    def remove_anchor(self,anchor:anchor) -> None:
        self.anchors.remove(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
    def calc(self,arr:np.array)->np.array:
        return np.dot(self.coef,arr)
    def print_anchors(self) -> None:
        for i in self.anchors:
            print("x={0},y={1},poles:{2}".format(i.x,i.y,i.poles))
    def print_coef(self) -> None:
        print(self.coef)