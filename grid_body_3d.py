## grid_body_3d.py
import numpy as np
import pole3d as pl3
class Load:
    def __init__(self,m,x=0,y=0,z=0,theta=0,phi=0,loadtype='F'):
        self.x=x
        self.y=y
        self.z=z
        self.m=m
        self.theta=theta
        self.phi=phi
        self.loadtype=loadtype
        if loadtype=='F':
            self.Fx=self.m*np.cos(self.phi)*np.cos(self.theta)
            self.Fy=self.m*np.cos(self.phi)*np.sin(self.theta)
            self.Fx=self.m*np.sin(self.phi)
            self.Mx,self.My,self.Mz=0,0,0
        elif loadtype=='M':
            self.Fx,self.Fy,self.Fz=0,0,0
            self.Mx=self.m*np.cos(self.phi)*np.cos(self.theta)
            self.My=self.m*np.cos(self.phi)*np.sin(self.theta)
            self.Mz=self.m*np.sin(self.phi)
        else:
            raise Exception('type accept F or M')
        return
# 定义二力杆与刚体的约束
class Anchor:
    
    # Brief: 反映二力杆的物理参数以及二力杆对于刚体的影响的矩阵列表
    # K : 二维弹力矩阵 F=K*x
    # G :反映二力杆和刚体的几何联系的矩阵
    class Matrix:
        mat_K=np.array([[0,0],[0,0]],dtype=np.float64)
        mat_G=np.zeros((3,2),dtype=np.float64)
        def __init__(self) -> None:
            pass
    matrix=Matrix()
    def __init__(self,x:float,y:float,z:float,poles:list):
        self.x=x
        self.y=y
        self.z=z
        self.poles=poles
        self.matrix.mat_K=pl3.GenMatrix(self.poles)

class Grid:
    anchors=list()
    loadlist=list()
    center=[0,0,0]
    m=0
    J=np.zeros((1,3))
    mode='auto'
    coef_matrix=list()
    def __init__(self,anchors:list):
        self.anchors=anchors
        self.__calcCenter()
        self.__clacCoefMatrix()
    def __init__(self):
        pass
    def __calcCenter(self):#仅在mode='auto'有效
        if self.mode=='auto':
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
            if type(anchor) is not Anchor:
                raise Exception('anchor : type error')
            anchor.matrix.mat_G=\
                np.array([[1,0,0],                                                #Fx
                        [0,1,0],                                                  #Fy
                        [0,0,1],                                                  #Fz
                        [0,-anchor.z+self.center[2],anchor.y-self.center[1]],     #Mcx
                        [-anchor.z+self.center[2],0,anchor.x-self.center[0]],     #Mcy
                        [-anchor.y+self.center[1],anchor.x-self.center[0],0]],    #Mcz
                        dtype=np.float64)
        self.coef_matrix.append(np.dot(anchor.matrix.mat_G,np.dot(anchor.matrix.mat_K,anchor.matrix.mat_G.T)))
        self.coef=np.sum(self.coef_matrix,axis=0)
    def add_anchor(self,anchor:Anchor):
        self.anchors.append(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
    def remove_anchor(self,anchor:Anchor):
        self.anchors.remove(anchor)
        self.__calcCenter()
        self.__clacCoefMatrix()
        
    def add_load(self,load:Load):
        self.loadlist.append(load)   
    def remove_load(self,load:Load):
        self.loadlist.remove(load) 
        
######        
    def __setCenter(self,x:float,y:float,z:float):
        self.center[0]=x
        self.center[1]=y
        self.center[2]=z
        self.mode='manual'
        self.__clacCoefMatrix()
    def __setJ(self,Jyz:float,Jzx:float,Jxy:float):
        self.J=np.array([Jyz,Jzx,Jxy],dtype=np.float64)
        self.__clacCoefMatrix()
    def __setM(self,M:float):
        self.m=M
    def set(self,M:float,J:list,center:list):
        self.__setM(M)
        self.__setJ(J[0],J[1],J[2])
        self.__setCenter(center[0],center[1])
        self.mode='manual'
    def calc(self,arr:np.array):
        return np.dot(self.coef,arr)
    def calc_loadMat(self): #计算负载向量,[Fx,Fy,Mz]
        loadMat=np.zeros((6,1),dtype=np.float64)
        for load in self.loadlist:
            if type(load) is not Load:
                raise Exception('grid_body_3d:load : type error')
            if load.loadtype=='F':
                loadMat[0]=loadMat[0]+load.Fx
                loadMat[1]=loadMat[1]+load.Fy
                loadMat[2]=loadMat[2]+load.Fz
                loadMat[3]=loadMat[3]+load.Fz*(load.y-self.center[1])-load.Fx*(load.z-self.center[2])
                loadMat[4]=loadMat[4]+load.Fx*(load.z-self.center[2])-load.Fz*(load.x-self.center[0])
                loadMat[5]=loadMat[5]+load.Fy*(load.x-self.center[0])-load.Fx*(load.y-self.center[1])
            elif load.loadtype=='M':
                loadMat[0],loadMat[1],loadMat[2]=0,0,0
                loadMat[3]=loadMat[3]+load.Mx
                loadMat[4]=loadMat[4]+load.My
                loadMat[5]=loadMat[5]+load.Mz
        return loadMat
     
    def print_anchors(self):
        for i in self.anchors:
            print("x={0},y={1},poles:{2}".format(i.x,i.y,i.poles))
    def print_coef(self):
        print(self.coef)