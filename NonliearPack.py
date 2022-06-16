'''
Brief: 非线性计算工具包,采用广义坐标的能量法对其进行计算,提供更加精确地解。
同时利用能量法特性提供对刚体进行动力学分析的工具。
calcV_2d(或_3d):计算某一广义坐标下的包含二力杆弹性势能和负载(等效)势能的值
gradV_2d(或_3d):计算某一广义坐标下的包含二力杆弹性势能和负载(等效)势能的梯度
next_param_2d(或_3d):根据梯度法计算下一个广义坐标,其中考虑了广义速度
'''
import grid_body_2d as gb2d
import grid_body_3d as gb3d
import numpy as np
from numpy import array, cos,sin
#刚体位移广义坐标
def calcV_2d(Grid:gb2d.Grid,x,y,theta,load_type='static'):
    V=0.0
    if load_type == 'static':
        k=0.5
    elif load_type == 'dynamic':
        k=0.5
    else:
        raise Exception("load_type must be 'static' or 'dynamic'")
    for anchor in Grid.anchors:
        if type(anchor) is not gb2d.Anchor:
            raise Exception("anchor must be a Anchor object")
        K=anchor.matrix.mat_K
        rotate_mat=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
        old=np.array([[anchor.x-Grid.center[0]],[anchor.y-Grid.center[1]]],dtype=np.float64)
        new=np.dot(rotate_mat,old)+np.array([[x],[y]],dtype=np.float64)
        delta=new-old
        V=V+np.dot(delta.T,np.dot(K,delta))*k#弹性势能
    grid_loadMat=Grid.calc_loadMat()
    V_load=-np.dot(grid_loadMat.T,np.array([[x],[y],[theta]]))#等效势能
    V=V+V_load
    return V

#刚体位移广义坐标,3d 静加载
def calcV_3d(Grid:gb3d.Grid,x,y,z,alpha,beta,gamma,load_type='static'):
    V=0
    if load_type == 'static':
        k=0.5
    elif load_type == 'dynamic':
        k=1
    else:
        raise Exception("load_type must be 'static' or 'dynamic'")
    for anchor in Grid.anchors:
        if type(anchor) is not gb3d.Anchor:
            raise Exception("anchor must be a Anchor object") 
        K=anchor.matrix.mat_K#锚点系数矩阵
        '''
        rotate_mat=np.array(\
            [[cos(gamma)*cos(alpha)-sin(alpha)*sin(beta)*sin(gamma),\
                -cos(gamma)*sin(alpha)-sin(gamma)*sin(beta)*cos(alpha),\
                    -sin(gamma)*cos(beta)],\
            [cos(beta)*sin(alpha),\
                cos(beta)*cos(alpha),\
                    -sin(beta)],\
            [sin(gamma)*cos(alpha)+sin(alpha)*sin(beta)*cos(gamma),\
                -sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)*sin(beta),\
                    cos(gamma)*cos(beta)]],dtype=np.float64)
                    '''
        rotate_mat_x=np.array([[1,0,0],[0,cos(alpha),-sin(alpha)],[0,sin(alpha),cos(alpha)]],dtype=np.float64)
        rotate_mat_y=np.array([[cos(beta),0,sin(beta)],[0,1,0],[-sin(beta),0,cos(beta)]],dtype=np.float64)
        rotate_mat_z=np.array([[cos(gamma),-sin(gamma),0],[sin(gamma),cos(gamma),0],[0,0,1]],dtype=np.float64)
        rotate_mat=np.dot(rotate_mat_z,np.dot(rotate_mat_y,rotate_mat_x))
        old=np.array([[anchor.x-Grid.center[0]],[anchor.y-Grid.center[1]]\
            ,[anchor.z-Grid.center[2]]],dtype=np.float64)
        new=np.dot(rotate_mat,old)+np.array([[x],[y],[z]],dtype=np.float64)
        delta=new-old
        V=V+np.dot(delta.T,np.dot(K,delta))*k
    grid_loadMat=Grid.calc_loadMat()
    V_load=-np.dot(grid_loadMat.T,\
            np.array([[x],[y],[z],[alpha],[beta],[gamma]]))
    V=V+V_load
    return V

#梯度,2d
def gradV_2d(Grid:gb2d,x=0,y=0,theta=0,params=None,step=0.001,load_type="static"):
    if params is None:
        param=np.array([x,y,theta])
    else:
        if np.shape(params)[0]!=3:
            raise Exception("params must be a list with 3 elements")
        param=params
    V0=calcV_2d(Grid,param[0],param[1],param[2],load_type)
    gradV=np.zeros((3))
    param1=np.zeros((3))
    for i in range(3):
        param1=param
        param1[i]=param[i]+step
        gradV[i]=((calcV_2d(Grid,param1[0],param1[1],param1[2],load_type)-V0)/step)
    return gradV
#计算梯度
def gradV_3d(Grid:gb3d,x=0,y=0,z=0,alpha=0,beta=0,gamma=0,params=None,step=0.001,load_type='static'):
    if params is None:
        param=np.array([x,y,z,alpha,beta,gamma])
    else:
        if np.shape(params)[0]!=6:
            raise Exception("params must be a list with 6 elements")
        param=params
    V0=calcV_3d(Grid,params[0],params[1],params[2],params[3],params[4],params[5],load_type)
    gradV=np.zeros((6))
    for i in range(6):
        param1=param
        param1[i]=param[i]+step
        gradV[i]=((calcV_3d(Grid,param1[0],param1[1],param1[2],param1[3],param1[4],param1[5],load_type)-V0)/step)
    return gradV    
        
#根据初始状态计算下一个时间点的状态 init_frame=[x,y,z,vx,vy,vz]
'''
brief: 根据动力学方程初始状态计算下一个时间点的状态  2d
init_frame=[x,y,z,vx,vy,vz] :起始帧
t:计算的下一帧到起始帧的时间
b :衰减因子,使得速度为不衰减的exp(-bt)倍,默认b=0
gead_step: 计算梯度时的步长
'''    
def next_frame_2d(Grid:gb2d.Grid,init_frame:np.array,t=0.01,b=0,grad_step=0.001):
    if Grid.mode=='auto':
        raise Exception('you must use "setCenter()" and "setJ"to set the center of the grid')
    if np.shape(init_frame)[0]!=6:
        raise Exception('init_frame must be a list with 6 elements')
    next=np.zeros((6))
    next[:]=init_frame[:]
    F=-gradV_2d(Grid,params=init_frame[0:3],step=grad_step,load_type='static')
    a=[0,0,0]
    a[0]=F[0]/Grid.m
    a[1]=F[1]/Grid.m
    a[2]=F[2]/Grid.J
    for i in range(3):
        next[i+3]=next[i+3]+a[i]*t  #更新速度 v=v+a*t
        next[i+3]=next[i+3]*np.exp(-b*t)
        next[i]=next[i]+(next[i+3]+init_frame[i+3])*t/2    #更新位移s=s0+(v0+v1)/2*t
        #该部分改为next[i]=init_frame[i]+(next[i+3]+init_frame[i+3])*t/2时存在较大差异
    return next
#进行动态计算，3d
'''
brief: 根据动力学方程初始状态计算下一个时间点的状态  2d
init_frame=[x,y,z,vx,vy,vz] :起始帧
t:计算的下一帧到起始帧的时间
b :衰减因子,使得速度为不衰减的exp(-bt)倍,默认b=0
gead_step: 计算梯度时的步长
'''    
def next_frame_3d(Grid:gb3d.Grid,init_frame:np.array,t=0.01,b=0,grad_step=0.001):
    if Grid.mode=='auto':
        raise Exception('you must use "setCenter()" and "setJ"to set the center of the grid')
    if np.shape(init_frame)[0]!=12:
        raise Exception('init_frame must be a list with 12 elements')
    next=np.zeros((12))
    next[:]=init_frame[:]
    F=-gradV_3d(Grid,params=init_frame[0:6],step=grad_step,load_type='dynamic')
    a=[0,0,0,0,0,0]
    a[0]=F[0]/Grid.m
    a[1]=F[1]/Grid.m
    a[2]=F[2]/Grid.m
    a[3]=F[3]/Grid.J[0]
    a[4]=F[4]/Grid.J[1]
    a[5]=F[5]/Grid.J[2]
    for i in range(6):
        next[i+6]=next[i+6]+a[i]*t  #更新速度 v=v+a*t
        next[i+6]=next[i+6]*np.exp(-b*t)
        next[i]=next[i]+(next[i+6]+init_frame[i+6])*t/2     #更新位移s=s0+(v0+v1)/2*t
    return next