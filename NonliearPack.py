import grid_body_2d as gb2d
import grid_body_3d as gb3d
import numpy as np
from numpy import cos,sin
#刚体位移广义坐标,静加载
def calcV_2d(Grid:gb2d.Grid,x,y,theta):
    V=0.0
    for anchor in Grid.anchors:
        K=anchor.matrix.mat_K
        rotate_mat=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        old=np.array([[anchor.x],[anchor.y]],dtype=np.float64)
        new=np.dot(rotate_mat,old)+np.array([[x],[y]],dtype=np.float64)
        delta=new-old
        grid_loadMat=Grid.calc_loadMat()
        V_load=-np.dot(grid_loadMat.T,np.array([[x],[y],[theta]]))/2
        V=V+np.dot(delta.T,np.dot(K,delta))/2+V_load
    return V
#刚体位移广义坐标,3d
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
            np.array([[x],[y],[z],[alpha],[beta],[gamma]]))/2
        V=V+np.dot(delta.T,np.dot(K,delta))/2+V_load
    return V
#梯度,2d
def gradV_2d(Grid:gb2d,x=0,y=0,theta=0,params=None,step=0.001):
    if params is None:
        param=np.array([x,y,theta])
    else:
        if len(params)!=3:
            raise Exception("params must be a list with 3 elements")
        param=params
    V0=calcV_2d(Grid,x,y,theta)
    gradV=np.zeros((1,3))
    for i in len(param):
        param1=param
        param1[i]=param+step
        gradV[i]=((calcV_2d(Grid,param1[0],param1[1],param1[2])-V0)/step)
    return gradV
#计算梯度
def gradV_3d(Grid:gb2d,x=0,y=0,z=0,alpha=0,beta=0,gamma=0,params=None,step=0.001):
    if params is None:
        param=np.array([x,y,z,alpha,beta,gamma])
    else:
        if len(params)!=6:
            raise Exception("params must be a list with 6 elements")
        param=params
    V0=calcV_3d(Grid,x,y,z,alpha,beta,gamma)
    gradV=np.zeros((1,6))
    for i in len(param):
        param1=param
        param1[i]=param+step
        gradV[i]=((calcV_2d(Grid,param1[0],param1[1],param1[2])-V0)/step)
    return gradV    
        
#根据初始状态计算下一个时间点的状态 init_frame=[x,y,z,vx,vy,vz]    
def next_frame_2d(Grid:gb2d.Grid,init_frame:np.array,t=0.01,step=0.001):
    if Grid.mode=='auto':
        raise Exception('you must use "setCenter()" and "setJ"to set the center of the grid')
    if len(init_frame)!=6:
        raise Exception('init_frame must be a list with 6 elements')
    next=init_frame
    F=gradV_2d(Grid,param=init_frame,step=step)
    a=[0,0,0]
    a[0]=F[0]/Grid.m
    a[1]=F[1]/Grid.m
    a[2]=F[2]/Grid.J
    for i in range(3):
        next[i+3]=next[i+3]+a[i]*t  #更新速度 v=v+a*t
        next[i]=next[1]+(next[i+3]+init_frame[i+3])*t/2     #更新位移s=s0+(v0+v1)/2*t
    return next
#进行动态计算，3d
def next_frame_3d(Grid:gb3d.Grid,init_frame:np.array,t=0.01,step=0.001):
    if Grid.mode=='auto':
        raise Exception('you must use "setCenter()" and "setJ"to set the center of the grid')
    if len(init_frame)!=12:
        raise Exception('init_frame must be a list with 12 elements')
    next=init_frame
    F=gradV_2d(Grid,param=init_frame[0,5],step=step)
    a=[0,0,0,0,0,0]
    a[0]=F[0]/Grid.m
    a[1]=F[1]/Grid.m
    a[2]=F[2]/Grid.m
    a[3]=F[3]/Grid.J[0]
    a[4]=F[4]/Grid.J[1]
    a[5]=F[5]/Grid.J[2]
    for i in range(3):
        next[i+3]=next[i+3]+a[i]*t  #更新速度 v=v+a*t
        next[i]=next[1]+(next[i+3]+init_frame[i+3])*t/2     #更新位移s=s0+(v0+v1)/2*t
    return next