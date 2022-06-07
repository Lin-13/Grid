import scipy.stats as st
import numpy as np
import grid_body_2d as gb2d
import grid_body_3d as gb3d
def calc_pole_force_2d(grid:gb2d.Grid,x:float,y:float,theta:float):
    r=np.array([[x],[y],[theta]])
    for anchor in grid.anchors:
        if type(anchor) is not gb2d.Anchor:
            raise Exception('grid.anchors must be a list of gb2d.Anchor')
        r_k=