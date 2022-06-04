import pole2d as pl2
import grid_body_2d as gb2d
class Grid_nonlear(gb2d.Grid):
    V=0
    def __init__(self,grid:gb2d.Grid):
        self=grid