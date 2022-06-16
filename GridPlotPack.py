import matplotlib.pyplot as plt
import numpy as np
import grid_body_2d as gb2d
def plot_grid_2d(grid:gb2d.Grid,title='grid'):
    plt.figure(title)
    plt.plot(grid.center[0],grid.center[1],'o')
    i,j=0,0
    for anchor in grid.anchors:
        if type(anchor) is not gb2d.Anchor:
            raise Exception('anchor type error')
        plt.plot(anchor.x,anchor.y,'o')
        plt.plot([grid.center[0],anchor.x],[grid.center[1],anchor.y],color='r')
        i=i+1
        j=0
        for pole in anchor.poles:
            j=j+1
            if type(pole) is not gb2d.pl2.Pole:
                raise Exception('pole type error')
            k=0.2
            line_end=k*pole.L*np.array([np.cos(pole.theta),np.sin(pole.theta)])+\
                np.array([anchor.x,anchor.y])
            plt.plot([anchor.x,line_end[0]],\
                [anchor.y,line_end[1]],'-')
            plt.text(line_end[0],line_end[1],str(i)+','+str(j))
    plt.axis('equal')
    plt.show()