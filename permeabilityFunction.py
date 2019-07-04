#!/usr/bin/env python
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##################################################################################################
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
# generate mesh.npy and load it to calculate permeability

def permeability(mesh):

    size = np.shape(mesh)
    # variable delaration
    nx = size[1]  # 800
    ny = size[0]  # 302

    nt = 1000
    unit =  1E-6   #1um
    xmin = 0
    xmax = size[1]*unit

    ymin = 0.0
    ymax = size[0]*unit

    dx = ((xmax-xmin)/(nx))  # 1*unit
    dy = ((ymax-ymin)/(ny))  # 1*unit

    const = 1  # dp/dz/mu

    # initializing all matrices
    u = np.zeros((ny, nx))
    un = np.zeros((ny, nx))
    # b = np.ones((nx,ny))
    b = np.ones((ny, nx))

    # the scheme is copied from step 9
    for n in range(nt):
        un = u.copy()
        u[1:-1,1:-1] = np.where(mesh[1:-1,1:-1]==0,(dx**2*(un[1:-1,0:-2]+un[1:-1,2:])+ dy**2*(un[0:-2,1:-1]+un[2:,1:-1])-b[1:-1,1:-1]*dx**2*dy**2)/(2*dx**2+2*dy**2),0)

        u[0,:]=u[-1,:]=0
    # # set boundary condition
        for i in range(0, ny):  #(0,302)
            j=0
            if mesh[i][j] == 0:
                u[i][j] = (dx**2*(un[i,(j+1)%nx]+un[i,(j-1)%nx]) + dy**2*(un[i+1,j] + un[i-1,j])-b[i,j]*dx**2*dy**2)/(2*dx**2+2*dy**2)
            k=-1
            if mesh[i][k] == 0:
                u[i][k] = (dx**2*(un[i,(k+1)%nx]+un[i,(k-1)%nx]) + dy**2*(un[i+1,k] + un[i-1,k])-b[i,k]*dx**2*dy**2)/(2*dx**2+2*dy**2)
    #calculate permeability
    Q_flowrate = u.sum()*dx*dy
    solid = mesh.sum()*dx*dy
    flux =  Q_flowrate/(xmax*ymax-solid)
    fpermeability = flux / (-const) #m^2
    #change the unit of permeability
    fpermeability = fpermeability*1013249965828.1448 #darcy
    return fpermeability

#run in this file 
if __name__ == "__main__":
    mesh = np.load('mesh.npy')	#can be run in linux; the filedirectory is default in linux(home/)
    permeability = permeability(mesh)
    print("the permeability of this grain packing is " + str(permeability))
    plt.figure(1)
    plt.imshow(mesh)

    # plt.figure(4)
    # # plt.tripcolor(x,y,u[:],cmap=cm.jet)
    # # plt.imshow(u,cmap=cm.jet)
    # plt.colorbar()
    plt.show()
