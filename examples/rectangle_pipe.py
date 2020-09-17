from pylab import *

from sys import path as sys_path
sys_path.append('./../')
from Poiseuille.Poiseuille import rectangle

##############################
#
# PARAMETERS
#
##############################

D = 0.3
W = 1

a = W/2
b = D

##############################
#
# MESH
#
##############################

y = linspace( -a, a, 30 )
z = linspace( -b, 0, int( len( y )*b/a ) )
y, z = meshgrid( y, z )

##############################
#
# plot
#
##############################

pipe = rectangle( a, b )

print( pipe.discharge() )

contourf( y, z, pipe.velocity( y, z ) )

y = linspace( -a, a, 15 )
z = linspace( -b, 0, int( len( y )*b/a ) )
y, z = meshgrid( y, z )

quiver( y, z, pipe.dyu( y, z ), pipe.dzu( y, z ), color = 'w' )

axis('equal')

show()
