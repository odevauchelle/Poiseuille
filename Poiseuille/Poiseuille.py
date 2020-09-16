
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Olivier Devauchelle

from numpy import cosh, cos, pi

#############################
#
# series expansion
#
#############################

bibtex_reference = '''
@book{white1991viscous,
  title={Viscous fluid flow},
  author={White, Frank M.},
  year={1991},
  ISBN={0-07-069712-4},
  publisher={McGraw-Hill, Inc.},
  series={Mechanical Engineering},
  page={120},
  edition={2},
}
'''

def expansion_term_rectangle_pipe( i, y, z, a, b ) :
    return ( -1 )**( ( i - 1 )/2 )*( 1 - cosh( i*pi*z/( 2*a ) )/cosh( i*pi*b/( 2*a ) ) )*cos( i*pi*y/( 2*a ) )/i**3

def flow_in_rectangle_pipe( y, z, a, b, imax = 30 ) :

    '''
    Poiseuille flow in a rectangular pipe.

    velocity = expansion( y, z, a, b, imax = 30 )

    Parameters :
        y, z (float or array) : coordinates
        a, b (float) : dimensions
        imax (int) : number of terms in series (only odd terms matter)

    Output :
        velocity (float or array): Poiseuille velocity
    '''

    result = 0*y

    for i in range( 1, imax + 1, 2 ) :
        result += expansion_term_rectangle_pipe( i, y, z, a, b )

    return result*16*a**2/pi**3


if __name__ == '__main__' :

    flow_in_rectangle_pipe( 1,1,1,1 )

    from pylab import *

    ##############################
    #
    # PARAMETERS
    #
    ##############################

    D = 0.1
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

    contourf( y, z, flow_in_rectangle_pipe( y, z, a, b ) )

    axis('equal')

    show()
