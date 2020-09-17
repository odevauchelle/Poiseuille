
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

from numpy import cosh, sinh, sin, cos, pi, tanh

#############################
#
# series expansion
#
#############################

class rectangle :

    def __init__( self, a, b ) :

        '''
        Poiseuille flow in a rectangular pipe.

        pipe = rectangle( a, b )

        Parameters :
            a, b (float) : dimensions
        '''

        self.a = a
        self.b = b

        self.reference = '''
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

    def velocity_expansion_term( self, y, z, i ) :
        return ( -1 )**( ( i - 1 )/2 )*( 1 - cosh( i*pi*z/( 2*self.a ) )/cosh( i*pi*self.b/( 2*self.a ) ) )*cos( i*pi*y/( 2*self.a ) )/i**3

    def velocity_sum( self, y, z, imax, expansion_term ) :

        result = 0*y

        for i in range( 1, imax + 1, 2 ) :
            result += expansion_term( y, z, i )

        return result*16*self.a**2/pi**3

    def velocity( self, y, z, imax = 30 ) :
        return self.velocity_sum( y, z, imax, self.velocity_expansion_term )

    def dzu_expansion_term( self, y, z, i ) :
        return ( -1 )**( ( i - 1 )/2 )*( - i*pi/( 2*self.a )*sinh( i*pi*z/( 2*self.a ) )/cosh( i*pi*self.b/( 2*self.a ) ) )*cos( i*pi*y/( 2*self.a ) )/i**3

    def dyu_expansion_term( self, y, z, i ) :
        return -( -1 )**( ( i - 1 )/2 )*( 1 - cosh( i*pi*z/( 2*self.a ) )/cosh( i*pi*self.b/( 2*self.a ) ) )*i*pi/( 2*self.a )*sin( i*pi*y/( 2*self.a ) )/i**3

    def dyu( self, y, z, imax = 30 ) :
        return self.velocity_sum( y, z, imax, self.dyu_expansion_term )

    def dzu( self, y, z, imax = 30 ) :
        return self.velocity_sum( y, z, imax, self.dzu_expansion_term )

    def discharge_expansion_term( self, i ) :
        return tanh( i*pi*self.b/( 2*self.a ) )/i**5

    def discharge( self, imax = 30 ) :

        result = 0

        for i in range( 1, imax + 1, 2 ) :
            result += self.discharge_expansion_term( i )

        result *= -192*self.a/( pi**5*self.b )
        result += 1
        result *= 4*self.b*self.a**3/3

        return result

if __name__ == '__main__' :

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

    pipe = rectangle( a, b )

    contourf( y, z, pipe.velocity( y, z ) )

    axis('equal')

    show()
