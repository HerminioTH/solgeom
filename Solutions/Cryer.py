#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 16:28:51 2018

@author: herminio
         Adapted from Ricieri's code.
"""


import numpy as np
import math;


# CLASS DEFINITION ==============================================================================

class Solution( object ):
    def __init__( self, radius, sigma, solid, fluid ):
        self.radius = radius;
        self.sigma = sigma;

        fluid_name = list(fluid.keys())[0]
        self.mu = fluid[fluid_name]["Viscosity"]
        self.c_f = fluid[fluid_name]["Compressibility"]

        solid_name = list(solid.keys())[0]
        self.c_s = solid[solid_name]["Compressibility"]
        self.permeability = solid[solid_name]["Permeability"]
        self.phi = solid[solid_name]["Porosity"]
        self.ni = solid[solid_name]["PoissonsRatio"]
        self.G = solid[solid_name]["ShearModulus"]

        self._calculate_E();
        self._calculate_K_s();
        self._calculate_K_f();
        self._calculate_K_phi();
        self._calculate_K();

        self._calculate_H();
        self._calculate_R();
        self._calculate_alpha();
        self._calculate_K_p();
        self._calculate_B();
        self._calculate_K_u();
##        self._calculate_ni();
        self._calculate_ni_u();
        self._calculate_K_ni_u();
        self._calculate_gama();
        self._calculate_c();

        self._calculate_S();
        self._calculate_eta();
#        self._calculate_cv();

        self._computeInitialPressure();


    # Internal functions --------------------------------------------------------------------

    def _calculate_E(self):
        self.E = 2*self.G*(1 + self.ni)

    def _calculate_K_s(self):
        if self.c_s == 0.0:
            self.K_s = 1.0e+100
        else:
            self.K_s = 1.0 / self.c_s

    def _calculate_K_f(self):
        if self.c_f == 0.0:
            self.K_f = 1.0e+100
        else:
            self.K_f = 1.0 / self.c_f

    def _calculate_K_phi(self):
        self.K_phi = self.K_s

    def _calculate_K(self):
        self.K = 2*self.G*( 1 + self.ni ) / ( 3 - 6*self.ni )

    def _calculate_H(self):
        self.H = 1.0 / ( ( 1.0 / self.K ) - ( 1.0 / self.K_s ) )

    def _calculate_R(self):
        self.R = 1.0 / ( ( 1.0 / self.H ) + self.phi * ( ( 1.0 / self.K_f ) - ( 1.0 / self.K_phi ) ) )

    def _calculate_alpha(self):
        self.alpha = 1.0 - ( self.K / self.K_s )

    def _calculate_K_p(self):
        self.K_p = self.phi * self.K / self.alpha

    def _calculate_B(self):
        self.B = self.R / self.H

    def _calculate_K_u(self):
        self.K_u = self.K / ( 1.0 - self.alpha * self.B )

    def _calculate_ni(self):
        self.ni = ( 3.0 * self.K - 2.0 * self.G ) / ( 2.0 * ( 3.0 * self.K + self.G ) )

    def _calculate_ni_u(self):
        self.ni_u = ( ( 3.0 * self.ni + self.alpha * self.B * ( 1.0 - 2.0 * self.ni ) ) /
                          ( 3.0 - self.alpha * self.B * ( 1.0 - 2.0 * self.ni ) ) )

    def _calculate_K_ni_u(self):
        self.K_ni_u = ( 3.0 * self.K_u * ( 1.0 - self.ni_u ) ) / ( 1.0 + self.ni_u )

    def _calculate_c_s(self):
        self.c_s = 1.0 / self.K_s

    def _calculate_c_f(self):
        self.c_f = 1.0 / self.K_f

    def _calculate_gama(self):
        self.gama = ( self.B * ( 1.0 + self.ni_u ) ) / ( 3.0 * ( 1.0 - self.ni_u ) )

    def _calculate_c(self):
        self.c = ( ( 2.0 * self.permeability * self.G * ( 1.0 - self.ni ) * ( self.ni_u - self.ni ) ) /
                   ( self.mu * ( self.alpha ** 2.0 ) * ( 1.0 - self.ni_u ) * ( ( 1.0 - 2.0 * self.ni ) ** 2.0 ) ) )

    def _calculate_eta(self):
        self.eta = (self.K + 4*self.G/3.)/(2*self.G)
        self.eta *= (1 + self.K*self.S/self.alpha**2)

    def _calculate_S(self):
        self.S = self.phi*self.c_f + (self.alpha - self.phi)*self.c_s



    def _computeInitialPressure(self):
        self.initialPressure = self.sigma*self.alpha/(self.alpha**2 + self.K*self.S)

    def _calculateFunctionToFindTheRoots(self, x):
        return (1 - self.eta*x*x/2)*math.tan(x) - x

    def _calculateRoots(self, numberOfRoots, maxIterations=100, maxResidue=1.0e-12):
        roots = []
        for i in range( 1, numberOfRoots+1 ):
            x_A = i * np.pi - np.pi/3.
            x_B = i * np.pi + np.pi/3.
            x_C = ( x_A + x_B ) / 2
            y_C = self._calculateFunctionToFindTheRoots( x_C )
            iteration = 0
            residue = 1.0
            while iteration < maxIterations and residue > maxResidue and y_C != 0:
                y_A = self._calculateFunctionToFindTheRoots( x_A )
                y_C = self._calculateFunctionToFindTheRoots( x_C )
                if y_A * y_C > 0.0:
                        x_A = x_C
                else:
                        x_B = x_C
                x_C = ( x_A + x_B ) / 2
                residue = x_B - x_A
                iteration += 1
            roots.append( x_C )
        return roots




    # Class interface -----------------------------------------------------------------------
    def getNonDimensionalTimeValues(self, timeValues):
        return (self.c/self.radius**2)*timeValues

    def getPressureValue(self, time, numberOfSummationTerms=200, roots=[]):
        if len(roots) == 0:
            roots = self._calculateRoots( numberOfSummationTerms )
        summationResult = 0.0
        for i in range(numberOfSummationTerms):
            root = roots[i]
            num = (math.sin(root) - root) * math.exp(-root**2*self.c*time/self.radius**2)
            den = (self.eta - 1)*math.sin(root) + self.eta*root*math.cos(root)/2.
            summationResult += num/den
        pressureValue = self.initialPressure*self.eta*summationResult
        return pressureValue

    def getPressureValues(self, timeValues, numberOfSummationTerms=200):
        roots = self._calculateRoots(numberOfSummationTerms)
        pressureValues = np.array([self.getPressureValue(time, numberOfSummationTerms, roots) for time in timeValues])
        return pressureValues

    def getNormalizedPressureValues(self, timeValues, numberOfSummationTerms=200):
        roots = self._calculateRoots(numberOfSummationTerms)
        pressureValues = []
        for i in range( 0, timeValues.size ):
            pressureValue = self.getPressureValue(timeValues[i], numberOfSummationTerms, roots)
            pressureValues.append(pressureValue/self.initialPressure)
        return pressureValues




if __name__ == '__main__':
    import pylab as pl
    import sys, os
    os.chdir("..")
    current_dir = os.getcwd()
    sys.path.append(current_dir+'/PhysicalPropertyTools')
    from PropertyParser import Properties

    radius = 1.0
    sigma = 1000.0
    rock = Properties( current_dir + "/PhysicalPropertyTools/Json_Files/solid.json" )
    fluid = Properties( current_dir + "/PhysicalPropertyTools/Json_Files/fluid.json" )

    c = Solution( radius, sigma, rock, fluid )

    t = np.logspace(-1, 4, 100)
    tn = c.getNormalizedTimeValues(t)
    pn2 = c.getNormalizedPressureAtCentre(t, 250)

    pl.semilogx(tn, pn2, 'b-')
    pl.grid(True)
    pl.show()
