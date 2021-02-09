import numpy as np
import pylab as pl
import math



class Solution( object ):
    def __init__(self, length, height, force, solid, fluid):
        self.length = length
        self.height = height
        self.F = force

        fluid_name = list(fluid.keys())[0]
        self.mu = fluid[fluid_name]["Viscosity"]
        self.c_f = fluid[fluid_name]["Compressibility"]

        solid_name = list(solid.keys())[0]
        self.c_s = solid[solid_name]["Compressibility"]
        self.permeability = solid[solid_name]["Permeability"]
        self.phi = solid[solid_name]["Porosity"]
        self.nu = solid[solid_name]["PoissonsRatio"]
        self.G = solid[solid_name]["ShearModulus"]


        self._calculate_K_s()
        self._calculate_K_f()
        self._calculate_K_phi()
        self._calculate_K()

        self._calculate_H()
        self._calculate_R()
        self._calculate_alpha()
        self._calculate_K_p()
        self._calculate_B()
        self._calculate_K_u()
        self._calculate_ni_u()
        self._calculate_c()


    # Internal functions --------------------------------------------------------------------

    def _calculate_K_s( self ):
        if self.c_s == 0.0:
            self.K_s = 1.0e+100
        else:
            self.K_s = 1.0 / self.c_s

    def _calculate_K_f( self ):
        if self.c_f == 0.0:
            self.K_f = 1.0e+100
        else:
            self.K_f = 1.0 / self.c_f

    def _calculate_K_phi( self ):
        self.K_phi = self.K_s

    def _calculate_K( self ):
        self.K = 2*self.G*( 1 + self.nu ) / ( 3 - 6*self.nu )

    def _calculate_H( self ):
        self.H = 1.0 / ( ( 1.0 / self.K ) - ( 1.0 / self.K_s ) )

    def _calculate_R( self ):
        self.R = 1.0 / ( ( 1.0 / self.H ) + self.phi * ( ( 1.0 / self.K_f ) - ( 1.0 / self.K_phi ) ) )

    def _calculate_alpha( self ):
        self.alpha = 1.0 - self.c_s*self.K

    def _calculate_K_p( self ):
        self.K_p = self.phi * self.K / self.alpha

    def _calculate_B( self ):
        self.B = self.R / self.H

    def _calculate_K_u( self ):
        self.K_u = self.K / ( 1.0 - self.alpha * self.B )

    def _calculate_ni( self ):
        self.nu = ( 3.0 * self.K - 2.0 * self.G ) / ( 2.0 * ( 3.0 * self.K + self.G ) )

    def _calculate_ni_u( self ):
        self.nu_u = ( ( 3.0 * self.nu + self.alpha * self.B * ( 1.0 - 2.0 * self.nu ) ) /
                          ( 3.0 - self.alpha * self.B * ( 1.0 - 2.0 * self.nu ) ) )

    def _calculate_c( self ):
        self.c = ( ( 2.0 * self.permeability * self.G * ( 1.0 - self.nu ) * ( self.nu_u - self.nu ) ) /
                   ( self.mu * ( self.alpha ** 2.0 ) * ( 1.0 - self.nu_u ) * ( ( 1.0 - 2.0 * self.nu ) ** 2.0 ) ) )


    def _calculateFunctionToFindTheRoots( self, x ):
            y = math.tan( x ) - ( ( 1 - self.nu ) / ( self.nu_u - self.nu ) ) * x
            return y


    def _calculateRoots(self, numberOfRoots, maxIterations=100, maxResidue=1.0e-12):
        roots = [ ]
        for i in range( 0, numberOfRoots ):
            x_A = i * math.pi + maxResidue
            x_B = i * math.pi + ( math.pi / 2 ) - maxResidue
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


    def __getPositionValuesAndSize(self, ny, funcPosition):
        if type( ny ) == int:
            positionValues = funcPosition( ny )
            size = ny
        elif type( ny ) == np.ndarray:
            positionValues = ny
            size = ny.size
        elif type( ny ) == list:
            positionValues = ny
            size = len( ny )
        return positionValues, size

    def __getPositionValues(self, D, n=200):
        '''If n is an interger, then a list of equally spaced vertical position will be created. However, n can be
            a list or an numpy array with specified vertical positions.'''
        if type(n) == int:
            return np.linspace(0, D, n)
        elif type(n)==np.ndarray or type(n) == list:
            return n
        else:
            raise


    def __getFieldValuesConstTime(self, fieldFunction, time, nx, numberOfSummationTerms):
        positionValues = self.getXPositionValues(nx)
        roots = self._calculateRoots(numberOfSummationTerms)
        fieldValues = np.zeros(len(positionValues))
        for i, pos in enumerate(positionValues):
            fieldValues[i] = fieldFunction(pos, time, numberOfSummationTerms, roots)
        return fieldValues


    def __getFieldValuesAtTime(self, fieldFunction, position, timeList, numberOfSummationTerms=200, timePoints=200):
        roots = self._calculateRoots( numberOfSummationTerms )
        fieldValues = np.zeros(len(timeList))
        for i, time in enumerate(timeList):
            fieldValues[i] = fieldFunction(position, time, numberOfSummationTerms, roots)
        return fieldValues



    # Class interface -----------------------------------------------------------------------
    def getXPositionValues(self, n=200):
        return self.__getPositionValues(self.length, n)

    def getYPositionValues(self, n=200):
        return self.__getPositionValues(self.height, n)


    def getPressureValue(self, xPosition, time, numberOfSummationTerms=200, roots=[]):
        if time == 0.0:
            pressureValue = self.F*self.B*(1 + self.nu_u)/(self.length * 3)
            return pressureValue
        else:
            if len(roots) == 0:
                roots = self._calculateRoots(numberOfSummationTerms)
            summationResult = 0.0
            for i in range(numberOfSummationTerms):
                term_1 = math.sin( roots[i] ) / ( roots[i] - math.sin( roots[i] ) * math.cos( roots[i] ) )
                term_2 = math.cos( roots[i] * xPosition / self.length ) - math.cos( roots[i] )
                term_3 = math.exp( - ( ( self.c * time * ( roots[i] ** 2 ) ) / ( self.length ** 2 ) ) )
                summationResult += term_1 * term_2 * term_3
            pressureValue = ( ( 2 * self.F * self.B * ( 1 + self.nu_u ) ) / ( 3 * self.length ) ) * summationResult
            return pressureValue


    def getVertStressValue(self, xPosition, time, numberOfSummationTerms=200, roots=[]):
        if time == 0.0:
            vertStressValue = - self.F / self.length
            return vertStressValue
        else:
            if len(roots) == 0:
                    roots = self._calculateRoots(numberOfSummationTerms)
            summationResult_1 = 0.0
            summationResult_2 = 0.0
            for i in range(numberOfSummationTerms):
                term_1 = math.sin( roots[i] ) / ( roots[i] - math.sin( roots[i] ) * math.cos( roots[i] ) )
                term_2 = math.cos( roots[i] * xPosition / self.length )
                term_3 = math.exp( - ( ( self.c * time * ( roots[i] ** 2 ) ) / ( self.length ** 2 ) ) )
                summationResult_1 += term_1 * term_2 * term_3
            for i in range(numberOfSummationTerms):
                term_1 = ( math.sin( roots[i] ) * math.cos( roots[i] ) ) / ( roots[i] - math.sin( roots[i] ) * math.cos( roots[i] ) )
                term_2 = math.exp( - ( ( self.c * time * ( roots[i] ** 2 ) ) / ( self.length ** 2 ) ) )
                summationResult_2 += term_1 * term_2
            vertStressValue = ( - ( 2 * self.F * ( self.nu_u - self.nu ) ) / ( self.length * ( 1 - self.nu ) ) * summationResult_1 -
                          self.F / self.length +
                          2 * self.F / self.length * summationResult_2 )
            return vertStressValue


    def getHorDisplacementValue(self, xPosition, time, numberOfSummationTerms=200, roots=[]):
        if len(roots) == 0:
            roots = self._calculateRoots(numberOfSummationTerms)
        summationResult_1 = 0.0
        summationResult_2 = 0.0
        for i in range(numberOfSummationTerms):
            term_1 = ( math.sin( roots[i] ) * math.cos( roots[i] ) ) / ( roots[i] - math.sin( roots[i] ) * math.cos( roots[i] ) )
            term_2 = math.exp( - ( ( self.c * time * ( roots[i] ** 2 ) ) / ( self.length ** 2 ) ) )
            summationResult_1 += term_1 * term_2
        for i in range(numberOfSummationTerms):
            term_1 = math.cos( roots[i] ) / ( roots[i] - math.sin( roots[i] ) * math.cos( roots[i] ) )
            term_2 = math.sin( roots[i] * xPosition / self.length )
            term_3 = math.exp( - ( ( self.c * time * ( roots[i] ** 2 ) ) / ( self.length ** 2 ) ) )
            summationResult_2 += term_1 * term_2 * term_3
        firstTerm = ( ( self.F * self.nu ) / ( 2 * self.G * self.length ) - ( self.F * self.nu_u ) / ( self.G * self.length ) * summationResult_1 ) * xPosition
        secondTerm = self.F / self.G * summationResult_2
        horDisplacementValue = firstTerm + secondTerm
        return horDisplacementValue


    def getVertDisplacementValue(self, yPosition, time, numberOfSummationTerms=200, roots=[]):
        if len(roots) == 0:
            roots = self._calculateRoots(numberOfSummationTerms)
        summationResult = 0.0
        for i in range(numberOfSummationTerms):
            term_1 = ( math.sin( roots[i] ) * math.cos( roots[i] ) ) / ( roots[i] - math.sin( roots[i] ) * math.cos( roots[i] ) )
            term_2 = math.exp( - ( ( self.c * time * ( roots[i] ** 2 ) ) / ( self.length ** 2 ) ) )
            summationResult += term_1 * term_2
        vertDisplacementValue = ( ( self.F * ( 1.0 - self.nu_u ) ) / ( self.G * self.length ) * summationResult -
                                  ( self.F * ( 1.0 - self.nu ) ) / ( 2.0 * self.G * self.length ) ) * yPosition
        return vertDisplacementValue




    def getPressureValuesConstTime(self, time, nx=200, numberOfSummationTerms=200):
        return self.__getFieldValuesConstTime(self.getPressureValue, time, nx, numberOfSummationTerms)

    def getVertStressValuesConstTime(self, time, nx=200, numberOfSummationTerms=200):
        return self.__getFieldValuesConstTime(self.getVertStressValue, time, nx, numberOfSummationTerms)

    def getHorDisplacementValuesConstTime(self, time, nx=200, numberOfSummationTerms=200):
        return self.__getFieldValuesConstTime(self.getHorDisplacementValue, time, nx, numberOfSummationTerms)

    def getVertDisplacementValuesConstTime(self, time, ny=200, numberOfSummationTerms=200):
        return self.__getFieldValuesConstTime(self.getVertDisplacementValue, time, ny, numberOfSummationTerms)


    def getPressureValuesAtPosition(self, xPosition, timeList, numberOfSummationTerms=200, timePoints=200):
        return self.__getFieldValuesAtTime(self.getPressureValue, xPosition, timeList, numberOfSummationTerms)

    def getHorDisplacementValuesAtPosition(self, xPosition, timeList, numberOfSummationTerms=200, timePoints=200):
        return self.__getFieldValuesAtTime(self.getHorDisplacementValue, xPosition, timeList, numberOfSummationTerms)

    def getVertDisplacementValuesAtPosition(self, yPosition, timeList, numberOfSummationTerms=200, timePoints=200):
        return self.__getFieldValuesAtTime(self.getVertDisplacementValue, yPosition, timeList, numberOfSummationTerms)

    def getVertStressValuesAtPosition(self, xPosition, totalTimeInterval, numberOfSummationTerms=200, timePoints=200):
        return self.__getFieldValuesAtTime(self.getVertStressValue, xPosition, timeList, numberOfSummationTerms)







if __name__ == '__main__':
    import pylab as pl
    import sys, os
    os.chdir("..")
    current_dir = os.getcwd()
    sys.path.append(current_dir+'/PhysicalPropertyTools')
    from PropertyParser import Properties

    length = 5.0
    height = 1.0
    Force = 5.0e+7
    rock = Properties( current_dir + "/PhysicalPropertyTools/Json_Files/solid.json" )
    fluid = Properties( current_dir + "/PhysicalPropertyTools/Json_Files/fluid.json" )

    mandel = Solution( length, height, Force, rock, fluid )
    x = mandel.getXPositionValues()
    y = mandel.getYPositionValues()

    time = 5000.
    p = mandel.getPressureValuesConstTime(time, x)

    pl.plot(x, p)
    pl.grid(True)
    pl.show()
