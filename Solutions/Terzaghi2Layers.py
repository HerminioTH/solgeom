# ABREVIATIONS USED =============================================================================

# PROPERTIES OF THE UPPER PART

# height_1 - Height
# permeability_1 - Reference permeability
# phi_1 - Reference porosity
# K_1 - Bulk modulus
# K_s_1 - Solid bulk modulus
# G_1 - Sheer modulus
# K_phi_1 - Unjacket pore compressibility

# H_1 - Poroelastic expansion coefficient
# R_1 - Biot modulus
# alpha_1 - Biot-Willis coefficient
# K_p_1 - Drained pore compressibility
# K_u_1 - Undrained bulk modulus
# B_1 - Skempton"s coefficient
# nu_1 - Drained Poisson"s ratio
# nu_u_1 - Undrained Poisson"s ratio
# K_nu_u_1 - Uniaxial undrained bulk modulus
# c_s_1 - Solid phase compressibility
# m_1 - Confined compressibility
# gama_1 - Loading efficiency
# c_1 - Consolidation coefficient / Hydraulic diffusivity
# t_1 - Consolidation time

# PROPERTIES OF THE LOWER PART

# height_2 - Height
# permeability_2 - Reference permeability
# phi_2 - Reference porosity
# K_2 - Bulk modulus
# K_s_2 - Solid bulk modulus
# G_2 - Sheer modulus
# K_phi_2 - Unjacket pore compressibility

# H_2 - Poroelastic expansion coefficient
# R_2 - Biot modulus
# alpha_2 - Biot-Willis coefficient
# K_p_2 - Drained pore compressibility
# K_u_2 - Undrained bulk modulus
# B_2 - Skempton"s coefficient
# nu_2 - Drained Poisson"s ratio
# nu_u_2 - Undrained Poisson"s ratio
# K_nu_u_2 - Uniaxial undrained bulk modulus
# c_s_2 - Solid phase compressibility
# m_2 - Confined compressibility
# gama_2 - Loading efficiency
# c_2 - Consolidation coefficient / Hydraulic diffusivity
# t_2 - Consolidation time

# PROPERTIES OF THE FLUID

# mi - Viscosity
# K_f - Fluid bulk modulus

# c_f - Fluid phase compressibility
# p_0_1 - Initial pore pressure in the upper part
# p_0_2 - Initial pore pressure in the lower part

# GENERAL DATA

# tao_0 - Normal stress
# beta - Ratio between the consolidation times
# w - Ratio between the compressibilities and the permeabilities


import math, sys
import numpy as np


# CLASS DEFINITION ==============================================================================

class Solution(object):
    def __init__(self, height_1, height_2, tao_0, solid, fluid, p_0_1=0, p_0_2=0):
        self.tao_0 = tao_0;
        self.height_1 = height_1
        self.height_2 = height_2
        self.height = height_1 + height_2

        fluid_name = list(fluid.keys())[0]
        self.mu = fluid[fluid_name]["Viscosity"]
        self.c_f = fluid[fluid_name]["Compressibility"]
        self.rho_f = fluid[fluid_name]["Density"]
        self.K_f = self._calculate_K_f(self.c_f)

        solid_upper_layer = list(solid.keys())[0]
        self.c_s_1 = solid[solid_upper_layer]["Compressibility"]
        self.permeability_1 = solid[solid_upper_layer]["Permeability"]
        self.phi_1 = solid[solid_upper_layer]["Porosity"]
        self.nu_1 = solid[solid_upper_layer]["PoissonsRatio"]
        self.G_1 = solid[solid_upper_layer]["ShearModulus"]

        self.K_s_1 = self._calculate_K_s(self.c_s_1)
        self.K_phi_1 = self._calculate_K_phi(self.K_s_1)
        self.K_1 = self._calculate_K(self.G_1, self.nu_1)
        self.alpha_1 = self._calculate_alpha(self.c_s_1, self.K_1)
        self.H_1 = self._calculate_H(self.K_1, self.K_s_1)
        self.R_1 = self._calculate_R(self.H_1, self.phi_1, self.K_f, self.K_phi_1)
        self.K_p_1 = self._calculate_K_p(self.phi_1, self.K_1, self.alpha_1)
        self.B_1 = self._calculate_B(self.R_1, self.H_1)
        self.K_u_1 = self._calculate_K_u(self.alpha_1, self.B_1, self.K_1)
        self.nu_u_1 = self._calculate_nu_u(self.nu_1, self.alpha_1, self.B_1)
        self.K_nu_u_1 = self._calculate_K_nu_u(self.K_u_1, self.nu_u_1)
        self.m_1 = self._calculate_m(self.K_1, self.G_1)
        self.gama_1 = self._calculate_gama(self.B_1, self.nu_u_1)
        self.c_1 = self._calculate_c(self.permeability_1, self.G_1, self.nu_1, self.nu_u_1, self.mu, self.alpha_1)
        self.t_1 = self._calculate_t(self.height_1, self.c_1)



        solid_lower_layer = list(solid.keys())[1]
        self.c_s_2 = solid[solid_lower_layer]["Compressibility"]
        self.permeability_2 = solid[solid_lower_layer]["Permeability"]
        self.phi_2 = solid[solid_lower_layer]["Porosity"]
        self.nu_2 = solid[solid_lower_layer]["PoissonsRatio"]
        self.G_2 = solid[solid_lower_layer]["ShearModulus"]

        self.K_s_2 = self._calculate_K_s(self.c_s_2)
        self.K_phi_2 = self._calculate_K_phi(self.K_s_2)
        self.K_2 = self._calculate_K(self.G_2, self.nu_2)
        self.alpha_2 = self._calculate_alpha(self.c_s_2, self.K_2)
        self.H_2 = self._calculate_H(self.K_2, self.K_s_2)
        self.R_2 = self._calculate_R(self.H_2, self.phi_2, self.K_f, self.K_phi_2)
        self.K_p_2 = self._calculate_K_p(self.phi_2, self.K_2, self.alpha_2)
        self.B_2 = self._calculate_B(self.R_2, self.H_2)
        self.K_u_2 = self._calculate_K_u(self.alpha_2, self.B_2, self.K_2)
        self.nu_u_2 = self._calculate_nu_u(self.nu_2, self.alpha_2, self.B_2)
        self.K_nu_u_2 = self._calculate_K_nu_u(self.K_u_2, self.nu_u_2)
        self.m_2 = self._calculate_m(self.K_2, self.G_2)
        self.gama_2 = self._calculate_gama(self.B_2, self.nu_u_2)
        self.c_2 = self._calculate_c(self.permeability_2, self.G_2, self.nu_2, self.nu_u_2, self.mu, self.alpha_2)
        self.t_2 = self._calculate_t(self.height_2, self.c_2)

        if p_0_1 == 0:
            self.p_0_1 = self._calculate_p_0(self.gama_1, self.tao_0)
            self.p_0_2 = self._calculate_p_0(self.gama_2, self.tao_0)
        else:
            self.p_0_1 = p_0_1
            self.p_0_2 = p_0_2

        self._calculate_beta()
        self._calculate_w()

    def _calculate_K_s(self, c_s):
        if c_s == 0.0:
            return 1.0e+100
        else:
            return 1.0 / c_s

    def _calculate_K_f(self, c_f):
        if c_f == 0.0:
            return 1.0e+100
        else:
            return 1.0 / c_f

    def _calculate_K_phi(self, K_s):
        return K_s

    def _calculate_K(self, G, ni):
        return 2*G*(1 + ni) / (3 - 6*ni)

    def _calculate_alpha(self, c_s, K):
        return 1.0 - c_s*K

    def _calculate_H(self, K, K_s):
        return 1.0 / ((1.0 / K) - (1.0 / K_s))

    def _calculate_R(self, H, phi, K_f, K_phi):
        return 1.0 / ((1.0 / H) + phi * ((1.0 / K_f) - (1.0 / K_phi)))

    def _calculate_K_p(self, phi, K, alpha):
        return phi * K / alpha

    def _calculate_B(self, R, H):
        return R / H

    def _calculate_K_u(self, alpha, B, K):
        if alpha * B != 1.0:
            return K / (1.0 - alpha * B)
        else:
            return 1.0e+100

    def _calculate_nu_u(self, nu, alpha, B):
        value = ((3.0 * nu + alpha * B * (1.0 - 2.0 * nu)) /
                  (3.0 - alpha * B * (1.0 - 2.0 * nu)))
        return value

    def _calculate_K_nu_u(self, K_u, nu_u):
        return (3.0 * K_u * (1.0 - nu_u)) / (1.0 + nu_u)

    def _calculate_gama(self, B, nu_u):
        return (B * (1.0 + nu_u)) / (3.0 * (1.0 - nu_u))

    def _calculate_c(self, permeability, G, nu, nu_u, mu, alpha):
        c = ((2.0 * permeability * G * (1.0 - nu) * (nu_u - nu)) /
                (mu * (alpha ** 2.0) * (1.0 - nu_u) * ((1.0 - 2.0 * nu) ** 2.0)))
        return c

    def _calculate_K_v_u(self, K_u, nu_u):
        return 3*K_u*(1 - nu_u) / (1 + nu_u)

    def _calculate_c_m(self, alpha, nu, K):
        return (alpha * (1.0 + nu)) / (3.0 * K * (1.0 - nu))

    def _calculate_m(self, K, G):
        return 1.0 / (K + G * (4.0 / 3.0))

    def _calculate_t(self, height, c):
        return (height ** 2) / c

    def _calculate_beta(self):
        self.beta = math.sqrt(self.t_1 / self.t_2)

    def _calculate_w(self):
        self.w = math.sqrt(self.permeability_2 * self.m_2) / math.sqrt(self.permeability_1 * self.m_1)

    def _calculate_p_0(self, gama, tao_0):
        return gama * tao_0

    def _calculateFunctionToFindTheRoots(self, x):
        y = - self.w * math.sin(self.beta * x) * math.sin(x) + math.cos(self.beta * x) * math.cos(x)
        return y

    def _calculateRoots(self, numberOfRoots, maxIterations = 100, increment = 0.001):
        roots = []
        while len(roots) < numberOfRoots:
            if len(roots) == 0:
                A = 0
            else:
                if self._calculateFunctionToFindTheRoots(B) != 0:
                    A = B
                else:
                    A = B + increment
            B = A + increment
            y_A = self._calculateFunctionToFindTheRoots(A)
            y_B = self._calculateFunctionToFindTheRoots(B)
            while y_A * y_B > 0:
                A = B
                B = A + increment
                y_A = self._calculateFunctionToFindTheRoots(A)
                y_B = self._calculateFunctionToFindTheRoots(B)
            if y_B == 0:
                roots.append(B)
            else:
                C = A
                D = B
                E = (C + D) / 2
                y_E = self._calculateFunctionToFindTheRoots(E)
                iteration = 0
                while iteration < maxIterations and y_E != 0:
                    y_C = self._calculateFunctionToFindTheRoots(C)
                    y_E = self._calculateFunctionToFindTheRoots(E)
                    if y_C * y_E < 0:
                            D = E
                    else:
                            C = E
                    E = (C + D) / 2
                    iteration += 1
                roots.append(E)
        return roots

    def _calculateCommonDivisors(self, roots):
        commonDivisors = []
        for i in range(len(roots)):
            term_1 = (1.0 + self.w * self.beta) * math.cos(self.beta * roots[i]) * math.sin(roots[i])
            term_2 = (self.w + self.beta) * math.sin(self.beta * roots[i]) * math.cos(roots[i])
            commonDivisor = roots[i] * (term_1 + term_2)
            commonDivisors.append(commonDivisor)
        return commonDivisors


    # Class interface -----------------------------------------------------------------------
    def getPositionValues(self, ny=200):
        '''If ny is an interger, then a list of equally spaced vertical position will be created. However, ny can be
            a list or an numpy array with specified vertical positions.'''
        if type(ny) == int:
            return np.linspace(0, self.height, ny)
        elif type(ny)==np.ndarray or type(ny) == list:
            return ny
        else:
            raise


    def getPressureValue(self, yPosition, time, numberOfSummationTerms=200, roots=[], commonDivisors=[]):
        position = yPosition - self.height_2
        if len(commonDivisors) == 0:
            roots = self._calculateRoots(numberOfSummationTerms)
            commonDivisors = self._calculateCommonDivisors(roots)
        summationResult = 0
        if yPosition <= self.height_2:
            for j in range(numberOfSummationTerms):
                term_1 = math.cos(roots[j]) * math.cos(roots[j] * position / self.height_2)
                term_2 = math.sin(roots[j]) * math.sin(roots[j] * position / self.height_2)
                term_3 = math.exp(- (roots[j]**2) * time / self.t_2)
                summationResult += term_3 * (term_1 - term_2) / commonDivisors[j]
        else:
            for j in range(numberOfSummationTerms):
                term_1 = math.cos(roots[j]) * math.cos(self.beta * roots[j] * position / self.height_1)
                term_2 = math.sin(roots[j]) * math.sin(self.beta * roots[j] * position / self.height_1)
                term_3 = math.exp(- (roots[j]**2) * time / self.t_2)
                summationResult += term_3 * (term_1 - term_2*self.w) / commonDivisors[j]
        pressureValue = 2.0 * self.p_0_2 * summationResult
        return pressureValue

    def getPressureValuesConstTime(self, time, numberOfSummationTerms=200, ny=200):
        numberOfSummationTerms = 6000 if time==0.0 else numberOfSummationTerms
        positionValues = self.getPositionValues(ny)
        roots = self._calculateRoots(numberOfSummationTerms)
        commonDivisors = self._calculateCommonDivisors(roots)
        pressureValues = np.zeros(len(positionValues))
        for i, pos in enumerate(positionValues):
            pressureValues[i] = self.getPressureValue(pos, time, numberOfSummationTerms, roots, commonDivisors)
        return pressureValues

    def getPressureValuesAtPosition(self, position, timeValues, numberOfSummationTerms=200, timePoints=200):
        roots = self._calculateRoots(numberOfSummationTerms)
        commonDivisors = self._calculateCommonDivisors(roots)
        pressureValues = np.zeros(len(timeValues))
        for i, time in enumerate(timeValues):
            pressureValues[i] = self.getPressureValue(position, time, numberOfSummationTerms, roots, commonDivisors)
        return pressureValues
