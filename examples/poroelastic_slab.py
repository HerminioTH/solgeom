import matplotlib.pyplot as plt
import numpy as np
import json
import sys
sys.path.append("../Solutions")
from Mandel import Solution

mm = 1000.
kPa = 1/1000.


H = 1.0
L = 2.0
force = 1.0e+4

solid = json.load(open("solid.json", "r"))
fluid = json.load(open("fluid.json", "r"))

mandel = Solution(L, H, force, solid, fluid)
x = mandel.getXPositionValues()
y = mandel.getYPositionValues()



fig, ax = plt.subplots(2, 4, figsize=(16,7))
fig.subplots_adjust(left=0.040, right=0.990, top=0.970, bottom=0.085, hspace=0.235, wspace=0.300)
fonts = {'fontname': 'serif'}

# Plot pressure profiles ---------------------------------------------------
times_1 = [0.01, 1., 10., 30., 80.]
for t in times_1:
	p = mandel.getPressureValuesConstTime(t)
	ax[0][0].plot(x, p*kPa)
ax[0][0].set_xlabel("Pressure (kPa)", size=12, **fonts)
ax[0][0].set_ylabel("Height (m)", size=12, **fonts)
ax[0][0].grid(True)
# --------------------------------------------------------------------------

# Plot pressure over time --------------------------------------------------
times_2 = np.logspace(-2, 3., 100)
p = mandel.getPressureValuesAtPosition(0.0, times_2)
ax[1][0].semilogx(times_2, p*kPa)
ax[1][0].set_xlabel("Time (s)", size=12, **fonts)
ax[1][0].set_ylabel("Pressure at %s"%(r"$x=0$"), size=12, **fonts)
ax[1][0].grid(True)
# --------------------------------------------------------------------------

# Plot horizontal displacement profiles ------------------------------------
for t in times_1:
	u = mandel.getHorDisplacementValuesConstTime(t)
	ax[0][1].plot(x, u*mm)
ax[0][1].set_xlabel("x (m)", size=12, **fonts)
ax[0][1].set_ylabel("Hor. displacement (mm)", size=12, **fonts)
ax[0][1].grid(True)
# --------------------------------------------------------------------------

# Plot horizontal displacement at x=L --------------------------------------
u = mandel.getHorDisplacementValuesAtPosition(L, times_2)
ax[1][1].semilogx(times_2, u*mm)
ax[1][1].set_xlabel("Time (s)", size=12, **fonts)
ax[1][1].set_ylabel("Hor. displacement at %s"%(r"$x=L$"), size=12, **fonts)
ax[1][1].grid(True)
# --------------------------------------------------------------------------

# Plot vertical displacement profiles --------------------------------------
for t in times_1:
	v = mandel.getVertDisplacementValuesConstTime(t)
	ax[0][2].plot(v*mm, y)
ax[0][2].set_xlabel("Vert. displacement (mm)", size=12, **fonts)
ax[0][2].set_ylabel("y (m)", size=12, **fonts)
ax[0][2].grid(True)
# --------------------------------------------------------------------------

# Plot vertical displacement at y=H ----------------------------------------
v = mandel.getVertDisplacementValuesAtPosition(H, times_2)
ax[1][2].semilogx(times_2, v*mm)
ax[1][2].set_xlabel("Time (s)", size=12, **fonts)
ax[1][2].set_ylabel("Hor. displacement at %s"%(r"$x=L$"), size=12, **fonts)
ax[1][2].grid(True)
# --------------------------------------------------------------------------

# Plot vertical stress profiles --------------------------------------------
for t in times_1:
	sigma_y = mandel.getVertStressValuesConstTime(t)
	ax[0][3].plot(sigma_y*kPa, y)
ax[0][3].set_xlabel("Vert. stress (kPa)", size=12, **fonts)
ax[0][3].set_ylabel("y (m)", size=12, **fonts)
ax[0][3].grid(True)
# --------------------------------------------------------------------------

# # Plot vertical displacement at y=H ----------------------------------------
# v = mandel.getVertDisplacementValuesAtPosition(H, times_2)
# ax[1][2].semilogx(times_2, v*mm)
# ax[1][2].set_xlabel("Time (s)", size=12, **fonts)
# ax[1][2].set_ylabel("Hor. displacement at %s"%(r"$x=L$"), size=12, **fonts)
# ax[1][2].grid(True)
# # --------------------------------------------------------------------------






plt.show()
