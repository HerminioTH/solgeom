import matplotlib.pyplot as plt
import numpy as np
import json
import sys
sys.path.append("../Solutions")
from Terzaghi2Layers import Solution

mm = 1000.
kPa = 1/1000.


# 		     | load
# 		     |
# 		 +---V---+  ---
# 		 |       |   |
# 		 |       |   |
# 		 |       |   |
# 		 | UPPER |   | h_upper
# 		 |       |   |
# 		 |       |   |
# 		 +-------+  ---
# 		 |       |   |
# 		 |       |   | h_lower
#  z ^	 | LOWER |   |
# 	 |	 |       |   |
# 	_|_  |_______|  _|_



h_upper = 6.0
h_lower = 4.0
load = 1.0e+4

solid = json.load(open("solid.json", "r"))
fluid = json.load(open("fluid.json", "r"))

terza = Solution(h_upper, h_lower, load, solid, fluid)
z = terza.getPositionValues()


# Define axes for plot -----------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,3))
fig.subplots_adjust(left=0.070, right=0.975, top=0.970, bottom=0.160, hspace=0.235, wspace=0.300)
fonts = {'fontname': 'serif'}
# --------------------------------------------------------------------------

# Plot pressure profiles ---------------------------------------------------
times_1 = [1., 10., 80., 180., 500.]
for t in times_1:
	p = terza.getPressureValuesConstTime(t, numberOfSummationTerms=50)
	ax1.plot(p*kPa, z)
ax1.set_xlabel("Pressure (kPa)", size=12, **fonts)
ax1.set_ylabel("Height (m)", size=12, **fonts)
ax1.grid(True)
# --------------------------------------------------------------------------

# Plot bottom and mid pressure over time -----------------------------------
times_2 = np.logspace(1, 3.5, 1000)
p_bottom = terza.getPressureValuesAtPosition(0.0, times_2)
p_mid = terza.getPressureValuesAtPosition(h_lower, times_2)
ax2.semilogx(times_2, p_mid*kPa, label="z=4.0")
ax2.semilogx(times_2, p_bottom*kPa, label="z=0.0")
ax2.set_xlabel("Time (s)", size=12, **fonts)
ax2.set_ylabel("Pressure (kPa)", size=12, **fonts)
ax2.grid(True)
ax2.legend(loc=0, fancybox=True, shadow=True)
# --------------------------------------------------------------------------




plt.show()
