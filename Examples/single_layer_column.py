import matplotlib.pyplot as plt
import numpy as np
import json
import sys
sys.path.append("../Solutions")
from Terzaghi import Solution

mm = 1000.
kPa = 1/1000.



height = 10.0
load = 1.0e+4
gravity = 0.0

rock = json.load(open("solid.json", "r"))
fluid = json.load(open("fluid.json", "r"))



terza = Solution(height, load, rock, fluid , gravity)
z = terza.getPositionValues()



fig, ax = plt.subplots(2, 2, figsize=(8,7))
fig.subplots_adjust(left=0.070, right=0.975, top=0.970, bottom=0.065, hspace=0.235, wspace=0.300)
fonts = {'fontname': 'serif'}


times_1 = [0., 1., 10., 30., 80.]
for t in times_1:
	p = terza.getPressureValuesConstTime(t)
	ax[0][0].plot(p*kPa, z)
ax[0][0].set_xlabel("Pressure (kPa)", size=12, **fonts)
ax[0][0].set_ylabel("Height (m)", size=12, **fonts)
ax[0][0].grid(True)

for t in times_1:
	w = terza.getDisplacementValuesConstTime(t)
	ax[0][1].plot(w*mm, z)
ax[0][1].set_xlabel("Displacement (mm)", size=12, **fonts)
ax[0][1].set_ylabel("Height (m)", size=12, **fonts)
ax[0][1].grid(True)


times_2 = np.linspace(0, 200., 100)
# for t in times_2:
p = terza.getPressureValuesAtPosition(0.0, times_2) # Bottom pressure (z=0.0)
ax[1][0].plot(times_2, p*kPa)
ax[1][0].set_xlabel("Time (s)", size=12, **fonts)
ax[1][0].set_ylabel("Bottom Pressure (kPa)", size=12, **fonts)
ax[1][0].grid(True)


w = terza.getDisplacementValuesAtPosition(height, times_2) # Top displacement (z=height)
ax[1][1].plot(times_2, w*mm)
ax[1][1].set_xlabel("Time (s)", size=12, **fonts)
ax[1][1].set_ylabel("Top Displacement", size=12, **fonts)
ax[1][1].grid(True)



plt.show()
