import matplotlib.pyplot as plt
import numpy as np
import json
import sys
sys.path.append("../Solutions")
from Cryer import Solution


radius = 1.0
load = 1.0e+4

rock = json.load(open("solid.json", "r"))
fluid = json.load(open("fluid.json", "r"))

cryer = Solution(radius, load, rock, fluid)
times = np.logspace(-4, 1, 100)
p = cryer.getPressureValues(times)


fig, ax = plt.subplots(1, figsize=(5,4))
fig.subplots_adjust(left=0.120, right=0.975, top=0.970, bottom=0.135, hspace=0.235, wspace=0.300)
fonts = {'fontname': 'serif'}

ax.semilogx(times, p/cryer.initialPressure)
ax.set_xlabel("Time (s)", size=12, **fonts)
ax.set_ylabel("Normalized pressure %s"%(r"$(p/p_0)$"), size=12, **fonts)
ax.grid(True)

plt.show()
