from VMPlot import VMPlot
from VonMises import VonMises

center = 0+0j
singularities = [1+0j, -1+0j]

v_inf = 1
alpha_rad = 0

vm1 = VonMises(center, v_inf, alpha_rad, singularities)

vmPlot = VMPlot(vm1)

### Outstanding Items ###
# Make X and Y Limits Scale Automatically
# Try to optimize
# Can I have an option to compute drag with boundary layer equations?
# Run online with PyScript or make an compile an executable. Executable might be better
# Change File Write commands to Print in pyscript version

### Notes ###
# Sometimes weird artifacts show up in pressure dist when running in XFOIL
# PANE command before running seems to remove these

# To Compile: pyinstaller VonMisesDesign.py --onefile