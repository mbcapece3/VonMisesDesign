from VMPlot import VMPlot
from VonMises import VonMises

center = 0+0j
singularities = [1+0j, -1+0j]

v_inf = 1
alpha_rad = 0

vm1 = VonMises(center, v_inf, alpha_rad, singularities)

vmPlot = VMPlot(vm1)

### Outstanding Items ###
# Try to optimize
# Rotate to ensure chord parallel to real axis. Is aoa inherently relative to chord in the match?
# Calculate Cm
# Calculate Cd with boundary layer equations?
# Test EXE on windows 11 CPU

### Notes ###
# Sometimes artifacts show up in pressure dist when running in XFOIL. PANE command before running seems to remove these
# To Compile: cxfreeze -c VonMisesDesign.py 