from VMPlot import VMPlot
from VonMises import VonMises

center = 0+0j
singularities = [1+0j, -1+0j]

v_inf = 1
alpha_rad = 0

vm1 = VonMises(center, v_inf, alpha_rad, singularities)

vmPlot = VMPlot(vm1)

### Current Outstanding Items ###
# Test EXE on windows 11 CPU
# Normalize airfoil x axis from 0 to 1.
# Calculate Cm
# Add title/credits when exe runs
# Update Readme. Provide download link for compiled version. Show example image. 
# Make Public

### Future Outstanding Items ###
# Rotate airfoil to ensure chord parallel to real axis by multiplying points by e^(-1j*theta) and fix TE at (1+0j). Assume LE is leftmost pt. Rotation angle theta can be determined by line btwn LE and TE. AOA must be adjusted to compensate for rotation.
    # XFOIL shows that aoa is relative to x axis, not chord.
# Try to optimize performance
# Calculate Cd with boundary layer equations?

### Notes ###
# Sometimes artifacts show up in pressure dist when running in XFOIL. PANE command before running seems to remove these
# To Compile: cxfreeze -c VonMisesDesign.py 