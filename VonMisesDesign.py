from VMPlot import VMPlot
from VonMises import VonMises

#singularities = [1+0j, -1+0j]
#singularities = [1+0j, -1+0j, .37-.471j, -.37+.471j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j]
singularities = [1+0j, -1+0j, .37-.471j, -.37+.471j]

v_inf = 1
alpha_rad = 0

vm1 = VonMises(-.072+.29j, v_inf, alpha_rad, singularities)

vmPlot = VMPlot(vm1)

### Outstanding Items ###
# Make X and Y Limits Scale Automatically
# Start with joukowski and have ability to Add/Remove points
# Toggle Pressure Dist instead of velocity dist
# Add Import/Export Options when you run the program
# Display exact singularity/center values. Maybe have sliders/manual entry box
# Compute lift from circulation
# Can I have an option to compute drag with boundary layer equations?
# Run online with PyScript or make an compile an executable. Executable might be better