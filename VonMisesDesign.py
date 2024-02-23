from VMPlot import VMPlot
from VonMises import VonMises

center = 0+0j
singularities = [1+0j, -1+0j]
#singularities = [1+0j, -1+0j, .37-.471j, -.37+.471j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j]
#singularities = [1+0j, -1+0j, .37-.471j, -.37+.471j]

v_inf = 1
alpha_rad = 0

#vm1 = VonMises(-.072+.29j, v_inf, alpha_rad, singularities)
vm1 = VonMises(center, v_inf, alpha_rad, singularities)

vmPlot = VMPlot(vm1)

### Outstanding Items ###
# Make X and Y Limits Scale Automatically
# Add Import/Export Options when you run the program
# Display exact singularity/center values. 
# Can I have an option to compute drag with boundary layer equations?
# Run online with PyScript or make an compile an executable. Executable might be better