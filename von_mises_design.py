from InteractivePlot import VMPlot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
from timeit import default_timer as timer

class VonMises:
     
    def __init__(self, center, v_inf, alpha_rad, singularities=[1+0j, -1+0j, 0, 0], num_pts=200):
        self.circle_TE = 1 + 0j  # fix trailing edge at (1,0)
        self.center = center # complex number
        self.singularities = singularities # list of singularities as complex points
        self.v_inf = v_inf
        self.alpha_rad = alpha_rad
        self.radius = abs(self.circle_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)
        self.num_pts = num_pts  # number of discreet points in circle
        self.circle_LE = self.radius * np.exp(1j * (np.pi + self.beta)) + self.center
        
        circle_data = {
            "angles": np.linspace(-1*self.beta, -1*self.beta + 2*np.pi, self.num_pts),
        }

        self.circle_df = pd.DataFrame(circle_data)

        self.circle_df["circle_pts"] = self.radius * np.exp(self.circle_df.angles*1j)+self.center

        self.C1 = self.compute_vonmises_coeff(1)
        self.C2 = self.compute_vonmises_coeff(2)
        self.C3 = self.compute_vonmises_coeff(3)


    def compute_vonmises_coeff(self, coeff_num):
            # Compute an Arbitrary Von Mises Transform Coefficient

            num_terms = coeff_num + 1
            k = len(self.singularities)
            sum = 0
        
            for indices in itertools.combinations(range(k), num_terms):
                product = 1
                for idx in indices:
                    product *= self.singularities[idx]
                sum += product

            if coeff_num % 2:
                return -1 / coeff_num * sum
            else:
                return 1 / coeff_num * sum 

    def updateCircle(self):
        self.radius = abs(self.circle_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)
        self.circle_LE = self.radius * np.exp(1j * (np.pi + self.beta)) + self.center

        self.circle_df["angles"] = np.linspace(-1*self.beta, -1*self.beta + 2*np.pi, self.num_pts)
        self.circle_df["circle_pts"] = self.radius * np.exp(self.circle_df.angles*1j)+self.center

    def updateCoefficients(self):
        self.C1 = self.compute_vonmises_coeff(1)
        self.C2 = self.compute_vonmises_coeff(2)
        self.C3 = self.compute_vonmises_coeff(3)

    def setV_inf(self, v_inf):
        self.v_inf = v_inf

    def setAlpha_rad(self, alpha_rad):
        self.alpha_rad = alpha_rad

    def conformalMap(self):

        # Map Airfoil Geometry
        self.circle_df["mapped_pts"] = self.circle_df.circle_pts + self.C1/self.circle_df.circle_pts + self.C2/(self.circle_df.circle_pts**2) + self.C3/(self.circle_df.circle_pts**3)

        # Scale Airfoil Geometry to X/C
        self.airfoil_TE = self.circle_TE + self.C1/self.circle_TE + self.C2/self.circle_TE**2 + self.C3/self.circle_TE**3
        self.airfoil_LE = self.circle_LE + self.C1/self.circle_LE + self.C2/self.circle_LE**2 + self.C3/self.circle_LE**3
        self.chord = abs(self.airfoil_TE - self.airfoil_LE)
        self.circle_df["airfoil_pts"] = (np.real(self.circle_df.mapped_pts) - np.real(self.airfoil_LE)) / self.chord + (np.imag(self.circle_df.mapped_pts) / self.chord)*1j
        
        # Map Airfoil Surface Velocities 
        self.circle_df["circle_vels"] = 2*v_inf*(np.sin(alpha_rad + self.beta) - np.sin(alpha_rad - self.circle_df.angles))
        self.circle_df["map_func_derivatives"] = 1 - self.C1/(self.circle_df.circle_pts**2) - 2*self.C2/(self.circle_df.circle_pts**3) - 3*self.C3/(self.circle_df.circle_pts**4)  # first derivative of airfoil pts wrt circle pts
        self.circle_df["airfoil_vels"] = abs(self.circle_df.circle_vels) / abs(self.circle_df.map_func_derivatives)

        # Manually Set TE_Vel with L'hopitals rule because indeterminate
        TE_map_func_second_deriv = 2*self.C1/(self.circle_TE**3) + 6*self.C2/(self.circle_TE**4) + 12*self.C3/(self.circle_TE**5)  # second derivative of airfoil pts wrt circle pts evaluated at trailing edge
        TE_vel = abs(2*v_inf*np.cos(alpha_rad+self.beta)/self.radius) / abs(TE_map_func_second_deriv)  # velocity at trailing edge
        self.circle_df.at[0, 'airfoil_vels'] = TE_vel
        self.circle_df.at[self.num_pts-1, 'airfoil_vels'] = TE_vel

        #print(self.circle_df)

    def plotMapping(self):
        fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(15,5), gridspec_kw={'height_ratios': [1]})
        ax1.plot(np.real(self.circle_df.circle_pts), np.imag(self.circle_df.circle_pts))
        ax1.scatter(np.real(self.center), np.imag(self.center))
        ax1.scatter(np.real(self.singularities),np.imag(self.singularities))
        ax1.set_aspect('equal')
        ax1.set_xlim(-2,2)
        ax1.set_ylim(-2,2)
        ax1.grid()
        ax1.vlines(0,-2,2,color='black')
        ax1.hlines(0,-2,2,color='black')

        ax2.plot(np.real(self.circle_df.airfoil_pts), np.imag(self.circle_df.airfoil_pts))
        ax2.set_aspect('equal')
        ax2.set_xlim(-.1,1.1)
        ax2.set_ylim(-.6,.6)
        ax2.grid()
        ax2.vlines(0,-2,2,color='black')
        ax2.hlines(0,-2,2,color='black')

        ax3.plot(np.real(self.circle_df.airfoil_pts), self.circle_df.airfoil_vels/self.v_inf)
        ax3.set_aspect(.5)
        ax3.set_xlim(0,1.1)
        ax3.set_ylim(0,2.2)
        ax3.grid()
        ax3.vlines(0,-2,2,color='black')
        ax3.hlines(0,-2,2,color='black')



        
#singularities = [1+0j, -1+0j]
#singularities = [1+0j, -1+0j, .37-.471j, -.37+.471j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j, 0+0j]
singularities = [1+0j, -1+0j, .37-.471j, -.37+.471j]
#singularities = [1+0j, .4+.25j, -.95+.15j, -.45-.4j]
#singularities = [1+0j, .4+.65j, -.95-.05j, -.45-.6j]

v_inf = 1
alpha_rad = 0
 
vm1 = VonMises(-.072+.29j, v_inf, alpha_rad, singularities)
vm1.conformalMap()

#vm1.plotMapping() # still image

epsilon = 5 #max pixel distance

fig1,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,5))

vmPlot = VMPlot(vm1, fig1, ax1, ax2, ax3, epsilon)

ax1.scatter (np.real(vm1.center), np.imag(vm1.center),color='r',marker='o')
ax1.scatter (np.real(vm1.singularities),np.imag(vm1.singularities),color='k',marker='o')
ax1.plot(np.real(vm1.circle_df.circle_pts), np.imag(vm1.circle_df.circle_pts))

ax1.set_aspect('equal')
ax1.set_xlim(-2,2)
ax1.set_ylim(-2,2)
ax1.grid()
ax1.vlines(0,-2,2,color='black')
ax1.hlines(0,-2,2,color='black')

ax2.plot(np.real(vm1.circle_df.airfoil_pts), np.imag(vm1.circle_df.airfoil_pts))
ax2.set_aspect('equal')
ax2.set_xlim(-.1,1.1)
ax2.set_ylim(-.6,.6)
ax2.grid()
ax2.vlines(0,-2,2,color='black')
ax2.hlines(0,-2,2,color='black')
ax2.set_xlabel('X/C')
ax2.set_ylabel('Y/C')

ax3.plot(np.real(vm1.circle_df.airfoil_pts), vm1.circle_df.airfoil_vels/vm1.v_inf)
ax3.set_aspect(.5)
ax3.set_xlim(0,1.1)
ax3.set_ylim(0,2.2)
ax3.grid()
ax3.vlines(0,-2,2,color='black')
ax3.hlines(0,-2,2,color='black')
ax3.set_xlabel('X/C')
ax3.set_ylabel('V/Vinf')

# Bind Event Functions
fig1.canvas.mpl_connect('button_press_event', vmPlot.button_press_callback)
fig1.canvas.mpl_connect('button_release_event', vmPlot.button_release_callback)
fig1.canvas.mpl_connect('motion_notify_event', vmPlot.motion_detect_callback)

plt.show()

### Outstanding Items ###
# Ensure correct number of coefficients used for singularities. Also, mapping functions only go up 1 3 coefficients
# Check results against examples
# Put high level code in a main file and pull in both classes
# Least squares?
# Make X and Y Limits Scale Automatically
# Start with joukowski and have ability to Add/Remove points
# Add Import/Export Options when you run the program
# Display exact singularity/center values. Maybe have sliders/manual entry box
# Compute lift from circulation
# Can I have an option to compute drag with boundary layer equations?
# Run online with PyScript or make an compile an executable


