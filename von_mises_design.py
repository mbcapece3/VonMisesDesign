import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
import pandas as pd
import itertools

class VonMises:
    
            
    def __init__(self,center,singularities=[1+0j, -1+0j, 0, 0],num_pts=200):
        self.zeta_TE = 1 + 0j  # fix trailing edge at (1,0)
        self.center = center # complex number
        self.singularities = singularities # list of singularities as complex points
        self.radius = abs(self.zeta_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)
        self.num_pts = num_pts  # number of discreet points in circle 

        circle_data = {
            "angles": np.arange(0, 2*np.pi, 2*np.pi / self.num_pts),
        }
        self.circle_df = pd.DataFrame(circle_data)

        self.circle_df["circle_pts"] = self.radius * np.exp(self.circle_df.angles*1j)+self.center

        s1 = self.singularities[0]
        s2 = self.singularities[1]
        s3 = self.singularities[2]
        s4 = self.singularities[3]

        self.C1 = self.compute_vonmises_coeff(1)
        self.C2 = self.compute_vonmises_coeff(2)
        self.C3 = self.compute_vonmises_coeff(3)

        # self.C1 = s1**2 - s2*(s3+s4) - s3*s4
        # self.C2 = 1/2*(s1*(s2*(s3+s4) + s3*s4) + s2*s3*s4)
        # self.C3 = -1/3*(s1*s2*s3*s4)

    def compute_vonmises_coeff(self, coeff_num):
            num_terms = coeff_num + 1
            k = len(self.singularities)
            sum = 0
        
            for indices in itertools.combinations(range(k), num_terms):
                product = 1
                for idx in indices:
                    product *= self.singularities[idx]
                sum += product

            if (coeff_num % 2) == 0:
                return 1 / coeff_num * sum
            else:
                return -1 / coeff_num * sum 

    def setCenter(self,new_center):
        self.center = new_center # complex number
        self.radius = abs(self.zeta_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)

        self.circle_df["circle_pts"] = self.radius * np.exp(self.circle_df.angles*1j)+self.center


    def setSingularities(self, newSingularities):
        self.singularities = newSingularities # list of singularities as complex points

        s1 = self.singularities[0]
        s2 = self.singularities[1]
        s3 = self.singularities[2]
        s4 = self.singularities[3]

        self.C1 = VonMises.compute_vonmises_coeff(test.singularities,1)
        self.C2 = VonMises.compute_vonmises_coeff(test.singularities,2)
        self.C3 = VonMises.compute_vonmises_coeff(test.singularities,3)
        # self.C1 = s1**2 - s2*(s3+s4) - s3*s4
        # self.C2 = 1/2*(s1*(s2*(s3+s4) + s3*s4) + s2*s3*s4)
        # self.C3 = -1/3*(s1*s2*s3*s4)

    def conformalMap(self, v_inf, alpha_rad):
        self.circle_df["airfoil_pts"] = self.circle_df.circle_pts + self.C1/self.circle_df.circle_pts + self.C2/(self.circle_df.circle_pts**2) + self.C3/(self.circle_df.circle_pts**3)

        # SCALE AIRFOIL PTS TO X/C
        #####################################################################
        #####################################################################
        #####################################################################
        #####################################################################
        #####################################################################
        #####################################################################
        
        self.circle_df["circle_vels"] = 2*v_inf*(np.sin(alpha_rad + self.beta) - np.sin(alpha_rad - self.circle_df.angles))
        self.circle_df["map_func_derivatives"] = 1 - self.C1/(self.circle_df.circle_pts**2) - self.C2/(self.circle_df.circle_pts**3) - self.C3/(self.circle_df.circle_pts**4)  # first derivative of airfoil pts wrt circle pts
        self.circle_df["airfoil_vels"] = abs(self.circle_df.circle_vels) / abs(self.circle_df.map_func_derivatives)

        
        TE_map_func_second_deriv = 2*self.C1/(self.zeta_TE**3) + 6*self.C2/(self.zeta_TE**4) + 12*self.C3/(self.zeta_TE**5)  # second derivative of airfoil pts wrt circle pts evaluated at trailing edge
        TE_vel = v_inf*(2*abs(np.cos(alpha_rad+self.beta))*np.cos(self.beta)) / ((1+np.real(self.center)) * TE_map_func_second_deriv)  # velocity at trailing edge

        #Need to find 2 locations of TE index so that they can be replaced with TE_vel
        ############################################################################
        ############################################################################
        ############################################################################
        ############################################################################
        ############################################################################
        ############################################################################
        ############################################################################
        ############################################################################
        ############################################################################


    def plotCircle(self):
        fig, ax = plt.subplots()
        #ax.plot(self.circle_df.x_pts, self.circle_df.y_pts)
        ax.plot(np.real(self.circle_df.circle_pts), np.imag(self.circle_df.circle_pts))
        ax.scatter(np.real(self.center), np.imag(self.center))
        ax.scatter(np.real(self.singularities),np.imag(self.singularities))
        ax.set_aspect('equal')
        ax.set_xlim(-2,2)
        ax.set_ylim(-2,2)
        ax.grid()
        ax.vlines(0,-2,2,color='black')
        ax.hlines(0,-2,2,color='black')

        fig, ax2 = plt.subplots()
        #ax.plot(self.circle_df.x_pts, self.circle_df.y_pts)
        ax2.plot(np.real(self.circle_df.airfoil_pts), np.imag(self.circle_df.airfoil_pts))
        #ax2.scatter(np.real(self.center), np.imag(self.center))
        #ax2.scatter(np.real(self.singularities),np.imag(self.singularities))
        ax2.set_aspect('equal')
        ax2.set_xlim(-2,2)
        ax2.set_ylim(-2,2)
        ax2.grid()
        ax2.vlines(0,-2,2,color='black')
        ax2.hlines(0,-2,2,color='black')

        
singularities = [1+0j, -1+0j, 0, 0]
singularities2 = [1+0j, -1+0j, .37-.471j, -.37+.471j]
 
test = VonMises(-.072+.29j,singularities2)
test.conformalMap(1,0)
#print(test.circle_df)
test.plotCircle()
#test.setCenter(-.15+.5j)
#test.plotCircle()
#test.setSingularities(singularities2)
#test.conformalMap(1,0)
#test.plotCircle()

plt.show()



