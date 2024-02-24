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

        # Generate Initial Von Mises Transform Coefficients
        self.coefficients = [self.compute_vonmises_coeff(c) for c in range(1,len(self.singularities))]

        self.conformalMap()



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
        # Updates circle parameters for center location changes
        self.radius = abs(self.circle_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)
        self.circle_LE = self.radius * np.exp(1j * (np.pi + self.beta)) + self.center

        self.circle_df["angles"] = np.linspace(-1*self.beta, -1*self.beta + 2*np.pi, self.num_pts)
        self.circle_df["circle_pts"] = self.radius * np.exp(self.circle_df.angles*1j)+self.center

    def updateCoefficients(self):
        # Updates Von Mises Coefficients for singularity location changes
        self.coefficients = [self.compute_vonmises_coeff(c) for c in range(1,len(self.singularities))]

        self.C1 = self.compute_vonmises_coeff(1)
        self.C2 = self.compute_vonmises_coeff(2)
        self.C3 = self.compute_vonmises_coeff(3)

    def setV_inf(self, v_inf):
        # Set freestream velocity
        self.v_inf = v_inf

    def setAlpha_rad(self, alpha_rad):
        # Set Angle of Attach
        self.alpha_rad = alpha_rad

    def vmTransform(self,value,deriv=0):
        # Exaluates Von Mises Transform and its derivatives
        if deriv == 0:
            trans_value = value
            for term_idx in range(len(self.coefficients)):
                trans_value = trans_value + self.coefficients[term_idx] / value**(term_idx+1)
        elif deriv == 1:
            trans_value = 1
            for term_idx in range(len(self.coefficients)):
                trans_value = trans_value - (term_idx+1) * self.coefficients[term_idx] / value**(term_idx+2)
        elif deriv == 2:
            trans_value = 0
            for term_idx in range(len(self.coefficients)):
                trans_value = trans_value + (term_idx**2 + 3*term_idx + 2) * self.coefficients[term_idx] / value**(term_idx+3)

        return trans_value

    def conformalMap(self):
        # Map Airfoil Geometry
        self.circle_df["mapped_pts"] = self.vmTransform(self.circle_df.circle_pts.copy())

        # Scale Airfoil Geometry to X/C
        self.airfoil_TE = self.vmTransform(self.circle_TE)
        self.airfoil_LE = self.vmTransform(self.circle_LE)
        self.chord = abs(self.airfoil_TE - self.airfoil_LE)
        self.circle_df["airfoil_pts"] = (np.real(self.circle_df.mapped_pts) - np.real(self.airfoil_LE)) / self.chord + (np.imag(self.circle_df.mapped_pts) / self.chord)*1j
        
        # Map Airfoil Surface Velocities 
        self.circle_df["circle_vels"] = 2*self.v_inf*(np.sin(self.alpha_rad + self.beta) - np.sin(self.alpha_rad - self.circle_df.angles))
        self.circle_df["map_func_derivatives"] = self.vmTransform(self.circle_df.circle_pts.copy(), deriv=1)# first derivative of airfoil pts wrt circle pts
        self.circle_df["airfoil_vels"] = abs(self.circle_df.circle_vels) / abs(self.circle_df.map_func_derivatives)

        # Manually Set TE_Vel with L'hopitals rule because indeterminate
        TE_map_func_second_deriv = self.vmTransform(self.circle_TE, deriv=2) # second derivative of airfoil pts wrt circle pts evaluated at trailing edge
        TE_vel = abs(2*self.v_inf*np.cos(self.alpha_rad+self.beta)/self.radius) / abs(TE_map_func_second_deriv)  # velocity at trailing edge
        self.circle_df.at[0, 'airfoil_vels'] = TE_vel
        self.circle_df.at[self.num_pts-1, 'airfoil_vels'] = TE_vel

        # Calculate Pressure Coefficients
        self.circle_df["airfoil_Cp"] = 1 - (self.circle_df.airfoil_vels / self.v_inf)**2
        
        # Calculate Lift Coefficient
        self.Cl = 8 * np.pi * (self.radius / self.chord) * np.sin(self.alpha_rad + self.beta)


    def plotMapping(self):
        # Plots a still image of the Von Mises Transform
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

    def saveFile(self):
        # Saves Circle Data to file
        file_name = input("Enter File Name: ")
        f = open(file_name + ".txt", "w")
        f.write(f"Von Mises Airfoil Data\n\nCenter\n{self.center}\n\nPoles\n")
        [f.write(f"{self.singularities[idx]}\n") for idx in range(len(self.singularities))]

    def loadFile(self):
        # Load Circle Data from file
        f=None

        file_name = input("Enter File Name: ")

        try:
            f = open(file_name, "r")
        except:
            print("Cannot Find File")

        getCenter = False
        getPoles = False
        center = []
        poles = []
        
        if f:
            try:

                for line in f:
                    if line.startswith('Center'):
                        getCenter = True
                        getPoles = False
                    elif line.startswith('Poles'):
                        getCenter = False
                        getPoles = True

                    if getCenter:
                        center.append(line.strip())

                    if getPoles:
                        poles.append(line.strip())

                center = complex(center[1])
                poles = [complex(poles[idx]) for idx in range(1,len(poles)) if poles[idx]]

                if center and poles:
                    self.center = center
                    self.singularities = poles
                    self.updateCircle()
                    self.updateCoefficients()
                else :
                    print('File Missing Data')

            except:
                print('Invalid File Format')


        
    def exportDAT(self):
        # Saves Circle Data to file
        airfoil_name = input("Enter Airfoil Name: ")
        file_name = input("Enter DAT File Name: ")
        f = open(file_name + ".dat", "w")
        f.write(f"{airfoil_name}\n")
        [f.write(f"  {np.real(self.circle_df.airfoil_pts[idx]):.6f}  {np.imag(self.circle_df.airfoil_pts[idx]):.6f}\n") for idx in range(len(self.circle_df.airfoil_pts))]