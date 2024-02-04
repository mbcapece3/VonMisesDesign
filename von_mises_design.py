import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
import pandas as pd

class Circle:
    def __init__(self,center,singularities=[1+0j, -1+0j],num_pts=200):
        self.zeta_TE = 1 + 0j  # fix trailing edge at (1,0)
        self.center = center # complex number
        self.singularities = singularities # list of singularities as complex points
        self.radius = abs(self.zeta_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)
        self.num_pts = num_pts  # number of discreet points in circle 
        self.angles = np.arange(0, 2*np.pi, 2*np.pi / self.num_pts)

        circle_x_vals = self.radius * np.cos(self.angles) + np.real(self.center)
        circle_y_vals = self.radius * np.sin(self.angles) + np.imag(self.center)

        circle_data = {
            "angle": self.angles,
            "x_pts": circle_x_vals,
            "y_pts": circle_y_vals,
        }

        self.circle_df = pd.DataFrame(circle_data)
        
        
    def setCenter(self,new_center):
        self.center = new_center # complex number
        self.radius = abs(self.zeta_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)

        circle_x_vals = self.radius * np.cos(self.angles) + np.real(self.center)
        circle_y_vals = self.radius * np.sin(self.angles) + np.imag(self.center)

        self.circle_df["x_pts"] = circle_x_vals
        self.circle_df["y_pts"] = circle_y_vals


    def setSingularities(self, newSingularities):
        self.singularities = newSingularities # list of singularities as complex points

    def plotCircle(self):
        fig, ax = plt.subplots()
        ax.plot(self.circle_df.x_pts, self.circle_df.y_pts)
        ax.scatter(np.real(self.center), np.imag(self.center))
        ax.scatter(np.real(self.singularities),np.imag(self.singularities))
        ax.set_aspect('equal')
        ax.set_xlim(-2,2)
        ax.set_ylim(-2,2)
        ax.grid()
        ax.vlines(0,-2,2,color='black')
        ax.hlines(0,-2,2,color='black')

#class ConformalMap:


        

# class Airfoil:
#     def __init__(self, x_pts, y_pts):
#         self.x_pts = x_pts
#         self.y_pts = y_pts

#     def drawAirfoil():
#         print('test')

singularities = [1+0j, -1+0j]
singularities2 = [1+0j, -1+0j, 0+1j, 0-1j]

test = Circle(-.072+.29j,singularities)
test.plotCircle()
test.setCenter(-.15+.5j)
test.plotCircle()
test.setSingularities(singularities2)
test.plotCircle()
plt.show()