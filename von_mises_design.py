import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
import pandas as pd

zeta_TE = 1 + 0j  # fix trailing edge at (1,0)


class Circle:
    def __init__(self,center,num_pts=1000):
        self.center = center # complex number
        self.radius = abs(zeta_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)
        self.num_pts = num_pts  # number of discreet points in circle
        #print(self.beta*180/np.pi)
        
    def setCenter(self,new_center):
        self.center = new_center # complex number
        self.radius = abs(zeta_TE - self.center)
        self.beta = np.arcsin(np.imag(self.center)/self.radius)

    def plotCircle(self):
        angles = np.arange(0, 2*np.pi, 2*np.pi / self.num_pts)
        circle_x_vals = self.radius * np.cos(angles) + np.real(self.center)
        circle_y_vals = self.radius * np.sin(angles) + np.imag(self.center)
        circle_data = {
            "angle": angles,
            "x_pt": circle_x_vals,
            "y_pt": circle_y_vals
        }
        circle_df = pd.DataFrame(circle_data)
        
        fig, ax = plt.subplots()
        ax.plot(circle_df.x_pt, circle_df.y_pt)
        ax.scatter(np.real(self.center), np.imag(self.center))
        ax.set_aspect('equal')
        ax.set_xlim(-2,2)
        ax.set_ylim(-2,2)
        ax.grid()

        plt.show()
        

# class Airfoil:
#     def __init__(self, x_pts, y_pts):
#         self.x_pts = x_pts
#         self.y_pts = y_pts

#     def drawAirfoil():
#         print('test')
        
test = Circle(-.072+.29j)
test.plotCircle()