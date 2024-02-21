from matplotlib.widgets import Slider, Button
from matplotlib import pyplot as plt
import numpy as np

class VMPlot:
    def __init__(self, data_obj, epsilon=5):
        self.data_obj = data_obj # VonMises style object
        #self.fig = fig
        #self.ax1 = ax1
        #self.ax2 = ax2
        #self.ax3 = ax3
        self.p_idx = None #active point
        self.epsilon = epsilon #max pixel distance

        self.fig,(self.ax1,self.ax2,self.ax3) = plt.subplots(1,3,figsize=(18,5))  # Create figure
        self.data_obj.conformalMap() # Compute initial transform
        self.update_plots()  # Initial Plot

        # Bind interactive actions
        self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.fig.canvas.mpl_connect('motion_notify_event', self.motion_detect_callback)

        plt.show() # show plot

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        self.p_idx = self.get_cursor_idx(event)  

    def button_release_callback(self,event):
        'whenever a mouse button is released'
        
        if event.button != 1:
            return
        self.p_idx = None

    def get_cursor_idx(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        x_vals = np.concatenate((np.array(np.real(self.data_obj.center)), np.real(self.data_obj.singularities)),axis=None) #center is idx 0
        y_vals = np.concatenate((np.array(np.imag(self.data_obj.center)), np.imag(self.data_obj.singularities)),axis=None) #center is idx 0

        tinv = self.ax1.transData
        xr = np.reshape(x_vals,(np.shape(x_vals)[0],1))
        yr = np.reshape(y_vals,(np.shape(y_vals)[0],1))
        xy_vals = np.append(xr,yr,1)
        xyt = tinv.transform(xy_vals) #transform to pixels
        xt, yt = xyt[:, 0], xyt[:, 1] # get initial points x and y pixels
        d = np.hypot(xt - event.x, yt - event.y) # calc distance from cursor to pt
        indseq, = np.nonzero(d == d.min()) # get index of point that is closest to cursor
        idx = indseq[0] # get only index in array
        
        if d[idx] >= self.epsilon:
            idx = None # if cursor is within epsilon tolerance of pt, set index to that pt
        
        return idx
    
    def motion_detect_callback(self, event):
        'on mouse movement'

        if self.p_idx is None:
            return
        if len(self.data_obj.singularities)==2 and self.p_idx!=0: # case of joukowski airfoil
            return
        if self.p_idx == 1:  # Case where trailing edge point is selected
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
            
        #update data
        if self.p_idx == 0:
            self.data_obj.center = event.xdata + event.ydata * 1j  # set center location
            self.data_obj.updateCircle() # update circle parameters
        else:
            self.old_sing_val = self.data_obj.singularities[self.p_idx-1] # store previous singularity value
            self.data_obj.singularities[self.p_idx-1] = event.xdata + event.ydata * 1j # set singularity location

            #Adjust singularity values to ensure all singularities sum to zero
            compensation = -1 * (self.data_obj.singularities[self.p_idx-1] - self.old_sing_val) / (len(self.data_obj.singularities) - 2) #amount each other singularity (excluding TE) moves to keep sum  of pts to 0

            for sing_idx in range(1,len(self.data_obj.singularities)):
                # Excludes idx 1 because trailing edge
                # Excludes self.p_idx-1 because it is the pt that was moved intentionally
                if sing_idx != (self.p_idx-1):
                    self.data_obj.singularities[sing_idx] = self.data_obj.singularities[sing_idx] + compensation
            
            self.data_obj.updateCoefficients() # update von mises coefficients

        #print(np.sum(self.data_obj.singularities))
        self.data_obj.conformalMap() #update conformal mapping
        
        #update plot
        self.update_plots()
        self.fig.canvas.draw_idle()

    def update_plots(self):
        'update plots'

        #Clear Previous Plot
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        
        # Circle Plot
        self.ax1.scatter (np.real(self.data_obj.center), np.imag(self.data_obj.center),color='r',marker='o')
        self.ax1.scatter (np.real(self.data_obj.singularities),np.imag(self.data_obj.singularities),color='k',marker='o')
        self.ax1.plot(np.real(self.data_obj.circle_df.circle_pts), np.imag(self.data_obj.circle_df.circle_pts))
        self.ax1.set_aspect('equal')
        self.ax1.set_xlim(-2,2)
        self.ax1.set_ylim(-2,2)
        self.ax1.grid()
        self.ax1.vlines(0,-2,2,color='black')
        self.ax1.hlines(0,-2,2,color='black')

        # Airfoil Plot
        self.ax2.plot(np.real(self.data_obj.circle_df.airfoil_pts), np.imag(self.data_obj.circle_df.airfoil_pts))
        self.ax2.set_aspect('equal')
        self.ax2.set_xlim(-.1,1.1)
        self.ax2.set_ylim(-.6,.6)
        self.ax2.grid()
        self.ax2.vlines(0,-2,2,color='black')
        self.ax2.hlines(0,-2,2,color='black')
        self.ax2.set_xlabel('X/C')
        self.ax2.set_ylabel('Y/C')

        # Velocity Dist Plot
        self.ax3.plot(np.real(self.data_obj.circle_df.airfoil_pts), self.data_obj.circle_df.airfoil_vels/self.data_obj.v_inf)
        self.ax3.set_aspect(.5)
        self.ax3.set_xlim(0,1.1)
        self.ax3.set_ylim(0,2.2)
        self.ax3.grid()
        self.ax3.vlines(0,-2,2,color='black')
        self.ax3.hlines(0,-2,2,color='black')
        self.ax3.set_xlabel('X/C')
        self.ax3.set_ylabel('V/Vinf')

