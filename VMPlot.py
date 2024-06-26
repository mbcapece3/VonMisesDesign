from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import pyplot as plt
import numpy as np

class VMPlot:
    def __init__(self, data_obj, epsilon=5):
        #Initialize Values
        self.data_obj = data_obj # VonMises style object
        self.p_idx = None #active point
        self.epsilon = epsilon #max pixel distance
        self.toggle_state = 'Velocity' # default distribution plot

        # Create Figure
        self.fig,(self.ax1,self.ax2,self.ax3) = plt.subplots(1,3,figsize=(18,5))  # Create figure
        self.fig.subplots_adjust(left=0.13, bottom=0.15) #Adjust to "push" subplots up and to the right
        
        # Compute initial transform
        self.data_obj.conformalMap()

        # Initialize Figure Textboxes
        self.alpha_text = self.fig.text(.9,.80,f'α  = {self.data_obj.alpha_rad * 180 / np.pi:.2f}',fontsize=12) # Display angle of attack
        self.Cl_text = self.fig.text(.9,.85,f'Cl = {self.data_obj.Cl:.2f}',fontsize=12) # Display lift coefficient

        # Initial Plot
        self.update_plots()  

        # Bind interactive actions
        self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.fig.canvas.mpl_connect('motion_notify_event', self.motion_detect_callback)

        # Define Angle of Attack Slider
        aoa = self.fig.add_axes([0.19,0.02,0.65,0.03])
        self.slider_aoa = Slider(aoa, 'Angle of Attack',-45,45, valinit=self.data_obj.alpha_rad * 180 / np.pi, valstep = np.linspace(-45,45,181))
        self.slider_aoa.on_changed(self.update_aoa)

        # Define Add Point Button
        add_pt = plt.axes([0.015,0.12,0.07,0.05])
        self.button_add_pt = Button(add_pt, 'Add Pole', color='.75', hovercolor='.9')
        self.button_add_pt.on_clicked(self.add_singularity)

        # Define Remove Point Button
        remove_pt = plt.axes([0.015,0.20,0.07,0.05])
        self.button_remove_pt = Button(remove_pt, 'Remove Pole', color='.75', hovercolor='.9')
        self.button_remove_pt.on_clicked(self.remove_singularity)

        # Define Distribution Type Radio Button
        radio = plt.axes([0.9, 0.15, 0.06, 0.15])
        self.radio_select = RadioButtons(radio, ('Velocity', 'Pressure'), active=0)
        self.radio_select.on_clicked(self.toggle_radio)

        # Define Save File Button
        save_file = plt.axes([0.015,0.8,0.07,0.05])
        self.button_save_file = Button(save_file, 'Create Save File', color='.75', hovercolor='.9')
        self.button_save_file.on_clicked(self.save_file_callback)

        # Define Load File Button
        load_file = plt.axes([0.015,0.72,0.07,0.05])
        self.button_load_file = Button(load_file, 'Load Save File', color='.75', hovercolor='.9')
        self.button_load_file.on_clicked(self.load_file_callback)

        # Define Export DAT Button
        export_dat = plt.axes([0.015,0.64,0.07,0.05])
        self.button_export_dat = Button(export_dat, 'Export DAT File', color='.75', hovercolor='.9')
        self.button_export_dat.on_clicked(self.export_DAT_callback)

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

    def update_aoa(self, event):
        'update angle of attack'
        self.data_obj.setAlpha_rad(self.slider_aoa.val * np.pi / 180) # get slider value in degrees and input in radians
        self.data_obj.conformalMap() # update conformal mapping
        self.update_plots() 
        self.fig.canvas.draw_idle()

    def add_singularity(self, event):
        'add singularity at 0'
        self.data_obj.singularities = self.data_obj.singularities + [0.0+0.0j]  # Add singularity at 0 as to not affect other values
        self.data_obj.updateCoefficients() # update von mises coefficients
        self.data_obj.conformalMap() # update conformal mapping
        self.update_plots() 
        self.fig.canvas.draw_idle()

    def remove_singularity(self, event):
        'remove last singularity'
        num_sing = len(self.data_obj.singularities)
        if num_sing > 2:
            xvals = np.linspace(1,-1,num_sing-1) # X values linearly spaces so sum is 0. Idx 0 must be trailing edge becuase idx 0 is fixed in place
            self.data_obj.singularities = [(val + 0j) for val in xvals] 
            self.data_obj.updateCoefficients() # update von mises coefficients
            self.data_obj.conformalMap() # update conformal mapping
            self.update_plots() 
            self.fig.canvas.draw_idle()

    def toggle_radio(self, event):
        'switch distribution plot'
        self.toggle_state = event
        self.update_plots() 
        self.fig.canvas.draw_idle()

    def save_file_callback(self, event):
        'save airfoil info to file'
        self.data_obj.saveFile()

    def load_file_callback(self, event):
        'load airfoil info from file'
        self.data_obj.loadFile()
        self.data_obj.conformalMap() # update conformal mapping
        self.update_plots() 
        self.fig.canvas.draw_idle()

    def export_DAT_callback(self, event):
        'save airfoil coordinates to file'
        self.data_obj.exportDAT()

    def update_plots(self):
        'update plots'
        #Clear Previous Plot
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.alpha_text.remove()
        self.Cl_text.remove()

        # Display Textboxes
        self.alpha_text = self.fig.text(.9,.80,f'α  = {self.data_obj.alpha_rad * 180 / np.pi:.2f}',fontsize=12) # Display angle of attack
        self.Cl_text = self.fig.text(.9,.85,f'Cl = {self.data_obj.Cl:.2f}',fontsize=12) # Display lift coefficient

        # Circle Plot
        self.ax1.scatter (np.real(self.data_obj.center), np.imag(self.data_obj.center),color='r',marker='o')
        self.ax1.scatter (np.real(self.data_obj.singularities),np.imag(self.data_obj.singularities),color='k',marker='o')
        self.ax1.plot(np.real(self.data_obj.circle_df.circle_pts), np.imag(self.data_obj.circle_df.circle_pts))
        #self.ax1.plot(np.real([self.data_obj.circle_TE,self.data_obj.circle_LE]),np.imag([self.data_obj.circle_TE,self.data_obj.circle_LE]), linewidth=3, color='red')
        self.ax1.set_aspect('equal')
        self.ax1.set_title('Circle Plane')
        self.ax1.set_xlim(-2,2)
        self.ax1.set_ylim(-2,2)
        self.ax1.grid()
        self.ax1.vlines(0,-2,2,color='black')
        self.ax1.hlines(0,-2,2,color='black')

        # Airfoil Plot
        self.ax2.plot(np.real(self.data_obj.circle_df.airfoil_pts)/np.real(self.data_obj.circle_df.airfoil_pts[0]), np.imag(self.data_obj.circle_df.airfoil_pts)) # Scaled from 0 to 1
        #self.ax2.plot(np.real(self.data_obj.circle_df.airfoil_pts), np.imag(self.data_obj.circle_df.airfoil_pts))
        #self.ax2.plot(np.real(self.data_obj.scaleAirfoilPts([self.data_obj.airfoil_TE,self.data_obj.airfoil_LE])),np.imag(self.data_obj.scaleAirfoilPts([self.data_obj.airfoil_TE,self.data_obj.airfoil_LE])), linewidth=3, color='red')
        self.ax2.set_aspect('equal')
        self.ax2.set_title('Airfoil Plane')
        self.ax2.set_xlim(-.1,1.1)
        self.ax2.set_ylim(-.6,.6)
        self.ax2.grid()
        self.ax2.vlines(0,-2,2,color='black')
        self.ax2.hlines(0,-2,2,color='black')
        self.ax2.set_xlabel('X/C')
        self.ax2.set_ylabel('Y/C')

        # Distribution Plot
        if self.toggle_state == 'Velocity':
            # Velocity Dist Plot
            self.ax3.plot(np.real(self.data_obj.circle_df.airfoil_pts)/np.real(self.data_obj.circle_df.airfoil_pts[0]), self.data_obj.circle_df.airfoil_vels/self.data_obj.v_inf) # Scaled from 0 to 1
            #self.ax3.plot(np.real(self.data_obj.circle_df.airfoil_pts), self.data_obj.circle_df.airfoil_vels/self.data_obj.v_inf)
            self.ax3.set_aspect(.5)
            self.ax3.set_title('Velocity Distribution')
            self.ax3.set_xlim(0,1.1)
            self.ax3.set_ylim(0,2.2)
            self.ax3.grid()
            self.ax3.vlines(0,-2,2,color='black')
            self.ax3.hlines(0,-2,2,color='black')
            self.ax3.set_xlabel('X/C')
            self.ax3.set_ylabel('V/Vinf')

        elif self.toggle_state == 'Pressure':
            # Presure Dist Plot
            self.ax3.plot(np.real(self.data_obj.circle_df.airfoil_pts)/np.real(self.data_obj.circle_df.airfoil_pts[0]), self.data_obj.circle_df.airfoil_Cp) # Scaled from 0 to 1
            #self.ax3.plot(np.real(self.data_obj.circle_df.airfoil_pts), self.data_obj.circle_df.airfoil_Cp)
            self.ax3.set_aspect(.25)
            self.ax3.set_title('Pressure Distribution')
            self.ax3.set_xlim(0,1.1)
            self.ax3.set_ylim(1,-3.4)
            self.ax3.grid()
            self.ax3.vlines(0,-2,2,color='black')
            self.ax3.hlines(0,-2,2,color='black')
            self.ax3.set_xlabel('X/C')
            self.ax3.set_ylabel('Cp')



