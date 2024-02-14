from matplotlib.widgets import Slider, Button
from matplotlib import pyplot as plt
import numpy as np

class InteractivePlot:
    def __init__(self, ax1, N, xmin, x_vals, y_vals, epsilon):
        self.ax1 = ax1
        self.N=N # number of points
        self.xmin = xmin
        self.xmax = xmax 
        self.x_vals = x_vals  # starting x values
        self.y_vals = y_vals  # starting y values
        self.p_idx = None #active point
        self.epsilon = epsilon #max pixel distance

    def update(self):
        'update plots'
        self.ax1.clear()
        self.ax1.scatter (self.x_vals,self.y_vals,color='k',marker='o')
        self.ax1.grid(True)

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        self.p_idx = self.get_ind_under_point(event)  

    def button_release_callback(self,event):
        'whenever a mouse button is released'
        
        if event.button != 1:
            return
        self.p_idx = None

    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'

        t = self.ax1.transData.inverted()
        tinv = ax1.transData 
        xy = t.transform([event.x,event.y])
        xr = np.reshape(self.x_vals,(np.shape(self.x_vals)[0],1))
        yr = np.reshape(self.y_vals,(np.shape(self.y_vals)[0],1))
        xy_vals = np.append(xr,yr,1)
        xyt = tinv.transform(xy_vals)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]

        if d[ind] >= self.epsilon:
            ind = None
        
        return ind
    
    def motion_notify_callback(self, event):
        'on mouse movement'

        if self.p_idx is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
            
        #update point loctions
        self.x_vals[self.p_idx] = event.xdata 
        self.y_vals[self.p_idx] = event.ydata 
        
        #update plot
        self.update()
        fig.canvas.draw_idle()

##################################################################################
##################################################################################

N=10 # number of points
xmin = 0 
xmax = 10 
x_vals = np.linspace(xmin,xmax,N)  # starting x values
y_vals = np.linspace(xmin,xmax,N)  # starting y values
p_idx = None #active point
epsilon = 5 #max pixel distance


fig,ax1 = plt.subplots(1,1,figsize=(9.0,8.0))

test = InteractivePlot(ax1, N, xmin, x_vals, y_vals, epsilon)

ax1.scatter (x_vals,y_vals,color='k',marker='o')

ax1.set_yscale('linear')
ax1.set_xlim(0, xmax)
ax1.set_ylim(0,xmax)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.yaxis.grid(True,which='minor',linestyle='--')
ax1.legend(loc=2,prop={'size':22})

# Find Event Functions
fig.canvas.mpl_connect('button_press_event', test.button_press_callback)
fig.canvas.mpl_connect('button_release_event', test.button_release_callback)
fig.canvas.mpl_connect('motion_notify_event', test.motion_notify_callback)

plt.show()