import numpy as np
import luch_basics as lb
import luch_elements as le
import luch_graphics as lg

from bokeh.plotting import figure, output_file, show
from bokeh.models import PrintfTickFormatter, HoverTool, Range1d, LabelSet, Label, WheelZoomTool

import streamlit as st

start1=np.array([[0.002],[0.00],[0.002],[0.00],[0],[0]])
start2=np.array([[0.000],[0.02],[0.00],[0.02],[0],[0]])

proton1=lb.Particle(1, 30, 1, start1)
proton2=lb.Particle(1, 30, 1, start2)

#------------------------- beamline definition---------------------------

#drift_length_1=st.sidebar.slider('first drift length [m]', min_value=0.0, max_value=3.0, value=0.9, step=0.05, format=None, key=None)
#drift_length_2=st.sidebar.slider('second drift length [m]', min_value=0.0, max_value=3.0, value=1.2, step=0.1, format=None, key=None)
#drift_length_3=st.sidebar.slider('third drift length [m]', min_value=0.0, max_value=5.0, value=4.3, step=0.1, format=None, key=None)
#drift_length_4=st.sidebar.slider('fourth drift length [m]', min_value=0.0, max_value=1.0, value=0.4, step=0.1, format=None, key=None)
drift_length_1=0.95
drift_length_2=1.2
drift_length_3=4.3
drift_length_4=0.4

gap1=le.driftspace(length=drift_length_1)
gap2=le.driftspace(length=drift_length_2)
gap3=le.driftspace(length=drift_length_3)
gap4=le.driftspace(length=drift_length_4)

dipole1=le.dipole(radius=0.33,angle=90, gap=0.055, index=0, pfa1=30.6, pfa2=30.6, K=0.5)




dr=st.sidebar.slider('radius step [mm]', min_value=-1.0, max_value=2.0, value=0.0, step=0.1, format=None, key=None)
Radius = np.ones(34) * 0.015
for i in range(len(Radius)):
    Radius[i] = Radius[i] + i * dr / 1000

tube_lens_voltage=st.sidebar.slider('tube lens voltage [kV]', min_value=-40.0, max_value=0.0, value=-16.5, step=0.5, format=None, key=None)
terminal_voltage=-1*st.sidebar.slider('tube lens voltage [kV]', min_value=0, max_value=460, value=460, step=1, format=None, key=None)

Voltages=np.linspace(0,terminal_voltage,34)
Voltages[1]=tube_lens_voltage
Voltages[2]=0
#Voltages=[0,tube_lens_voltage,  0,-70,-100,-130,-150,-180,-210, -240]

gap=0.025

acc_tube=le.tube(r=Radius,V=Voltages, gap=gap)

quad_voltage=st.sidebar.slider('Central quad voltage [V]', min_value=0, max_value=20000, value=1000, step=100, format=None, key=None)
quad_atigmatism=st.sidebar.slider('quad astigmatism [%]', min_value=-50.0, max_value=50.0, value=0.0, step=0.2, format=None, key=None)

quad1=le.quad(r=0.021,V=quad_voltage*(1+quad_atigmatism/100), L=0.102, dir=1)
quad2=le.quad(r=0.021,V=quad_voltage, L=0.175, dir=0)
quad3=le.quad(r=0.021,V=quad_voltage*(1+quad_atigmatism/100), L=0.102, dir=1)

#lens4=le.two_app_lens(r1=0.01, r2=0.01, gap=0.0254, V1=0, V2=20)


#-------------------------beamline assembly-------------------

beamline2=[gap1, dipole1, gap2, acc_tube, gap4, quad1,quad2, quad3, gap3]
#beamline2=[gap1, acc_tube, gap3]
beamline=lb.BeamLine(beamline2)

#print(beamline)
#------------------------tracing--------------------------

ray1=beamline.trace(proton1)
ray2=beamline.trace(proton2)

plot_x=np.array(ray1.trajectory[:,6]).reshape(-1,)
plot_y=np.array(ray1.trajectory[:,0]).reshape(-1,)
plot_y1=-1*np.array(ray1.trajectory[:,2]).reshape(-1,)
plot_y2=np.array(ray2.trajectory[:,0]).reshape(-1,)
plot_y3=-1*np.array(ray2.trajectory[:,2]).reshape(-1,)

#plot_x=lens4.z
#plot_y=lens4.A-lens4.B

#plot_y1=lens4.B

# ------------------------ graphics ---------------------

plot=lg.default_graph()

# generate the figure
plot.line(plot_x, plot_y, line_width=2, muted_alpha=0.2, color='blue',  legend_label='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y1, line_width=2, muted_alpha=0.2, color='blue',  legend_label='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y2, line_width=2, muted_alpha=0.2, color='red',  legend_label='2nd axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y3, line_width=2, muted_alpha=0.2, color='red',  legend_label='2nd axial')  # add curve as a line to the figure

st.bokeh_chart(plot)
