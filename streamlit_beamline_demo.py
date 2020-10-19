import numpy as np
import luch_basics as lb
import luch_elements as le
import luch_graphics as lg
from PIL import Image  # to import image into streamlit
import streamlit as st  # to create web-based GUI

# --------------------- particle definition---------------------------------
start1 = np.array([[0.000], [0.03], [0.000], [0.03], [0], [0]])  # create start position for 1-st particle
start2 = np.array([[0.003], [0.00], [0.003], [0.00], [0], [0]])  # create start position for 2-nd particle

proton1 = lb.Particle(1, 30, 1, start1)  # create first particle
proton2 = lb.Particle(1, 30, 1, start2)  # create second particle

# ------------------------- beamline definition---------------------------

# ------------------------- define drifts---------------------------

drift_length_1 = st.sidebar.slider('drift before magnet [m]', min_value=0.0, max_value=3.0, value=0.9, step=0.05,
                                   format=None, key=None)
drift_length_2 = st.sidebar.slider('drift before accelerator [m]', min_value=0.0, max_value=3.0, value=1.1, step=0.05,
                                   format=None, key=None)
drift_length_3 = st.sidebar.slider('drift before quad [m]', min_value=0.0, max_value=1.0, value=0.35, step=0.05,
                                   format=None, key=None)
drift_length_4 = 4.3 # drift after last optical element

drift_before_magnet = le.driftspace(length=drift_length_1)
drift_before_accelerator = le.driftspace(length=drift_length_2)
drift_before_quad = le.driftspace(length=drift_length_3)
drift_before_end_chamber = le.driftspace(length=drift_length_4)

#  -------------------------create dipole magnet-----------------------------
dipole = le.dipole(radius=0.33, angle=90, gap=0.055, index=0, pfa1=30.615, pfa2=30.615, K=0.6)


#  -------------------------generate accelerator tube-------------------

# create variation of electrode opening
dr = st.sidebar.slider('accelerator electrode radius step [mm]', min_value=-1.0, max_value=2.0, value=-0.5, step=0.1,
                       format=None, key=None)
Radius = np.ones(40) * 0.030  # create a numpy array of equal 30 mm openings
Radius = np.array([Radius[i] + i * dr / 1000 for i in range(len(Radius))])

tube_lens_voltage = st.sidebar.slider('tube lens voltage [kV]', min_value=-60.0, max_value=0.0, value=-56.0, step=0.5,
                                      format=None, key=None)
terminal_voltage = -1 * st.sidebar.slider('terminal voltage [kV]', min_value=0, max_value=460, value=460, step=1,
                                          format=None, key=None)

Voltages_tube = np.linspace(0, terminal_voltage, 38)
Voltages_lens = np.array([0, tube_lens_voltage])
Voltages = np.concatenate((Voltages_lens, Voltages_tube), axis=None)

acc_electrode_gap = 0.025  # gap between electrodes in the accelerator tube

acc_tube = le.tube(r=Radius, V=Voltages, gap=acc_electrode_gap)  # generate accelerator tube


#  -------------------------generate quad triplet-------------------

quad_voltage = st.sidebar.slider('Central quad voltage [V]', min_value=0, max_value=10000, value=5275, step=25,
                                 format=None, key=None)
quad_astigmatism = st.sidebar.slider('quad astigmatism [%]', min_value=-50.0, max_value=50.0, value=-9.6, step=0.2,
                                     format=None, key=None)
quad1 = le.quad(r=0.025, V=quad_voltage * (1 + quad_astigmatism / 100), L=0.102, dir=1)
quad2 = le.quad(r=0.025, V=quad_voltage, L=0.175, dir=0)
quad3 = le.quad(r=0.025, V=quad_voltage * (1 + quad_astigmatism / 100), L=0.102, dir=1)

# -------------------------beamline assembly-------------------

beamline_list = [drift_before_magnet, dipole, drift_before_accelerator, acc_tube, drift_before_quad, quad1, quad2,
                 quad3, drift_before_end_chamber]

beamline = lb.BeamLine(beamline_list)

# ------------------------tracing--------------------------

ray1 = beamline.trace(proton1)
ray2 = beamline.trace(proton2)

# ------------------------plotting--------------------------
plot = lg.plot_trajectory([ray1, ray2])
st.bokeh_chart(plot)
image = Image.open('images/optics.png')
st.image(image, caption='Example optics to reproduce', use_column_width=True)
