import matplotlib.pyplot as plt
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, Range1d, WheelZoomTool
import numpy as np
import math as m
import luch_graphics as lg
import luch_basics as lb
import luch_elements as le


class tube:
    def __init__(self, **kwargs):
        if (len(kwargs.get('r'))!= kwargs.get('V')):
            print('radii or voltages specified incorrectly')
        print(len(kwargs.get('r')))
        print(len(kwargs.get('V')))
        self.r = kwargs.get('r')  # length in meters
        self.V = kwargs.get('V')
        self.steps = 20
        self.gap = kwargs.get('gap')
        self.edges=2 # how many gap sizes as drift spaces on each side
        self.z=np.linspace(0,(len(self.r)-1+2*self.edges)*self.gap, (len(self.r)-1+2*self.edges)*self.steps)
        self.z_step = self.gap/self.steps

        #create array with gaps
        self.gaps=np.concatenate(([self.gap * self.edges],np.ones(len(self.r) - 1) * self.gap, [self.gap * self.edges]))

        # provide fake drifts on both eds to account for voltage tails
        self.r=np.concatenate(([0.000001],self.r, [0.000001]))
        self.V=np.concatenate (([self.V[0]],self.V,[self.V[-1]]))
        self.potential=np.zeros(len(self.z))
        self.voltagetable=np.zeros((len(self.z),len(self.r)-1))

        for i in range(len(self.r)-1):
            """
            phi=1/pi 1/g (A-B)
            A=2*(z+g/2)arctan(z+g/2)/R1-2R1
            B=2*(z-g/2)arctan(z-g/2)/R2-2R2
            """
            A=2*(self.z-np.sum(self.gaps[:i])+self.gaps[i]/2)*np.arctan((self.z-np.sum(self.gaps[:i])+self.gaps[i]/2)/self.r[i])-2*self.r[i]
            B=2*(self.z-np.sum(self.gaps[:i])-self.gaps[i]/2)*np.arctan((self.z-np.sum(self.gaps[:i])-self.gaps[i]/2)/self.r[i+1])-2*self.r[i+1]
            phi=1/(m.pi*self.gaps[i])*(A-B)
            voltage=(self.V[i]+self.V[i+1])/2+(self.V[i+1]-self.V[i])/2*phi
            self.voltagetable[:,i]=voltage-self.V[i]

            self.potential=self.potential+voltage-self.V[i]

    def trace(self, particle):
        self.V_eff=self.potential-(particle.energy/1.6E-16)

        for index,Vz in enumerate(self.V_eff):
            if index==0:
                V_prev=self.V_eff[0]

            else: V_prev=self.V_eff[index-1]

            if index==len(self.V_eff)-1:
                V_next=self.V_eff[index]


            else: V_next=self.V_eff[index+1]

            if V_prev!=0:
                eta_minus=self.V_eff[index]/V_prev
            else: eta_minus=1

            if self.V_eff[index] != 0:
                eta_plus =V_next/ self.V_eff[index]
            else: eta_plus=1

            # ------ acc field matrix --------------
            if eta_minus<=0:
                eta_minus=1
            if eta_plus<=0:
                eta_plus=1

            R12=2*self.z_step/(1+eta_minus**0.5)
            R34=R12
            if eta_minus!=0:
                R22 = 1 / (eta_minus ** 0.5)
            else: R22=1

            R44=R22
            R66=R22
            R56=self.z_step/particle.gamma()**2

            R_acc = np.matrix([[1, R12, 0, 0, 0, 0],
                               [0, R22, 0, 0, 0, 0],
                               [0, 0, 1, R34, 0, 0],
                               [0, 0, 0, R44, 0, 0],
                               [0, 0, 0, 0, 1, R56],
                               [0, 0, 0, 0, 0, R66]])

            # ----- lens matrix -------------------

            R21=-1*(eta_minus*eta_plus-2*eta_minus+1)/(4*eta_minus*self.z_step)
            R43=R21

            R_len = np.matrix([[1, 0, 0, 0, 0, 0],
                               [R21, 1, 0, 0, 0, 0],
                               [0, 0, 1, 0, 0, 0],
                               [0, 0, R43, 1, 0, 0],
                               [0, 0, 0, 0, 1, 0],
                               [0, 0, 0, 0, 0, 1]])


            outvector1 = np.dot(R_acc, particle.position)
            outvector2=np.dot(R_len, outvector1)
            particle.position = outvector2
            particle.energy=particle.energy+(self.V_eff[index]-V_prev)*1.6E-16

            z_energy_charge = np.array(
                [[particle.trajectory[-1, 6] + self.z_step], [particle.energy / 1.6E-16], [particle.charge]])
            particle.trajectory = np.concatenate(
                (particle.trajectory, np.concatenate((outvector2.T, z_energy_charge.T), axis=1)), axis=0)

        return (particle)


Radius = np.ones(15) * 0.015
for i in range(len(Radius) - 1):
    Radius[i] = Radius[i] + i * 0.001

Voltages = np.zeros(15)
for i in range(len(Voltages)):
    Voltages[i] = -i * 20

Voltages[0] = 0
Voltages[1] = 10
gap = 0.025

acc_tube = le.tube(r=Radius, V=Voltages, gap=gap)



print(len(Radius))
print(len(Voltages))




start=np.array([[0.001],[0.005],[0.001],[0.00],[0],[0]])
proton=lb.Particle(1, 20, 1, start)



print(acc_tube.r)
print(acc_tube.V)
print(acc_tube.gaps)
#print(acc_tube.z)


plot_x=acc_tube.z
plot_y=acc_tube.potential
plot_y1=acc_tube.voltagetable[:,0]
plot_y2=acc_tube.voltagetable[:,1]
plot_y3=acc_tube.voltagetable[:,2]
plot_y4=acc_tube.voltagetable[:,5]

# ------------------------ graphics ---------------------

plot=lg.default_graph()

# generate the figure
plot.line(plot_x, plot_y, line_width=2, muted_alpha=0.2, color='blue',  legend='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y1, line_width=2, muted_alpha=0.2, color='red',  legend='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y2, line_width=2, muted_alpha=0.2, color='green',  legend='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y3, line_width=2, muted_alpha=0.2, color='purple',  legend='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y4, line_width=2, muted_alpha=0.2, color='black',  legend='1st axial')  # add curve as a line to the figure

show(plot)
"""

ray=acc_tube.trace(proton)

plot_x=np.array(ray.trajectory[:,6]).reshape(-1,)
plot_y=np.array(ray.trajectory[:,0]).reshape(-1,)
plot_y1=-1*np.array(ray.trajectory[:,2]).reshape(-1,)

plot=lg.default_graph()

# generate the figure
plot.line(plot_x, plot_y, line_width=2, muted_alpha=0.2, color='blue',  legend='1st axial')  # add curve as a line to the figure
plot.line(plot_x, plot_y1, line_width=2, muted_alpha=0.2, color='red',  legend='1st axial')  # add curve as a line to the figure

show(plot)

"""