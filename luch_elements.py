import numpy as np
import luch_basics as lb
import math as m

class driftspace:
    def __init__(self,  steps=100, **kwargs):
        self.length = kwargs.get('length') # length in meters


        self.steps=steps
        dz=self.length/self.steps
        drift_list=[]
        for i in range(steps):
            drift_list=drift_list+[lb.Drift(dz)]
        self.list=drift_list

    def trace(self, particle):
        beamline=lb.BeamLine(self.list)
        particle=beamline.trace(particle)
        return particle



class einzellens: # three tube einzel lens see George H. Gillespiea and Thomas A. Brownb [1]
    def __init__(self, steps=25, **kwargs):
        self.R = kwargs.get('radius') # length in meters
        self.steps=steps
        self.df=kwargs.get('length1')
        self.gap=kwargs.get('gap')
        self.a=kwargs.get('length2')/2
        self.V1=kwargs.get('V1')
        self.V2=kwargs.get('V2')
        omega=1.31835 # fitinng constant from G&B
        omegaprime=1.67 # fitting constant from G&B
        self.halflength=self.df+self.gap+self.a
        self.z_step=self.halflength/self.steps
        half_z=np.linspace(0, self.halflength, self.steps)
        A=(np.cosh(2*omega*half_z/self.R)+np.cosh( (omega*self.a/self.R)+(omegaprime*self.gap/self.R) ))
        B= (np.cosh(2*omega*half_z/self.R)+np.cosh( (omega*self.a/self.R)-(omegaprime*self.gap/self.R) ))
        self.phi=self.R/(omegaprime*self.gap)*np.log(A/B)
        V_right=self.V1+((self.V2-self.V1)/2)*self.phi # right half field

        V_left=np.flipud(V_right) # flipping the array to use even nature of phi function
        V_left=V_left[:-1]        # drop central element
        self.V=np.concatenate((V_left, V_right))  # concatenate arrays to the potential function in the entire lens

    def trace(self, particle):
        V_right = self.V1-(particle.energy/1.6E-16) + ((self.V2 - self.V1) / 2) * self.phi  # right half field

        V_left = np.flipud(V_right)  # flipping the array to use even nature of phi function
        V_left = V_left[:-1]  # drop central element
        self.V = np.concatenate((V_left, V_right))  # concatenate arrays to the potential function in the entire lens

        for index,Vz in enumerate(self.V):
            if index==0:
                V_prev=self.V[0]

            else: V_prev=self.V[index-1]

            if index==len(self.V)-1:
                V_next=self.V[index]


            else: V_next=self.V[index+1]

            if V_prev!=0:
                eta_minus=self.V[index]/V_prev
            else: eta_minus=1

            if self.V[index] != 0:
                eta_plus =V_next/ self.V[index]
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
            particle.energy=particle.energy+(-self.V[index]+V_prev)*1.6E-16

            z_energy_charge = np.array(
                [[particle.trajectory[-1, 6] + self.z_step], [particle.energy / 1.6E-16], [particle.charge]])
            particle.trajectory = np.concatenate(
                (particle.trajectory, np.concatenate((outvector2.T, z_energy_charge.T), axis=1)), axis=0)

        return (particle)

class dipole:
    def __init__(self, **kwargs):
        self.angle=kwargs.get('angle')
        self.n=kwargs.get('index')
        self.r = kwargs.get('radius')
        self.pfa1=kwargs.get('pfa1')*2*m.pi/360
        self.pfa2=kwargs.get('pfa2')*2*m.pi/360
        self.gap=kwargs.get('gap')
        self.K=kwargs.get('K')



    def trace(self, particle):
        steps = 100
        da=self.angle/360*2*m.pi/steps
        ds=self.r*da
        h=self.angle/(abs(self.r*abs(self.angle)))
        kx=((1-self.n)*h**2)**0.5
        ky=(self.n*h**2)**0.5
        Cx=m.cos(kx*ds)
        Sx=m.sin(kx*ds)
        Cy=m.cos(ky*ds)
        Sy = m.sin(ky * ds)
        if self.n!=0:
            R34=Sy/ky
        else:
            R34=0
        if kx!=0:
            R16=h*(1-Cx)/(kx**2)
            R26=h*Sx/kx
            R51=-1*Sx/kx
            R52=-1*(1-Cx)/(kx**2)
        else:
            R16=0
            R26=0
            R51=0
            R52=0

        phi1=self.K*self.gap*(1+(m.sin(self.pfa1)**2))/(self.r*m.cos(self.pfa1))
        phi2 = self.K * self.gap * (1 + (m.sin(self.pfa2) ** 2)) / (self.r * m.cos(self.pfa2))

        R_pfa1=np.matrix([[1, 0, 0, 0, 0, 0],
                           [m.tan(self.pfa1)/self.r, 1, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, -1*m.tan(self.pfa1-phi1)/self.r, 1, 0, 0],
                           [0, 0, 0, 0, 1, 0],
                           [0, 0, 0, 0, 0, 1]])
        R_pfa2 = np.matrix([[1, 0, 0, 0, 0, 0],
                            [m.tan(self.pfa2) / self.r, 1, 0, 0, 0, 0],
                            [0, 0, 1, 0, 0, 0],
                            [0, 0, -1 * m.tan(self.pfa2 - phi2) / self.r, 1, 0, 0],
                            [0, 0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 0, 1]])

        R_dip = np.matrix([[Cx, Sx/kx, 0, 0, 0, R16],
                           [-1*kx*Sx, Cx, 0, 0, 0, R26],
                           [0, 0, Cy, R34, 0, 0],
                           [0, 0, -1*ky*Sy, Cy, 0, 0],
                           [R51, R52, 0, 0, 1, 0],
                           [0, 0, 0, 0, 0, 1]])

        #------------------ front pfa lens----------------------------------------------------
        outvector = np.dot(R_pfa1, particle.position)
        particle.position = outvector
        z_energy_charge = np.array(
            [[particle.trajectory[-1, 6] ], [particle.energy / 1.6E-16], [particle.charge]])
        particle.trajectory = np.concatenate(
            (particle.trajectory, np.concatenate((outvector.T, z_energy_charge.T), axis=1)), axis=0)

        #--------------- main body ---------------------------------------------------------



        for i in range(steps):
            outvector = np.dot(R_dip, particle.position)

            particle.position = outvector

            z_energy_charge = np.array(
                [[particle.trajectory[-1, 6] + ds], [particle.energy / 1.6E-16], [particle.charge]])
            particle.trajectory = np.concatenate(
                (particle.trajectory, np.concatenate((outvector.T, z_energy_charge.T), axis=1)), axis=0)

        # ------------------ exit pfa lens----------------------------------------------------

        outvector = np.dot(R_pfa2, particle.position)
        particle.position = outvector
        z_energy_charge = np.array(
                [[particle.trajectory[-1, 6]], [particle.energy / 1.6E-16], [particle.charge]])
        particle.trajectory = np.concatenate(
                (particle.trajectory, np.concatenate((outvector.T, z_energy_charge.T), axis=1)), axis=0)



        return (particle)

"""
class two_app_lens:
    def __init__(self, **kwargs):
        self.r1 = kwargs.get('r1')  # length in meters
        self.r2 = kwargs.get('r2')
        self.steps = 40
        self.gap = kwargs.get('gap')
        self.V1 = kwargs.get('V1')
        self.V2 = kwargs.get('V2')
        z = np.linspace(-1*(self.gap/2+3*self.r1), self.gap/2+3*self.r2, self.steps)
        self.z=z

        A = ( (2*z+self.gap)*(np.arctan((z+self.gap/2)/self.r1))-2*self.r1  ) # generate A
        B = np.ones(len(z))*(2-2*self.r2)                                   # iitiate B with its value for the central point
        self.A=A

        B[z!=0] = ( (2*z-self.gap)*(np.arctan((z-self.gap/2)/self.r2))-2*self.r2   ) # conditioally update B values in all but central point
        self.B=B

        self.phi = (1/m.pi)*(1/self.gap)*(A-B)
        self.V=(self.V1+self.V2)/2+((self.V2-self.V1)/2)*self.phi

        print(self.phi)
        print (self.V)

    def trace(self, particle):
        ds=self.gap/self.steps

"""

class tube:
    def __init__(self, **kwargs):
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
        #self.V_eff=self.potential
        self.V_eff = self.potential - (particle.energy / 1.6E-16)


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
            particle.energy=particle.energy+(-self.V_eff[index]+V_prev)*1.6E-16

            z_energy_charge = np.array(
                [[particle.trajectory[-1, 6] + self.z_step], [particle.energy / 1.6E-16], [particle.charge]])
            particle.trajectory = np.concatenate(
                (particle.trajectory, np.concatenate((outvector2.T, z_energy_charge.T), axis=1)), axis=0)

        return (particle)

    def __str__(self):  # string representation of the object method for human-readable output
            selfdescription = 'An  accelerator tube with electrode holes ' + str(self.r) + ' and potential distribution of ' + str(
                self.V)
            return selfdescription


class quad:
    def __init__(self, **kwargs):
        self.r = kwargs.get('r')
        self.V = kwargs.get('V')
        self.L = kwargs.get('L')
        self.dir=kwargs.get('dir')
        self.steps = 100
        self.entrance=1 # how may radii for piecewise linear field ramp up
        self.length = self.L + self.r*self.entrance*2 # length of the quad with two fringe areas
        self.z_step = self.length / self.steps
        self.z = np.linspace(0, self.length, self.steps)
        condlist=[self.z<=self.r*self.entrance, (self.z>self.r*self.entrance)&(self.z<=self.r*self.entrance+self.L),self.z>=self.r*self.entrance+self.L ]
        funclist=[lambda x: x/self.r*self.entrance,lambda x: 1, lambda x: 1-(x-self.r*self.entrance-self.L)/self.r*self.entrance]
        self.phi = np.piecewise(self.z, condlist, funclist)

    def trace(self, particle):
        half_step_drift=lb.Drift(self.z_step / 2)
        m = particle.mass*1.672E-27 # get particle rest mass
        a = self.r # get quad aperture radius
        q=particle.charge*1.6E-19
        L=self.z_step
        c=299792458
        for phi in self.phi:
            particle=half_step_drift.trace(particle) # let particle drift over half of step using standard drift
            V0=phi*self.V # define local voltage based on full voltage and form function phi

            gamma = particle.gamma()  # calculate particle gamma on a specific step
            beta = particle.beta()  # calculate particle beta factor on the specific step
            R56 = 0  # matrix element same as in drift
            if V0!=0:

                k = abs( 2 * V0 * q / ((a ** 2) * m * gamma * (beta ** 2) * (c ** 2)))

                # calculate focusing submatrix
                F11 = np.cos(k ** 0.5 * L)
                F12= np.sin(k ** 0.5 * L)/k ** 0.5
                F21=-1*k ** 0.5*np.sin(k ** 0.5 * L)
                # calculate defocusing submatrix
                D11 = np.cosh(k ** 0.5 * L)
                D12 = np.sinh(k ** 0.5 * L) / k ** 0.5
                D21 = 1 * k ** 0.5 * np.sinh(k ** 0.5 * L)
                # assign focusing and defoccusing sub matrixes depending on quad polarity
                if self.dir==0:
                    R11,R22=F11, F11
                    R12=F12
                    R21=F21
                    R33,R44=D11, D11
                    R34=D12
                    R43=D21
                if self.dir==1:
                    R11,R22=D11, D11
                    R12=D12
                    R21=D21
                    R33,R44=F11, F11
                    R34=F12
                    R43=F21


            if V0==0:
                R11,R22, R33,R44 =1,1,1,1
                R12,R21,R34,R43=0,0,0,0


            R_len = np.matrix([[R11, R12, 0, 0, 0, 0],
                               [R21, R22, 0, 0, 0, 0],
                               [0, 0, R33, R34, 0, 0],
                               [0, 0, R43, R44, 0, 0],
                               [0, 0, 0, 0, 1, R56],
                               [0, 0, 0, 0, 0, 1]])


            outvector1 = np.dot(R_len, particle.position)

            particle.position = outvector1
            #particle.energy = particle.energy + (-self.V[index] + V_prev) * 1.6E-16

            z_energy_charge = np.array(
                [[particle.trajectory[-1, 6] ], [particle.energy / 1.6E-16], [particle.charge]])
            particle.trajectory = np.concatenate(
                (particle.trajectory, np.concatenate((outvector1.T, z_energy_charge.T), axis=1)), axis=0)

            particle = half_step_drift.trace(particle)  # let particle drift over half of step using standard drift

        return (particle)

    def __str__(self):  # string representation of the object method for human-readable output
        selfdescription = 'A quadrupole of '+str(self.L)+' meters' + ' Bore diamter of '+ str(self.r) +' meters'
        return selfdescription


class esa:
    def __init__(self, **kwargs):
        self.type=kwargs.get('type')




            




