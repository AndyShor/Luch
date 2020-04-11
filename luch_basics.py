"""
This file contains basic definition for matrix calculation of beam envelope
"""
# pylint: disable=too-few-public-methods
import numpy as np


class Particle:
    """
    fundamental particle class keeps information of particle present state and its history
    """

    def __init__(self, mass, energy, charge, origin):
        self.mass = mass  # mass in AMU
        self.energy = energy * 1.6E-16  # transform from keV
        # self.gamma=1+self.energy/(self.mass*1.660539E-27*(299792458)**2)      # Lorentz factor
        self.charge = charge  # charge in elementary charges
        # trajectory will be initiated at the particle origin coordinates
        self.position = origin
        # column vector for z coordiate along the trajectory, particle energy and charge
        z_energy_charge = np.array([[0], [self.energy / 1.6E-16], [self.charge]])
        self.trajectory = np.concatenate((np.asarray(origin).T,
                                          np.asarray(z_energy_charge).T), axis=1)

    def __str__(self):
        """
        string representation of the object method for human-readable output
        """
        selfdescription = 'A particle of ' + str(self.mass) + ' AMU and energy ' \
                          + str(self.energy / 1.6E-16) + ' keV'

        return selfdescription

    def gamma(self):
        """
        returns relativistic gamma factor of a particle in present state
        """
        gamma = 1 + self.energy / (self.mass * 1.660539E-27 * (299792458) ** 2)
        return gamma

    def beta(self):
        """
        returns relativistic beta factor of a particle in present state
        """
        beta = (1 - 1 / (self.gamma()) ** 2) ** 0.5
        return beta


class Drift:
    """
    fundamental class for an elementary drift
    """
    # pylint: disable=too-few-public-methods

    def __init__(self, length):
        self.length = length  # length in meters

    def __str__(self):
        """"
        string representation of the object method for human-readable output
        """
        selfdescription = 'A drift of ' + str(self.length) + ' meters'

        return selfdescription

    def trace(self, particle):
        """
        mandatory trace method for any beamline-related class
        """
        gamma = particle.gamma()
        r_matrix = np.array([[1, self.length, 0, 0, 0, 0],
                              [0, 1, 0, 0, 0, 0],
                              [0, 0, 1, self.length, 0, 0],
                              [0, 0, 0, 1, 0, 0],
                              [0, 0, 0, 0, 1, self.length / gamma ** 2],
                              [0, 0, 0, 0, 0, 1]])
        outvector = np.dot(r_matrix, particle.position)
        particle.position = outvector
        z_energy_charge = np.array([[particle.trajectory[-1, 6] + self.length],
                                    [particle.energy / 1.6E-16], [particle.charge]])
        particle.trajectory = np.concatenate(
            (particle.trajectory,
             np.concatenate((outvector.T, z_energy_charge.T), axis=1)),
            axis=0)
        return particle


class ThinLens:
    """
    fundamental class for a thin lens
    """
    # pylint: disable=too-few-public-methods

    def __init__(self, fx, fy, fz):
        self.fx = fx
        self.fy = fy
        self.fz = fz

    def trace(self, particle):
        """
        trace method for thin lens accepts particle,
        returns particle with lens affects in the trajectory
        """
        if self.fx == 0:
            rxx = 0
        else:
            rxx = -1 / self.fx
        if self.fy == 0:
            ryy = 0
        else:
            ryy = -1 / self.fy
        if self.fz == 0:
            rzz = 0
        else:
            rzz = -1 / self.fz
        r_matrix = np.array([[1, 0, 0, 0, 0, 0],
                              [rxx, 1, 0, 0, 0, 0],
                              [0, 0, 1, 0, 0, 0],
                              [0, 0, ryy, 1, 0, 0],
                              [0, 0, 0, 0, 1, 0],
                              [0, 0, 0, 0, rzz, 1]])
        outvector = np.dot(r_matrix, particle.position)
        particle.position = outvector
        z_energy_charge = np.array([[particle.trajectory[-1, 6]],
                                    [particle.energy / 1.6E-16],
                                    [particle.charge]])
        particle.trajectory = np.concatenate(
            (particle.trajectory, np.concatenate((outvector.T, z_energy_charge.T), axis=1)), axis=0)
        return particle

    def __str__(self):
        """
        string representation of the object method for human-readable output
        """
        selfdescription = 'A thin lens with x,y,z focal lengths of ' + str(self.fx) \
                          + ', ' + str(self.fy) + ', ' + str(self.fz) + ' meters'
        return selfdescription


class BeamLine:
    """
    special class to trace particle through a list of beamline elements having trace method
    """
    # pylint: disable=too-few-public-methods
    def __init__(self, beamline):
        self.beamline = beamline

    def trace(self, particle):
        """
        trace method for BeamLine object
        loops through elements of BeamLine list and traces particle through each of them
        """
        for element in self.beamline:
            particle = element.trace(particle)
        return particle

    def __str__(self):
        """
        string representation of the beamline object
        """
        selfdescription = 'A beamline containing ' \
                          + str([type(element) for element in self.beamline])
        return selfdescription

class Matrix:
    """
    basic class for tracing through a matrix element that has zero length and conserves energy and charge
    """
    # pylint: disable=too-few-public-methods
    def __init__(self, **kwargs):
        self.matrix = kwargs.get('matrix')
        self.dz = kwargs.get('dz',0)
        self.de = kwargs.get('de', 0)
        self.dq = kwargs.get('dq', 0)

    def trace(self, particle):
        outvector = np.dot(self.matrix, particle.position)
        particle.position = outvector

        particle.energy = particle.energy + self.de

        z_energy_charge = np.array(
            [[particle.trajectory[-1, 6] + self.dz], [particle.energy / 1.6E-16], [particle.charge+self.dq]])
        particle.trajectory = np.concatenate(
            (particle.trajectory, np.concatenate((outvector.T, z_energy_charge.T), axis=1)), axis=0)

        return particle