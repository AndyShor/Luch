"""
this script contains set of unit tests for Luch
unit test check basic functionality
"""
import pytest
import numpy as np
import luch_basics as lb
import luch_elements as le


def test_drift_space_basic():
    """test creation of drift space using basic class Drift and tracing through it
    it is tested whether central particle with 10 mrad angle
    after 1 meter drift has appropriate coordinates"""
    start = np.array([[0.000], [0.010], [0.000], [0.00], [0], [0]])
    proton = lb.Particle(1, 30, 1, start)
    driftspace = lb.Drift(length=1.0)
    beamline = lb.BeamLine([driftspace])
    ray = beamline.trace(proton)
    endpoint_z=np.array(ray.trajectory[:,6]).reshape(-1,)[-1]
    endpoint_x = np.array(ray.trajectory[:, 0]).reshape(-1, )[-1]
    endpoint_x_p = np.array(ray.trajectory[:, 1]).reshape(-1, )[-1]
    endpoint_y = np.array(ray.trajectory[:, 2]).reshape(-1, )[-1]
    endpoint_y_p = np.array(ray.trajectory[:, 3]).reshape(-1, )[-1]
    #print('endpoint is equal to', endpoint_z)

    assert ((endpoint_z == pytest.approx(1)) and (endpoint_x == pytest.approx(0.010)) and (endpoint_x_p == pytest.approx(0.010)) and (endpoint_y == pytest.approx(0.00)) and (endpoint_y_p == pytest.approx(0.00)))

def test_drift_space():
    """test creation of drift space and tracing through it
    it is tested whether central particle with 10 mrad angle
    after 1 meter drift has appropriate coordinates"""
    start = np.array([[0.000], [0.010], [0.000], [0.00], [0], [0]])
    proton = lb.Particle(1, 30, 1, start)
    driftspace = le.driftspace(length=1.0)
    beamline = lb.BeamLine([driftspace])
    ray = beamline.trace(proton)
    endpoint_z=np.array(ray.trajectory[:,6]).reshape(-1,)[-1]
    endpoint_x = np.array(ray.trajectory[:, 0]).reshape(-1, )[-1]
    endpoint_x_p = np.array(ray.trajectory[:, 1]).reshape(-1, )[-1]
    endpoint_y = np.array(ray.trajectory[:, 2]).reshape(-1, )[-1]
    endpoint_y_p = np.array(ray.trajectory[:, 3]).reshape(-1, )[-1]
    #print('endpoint is equal to', endpoint_z)

    assert ((endpoint_z == pytest.approx(1)) and (endpoint_x == pytest.approx(0.010)) and (endpoint_x_p == pytest.approx(0.010)) and (endpoint_y == pytest.approx(0.00)) and (endpoint_y_p == pytest.approx(0.00)))



def test_drift_space_with_thin_lens():
    """test creation of drift space with thin lens in the middle of it
    with the focus length of the half of the drift space.
    It is tested whether  particle on axis with 10 mrad angle
    turns into parallel after the lens"""
    start = np.array([[0.000], [0.010], [0.000], [0.00], [0], [0]])
    proton = lb.Particle(1, 30, 1, start)
    driftspace = le.driftspace(length=0.5)
    thin_lens=lb.ThinLens(0.5,0.5,0)
    beamline = lb.BeamLine([driftspace, thin_lens, driftspace])
    ray = beamline.trace(proton)
    endpoint_z=np.array(ray.trajectory[:,6]).reshape(-1,)[-1]
    endpoint_x = np.array(ray.trajectory[:, 0]).reshape(-1, )[-1]
    endpoint_x_p = np.array(ray.trajectory[:, 1]).reshape(-1, )[-1]
    endpoint_y = np.array(ray.trajectory[:, 2]).reshape(-1, )[-1]
    endpoint_y_p = np.array(ray.trajectory[:, 3]).reshape(-1, )[-1]
    #print('endpoint is equal to', endpoint_z)

    assert ((endpoint_z == pytest.approx(1)) and (endpoint_x == pytest.approx(0.005)) and (endpoint_x_p == pytest.approx(0.000)) and (endpoint_y == pytest.approx(0.00)) and (endpoint_y_p == pytest.approx(0.00)))