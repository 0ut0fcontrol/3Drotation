#!/usr/bin/env python

import numpy as np


# Functions for quaternion algebra, quaternion-rotation matrix conversion,
# Cartesian-spherical coordinate conversion, etc.

def qmirror(M, q):
    newq = np.copy(q)
    newq[1:] = - np.dot(M, newq[1:])
    return newq

def qequal(q0, q1):
    eps = 1.E-6
    return np.linalg.norm(q0-q1) < eps or np.linalg.norm(q0+q1) < eps

# q = a * b 
# Rotate by q equivalent to rotate by a, then rotate by b
def qmult(a, b):
    q = np.zeros(4)
    q[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
    q[1] = a[0]*b[1] + a[1]*b[0] - a[2]*b[3] + a[3]*b[2]
    q[2] = a[0]*b[2] + a[2]*b[0] + a[1]*b[3] - a[3]*b[1]
    q[3] = a[0]*b[3] + a[3]*b[0] - a[1]*b[2] + a[2]*b[1]
    return q

def qinv(q):
    q = np.copy(q)
    q[0] = -q[0]
    return q

# a = q * b
# Rotate by a equivalent to rotate by q, then rotate by b
def qdiv(a, b):
    q = qmult(a, qinv(b))
    return q

def q2R(q):
    q0, q1, q2, q3 = q
    R = np.array([[1-2*q2**2-2*q3**2, 2*(q1*q2-q3*q0), 2*(q1*q3+q2*q0)],
                  [2*(q1*q2+q3*q0), 1-2*q1**2-2*q3**2, 2*(q2*q3-q1*q0)],
                  [2*(q1*q3-q2*q0), 2*(q2*q3+q1*q0), 1-2*q1**2-2*q2**2]])
    return R

def R2q(R):
    q0 = ( R[0,0] + R[1,1] + R[2,2] + 1.) / 4.
    q1 = ( R[0,0] - R[1,1] - R[2,2] + 1.) / 4.
    q2 = (-R[0,0] + R[1,1] - R[2,2] + 1.) / 4.
    q3 = (-R[0,0] - R[1,1] + R[2,2] + 1.) / 4.
    if(q0 < 0.): q0 = 0.
    if(q1 < 0.): q1 = 0.
    if(q2 < 0.): q2 = 0.
    if(q3 < 0.): q3 = 0.
    q0 = np.sqrt(q0)
    q1 = np.sqrt(q1)
    q2 = np.sqrt(q2)
    q3 = np.sqrt(q3)
    if q0 >= q1 and q0 >= q2 and q0 >= q3: 
        q0 *= +1.
        q1 *= np.sign(R[2,1] - R[1,2])
        q2 *= np.sign(R[0,2] - R[2,0])
        q3 *= np.sign(R[1,0] - R[0,1])
    elif q1 >= q0 and q1 >= q2 and q1 >= q3:
        q0 *= np.sign(R[2,1] - R[1,2])
        q1 *= +1.
        q2 *= np.sign(R[1,0] + R[0,1])
        q3 *= np.sign(R[0,2] + R[2,0])
    elif q2 >= q0 and q2 >= q1 and q2 >= q3:
        q0 *= np.sign(R[0,2] - R[2,0])
        q1 *= np.sign(R[1,0] + R[0,1])
        q2 *= +1.
        q3 *= np.sign(R[2,1] + R[1,2])
    elif q3 >= q0 and q3 >= q1 and q3 >= q2:
        q0 *= np.sign(R[1,0] - R[0,1])
        q1 *= np.sign(R[2,0] + R[0,2])
        q2 *= np.sign(R[2,1] + R[1,2])
        q3 *= +1.
    else:
        printf("coding error\n")
    q = np.array([q0, q1, q2, q3])
    q /= np.linalg.norm(q)
    return q

def xyz2spherical(X):
    x, y, z = X
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(z, np.sqrt(x**2 + y**2))
    theta = np.arctan2(y, x)
    return r, phi, theta

def spherical2xyz(r, phi, theta):
    x = r * np.cos(phi) * np.cos(theta)
    y = r * np.cos(phi) * np.sin(theta)
    z = r * np.sin(phi)
    return np.array([x, y, z])

def q2spherical(q):
    if q[0] < 0: q = -q
    q0, q1, q2, q3 = q
    theta = np.arctan2(q2, q3)
    phi2 = np.arctan2(q1, np.sqrt(q2**2 + q3**2))
    phi1 = np.arctan2(q0, np.sqrt(q1**2 + q2**2 + q3**2))
    return phi1, phi2, theta

def spherical2q(phi1, phi2, theta):
    q0 = np.sin(phi1)
    q1 = np.cos(phi1) * np.sin(phi2)
    q2 = np.cos(phi1) * np.cos(phi2) * np.sin(theta)
    q3 = np.cos(phi1) * np.cos(phi2) * np.cos(theta)
    return np.array([q0, q1, q2, q3])

