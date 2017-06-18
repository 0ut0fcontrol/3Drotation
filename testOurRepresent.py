#!/usr/bin/env python2
import os,sys
import numpy as np
import mol
from tools import *


def Xq2traj(mol, Xq_list, filename="xq2traj.pdb"):
    """
    mol is mol object
    Xq list:
    [(X1, q1) 
     (X1, q1)
    ]

    """
    f = open(filename,'w')
    for i, coor in enumerate(Xq_list):
        if len(coor) == 2:pdbstr = mol.Xq2PDB(Xcom=coor[0],q=coor[1])
        if len(coor) == 3:pdbstr = mol.Xq2PDB(Xcom=coor)
        if len(coor) == 4:pdbstr = mol.Xq2PDB(q=coor)
        f.write("MODEL%8d\n"%(i))
        f.write(pdbstr)
        f.write("TER\nENDMDL\n")

mol = mol.molType('alc')

def tergetAxis2q(axis):
    origin = np.array((1.,0.,0.,))
    target = np.array(axis, dtype=np.float64)
    target /= np.linalg.norm(target)
    midAxis = (origin + target)/2.
    mid_norm = np.linalg.norm(midAxis)
    if mid_norm < 0.000001:
        midAxis = (0., 1., 0.)
    else:
        midAxis /= mid_norm
    return np.array((0.,midAxis[0],midAxis[1],midAxis[2]))

Axises = [(1,0,0), (0,1,0),(0,0,1)]
q_list = []
for axis in Axises:
    q = tergetAxis2q(axis)
    q_list.append(q)
Xq2traj(mol,q_list,"targetAxis.pdb")


def tergetAxisAng2q(axis,angle):
    origin = np.array((1.,0.,0.,))
    target = np.array(axis, dtype=np.float64)
    target /= np.linalg.norm(target)
    midAxis = (origin + target)/2.
    mid_norm = np.linalg.norm(midAxis)
    if mid_norm < 0.000001:
        midAxis = (0., 1., 0.)
    else:
        midAxis /= mid_norm
    AxisQ = np.array((0.,midAxis[0],midAxis[1],midAxis[2])) # rotation pi 
    angle = (angle + np.pi)/2.
    sinA = np.sin(angle)
    angQ = np.array((np.cos(angle),sinA*target[0],sinA*target[1], sinA*target[2]))
    return qmult(AxisQ, angQ)


Axises = [(1,0,0), (0,1,0),(0,0,1)]
angles = np.linspace(0,2*np.pi, 10)
q_list = []
for axis in Axises:
    for ang in angles:
        q = tergetAxisAng2q(axis,ang)
        q_list.append(q)
Xq2traj(mol,q_list,"targetAxisAng.pdb")


def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

#Axises = sample_spherical(100)
Axises = [(1,0,0), (1,1,1),(0,1,1),(-1,1,1),(-1,0,0),(-1,-1,-1),(0,-1,-1),(1,-1,-1),(1,0,0)]
angles = np.linspace(0,2*np.pi, 10)
q_list = []
#for axis in Axises.T:
for axis in Axises:
    for ang in angles:
        q = tergetAxisAng2q(axis,ang)
        q_list.append(q)
Xq2traj(mol,q_list,"targetAxisAng.pdb")
"""
supersph_list= [(0.,0.,i) for i in np.linspace(0., np.pi, 10)]
q_list = [spherical2q(r,t,p) for r, t, p in supersph_list]
Xq2traj(mol,q_list,"superSph3.pdb")

supersph_list= [(0.,i,0.) for i in np.linspace(0., np.pi, 10)]
q_list = [spherical2q(r,t,p) for r, t, p in supersph_list]
Xq2traj(mol,q_list,"superSph2.pdb")

supersph_list= [(0.,i,0.) for i in np.linspace(0., np.pi, 10)]
q_list = [spherical2q(r,t,p) for r, t, p in supersph_list]
Xq2traj(mol,q_list,"superSph1.pdb")

q_list = []
for r in np.linspace(0., np.pi/2, 10):
    for t in np.linspace(0., np.pi, 10):
        for p in np.linspace(-np.pi, np.pi, 10):
            q_list.append(spherical2q(r,t,p))
Xq2traj(mol,q_list,"superSph.pdb")
    
# hopf global sampling   
q_list = []
for a in np.linspace(0., np.pi, 10):
    for b in np.linspace(-np.pi, np.pi, 20):
        for g in np.linspace(-np.pi, np.pi, 20):
            q_list.append(hopf2q(a,b,g))
Xq2traj(mol,q_list,"hopf.pdb")

# hopf local sampling   
q_list = []
log = open("localhopf.pdb.log", 'w')
#for a in np.linspace(0., np.pi, 1a0):
for a in (np.pi/4., np.pi/6.):
    for b in (-np.pi/4., np.pi/4.):
        for g in np.linspace(-np.pi, np.pi, 10):
            q = hopf2q(a,b,g)
            q_list.append(q)
#            print((a,b,g),q)
Xq2traj(mol,q_list,"localhopf.pdb")


# Axis-Angle local sampling
q_list = []
Axises = [[1,1,i] for i in range(10)]
Axises = [[1,0,0], [-1,0,0], [0,-1,0]] + Axises
angles = np.linspace(-np.pi, np.pi, 10)
for axis in Axises:
    for ang in angles:
        q_list.append(axisAng2q(axis, ang))
Xq2traj(mol,q_list,"localAxisAng.pdb")
"""
