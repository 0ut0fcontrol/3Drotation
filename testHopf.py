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
