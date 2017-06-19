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

def XR2traj(mol, XR_list, filename="xr2traj.pdb"):
    """
    mol is mol object
    Xq list:
    [(X1, q1) 
     (X1, q1)
    ]

    """
    f = open(filename,'w')
    for i, coor in enumerate(XR_list):
        if len(coor) == 2:
            pdbstr = mol.XR2PDB(Xcom=coor[0],R=coor[1])
        else:
            pdbstr = mol.XR2PDB(R=coor)
        f.write("MODEL%8d\n"%(i))
        f.write(pdbstr)
        f.write("TER\nENDMDL\n")

mol = mol.molType('wtr')
Xq2traj(mol, [[1,0,0,0]], "ref.pdb")



def ang2M(y_ang,z_ang,x_ang):
    a = x_ang
    b = y_ang 
    g = z_ang
    sinA = np.sin(a)
    cosA = np.cos(a)
    sinB = np.sin(b)
    cosB = np.cos(b)
    sinG = np.sin(g)
    cosG = np.cos(g)
    M = [[ cosB*cosG,  sinA*sinB*cosG + cosA*sinG, -cosA*sinB*cosG + sinA*sinG ],
         [-cosB*sinG, -sinA*sinB*sinG + cosA*cosG,  cosA*sinB*sinG + sinA*cosG ],
         [ sinB,      -sinA*cosB,                   cosA*cosB                  ]]
    M = np.array(M)
    return M

def M2ang(M):
    a = np.arctan2(-M[2,1], M[2,2])
    b = np.arcsin(M[2,0])
    g = np.arctan2(-M[1,0], M[0,0])
    x_ang = a
    y_ang = b
    z_ang = g
    return y_ang,z_ang,x_ang

#y_angs = np.linspace(-np.pi/2, np.pi/2, 10)
y_angs = np.linspace(-np.pi/2, np.pi/2, 9)
z_angs = np.linspace(-np.pi, np.pi, 9)
axis_angs = np.linspace(-np.pi, np.pi, 9)
#axis_angs = np.linspace(0, 2*np.pi, 9)
R_list = []
for y_ang in y_angs:
    for z_ang in z_angs:
        for axis_ang in axis_angs:
            R = ang2M(y_ang,z_ang,axis_ang)
            (a,b,g) = M2ang(R)
            eps = 1.E-6
            if abs(z_ang - np.pi) < eps or abs(z_ang + np.pi) < eps:
                z_ang = -np.pi
                b = -np.pi
            if abs(axis_ang - np.pi) < eps or abs(axis_ang + np.pi) < eps:
                axis_ang = -np.pi
                g = -np.pi
            if not np.allclose((y_ang,z_ang,axis_ang),(a,b,g),1.E-6):
                print("orgin: %.3f,%3f,%3f"%(y_ang,z_ang,axis_ang))
                print(axisR0)
                print("recov: %.3f,%3f,%3f"%(a,b,g))
                print(axisR1)
                print("\n")
            R_list.append(R)
XR2traj(mol, R_list, "M2ang.pdb")
