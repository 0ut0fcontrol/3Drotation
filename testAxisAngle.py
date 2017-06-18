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


def ang2R(y_ang,z_ang,axis_ang):
    #y_ang: ang of xy face
    #z_ang: ang of z axis
    #axis_ang: ang of axis
    yR = axisAng2R(np.array((0.,-1.,0.)),y_ang)
    zR = axisAng2R(np.array((0.,0.,1.)),z_ang) 
    zyR = np.dot(zR,yR)
    axis = zyR[:,0]
    #axis = np.array((np.cos(y_ang)*np.cos(z_ang),np.cos(y_ang)*np.sin(z_ang),np.sin(y_ang)))
    axisR = axisAng2R(axis, axis_ang)
    #R = np.dot(axisR, zyR)
    R = np.dot(axisR, zyR)
    #z_ang2 = np.arctan2(zyR[1,0], zyR[0,0])
    #if np.linalg.norm(R[:,0]-zyR[:,0]) >0.01:
    z_ang2 = np.arctan2(R[1,0], R[0,0])
    eps = 1.E-4
    if abs(z_ang - np.pi) < eps or abs(z_ang + np.pi) < eps: z_ang = -np.pi
    if abs(z_ang2 - np.pi) < eps or abs(z_ang2 + np.pi) < eps:z_ang2 = -np.pi
    if z_ang2 - z_ang >0.001:
        print("cont not recovery z_ang: org:%.5f rev:%.5f"%(z_ang,z_ang2))
        print("zyR:")
        print(zyR)
        print("R:")
        print(R)
    return R,axisR

def R2ang(R):
    R = np.mat(R)
    y_ang = np.arcsin(R[2,0])
    z_ang = np.arctan2(R[1,0], R[0,0])
    eps = 1.E-6
    yR = axisAng2R(np.array((0.,-1.,0.)),y_ang)
    zR = axisAng2R(np.array((0.,0.,1.)),z_ang)
    zyR = np.dot(zR,yR)
    if abs(z_ang - np.pi) < eps or abs(z_ang + np.pi) < eps: z_ang = -np.pi
    axisR = np.dot(R, zyR.T)
    axis_ang = np.arccos((axisR[0,0]+axisR[1,1]+axisR[2,2]-1.)/2)
    return (y_ang, z_ang, axis_ang),axisR


#y_angs = np.linspace(-np.pi/2, np.pi/2, 10)
y_angs = np.linspace(-np.pi/2, np.pi/2, 9)
z_angs = np.linspace(-np.pi, np.pi, 9)
axis_angs = np.linspace(-np.pi, np.pi, 9)
#axis_angs = np.linspace(0, 2*np.pi, 9)
R_list = []
for y_ang in y_angs:
    for z_ang in z_angs:
        for axis_ang in axis_angs:
            R,axisR0 = ang2R(y_ang,z_ang,axis_ang)
            (a,b,g), axisR1 = R2ang(R)
            eps = 1.E-6
            if abs(z_ang - np.pi) < eps or abs(z_ang + np.pi) < eps:z_ang = -np.pi
            if not np.allclose((y_ang,z_ang,axis_ang),(a,b,g),1.E-6):
                print("orgin: %.3f,%3f,%3f"%(y_ang,z_ang,axis_ang))
                print(axisR0)
                print("recov: %.3f,%3f,%3f"%(a,b,g))
                print(axisR1)
                print("\n")
            R_list.append(R)
XR2traj(mol, R_list, "XR.pdb")
