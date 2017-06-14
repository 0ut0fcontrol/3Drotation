import numpy as np
import tools

MASS =    {'O':15.999,  'N':14.010,
           'C':12.010,  'H': 1.008,
           'F':19.000,
           'Na':22.99,  'NA':22.99,
           'P':30.970,  'S':32.060,
           'Cl':35.45,  'CL':35.45,
           'Br':79.90,  'BR':79.90,
           'CU':63.55,
           }

ZCharge = {'O':'8.0','N':'7.0',\
           'C':'6.0','H':'1.0',\
           'B':'5.0','F':'9.0',\
           'Na':'11.0','NA':'11.0',\
           'P':'15.0','S':'16.0',\
           'Cl':'17.0','CL':'17.0',\
           'Br':'35.0','BR':'35.0',\
           'CU':'29.0','RU':'44.0',\
           'AU':'79.0','Au':'79.0','Se':'34.0'}


class Mol:
    def __init__(self, mol, xyz, refl_axes=(), frg = 'FRG'):
        """ init a molecular type.
        mol is a list of element and coordinate, like:
            [['H', 0.0, 0.0, 0.1]
             ['O', 0.2, 0.2, 0.2]
             ['H', 0.1, 0.1, 0.1]
             ]
        
        xyz is the atom to define xyz.
        it has two type: 
        one x axis in Angle bisector:like wtr
        ['bisec', 0, 1, 2 ]
        one x axis is along 1st atom 2 2nd atom, like alc:
        ['along', 0, 4, 5 ]

        refl_axes is a list of indices of axes. Reflection along each of these
        axes corresponds to a mirror symmetry of the molecule
        """
        self.mol = mol
        self.xyz = xyz
        self.frg = frg
        self.n = len(mol)
        self.refl_axes = refl_axes
        self.mass = np.array([MASS[i[0]] for i in self.mol])
        refCoor   = np.array([i[1:] for i in self.mol])

        # The following code ensures that refCoor has COM at origin and orientation
        # aligned with the getR() method
        refCoor = refCoor - self.getCOM(refCoor)
        R = self.getR(refCoor) 
        refCoor = np.dot(refCoor, R)
        self.refCoor = refCoor
    # Calculate the rotation matrix R that relates coors to self.refCoor
    # vec_in_reference_frame = R \dot vec_in_body_frame 
    # R(coors) \dot refCoor = coors - COM(coors)
    # This function defines the orientation of self.refCoor
    # Need to be consistent with self.refl_axes
    def getR(self, coors):
        coors = np.copy(coors)
        i, j, k = self.xyz[1:]
        offset = coors[i]
        coors -= offset
        xvec = None
        if self.xyz[0] is 'bisec':
            xvec = coors[j] + coors[k]
        if self.xyz[0] is 'along':
            xvec = coors[j]
        zvec = np.cross(coors[j], coors[k])
        yvec = np.cross(zvec, xvec)
        xvec /= np.linalg.norm(xvec)
        yvec /= np.linalg.norm(yvec)
        zvec /= np.linalg.norm(zvec)
        R = np.array([xvec, yvec, zvec]).T 
        return R

    # Calculate the center of mass
    def getCOM(self, coors):
        return np.dot(self.mass, coors) / self.mass.sum()

    # Convert atomic coordinates to Xcom and q
    def atomic2Xq(self, coors):
        Xcom = self.getCOM(coors)
        R = self.getR(coors)
        q = tools.R2q(R)
        return Xcom, q

    # Given Xcom and q, rebuild the atomic coordinates
    def Xq2Atomic(self, Xcom, q):
        R = tools.q2R(q)
        coor = np.dot(self.refCoor, R.T)
        coor += Xcom
        return coor

    def coor2PDB(self, coor, frg_idx=0, occupancy=1.0, bfactor=0.0):
        pdb = ''
        for i in range(self.n):
            ele = self.mol[i][0]
            x,y,z = coor[i]
            pdb += 'ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%(
                    i, ele, self.frg, 'A',frg_idx, x, y, z, occupancy, bfactor)
        return pdb
        
    # Given Xcom and q, rebuild the atomic coordinates in PDB format
    def Xq2PDB(self, Xcom=(0.,0.,0), q=(1.,0.,0.,0.), frg_idx=0, occupancy=1.0, bfactor=0.0):
        coor = self.Xq2Atomic(Xcom, q)
        return self.coor2PDB(coor, frg_idx, occupancy, bfactor)
    
    def coor2INP(self, coor):
        inp = ''
        for i in range(self.n):
            ele = self.mol[i][0]
            x, y, z = coor[i]
            inp += "%2s%8s%15.10f%15.10f%15.10f\n"%(
                ele, ZCharge[ele.upper()], x, y, z)
        return inp

    def Xq2INP(self, Xcom=(0.,0.,0), q=(1.,0.,0.,0.)):
        coor = self.Xq2Atomic(Xcom, q)
        return self.coor2INP(coor)

# A class that holds information related to the atomic structure of a water
# molecule. It also includes several methods that carries out operations 
# related to the atomic coordinates.
class wtr(Mol):
    def __init__(self):
        mol = [['O', -0.06556939,   0.00000000,    0.00000000],
               ['H',  0.52035943,   0.76114632,    0.00000000],
               ['H',  0.52035943,  -0.76114632,    0.00000000]
               ]
        xyz = ['bisec', 0, 1, 2 ]
        refl_axes = (1,2)
        frg = "HOH"
        Mol.__init__(self, mol, xyz, refl_axes, frg)

class alc(Mol):
    def __init__(self):
        mol = [['C', -0.73090998,   -0.02638086,   -0.00000969],
               ['H', -1.03950734,   -1.06465416,   -0.00019176],
               ['H', -1.13674450,    0.45621966,    0.88767058],
               ['H', -1.13699250,    0.45684043,   -0.88719977],
               ['O',  0.69497026,   -0.02638086,   -0.00000969],
               ['H',  0.99194520,    0.88465732,   -0.00000969]
               ]
        xyz = ['along', 0, 4, 5 ]
        refl_axes = (2,)
        frg = "ALC"
        Mol.__init__(self, mol, xyz, refl_axes, frg)

def molType(mol_type):
    molDict = {'wtr': wtr(),
               'alc': alc()
               }
    return molDict[mol_type]
