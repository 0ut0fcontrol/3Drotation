import numpy as np

def rotate(vec, axis, theta):
        """
        According to "Euler-Rodrigues formula"
                http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
        There was a negative sign '-' in front of the right term in "a = -np.cos(theta*0.5)"
                but not in Wikipedia formula.
        It seems, without the '-' sign results in "right hand role".
        """
        axis = axis/np.sqrt(np.dot(axis,axis))
        a = np.cos(theta*0.5)
        b,c,d = axis*np.sin(theta*0.5)
        rotmat = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                           [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                           [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
        return np.dot(vec,rotmat)

def axisAng2R(axis,angle):
        axis /= np.linalg.norm(axis)
        a = np.cos(angle*0.5)
        b,c,d = axis*np.sin(angle*0.5)
        rotmat = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                           [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                           [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
        return rotmat
def R2axisAng(R):
    """Conversion Matrix to Axis Angle
    http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/
    """
    epsilon = 0.01 # margin to allow for rounding errors
    epsilon2 = 0.1 # margin to distinguish between 0 and 180 degrees
    if (abs(R[0,1]-R[1,0]<epsilon 
        and abs(R[0,2]-R[2,0]<epsilon
        and abs(R[1,2]-R[2,1]<epsilon):
        # singularity found
        # first check for identity matrix which must have +1 for all terms
        # in leading diagonaland zero in other terms
        if (abs(R[0,1]+R[1,0]< epsilon2
            and abs(R[0,2]+R[2,0]< epsilon2
            and abs(R[1,2]+R[2,1]< epsilon2
            and abs(R[0,0]+R[1,1]+R[2,2]-3.) < epsilon2):
            # this singularity is identity matrix so angle = 0
            return np.array((1., 0., 0.)), 0. # arbitrary axis, zero angle
        # otherwise this singularity is angle = 180
        angle = np.pi
        axis = [0.0]
        xx = (R[0,0]+1)/2.
        yy = (R[1,1]+1)/2.
        zz = (R[2,2]+1)/2.
        xy = (R[0,1]+R[1,0])/4.
        xz = (R[0,2]+R[2,0])/4.
        yz = (R[1,2]+R[2,1])/4.
        if ((xx > yy) and (xx > zz)):
            if xx< epsilon:
                x 

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

def q2hopf(q):
    """convert quaternion to hopf coordinate

    hopf coordinate = 2-spherical + circle
    2-spherical in spherical coordinate system:(r=1,theta[0, pi),phi[-pi, pi))
    2-spherical in cartesian coordinate system:x, y ,z
    circle in angele: psi[-pi,pi)
    # see Spherical coordinate system
    #(https://en.wikipedia.org/wiki/Spherical_coordinate_system)

    # theta: angle of z, 
    # phi: angle of x,
    # psi: angle of norm
    """
    q0, q1, q2, q3 = q
    psi = np.arctan2(q1, q0) # only a half of psi
    phi = np.arctan2(q3, q2) - psi
    theta = np.arccos(          (q0+q1)
                      / (np.cos(psi)+np.sin(psi)) # avoid the zero point
                     ) * 2

    psi = psi * 2
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return (np.array([x,y,z]), psi)

#def hopf2q(vector, angle):
#    # Hopf Coordinates for SO(3)
#    # ref: Generating Uniform Incremental Grids on SO(3) Using the Hopf Fibration
#    x, y, z = vector
#    r = np.sqrt(x**2 + y**2 + z**2)
#    theta = np.arccos(z/r)
#    phi = np.arctan2(y, x)
#    psi = angle 
#    cosTheta = np.cos(theta/2.0)
#    sinTheta = np.sin(theta/2.0)
#
#    q0 = cosTheta * np.cos(psi/2.0)
#    q1 = cosTheta * np.sin(psi/2.0)
#    q2 = sinTheta * np.cos(phi + psi/2.0)
#    q3 = sinTheta * np.sin(phi + psi/2.0) 
#    return np.array([q0, q1, q2, q3])
#
def hopf2q(a,b,g):
    # Hopf Coordinates for SO(3)
    # ref: http://marc-b-reynolds.github.io/quaternions/2017/05/12/HopfCoordConvert.html
    ha = a/2.
    hg = g/2.
    sinHa = np.sin(ha)
    cosHa = np.cos(ha)
    #q0 = sinHa*np.sin(hg-b)
    #q1 = sinHa*np.cos(hg-b)
    #q2 = cosHa*np.sin(hg)
    #q3 = cosHa*np.cos(hg)
    q3 = sinHa*np.sin(hg + b)
    q2 = sinHa*np.cos(hg + b)
    q1 = cosHa*np.sin(hg)
    q0 = cosHa*np.cos(hg)
    return np.array([q0, q1, q2, q3])

def axisAng2q(axis,angle):
    n = np.array(axis,dtype=np.float64)
    norm = np.linalg.norm(n)
    n /= np.linalg.norm(n)
    angle /= 2.
    sinA = np.sin(angle)
    return np.array([np.cos(angle), sinA*n[0], sinA*n[1], sinA*n[2]])

def vet2ang(x, y):
    """get the angle of 2 vector

    """
    lx = np.sqrt(np.dot(x,x))
    ly = np.sqrt(np.dot(y,y))
    cos_angle = np.dot(x,y)/(lx * ly)
    angle = np.arccos(cos_angle)
    return angle

def sphere_triang_area(OA,OB,OC, r = 1):
    """get area of spherical triangle from 3 vectors (O point to surface).

    """
    a = vet2ang(OB,OC)
    b = vet2ang(OA,OC)
    c = vet2ang(OA,OB)
    cosA = (np.cos(a) - np.cos(b)*np.cos(c))/(np.sin(b)*np.sin(c))
    cosB = (np.cos(b) - np.cos(a)*np.cos(c))/(np.sin(a)*np.sin(c))
    cosC = (np.cos(c) - np.cos(b)*np.cos(a))/(np.sin(b)*np.sin(a))
    E = np.arccos(cosA) + np.arccos(cosB) + np.arccos(cosC) - np.pi
    return (E * r**2)
    
