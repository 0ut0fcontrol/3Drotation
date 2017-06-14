#!/usr/bin/env python

import numpy as np

# A structure for a grid that discretizes six dimensions one by one.
# The six dimensions are in the sequence of R, PHI, THETA for the translational space,
# followed by PHI1, PHI2, THETA for the rotational space represented by the quaternion q.
class Node:
    def __init__(self, xs=[]):
        self.xs = xs # for leaves, xs=[]
        self.y = None # for leaves, y is a np-array
        self.next = [] # for leaves, next=[] 

    def count(self):
        if not self.next:
            return 1
        else:
            return sum(node.count() for node in self.next)

class Grid:
    # Set the parameters that control the density of the grid points
    # For the meaning of ang_params and ori_params, check the _get_dl() method
    def __init__(self):
        self.headNode = None
        self.n = None
        # wtr-wtr
        self.rs = np.concatenate((np.linspace(2.0, 3.5, 15), np.linspace(3.5, 6.0, 9)[1:], np.linspace(6.0, 12.0, 8)[1:]))
        
        #self.rs = np.concatenate((np.linspace(2.0, 3.8, 18), np.linspace(3.8, 4.7, 4)[1:], 
        #                          np.linspace(5.0, 7.0, 5),  np.linspace(7.0, 12.0, 6)[1:]))
        self.ang_params = (0.20, 0.45, 3.8, 0.4)
        self.ori_params = (0.23, 0.45, 3.8, 0.4)
    
    # Setup the grid structure without assigning the y values
    def setup(self):
        rs = self.rs
        node_r = Node(rs)
        for r in rs:
            phis = self._discretize_phi(r)
            node_phi = Node(phis)
            node_r.next.append(node_phi)
            for phi in phis:
                thetas = self._discretize_theta(r, phi)
                node_theta = Node(thetas)
                node_phi.next.append(node_theta)
                for theta in thetas:
                    ophi1s = self._discretize_ophi1(r)
                    node_ophi1 = Node(ophi1s)
                    node_theta.next.append(node_ophi1)
                    for ophi1 in ophi1s:
                        ophi2s = self._discretize_ophi2(r, ophi1)
                        node_ophi2 = Node(ophi2s)
                        node_ophi1.next.append(node_ophi2)
                        for ophi2 in ophi2s:
                            othetas = self._discretize_otheta(r, ophi1, ophi2)
                            node_otheta = Node(othetas)
                            node_ophi2.next.append(node_otheta)
                            for otheta in othetas:
                                leaf = Node()
                                node_otheta.next.append(leaf)
        self.headNode = node_r
        self.n = self._count()
        print "The grid structure is now set up. The grid consists of %d points." % self.n

    # carry out interpolation
    def interpolate(self, coor, order=2):
        return self._interpolate_help(coor, self.headNode, order)

    # recursive helper
    def _interpolate_help(self, coor, node, order):
        if not node.next:
            return node.y
        my_x = coor[0]
        coor = coor[1:]
        if order == 1:
            neighbors = self._find_neighbors2(node.xs, my_x)
        elif order == 2:
            neighbors = self._find_neighbors3(node.xs, my_x)
        elif order == 3:
            neighbors = self._find_neighbors4(node.xs, my_x)
        else:
            raise Exception("Invalid order for interpolant. Choose from 1, 2 and 3")
        xs = []
        ys = []
        for i in neighbors:
            xs.append(node.xs[i])
            ys.append(self._interpolate_help(coor, node.next[i], order))
        my_y = self._interp_1D(xs, ys, my_x)
        fmt = '[' + '%.2f,'*(len(xs)-1) + '%.2f]'
        return my_y 

    # fill the y values for the grid points with the supplied objective function f
    def fill(self, f):
        for leaf, x in self._gen_leaves_with_x():
            leaf.y = f(x)

    # save the grid parameters and y values on each grid point to a text file
    def save(self, filename):
        with open(filename, 'w') as file:
            fmt = 'ANGPRM\t%f\t%f\t%f\t%f\n'
            file.write(fmt % self.ang_params)
            fmt = 'ORIPRM\t%f\t%f\t%f\t%f\n'
            file.write(fmt % self.ori_params)
            fmt = 'RS\t' + '%f\t' * (len(self.rs)-1) + '%f\n'
            file.write(fmt % tuple(self.rs))
            fmt = '%f\t' * 6 + '%f\n'
            for leaf, x in self._gen_leaves_with_x():
                if leaf.y == None:
                    leaf.y = [0.0 for i in range(7)]
                file.write(fmt % tuple(leaf.y))

    # load the grid parameters from a text file, build the grid structure, then 
    # load y values for the grid points from the text file
    def load(self, filename):    
        print "Loading y values for each grid point from", filename
        with open(filename) as file:
            for i in range(3):
                line = file.readline().split()
                if line[0] == 'ANGPRM': 
                    self.ang_params = tuple(float(entry) for entry in line[1:])
                elif line[0] == 'ORIPRM':
                    self.ori_params = tuple(float(entry) for entry in line[1:])
                elif line[0] == 'RS':
                    self.rs = np.array([float(entry) for entry in line[1:]])
                else:
                    raise Exception('Wrong data file format!')
            self.setup()
            for leaf, x in self._gen_leaves_with_x():
                line = file.readline().split()
                if len(line) != 7:
                    raise Exception('Data file format error!')
                leaf.y = np.array([float(entry) for entry in line])
            line = file.readline().split()
            if line:
                raise Exception('Extra lines in data file!')

    # Generate x values for all points
    def gen_x(self):
        for leaf, x in self._gen_leaves_with_x_help(self.headNode, []):
            yield x

    # Generate node leaves together with their x values, used by fill()
    def _gen_leaves_with_x(self):
        for leaf, x in self._gen_leaves_with_x_help(self.headNode, []):
            yield leaf, x

    # recursive helper
    def _gen_leaves_with_x_help(self, node, pre_x):
        if not node.next:
            yield node, pre_x
            return
        for i in range(len(node.next)):
            for _ in self._gen_leaves_with_x_help(node.next[i], pre_x+[node.xs[i]]):
                yield _

    def gen_grid_x(self):
        for grid, x in self._gen_grids_with_x_help(self.headNode, []):
            yield x

    def _gen_grids_with_x(self):
        for grid, x in self._gen_grids_with_x_help(self.headNode, []):
            yield grid, x

    def _gen_grids_with_x_help(self, node, pre_x, level = 3):
        if level == 0:
            yield node, pre_x + [0., 0., 0.]
            return
        for i in range(len(node.next)):
            for _ in self._gen_grids_with_x_help(node.next[i], pre_x+[node.xs[i]], level - 1):
                yield _
    # Count the total number of points in the grid
    def _count(self):
        return self.headNode.count()

    # Use a 4-parameter sigmoidal function to calculate the approximate distance between two points 
    # on the 3D sphere for angular DOFs or 4D sphere for rotational DOFs.
    # The four parameters specify the lower limit, upper limit, transition position and transition 
    # steepness of the sigmoidal function.
    def _get_dl(self, r, params):
        dl0, dl1, r0, s = params
        return (dl1 - dl0) / (1 + np.exp((r0 - r) / s)) + dl0

    # The following methods discretizes each dimension.
    # The density depends only on r.
    # attention: the discretizing range for each dimension fits the reflection 
    # planes of water's refCoor, need to consider propers ways to generalize
    def _discretize_phi(self, r):
        dl = self._get_dl(r, self.ang_params)
        mini, maxi = 0., np.pi/2
        l = maxi - mini
        n = max(1, np.ceil(l/dl))
        if n == 1: 
            return np.linspace(mini, maxi, n)
        else:
            da = l / (n-1)
            # attention: for periodicity, leave out for the moment
            #return np.linspace(mini-da, maxi, n+1)
            return np.linspace(mini, maxi, n)

    def _discretize_theta(self, r, phi):
        dl = self._get_dl(r, self.ang_params)
        mini, maxi = 0., np.pi
        l = (maxi - mini) * np.cos(phi)
        n = max(1, np.ceil(l/dl))
        if n == 1: 
            return np.linspace(mini, maxi, n)
        else:
            da = l / (n-1)
            #return np.linspace(mini-da, maxi+da, n+2)
            return np.linspace(mini, maxi, n)

    def _discretize_ophi1(self, r):
        dl = self._get_dl(r, self.ori_params)
        mini, maxi = 0., np.pi/2
        l = maxi - mini
        n = max(1, np.ceil(l/dl))
        if n == 1: 
            return np.linspace(mini, maxi, n)
        else:
            da = l / (n-1)
            #return np.linspace(mini-da, maxi, n+1)
            return np.linspace(mini, maxi, n)

    def _discretize_ophi2(self, r, ophi1):
        dl = self._get_dl(r, self.ori_params)
        mini, maxi = 0., np.pi/2
        l = (maxi - mini) * np.cos(ophi1)
        n = max(1, np.ceil(l/dl))
        if n == 1: 
            return np.linspace(mini, maxi, n)
        else:
            da = l / (n-1)
            #return np.linspace(mini-da, maxi, n+1)
            return np.linspace(mini, maxi, n)

    def _discretize_otheta(self, r, ophi1, ophi2):
        dl = self._get_dl(r, self.ori_params)
        mini, maxi = -np.pi, np.pi
        l = (maxi - mini) * np.cos(ophi1) * np.cos(ophi2)
        n = max(1, np.ceil(l/dl))
        if n == 1: 
            return np.linspace(mini, maxi, n)
        else:
            da = l / (n-1)
            #return np.linspace(mini-da, maxi+da, n+2)
            return np.linspace(mini, maxi, n)

    # The following methods carries out 1D neighbor search for the query value
    # attention: need better range check

    # Two points for linear interpolation
    def _find_neighbors2(self, xs, my_x):
        n = len(xs)
        if n <= 2:
            return range(len(xs))
        i = 1
        while i < n:
            if xs[i] > my_x:
                break
            i += 1
        if i == n:
            print xs, my_x
            raise Exception("x value out of range!")
        return [i-1, i]

    # Three points for 2nd order interpolation
    def _find_neighbors3(self, xs, my_x):
        n = len(xs)
        if n <= 3:
            return range(len(xs))
        i = 1
        while i < n:
            if xs[i] > my_x:
                break
            i += 1
        if i == n:
            raise Exception("x value out of range!")
        if i < 2:
            return [i-1, i, i+1]
        if i > n - 2:
            return [i-2, i-1, i]
        return [i-1, i, i+1]

    # Four points for 3rd order interpolation
    def _find_neighbors4(self, xs, my_x):
        n = len(xs)
        if n < 4:
            return range(len(xs))
        i = 1
        while i < n:
            if xs[i] > my_x:
                break
            i += 1
        if i == n:
            raise Exception("x value out of range!")
        if i == 1:
            return [0, 1, 2]
        if i == n-1:
            return [n-3, n-2, n-1]
        return [i-2, i-1, i, i+1]

    # Carries out 1D Lagrange interpolation
    # choose the order of the interpolant based on the number of points
    def _interp_1D(self, xs, ys, my_x):
        if len(xs) == 1:
            return ys[0]
        if len(xs) == 2:
            x0, x1 = xs
            y0, y1 = ys
            return y0 + (my_x - x0) * (y1 - y0) / (x1 - x0)
        if len(xs) == 3:
            x0, x1, x2 = xs
            y0, y1, y2 = ys
            a1 = (y1 - y0) / (x1 - x0)
            a2 = (y2 - y0 - a1 * (x2 - x0)) / (x2 - x0) / (x2 - x1)
            return y0 + a1 * (my_x - x0) + a2 * (my_x - x0) * (my_x - x1)
        if len(xs) == 4:
            x0, x1, x2, x3 = xs
            y0, y1, y2, y3 = ys
            a1 = (y1 - y0) / (x1 - x0)
            a2 = (y2 - y0 - a1 * (x2 - x0)) / (x2 - x0) / (x2 - x1)
            a3 = (y3 - y0 - a1 * (x3 - x0) - a2 * (x3 - x0) * (x3 - x1)) / (x3 - x0) / (x3 - x1) / (x3 - x2)
            return y0 + a1 * (my_x - x0) + a2 * (my_x - x0) * (my_x - x1) + \
                    a3 * (my_x - x0) * (my_x - x1) * (my_x - x2)
            

    
if __name__ == "__main__":
    grid = Grid() 
    print "setting up"
    grid.setup()
    print "counting"
    print "Total number of points:", grid.n

