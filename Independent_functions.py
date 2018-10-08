import numpy as np
import math
from  scipy.spatial.distance import cdist
import random

class Molecule:
    def __init__(self, atids, xyz, bonds, labels, types):
        self.atids = np.array(atids)
        self.xyz = np.array(xyz)
        self.bonds = np.array(bonds)
        self.labels = np.array(labels)
        self.types = np.array(types)

    def remove_clashes(self, sys, bx, by, bz):
        np.random.seed(len(sys.xyz))
        random.seed(len(sys.xyz))

        d = np.min(cdist(self.xyz, sys.xyz))
        while d < 0.94: #0.94 is twice the standard bead radius of Martini
            self.xyz = center(self.xyz)

            ax = np.random.rand(3)
            alf = random.random()*2*math.pi
            self.xyz = rotate(self.xyz, ax, alf)

            shift = np.multiply(np.random.rand(3) - 0.5, [bx, by, bz])
            self.xyz = displace(self.xyz, shift)

            d = np.min(cdist(self.xyz, sys.xyz))

        return self.xyz

class System:
    def __init__(self, bx, by, bz):
        self.BX = bx
        self.BY = by
        self.BZ = bz
        self.molecules = []
        self.xyz = []
        self.labels = []

    def add_molecule(self, molecule):
        self.molecules.append(molecule)
        if len(self.xyz) == 0:
            self.xyz = molecule.xyz
            self.labels = molecule.labels
        else:
            self.xyz = np.vstack((self.xyz, molecule.xyz))
            self.labels = np.append(self.labels, molecule.labels)


def center(xyz):
    COM = np.average(xyz, axis=0)
    xyz = xyz - COM
    return xyz

def rotate(pts, vec, ang):
    #rotates the points along vector vec in an angle ang in radians
    vec = vec / np.linalg.norm(vec)
    ct = math.cos(ang)
    st = math.sin(ang)
    x = vec[0]
    y = vec[1]
    z = vec[2]
    rot = np.array([[ct + x**2*(1-ct), x*y*(1-ct)-z*st, x*z*(1-ct)+y*st],     [x*y*(1-ct)+z*st, ct+y**2*(1-ct), y*z*(1-ct)-x*st],    [x*z*(1-ct)-y*st, y*z*(1-ct)+x*st, ct+z**2*(1-ct)]])
    return np.dot(rot, pts.T).T

def displace(pts, shift):
    for i in range(len(pts)):
        pts[i] = pts[i] + shift
    return pts

def print_xyz(outname, coords, names):
    coords = coords * 10
    output = open(outname + ".xyz", "w")
    output.write(str(len(coords)) + "\n\n")
    for i in range(len(coords)):
        output.write(names[i] + '{:.3f}'.format(coords[i,0]).rjust(10) + "{:.3f}".format(coords[i,1]).rjust(10) + "{:.3f}".format(coords[i,2]).rjust(10) + "\n")
    output.close()
