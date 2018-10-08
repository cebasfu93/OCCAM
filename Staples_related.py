import numpy as np
import math
from  scipy.spatial.distance import cdist

#For future placement of ligands
def sunflower_pts(num_pts, rad):
    indices = np.arange(0, num_pts, dtype=float) + 0.5

    phi = np.arccos(1 - 2*indices/num_pts)
    theta = math.pi * (1 + 5**0.5) * indices

    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    xyz = rad*np.array([x,y,z]).T
    return xyz

def put_staples(shell, radius, num_ligs):
    S_atoms = []
    if num_ligs > 12:
        S_atoms = sunflower_pts(num_ligs, radius)
        #S_atoms = hollow_sphere(radius, lignum_opt)
        distances = cdist(S_atoms, shell)
        mins = np.argmin(distances, axis=1)
    elif num_ligs == 6:
        r = radius
        a = math.sqrt(2)/2
        S_atoms = np.array([[r,0,0],[-r,0,0],[0,r,0],[0,-r,0],[0,0,r],[0,0,-r]])
        #S_atoms = np.array([[a*r,0,-a*r],[-a*r,0,a*r],[0,r,0],[0,-r,0],[a*r,0,a*r],[-a*r,0,-a*r]])
        distances = cdist(S_atoms, shell)
        mins = np.argmin(distances, axis=1)
    else:
        print("There are not enough ligands to be considered a homogeneous distribution")
    for i in range(len(S_atoms)):
        norma=np.linalg.norm(shell[mins[i]])
        scaling = (norma + 0.47)/norma #0.47 is the standard bead-bead distance in Martini
        S_atoms[i,:] = scaling*shell[mins[i]]
    return S_atoms
