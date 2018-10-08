import numpy as np
from Independent_functions import center
def make_fcc_lattice(const, cellsX, cellsY, cellsZ):
    N_unit_cells = cellsX* cellsY* cellsZ
    N_beads = N_unit_cells * 14
    fcc_block = np.array([])

    for i in range(cellsX):
        for j in range(cellsY):
            for k in range(cellsZ):

                fcc_block = np.append(fcc_block, [i,j,k,i+1,j,k,i,j+1,k,i,j,k+1,i+1,j+1,k,i+1,j,k+1,i,j+1,k+1,i+1,j+1,k+1])
                fcc_block = np.append(fcc_block, [i,j+0.5,k+0.5,i+0.5,j,k+0.5,i+0.5,j+0.5,k,i+1,j+0.5,k+0.5,i+0.5,j+1,k+0.5,i+0.5,j+0.5,k+1])

    fcc_block = fcc_block * const
    fcc_block = fcc_block.reshape((N_beads,3))
    fcc_block = np.unique(fcc_block, axis=0)
    fcc_block = center(fcc_block)
    return fcc_block
