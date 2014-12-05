import numpy as np


class Field:
    def __init__(self,nx,ny):
        self.nx = nx
        self.ny = ny
        self.R = np.zeros(nx)
        self.z = np.zeros(ny)
        self.psi = np.zeros((nx,ny))
        self.p = np.zeros((nx,ny))
