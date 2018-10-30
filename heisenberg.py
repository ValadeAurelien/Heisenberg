#! /bin/env python

import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
import copy
from imp import reload

class Lattice:
    N = None
    h, k = None, None
    a = None
    lattice = None
    latcos, latsin, latinter, latenergy = None, None, None, None
    distmat = None
    M, Q, m, q, E = None, None, None, None, None

    def __init__(self, N, h, k, a):
        self.N = N
        self.h = h
        self.k = k
        self.a = a

    def init(self, initype="random", mat=None):
        if (initype=="random"):
            self.lattice = np.pi*rd.rand(self.N,self.N)
        elif (initype=="up"):
            self.lattice = np.ones((self.N, self.N))
        elif (initype=="down"):
            self.lattice = -np.ones((self.N, self.N))
        elif (initype=="null"):
            self.lattice = -np.zeros((self.N, self.N))
        elif (initype=="chessboard"):
            self.lattice = np.zeros((self.N, self.N))
            self.lattice[1::2,1::2] = 1
            self.lattice[::2,::2] = 1
        elif (initype=="copy"):
            if not mat:
                print('initialisation type is copy, yet there is no matrix to copy')
                return (-1)
            self.lattice = mat
        else:
            print("type not understood")
            return(-1)
        self.calc_latcossin()
        self.latinter = np.zeros((self.N, self.N))

    def calc_latcossin(self):
        self.latcos = np.cos(self.lattice)
        self.latsin = np.sin(self.lattice)

    def calc_dist(self, i, j, k, l):
        return ( np.power(np.sqrt((i-k)**2 + (j-l)**2) +
                          (k==i)*(l==j), -self.a) )

    def calc_distmat(self, i, j):
        self.distmat = np.fromfunction(
            lambda k, l: self.calc_dist(i, j, k, l),
            (self.N, self.N))
        self.distmat[i, j] = 0

    def calc_interactions_one_spin(self, i, j):
        self.calc_distmat(i, j)
        return( np.sum( self.distmat * ( self.latcos[i, j]*self.latcos +
                                         self.latsin[i, j]*self.latsin ) ) )

    def calc_latinter(self):
        for i in range(self.N):
            for j in range(self.N) :
                self.latinter[i, j] = self.calc_interactions_one_spin(i, j)

    def calc_latenergy(self):
        self.latenergy = self.k*self.latinter+self.h*self.latcos

    def calc_Q(self):
        self.Q = np.sum(self.latinter)

    def calc_M(self):
        self.M = np.sum(self.latcos)

    def calc_q(self):
        self.q = np.sum(self.latinter)/self.N**2

    def calc_m(self):
        self.m = np.sum(self.latcos)/self.N**2

    def calc_E(self):
        self.E = -(self.k*self.Q + self.h*self.M)

    def calc_all(self):
        self.calc_latcossin()
        self.calc_latinter()
        self.calc_M()
        self.calc_Q()
        self.calc_E()

class Heisenberg:
    N = None
    OldLat, NewLat = None, None
    initype, MCtype = None, None

    def __init__(self, N, h, k, a, initype="random", MCtype="Metropolis"):
        self.N = N
        self.initype = initype
        self.MCtype = MCtype
        self.OldLat = Lattice(N, h, k, a)

    def init(self):
        self.OldLat.init(self.initype)
        self.NewLat = copy.deepcopy(self.OldLat)
        self.OldLat.calc_all()

    def update(self):
        i = rd.randint(self.N)
        j = rd.randint(self.N)
        newangle = np.pi*rd.rand()
        acc = rd.rand()

        self.NewLat.lattice[i, j] = newangle
        self.NewLat.calc_all()

        # print(acc, self.OldLat.E, self.NewLat.E,
              # np.exp(-(self.NewLat.E-self.OldLat.E)))
        if (self.MCtype == "Metropolis"):
            if self.NewLat.E < self.OldLat.E:
                self.OldLat = copy.deepcopy(self.NewLat)
            elif acc < np.exp(-(self.NewLat.E-self.OldLat.E)):
                self.OldLat = copy.deepcopy(self.NewLat)
            else :
                self.NewLat = copy.deepcopy(self.OldLat)

        elif (self.MCtype == "Glauber"):
            if acc < 1/(1+np.exp(-(self.NewLat.E-self.OldLat.E))):
                self.OldLat = copy.deepcopy(self.NewLat)
            else:
                self.NewLat = copy.deepcopy(self.OldLat)

        else:
            print("MCtype not understood")
            return(-1)

    def many_updates(self, U):
        for i in range(U):
            self.update()

if __name__=="__main__":
    H = Heisenberg(16, 0, 1, 2, "random")
    H.init()
    im = plt.imshow(H.OldLat.lattice, vmin=0, vmax=np.pi)
    plt.show(block=False)
    for i in range(100):
        # np.savetxt("test_%02u.txt"%i, H.OldLat.lattice, delimiter="\t\t")
        print(H.OldLat.Q, H.OldLat.M, H.OldLat.E)
        H.many_updates(100)
        im.set_data(H.OldLat.lattice)
        plt.draw()
        plt.pause(0.001)
