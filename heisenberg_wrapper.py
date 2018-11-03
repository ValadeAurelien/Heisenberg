#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import sys
from scipy.signal import correlate2d
from threading import Thread

class HeisenbergLauncher:
    dircty, latf, scaf = None, None, None
    argfile = None
    N		= None
    h		= None
    k		= None
    neighb	= None
    a		= None
    t_max	= None
    t_lat	= None
    t_sca	= None

    def __init__(self, dircty, latf="lattice", intf="inter",
                 scaf="scalars", ext=".dat", argfile="args"):
        self.dircty = dircty + "/"
        self.latf = latf
        self.intf = intf
        self.scaf = scaf
        self.ext = ext
        self.argfile = argfile

    def set_args(self, N=None, h=None, k=None, neighb=None, a=None,
                 t_max=None, t_lat=None, t_sca=None):
        self.N = N           if not N==None      else self.N
        self.h = h           if not h==None      else self.h
        self.k = k           if not k==None      else self.k
        self.neighb	= neighb if not neighb==None else self.neighb
        self.a = a           if not a==None      else self.a
        self.t_max = t_max   if not t_max==None  else self.t_max
        self.t_lat = t_lat   if not t_lat==None  else self.t_lat
        self.t_sca = t_sca   if not t_sca==None  else self.t_sca

    def make_or_clear_dir(self, dircty=None):
        subprocess.run(["rm", self.dircty, "-r"], stdout=None)
        subprocess.run(["mkdir", "-p",  self.dircty], stdout=None)

    def run(self):
        subprocess.run(["./heisenberg", str(self.N),
                        str(self.h), str(self.k), str(self.neighb), str(self.a),
                        str(self.t_max), str(self.t_lat), str(self.t_sca),
                        self.dircty + self.latf,
                        self.dircty + self.intf, self.dircty + self.scaf,
                        self.ext, self.dircty + self.argfile],
                       stdout=sys.stdout)

    def run_manual(self, N=None, h=None, k=None, neighb=None, a=None,
                   t_max=None, t_lat=None, t_sca=None,
                   latfname=None, intfname=None, scafname=None, ext=None,
                   argfile=None):
        subprocess.run(["./heisenberg", str(N),
                        str(h), str(k), str(neighb), str(a),
                        str(t_max), str(t_lat), str(t_sca),
                        latfname, intfname, scafname, ext, argfile],
                       stdout=sys.stdout)

    def print(self):
        print('dircty, latf, intf, scaf       = ', self.dircty,
              self.latf, self.intf, self.scaf)
        print('N                              = ', self.N)
        print('h                              = ', self.h)
        print('k                              = ', self.k)
        print('neighb                         = ', self.neighb)
        print('a                              = ', self.a)
        print('t_max                          = ', self.t_max)
        print('t_lat                          = ', self.t_lat)
        print('t_sca                          = ', self.t_sca)

class HeisenbergLauncherThreaded(HeisenbergLauncher, Thread):

    def __init__(self,  dircty, latf="lattice", intf="inter",
                 scaf="scalars", ext=".dat"):
        Thread.__init__(self)
        HeisenbergLauncher.__init__(self, dircty, latf="lattice", intf="inter",
                 scaf="scalars", ext=".dat")

    def run(self):
        HeisenbergLauncher.run(self)

class HeisenbergSimpleSampling:
    rootdircty, latf, scaf = None, None, None
    N		= None
    neighb	= None
    a		= None
    n_trials = None

    def __init__(self, rootdircty, latf="lattice", intf="inter",
                 scaf="scalars", ext=".dat", n_trials=None):
        self.rootdircty = rootdircty + "/"
        self.latf = latf
        self.intf = intf
        self.scaf = scaf
        self.ext = ext
        self.n_trials = n_trials

    def set_args(self, N=None, neighb=None, a=None,
                 n_trials=None):
        self.N = N           if not N==None      else self.N
        self.neighb	= neighb if not neighb==None else self.neighb
        self.a = a           if not a==None      else self.a

    def run(self):
        for n in range(self.n_trials):
            HL = HeisenbergLauncher(self.rootdircty + "trial_" + str(n),
                                    self.latf, self.intf, self.scaf, self.ext)
            HL.make_or_clear_dir()
            HL.set_args(self.N, h=0, k=0, neighb=self.neighb, a=self.a,
                        t_max=0, t_lat=0, t_sca=0)
            HL.run()

class HeisenbergManyLaunches:
    rootdircty, latf, scaf = None, None, None
    N		= None
    h		= None
    k		= None
    neighb	= None
    a		= None
    t_max	= None
    t_lat	= None
    t_sca	= None
    parallel = None
    HLs = []

    def __init__(self, rootdircty, latf="lattice", intf="inter",
                 scaf="scalars", ext=".dat", parallel=False):
        self.rootdircty = rootdircty + "/"
        self.latf = latf
        self.intf = intf
        self.scaf = scaf
        self.ext = ext
        self.parallel = parallel

    def set_args(self, N=None, h=None, k=None, neighb=None, a=None,
                 t_max=None, t_lat=None, t_sca=None):
        self.N = N           if not N==None      else self.N
        self.h = h           if not h==None      else self.h
        self.k = k           if not k==None      else self.k
        self.neighb	= neighb if not neighb==None else self.neighb
        self.a = a           if not a==None      else self.a
        self.t_max = t_max   if not t_max==None  else self.t_max
        self.t_lat = t_lat   if not t_lat==None  else self.t_lat
        self.t_sca = t_sca   if not t_sca==None  else self.t_sca

    def set_trials(self, n_trials):
        for n in range(n_trials):
            if (self.parallel):
                self.HLs.append(HeisenbergLauncherThreaded(self.rootdircty + "trial_" + str(n),
                                                           self.latf, self.intf, self.scaf, self.ext))
            else:
                self.HLs.append(HeisenbergLauncher(self.rootdircty + "trial_" + str(n),
                                                   self.latf, self.intf, self.scaf, self.ext))
                self.HLs[-1].set_args(self.N, self.h, self.k, self.neighb,
                                      self.a, self.t_max, self.t_lat, self.t_sca)
            self.HLs[-1].make_or_clear_dir()

    def run(self):
        for HL in self.HLs:
            if (self.parallel):
                HL.start()
            else:
                HL.run()
        if (self.parallel):
            for HL in self.HLs:
                HL.join()

class HeisenbergSweeps:
    rootdircty, latf, scaf = None, None, None
    N		= None
    h		= None
    k		= None
    neighb	= None
    a		= None
    t_max	= None
    t_lat	= None
    t_sca	= None
    parallel = None
    HLs = []

    def __init__(self, rootdircty, latf="lattice", intf="inter",
                 scaf="scalars", ext=".dat", parallel=False):
        self.rootdircty = rootdircty + "/"
        self.latf = latf
        self.intf = intf
        self.scaf = scaf
        self.ext = ext
        self.parallel = parallel

    def set_args(self, N=None, h=None, k=None, neighb=None, a=None,
                 t_max=None, t_lat=None, t_sca=None):
        self.N = N           if not N==None      else self.N
        self.h = h           if not h==None      else self.h
        self.k = k           if not k==None      else self.k
        self.neighb	= neighb if not neighb==None else self.neighb
        self.a = a           if not a==None      else self.a
        self.t_max = t_max   if not t_max==None  else self.t_max
        self.t_lat = t_lat   if not t_lat==None  else self.t_lat
        self.t_sca = t_sca   if not t_sca==None  else self.t_sca

    def set_sweeps(self, minh=None, maxh=None, nh=None,
                   mink=None, maxk=None, nk=None):
        hs = np.linspace(minh, maxh, nh)
        ks = np.linspace(mink, maxk, nk)
        for h in range(nh):
            for k in range(nk):
                if (self.parallel):
                    self.HLs.append(HeisenbergLauncherThreaded(self.rootdircty +
                                                               "h%ik%i"%(h, k),
                                                               self.latf, self.intf,
                                                               self.scaf, self.ext))
                else:
                    self.HLs.append(HeisenbergLauncher(self.rootdircty +
                                                       "h%ik%i"%(h, k),
                                                       self.latf, self.intf,
                                                       self.scaf, self.ext))
                self.HLs[-1].set_args(self.N, hs[h], ks[k], self.neighb,
                                        self.a, self.t_max, self.t_lat, self.t_sca)
                self.HLs[-1].make_or_clear_dir()

    def run(self):
        for HL in self.HLs:
            if (self.parallel):
                HL.start()
            else:
                HL.run()
        if (self.parallel):
            for HL in self.HLs:
                HL.join()


def avg_2p2(X):
    if len(X)%2:
        X = X[1:]
    return(len(X)/2, (X[0::2] + X[1::2])/2)

def prod_shift(tab, i, j):
    return ( tab * np.roll(np.roll(tab, i, axis=0), j, axis=1) )

class HeisenbergPlot:
    rootdircty, rootlatf, rootintf, rootscaf = None, None, None, None
    dircty, latf, intf, scaf = None, None, None, None
    fig_sca1, fig_sca2 = None, None
    n_lat = None
    fig_lat, fig_int, fig_stat = None, None, None
    ax_sca1, ax_sca2 = None, None
    ax_lat, ax_lat, ax_stat = None, None, None
    plot_sca1, plot_sca2 = None, None
    plot_lat, plot_int, plot_stat = None, None, None
    show_corr = True
    scalars = None
    lattice, inter = None, None
    ft2d_lat = None
    header_scalars, header_lattice, header_inter = None, None, None
    many_scalars, many_stabs = None, None
    stats = None
    corr_x, corr_y = None, None
    stats_mean_Ms, stats_mean_Qs, stats_mean_Es = [], [], []
    stats_std_Ms, stats_std_Qs, stats_std_Es = [], [], []
    stab_crit = .005
    serr_avgs = None
    png = None
    argfile = None
    N		= None
    h		= None
    k		= None
    neighb	= None
    a		= None
    t_max	= None
    t_lat	= None
    t_sca	= None

    def __init__(self, rootdircty, rootlatf="lattice", rootintf="inter",
                 rootscaf="scalars", ext=".dat", argfile="args", png=None):
        self.rootdircty = rootdircty + "/"
        self.rootlatf = rootlatf
        self.rootintf = rootintf
        self.rootscaf = rootscaf
        self.ext = ext
        self.dircty = self.rootdircty
        self.latf =   self.rootlatf
        self.intf =   self.rootintf
        self.scaf =   self.rootscaf
        self.argfile = argfile
        self.png = png
        try:
            self.read_args()
        except:
            pass


    def read_args(self):
        fargs = open(self.dircty + self.argfile).readlines()
        args = [ line.split() for line in fargs ]
        self.N		= int(args[0][2])
        self.h		= float(args[1][2])
        self.k		= float(args[2][2])
        self.neighb	= int(args[3][2])
        self.a		= int(args[4][2])
        self.t_max	= int(args[5][2])
        self.t_lat	= int(args[6][2])
        self.t_sca	= int(args[7][2])

    def read_lattice(self):
        flat = open(self.dircty + self.latf + "%016u"%(self.n_lat) + self.ext)
        header = flat.readline().split()
        self.header = [int(header[1]), int(header[2]), int(header[3]),
                       float(header[4]), float(header[5]),
                       float(header[6]), float(header[7])];
        flat.close()
        self.lattice = np.loadtxt(self.dircty + self.latf
                                  + "%016u"%(self.n_lat) + self.ext).transpose()

    def plot_lattice(self):
        if (self.plot_lat == None ):
            self.fig_lat = plt.figure("lattice")
            self.ax_lat = self.fig_lat.add_subplot(111)
            self.plot_lat = self.ax_lat.imshow(self.lattice,
                                               interpolation="bicubic",
                                               vmin=0, vmax=np.pi)
            plt.show(block=False)
        else:
            self.plot_lat.set_data(self.lattice)

    def show_lattice(self):
        self.fig_lat.canvas.draw()
        self.fig_lat.canvas.flush_events()
        if (self.png):
            fig_lat.savefig("lattice_" + self.rootdircty +
                            str(self.n_lat) + ".png")

    def refresh_n_show_lattice(self):
        self.read_lattice()
        self.plot_lattice()
        self.ax_lat.set_title("n=%u"%(self.n_lat))
        self.show_lattice()

    def plot_lattices(self, begin=0, end=0, t_lat=0, stype=None):
        if (stype=="last"):
            begin = self.t_max - self.t_lat
            end = self.t_max - self.t_lat+1
            t_lat = 1
        if (stype=="all"):
            begin = 0
            end = self.t_max - self.t_lat
            t_lat = self.t_lat
        for i in range(begin, end, t_lat):
            self.n_lat = i
            self.refresh_n_show_lattice()
            plt.pause(.2)

    def read_inter(self, n):
        fint = open(self.dircty + self.intf + "%016u"%(n) + self.ext)
        header = fint.readline().split()
        self.header = [int(header[1]), int(header[2]), int(header[3]),
                       float(header[4]), float(header[5]),
                       float(header[6]), float(header[7])];
        fint.close()
        self.inter = np.loadtxt(self.dircty + self.intf
                                  + "%016u"%(n) + self.ext).transpose()

    def plot_inter(self):
        if (self.plot_int == None ):
            self.fig_int = plt.figure("inter")
            self.ax_int = self.fig_int.add_subplot(111)
            self.plot_int = self.ax_int.imshow(-self.inter, vmin=0, vmax=np.pi)
            plt.show(block=False)
        else:
            self.plot_int.set_data(-self.inter)

    def show_inter(self):
        self.fig_int.canvas.draw()
        self.fig_int.canvas.flush_events()

    def refresh_n_show_inter(self, n):
        self.read_inter(n)
        self.plot_inter()
        self.ax_int.set_title("n=%u"%(n))
        self.show_inter()

    def plot_inters(self, begin=0, end=0, t_int=0):
        for i in range(begin, end, t_int):
            self.refresh_n_show_inter(i)
            plt.pause(.2)

    def read_scalars(self):
        fsca = open(self.dircty + self.scaf + self.ext)
        header = fsca.readline().split()
        self.header = [int(header[1]), int(header[2]),
                       float(header[3]), float(header[4]),
                       self.dircty];
        fsca.close()
        self.scalars = np.loadtxt(self.dircty + self.scaf
                                  + self.ext).transpose()
        if (self.scalars.shape == (4,)):
            self.scalars = np.array([[x,x] for x in self.scalars])

        self.scalars[1:] /= int(self.header[0])**2
        try :
            self.stab = np.where(np.abs((self.scalars[3] - self.scalars[3][-1])/
                                (self.scalars[3][-1]-self.scalars[3][0]))<self.stab_crit)[0][0]
            if len(self.scalars[0]) - self.stab < 2:
                self.stab = 0
        except:
            self.stab = 0
        print("Read scalars :\n"
              "h, k   : %.2f %.2f \n"
              "Length : %i\n"
              "N grid : %i\n"
              "Stab crit : %f\n"
              "Stab   : %i"%(self.header[2], self.header[3],
                             len(self.scalars[0]), int(self.header[0]),
                             self.stab_crit, self.stab))

    def plot_scalars(self):
        if (self.plot_sca1 == None or self.plot_sca2 == None):
            self.fig_sca1, self.ax_sca1 = plt.subplots(3, 1)
            self.fig_sca1.suptitle("scalars")
            self.plot_sca1 = [[], [], []]
            self.plot_sca1[0] = self.ax_sca1[0].plot(self.scalars[0],
                                                     self.scalars[1], "b",
                                                     label="M")
            self.plot_sca1[1] = self.ax_sca1[1].plot(self.scalars[0],
                                                     self.scalars[2], "r",
                                                     label="Q")
            self.plot_sca1[2] = self.ax_sca1[2].plot(self.scalars[0],
                                                     self.scalars[3], "g",
                                                     label="E")
            # self.ax_sca1[2].plot([self.scalars[0, self.stab],
                                  # self.scalars[0, self.stab]],
                                 # [np.min(self.scalars[1]),
                                  # np.max(self.scalars[1])], "k")
            self.fig_sca1.legend()
            self.fig_sca2, self.ax_sca2 = plt.subplots(1, 1)
            self.fig_sca2.suptitle("scalars")
            self.ax_sca2.set_xlabel("M")
            self.ax_sca2.set_ylabel("Q")
            self.plot_sca2 = self.ax_sca2.scatter(self.scalars[1, self.stab:],
                                                  self.scalars[2, self.stab:],
                                                  label="M/Q")
            self.fig_sca2.legend()
        else:
            self.plot_sca1[0].set_data(self.scalars[0], self.scalars[1])
            self.plot_sca1[0].set_data(self.scalars[0], self.scalars[2])
            self.plot_sca1[0].set_data(self.scalars[0], self.scalars[3])

    def show_scalars(self):
        self.fig_sca1.canvas.draw()
        self.fig_sca1.canvas.flush_events()
        self.fig_sca2.canvas.draw()
        self.fig_sca2.canvas.flush_events()

        if (self.png):
            self.fig_sca1.savefig("sca1_" + self.rootdircty.replace("/", "_")[:-2]
                                  + ".png")
            self.fig_sca2.savefig("sca2_" + self.rootdircty.replace("/", "_")[:-2]
                                  + ".png")

    def refresh_n_show_scalars(self):
        self.read_scalars()
        self.plot_scalars()
        self.show_scalars()

    def calc_stat_err(self):
        serr_avgs = [[], [], []]
        for i in range(3):
            sca = self.scalars[1+i, self.stab:]
            ns, n = [], len(sca)
            while(len(sca)>1) :
                serr_avgs[i].append(np.std(sca)/np.sqrt(len(sca)))
                ns.append(n)
                n, sca = avg_2p2(sca)
        n0 = np.argmax(serr_avgs[0])
        n1 = np.argmax(serr_avgs[1])
        n2 = np.argmax(serr_avgs[2])
        n  = max((n0, n1, n2))
        serr = [savg[n] for savg in serr_avgs]
        fig_serr, ax_serr = plt.subplots(3, 1)
        ax_serr[0].plot(np.log2(ns), serr_avgs[0], "b", label="stat err M")
        ax_serr[1].plot(np.log2(ns), serr_avgs[1], "r", label="stat err Q")
        ax_serr[2].plot(np.log2(ns), serr_avgs[2], "g", label="stat err E")
        self.stat['serr'] = serr
        self.stat['serrn'] = 2**(np.log2(len(self.scalars[0])-self.stab) - n -1)

    def calc_stat(self):
        self.stat = {}
        self.stat['mean'] = np.mean(self.scalars[1:, self.stab:],
                                    axis=1)
        self.stat['std'] = np.std(self.scalars[1:, self.stab:],
                                  axis=1)
        self.calc_stat_err()
        print("Statistics : \n"
              "means  : \t%-8g \t\t%-8g \t\t%-8g \n"
              "std    : \t%-8g \t\t%-8g \t\t%-8g \n"
              "serr   : \t%-8g \t\t%-8g \t\t%-8g \n"
              "serrn  : \t%-8i \n"
              %(self.stat['mean'][0], self.stat['mean'][1], self.stat['mean'][2],
                self.stat['std'][0], self.stat['std'][1], self.stat['std'][2],
                self.stat['serr'][0], self.stat['serr'][1], self.stat['serr'][2],
                self.stat['serrn']))
        print(self.header)

#     def correlation_lattice(self):
        # self.ft2d_lat = np.fft.rfft2(self.lattice)
        # self.cor2d_lat = np.fft.irfft2(np.abs(self.ft2d_lat)**2)
        # plt.figure()
        # print(np.max(self.cor2d_lat), np.min(self.cor2d_lat), np.average(self.cor2d_lat), np.std(self.cor2d_lat))
        # plt.imshow(np.log(np.abs(self.cor2d_lat - np.average(self.cor2d_lat))))
        # plt.imshow(np.log(np.abs(np.fft.fftshift(np.fft.rfft2(self.cor2d_lat), axes=0))))
        # plt.show()

    # def ft_lattice(self):
        # self.ft2d_lat = np.log(np.abs(np.fft.rfft2(self.lattice)))
        # Y, X = np.meshgrid(np.arange(self.ft2d_lat.shape[1]),
                           # 1j*np.arange(self.ft2d_lat.shape[0]))
        # self.ft2d_r = np.abs( X + Y )
        # allrad = []
        # for i in range(int(self.ft2d_lat.shape[0]/2)):
            # _allrad = [ (self.ft2d_r[i, j],
                         # (self.ft2d_lat[i, j] + self.ft2d_lat[-i, j])/2.)
                       # for j in range(self.ft2d_lat.shape[1]) ]
            # allrad += _allrad
        # allrad_r = [ item[0] for item in allrad ]
        # allrad_f = [ item[1] for item in allrad ]
        # interp_r = np.linspace(0, self.ft2d_lat.shape[1],
                               # self.ft2d_lat.shape[1])
        # interp_f = np.interp(interp_r, allrad_r, allrad_f)
        # plt.figure()
        # # plt.plot(allrad_r, allrad_f)
        # # plt.plot(interp_r, interp_f)
        # plt.plot(np.fft.ifft(interp_f))
        # plt.show()

    def correlation_function(self):
        corr_x, corr_y = [], []
        for i in range(int(self.N/4)):
            print("\r %i %%"%(int(i/int(self.N/4)*100)), end="")
            for j in range(int(self.N/4)):
                corr_x.append(np.abs(np.sqrt(i*i + j*j)))
                corr_y.append(np.sum( prod_shift(self.lattice, i, j) +
                                     prod_shift(self.lattice, -i, j) +
                                     prod_shift(self.lattice, i, -j) +
                                     prod_shift(self.lattice, -i, -j) ) /
                              (4*self.N**2))
        asort = np.argsort(corr_x)
        corr_x = np.array(corr_x)[asort]
        corr_y = np.array(corr_y)[asort]
        corr_x_range = corr_x<10
        self.corr_x = corr_x[corr_x_range]
        self.corr_y = corr_y[corr_x_range]
        expfit, eres, _, _, _ = np.polyfit(corr_x, np.log(corr_y), 1, full=True)
        powfit, pres, _, _, _ = np.polyfit(corr_x, 1/corr_y, 1, full=True)
        print("expfit : %.2fx+%.2f, R^2 = %f\n"
              "powfit : %.2fx+%.2f, R^2 = %f\n"%
              (expfit[0], expfit[1], eres,
               powfit[0], powfit[1], pres))
        if (self.show_corr):
            fig_corr, ax_corr = plt.subplots(3, 1)
            ax_corr[0].plot(corr_x, corr_y, "b", label="correlation")
            ax_corr[1].plot(corr_x, np.log(corr_y), "g", label="log(corr)")
            ax_corr[1].plot(corr_x, expfit[0]*corr_x + expfit[1], "grey")
            ax_corr[2].plot(corr_x, 1/corr_y-1, "r", label="1/corr-1")
            ax_corr[2].plot(corr_x, powfit[0]*corr_x + powfit[1], "grey")
            fig_corr.legend()

            if (self.png):
                fig_corr.savefig("corr_" + self.rootdircty.replace("/", "_")[:-2]
                                 + ".png")

    def plot_hsweep_vs_h(self):
        ks = np.unique( np.sort( [ head[3] for head in self.stats ] ) )
        hs_x = np.unique( np.sort( [ head[2] for head in self.stats ] ) )
        hs_y_mean_Ms= [[0 for h in hs_x] for h in hs]
        hs_y_mean_Qs= [[0 for h in hs_x] for h in hs]
        hs_y_mean_Es= [[0 for h in hs_x] for h in hs]
        hs_y_std_Ms= [[0 for h in hs_x] for h in hs]
        hs_y_std_Qs= [[0 for h in hs_x] for h in hs]
        hs_y_std_Es= [[0 for h in hs_x] for h in hs]
        hs_y_serr_Ms= [[0 for h in hs_x] for h in hs]
        hs_y_serr_Qs= [[0 for h in hs_x] for h in hs]
        hs_y_serr_Es= [[0 for h in hs_x] for h in hs]
        for head, stat in self.stats.items() :
            i = np.where( ks==head[3] )[0][0]
            j = np.where( hs_x==head[2] )[0][0]
            hs_y_mean_Ms[i][j] = stat['mean'][0]
            hs_y_mean_Qs[i][j] = stat['mean'][1]
            hs_y_mean_Es[i][j] = stat['mean'][2]
            hs_y_std_Ms[i][j] = stat['std'][0]
            hs_y_std_Qs[i][j] = stat['std'][1]
            hs_y_std_Es[i][j] = stat['std'][2]
            hs_y_serr_Ms[i][j] = stat['serr'][0]
            hs_y_serr_Qs[i][j] = stat['serr'][1]
            hs_y_serr_Es[i][j] = stat['serr'][2]
        self.fig_hsweep, self.ax_hsweep = plt.subplots(3, 1)
        self.fig_hsweep.suptitle("h-sweep")
        self.plot_hsweep = [[], [], []]
        for i in range(len(hs)):
            print(hs_y_mean_Qs[i])
            # self.plot_hsweep[0].append(self.ax_hsweep[0].plot(hs_x,
                                                              # hs_y_mean_Ms[i],
                                                              # label="h = %.1f"%(hs[i])))
            # self.plot_hsweep[1].append(self.ax_hsweep[1].plot(hs_x,
                                                              # hs_y_mean_Qs[i],
                                                              # label="h = %.1f"%(hs[i])))
            # self.plot_hsweep[2].append(self.ax_hsweep[2].plot(hs_x,
                                                              # hs_y_mean_Es[i],
                                                              # hs_y_serr_Es[i],
                                                              # label="h = %.1f"%(hs[i])))
            self.plot_hsweep[0].append(self.ax_hsweep[0].errorbar(hs_x,
                                                                  hs_y_mean_Ms[i],
                                                                  hs_y_serr_Ms[i],
                                                                  label="h = %.1f"%(hs[i])))
            self.plot_hsweep[1].append(self.ax_hsweep[1].errorbar(hs_x,
                                                                  hs_y_mean_Qs[i],
                                                                  hs_y_serr_Qs[i],
                                                                  label="h = %.1f"%(hs[i])))
            self.plot_hsweep[2].append(self.ax_hsweep[2].errorbar(hs_x,
                                                                  hs_y_mean_Es[i],
                                                                  hs_y_serr_Es[i],
                                                                  label="h = %.1f"%(hs[i])))

            if (self.png):
                fig_hsweep.savefig("hsweep_" + self.rootdircty.replace("/", "_")[:-2]
                                   + ".png")

    def plot_ksweep_vs_h(self):
        hs = np.unique( np.sort( [ head[2] for head in self.stats ] ) )
        ks_x = np.unique( np.sort( [ head[3] for head in self.stats ] ) )
        ks_y_mean_Ms= [[0 for k in ks_x] for h in hs]
        ks_y_mean_Qs= [[0 for k in ks_x] for h in hs]
        ks_y_mean_Es= [[0 for k in ks_x] for h in hs]
        ks_y_std_Ms= [[0 for k in ks_x] for h in hs]
        ks_y_std_Qs= [[0 for k in ks_x] for h in hs]
        ks_y_std_Es= [[0 for k in ks_x] for h in hs]
        ks_y_serr_Ms= [[0 for k in ks_x] for h in hs]
        ks_y_serr_Qs= [[0 for k in ks_x] for h in hs]
        ks_y_serr_Es= [[0 for k in ks_x] for h in hs]
        for head, stat in self.stats.items() :
            i = np.where( hs==head[2] )[0][0]
            j = np.where( ks_x==head[3] )[0][0]
            ks_y_mean_Ms[i][j] = stat['mean'][0]
            ks_y_mean_Qs[i][j] = stat['mean'][1]
            ks_y_mean_Es[i][j] = stat['mean'][2]
            ks_y_std_Ms[i][j] = stat['std'][0]
            ks_y_std_Qs[i][j] = stat['std'][1]
            ks_y_std_Es[i][j] = stat['std'][2]
            ks_y_serr_Ms[i][j] = stat['serr'][0]
            ks_y_serr_Qs[i][j] = stat['serr'][1]
            ks_y_serr_Es[i][j] = stat['serr'][2]
        self.fig_ksweep, self.ax_ksweep = plt.subplots(3, 1)
        self.fig_ksweep.suptitle("k-sweep")
        self.plot_ksweep = [[], [], []]
        for i in range(len(hs)):
            print(ks_y_mean_Qs[i])
            # self.plot_ksweep[0].append(self.ax_ksweep[0].plot(ks_x,
                                                              # ks_y_mean_Ms[i],
                                                              # label="h = %.1f"%(hs[i])))
            # self.plot_ksweep[1].append(self.ax_ksweep[1].plot(ks_x,
                                                              # ks_y_mean_Qs[i],
                                                              # label="h = %.1f"%(hs[i])))
            # self.plot_ksweep[2].append(self.ax_ksweep[2].plot(ks_x,
                                                              # ks_y_mean_Es[i],
                                                              # ks_y_serr_Es[i],
                                                              # label="h = %.1f"%(hs[i])))
            self.plot_ksweep[0].append(self.ax_ksweep[0].errorbar(ks_x,
                                                                  ks_y_mean_Ms[i],
                                                                  ks_y_serr_Ms[i],
                                                                  label="h = %.1f"%(hs[i])))
            self.plot_ksweep[1].append(self.ax_ksweep[1].errorbar(ks_x,
                                                                  ks_y_mean_Qs[i],
                                                                  ks_y_serr_Qs[i],
                                                                  label="h = %.1f"%(hs[i])))
            self.plot_ksweep[2].append(self.ax_ksweep[2].errorbar(ks_x,
                                                                  ks_y_mean_Es[i],
                                                                  ks_y_serr_Es[i],
                                                                  label="h = %.1f"%(hs[i])))

            if (self.png):
                fig_ksweep.savefig("ksweep_" + self.rootdircty.replace("/", "_")[:-2]
                                   + ".png")




    def analyse(self):
        self.many_scalars, self.many_stabs, self.stats = [], [], {}
        subdirs = [x[0] for x in os.walk(self.rootdircty)]
        print(subdirs)
        for subdir in subdirs :
            if not len(subdir) - len(self.rootdircty) < 2:
                self.dircty =  subdir + "/"
                print("Processing : " + self.dircty)
                self.read_args()
                self.read_scalars()
                self.many_stabs.append(self.stab)
                self.many_scalars.append(self.scalars)
                self.calc_stat()
                self.stats[tuple(self.header)] = self.stat

    def plot_analyse(self):
        if (self.plot_stat == None ):
            self.fig_stat, self.ax_stat = plt.subplots(1,1)
            self.fig_stat.suptitle("statistics")
            self.plot_stat = [[], [], []]
            for sca, stab in zip(self.many_scalars, self.many_stabs):
                if self.stab != 0:
                    step=int(np.ceil(len(sca[1, stab:])/100))
                    self.plot_stat[0] = self.ax_stat.scatter(sca[1, stab::step],
                                                             sca[2, stab::step])
                    if (self.png):
                        fig.savefig("MQ_" + self.rootdircty.replace("/", "_")[:-2]
                                    + ".png")

if __name__ == "__main__":
    import argparse as ap
    parser = ap.ArgumentParser(description='Heisenberg analyse tools.')
    parser.add_argument('action', help='What to do ? launch (one simulation)'
                        ', sweep (many simulations), simple (simple sampling)'
                        'or plot results')
    parser.add_argument('datadir', help='Directory to write or read the data')
    parser.add_argument('-N', type=int,
                        help='size of the lattice')
    parser.add_argument('-hc', type=float,
                        help='magnetic field strength')
    parser.add_argument('-kc', type=float,
                        help='magnetic field strength')
    parser.add_argument('-n', type=int,
                        help='number of neighbours to consider')
    parser.add_argument('-a', type=int,
                        help='power of distance coupling (int)')
    parser.add_argument('-tmax', type=int,
                        help='nb of MC sweeps')
    parser.add_argument('-tlat', type=int,
                        help='save lattice every tlat sweeps')
    parser.add_argument('-tsca', type=int,
                        help='save scalars every tsca sweeps')
    parser.add_argument('-hsweep', type=int, nargs='*',
                        help='if launch : h sweep begin, end, nsteps\n')
    parser.add_argument('-ksweep', type=int, nargs='*',
                        help='if launch : k sweep begin, end, nsteps\n')
    parser.add_argument('-par', action='count',
                        help='parallel or sequential sweeps')
    parser.add_argument('-ntrials', type=int,
                        help='number of trials for simple sampling or'
                        'important samplings')
    parser.add_argument('-sca', action='count',
                        help='plot scalars')
    parser.add_argument('-lat', nargs='*',
                        help='plot lattices and begin, end, tlat')
    parser.add_argument('-int', nargs='*',
                        help='plot interaction matrixes and begin, end, tlat')
    parser.add_argument('-corr', action='count',
                        help='plot fourrier transform')
    parser.add_argument('-phsweep', action='count',
                        help='plot hsweep')
    parser.add_argument('-pksweep', action='count',
                        help='plot ksweep')
    parser.add_argument('-ana', action='count',
                        help='plot whole analyse')
    parser.add_argument('-png', action='count',
                        help='create pngs instead of plotting')

    args = vars(parser.parse_args())
    if (args['action'] == "launch"):
        if not args.get('ntrials'):
            HL = HeisenbergLauncher(args["datadir"],
                                    latf="lattice",
                                    scaf="scalars",
                                    ext=".dat")
            HL.set_args(N=args.get('N'), h=args.get('hc'), k=args.get('kc'),
                        neighb=args.get('n'), a=args.get('a'),
                        t_max=args.get('tmax'), t_lat=args.get('tlat'),
                        t_sca=args.get('tsca'))
            HL.make_or_clear_dir()
            HL.run()
        else :
            HML = HeisenbergManyLaunches(args["datadir"],
                                         latf="lattice",
                                         scaf="scalars",
                                         ext=".dat",
                                         parallel=args.get('par'))
            HML.set_args(N=args.get('N'), h=args.get('hc'), k=args.get('kc'),
                        neighb=args.get('n'), a=args.get('a'),
                        t_max=args.get('tmax'), t_lat=args.get('tlat'),
                        t_sca=args.get('tsca'))
            HML.set_trials(args.get('ntrials'))
            HML.run()
    elif(args['action'] == "sweep"):
        if not (args.get('hsweep') and args.get('ksweep')):
            print("Missing ksweep and/or hsweep args")
        HS = HeisenbergSweeps(args["datadir"],
                              latf="lattice",
                              scaf="scalars",
                              ext=".dat",
                              parallel=args.get('par'))
        HS.set_args(N=args.get('N'), h=args.get('hc'), k=args.get('kc'),
                    neighb=args.get('n'), a=args.get('a'),
                    t_max=args.get('tmax'), t_lat=args.get('tlat'),
                    t_sca=args.get('tsca'));

        if (args.get('hsweep') and len(args['hsweep'])==3):
            hs = args['hsweep']
        else :
            hs = [args.get('hc'), args.get('hc'), 1]
        if (args.get('ksweep') and len(args['ksweep'])==3):
            ks = args['ksweep']
        else :
            ks = [args.get('kc'), args.get('kc'), 1]
        HS.set_sweeps(minh=hs[0], maxh=hs[1], nh=hs[2],
                      mink=ks[0], maxk=ks[1], nk=ks[2])
        HS.run()
    elif(args['action'] == "simple"):
        HSS = HeisenbergSimpleSampling(args["datadir"],
                                       latf="lattice",
                                       scaf="scalars",
                                       ext=".dat",
                                       n_trials = args['ntrials'])
        HSS.set_args(N=args.get('N'), neighb=args.get('n'), a=args.get('a'))
        HSS.run()
    elif(args['action'] == "plot"):
        HP = HeisenbergPlot(args["datadir"],
                            rootlatf="lattice",
                            rootscaf="scalars",
                            ext=".dat",
                            png=args.get('png'))
        if (args.get('phsweep')):
            HP.analyse()
            HP.plot_hsweep_vs_k()
        if (args.get('pksweep')):
            HP.analyse()
            HP.plot_ksweep_vs_h()
        if (args.get('ana')):
            HP.analyse()
            HP.plot_analyse()
        else:
            if(args.get('sca')):
                HP.refresh_n_show_scalars()
                HP.calc_stat()
            if (args.get('lat')):
                if (len(args.get('lat')) == 3):
                    HP.plot_lattices(begin=int(args.get('lat')[0]),
                                     end=int(args.get('lat')[1]),
                                     t_lat=int(args.get('lat')[2]))
                else:
                    if (args.get('lat')[0] == "last"):
                        HP.plot_lattices(stype="last")
                    elif (args.get('lat')[0] == "all"):
                        HP.plot_lattices(stype="all")
                    else:
                        HP.plot_lattices(begin=int(args.get('lat')[0]),
                                         end=int(args.get('lat')[0])+1,
                                         t_lat=1)
            if (args.get('int')):
                if (len(args.get('int')) == 3):
                    HP.plot_inters(begin=args.get('int')[0],
                                    end=args.get('int')[1],
                                    t_int=args.get('int')[2])
                else:
                    HP.plot_inters(begin=args.get('int')[0],
                                    end=args.get('int')[0]+1,
                                    t_int=1)
            if(args.get('corr')):
                HP.correlation_function()
        plt.show()

