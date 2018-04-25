import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy as sp
from scipy.integrate import simps
from scipy import sin, cos, pi, exp, tanh, log
from scipy.special import erf, erfc
from scipy.optimize import curve_fit
from scipy import stats

import os
import os.path
import csv
import pickle

plt.rc('text', usetex=True)
plt.rc('font', family='Times')
params = {'text.usetex': True,
          'font.size': 8,
          'font.family': 'Times',
          'text.latex.unicode': True,
          'legend.fontsize': 4,
          'legend.handlelength': 3}
plt.rcParams.update(params)

def derivO4(F, h):
    n = len(F)
    dF = np.zeros(n)
    for ii in range(0, n):
        if ii == 0:
            dF[ii] = F[ii + 1] - F[ii]
            dF[ii] /= h
        elif ii == 1:
            dF[ii] = F[ii + 1] - F[ii - 1]
            dF[ii] /= 2. * h
        elif ii == (n - 2):
            dF[ii] = F[ii + 1] - F[ii - 1]
            dF[ii] /= 2. * h
        elif ii == (n - 1):
            dF[ii] = F[ii] - F[ii - 1]
            dF[ii] /= h
        else:
            dF[ii] = F[ii - 2] - 8. * F[ii - 1] + 8. * F[ii + 1] - F[ii + 2]
            dF[ii] /= 12. * h
    return dF

def smooth(F, box_pts):
    box = np.ones(box_pts) / box_pts
    F_smooth = np.convolve(F, box, mode='same')
    return F_smooth

class Statistics(object):
    def __init__(self, nphi, filename):
        self.nphi = nphi
        with open(filename + 'budget', 'r') as f:
            data = list(csv.reader(f, delimiter=' ', skipinitialspace=True,
                                   quoting=csv.QUOTE_NONNUMERIC))
        data = np.transpose(data)
        self.t_energy = data[0]
        self.ek = data[1]
        self.dek = data[2]
        self.ep = data[3]
        self.dep = data[4]
        del data
        # self.dek[0] = 0.
        # self.dep[0] = 0.
        # self.ek[0] = 0.
        n = self.t_energy.size
        self.Idek = np.zeros(n)
        self.Idep = np.zeros(n)
        for t in range(n):
            self.Idek[t] = np.trapz(self.dek[0:t], self.t_energy[0:t])
            self.Idep[t] = np.trapz(self.dep[0:t], self.t_energy[0:t])
        self.total = self.ek + self.ep + self.Idek + self.Idep

        with open(filename + 'statistics', 'r') as f:
            data = list(csv.reader(f, delimiter=' ', skipinitialspace=True,
                                   quoting=csv.QUOTE_NONNUMERIC))
        data = np.transpose(data)
        counter = 0
        self.t = data[counter]
        counter += 1
        self.mp = []
        for n in range(nphi):
            self.mp.append(data[counter])
            counter += 1
        del data

transp = False
formt = 'png'
qual = 800
siz1 = 2.25
siz2 = 2.25

#Lx=10
Ly=57.
#Lz=10
#nx=600
ny=193
#nz=600
simulation = {}
simulation["N500"] = Statistics(1,'./out/')
simulation["N500"].style = 'solid'
simulation["N500"].color = 'black'
simulation["N500"].width = 0.7

############
# Budget
############
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.xlabel(r'$t$')
plt.ylabel(r'$E_T$')
#plt.xlim(0., 40.)
#plt.ylim(0., 1.)
for key, data in simulation.items():
    et0 = data.total[0]
    plt.plot(data.t_energy,data.total/et0, label='total', linestyle='solid', color='black', linewidth=1.)
    plt.plot(data.t_energy,data.ek/et0, label='ek', linestyle='solid', color='blue', linewidth=1.)
    plt.plot(data.t_energy,data.ep/et0, label='ep', linestyle='solid', color='red', linewidth=1.)
    plt.plot(data.t_energy,(data.Idek+data.Idep)/et0, label='diss', linestyle='dotted', color='black', linewidth=1.)

#plt.axhline(y=1., xmin=0., xmax=1000, linewidth=0.5, color='black', linestyle='dotted')
plt.grid(True)
plt.minorticks_on()
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
plt.legend(loc='center right', shadow=False, fontsize=5.)
plt.savefig('energy_total' + '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close('all')

############
# Mass control
############
fig, ax1 = plt.subplots()
fig.set_size_inches(siz1, siz2)
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$m_p/m_{p0}$')
#ax1.set_ylim(0.0,1.05)
#ax1.set_xlim(0,40.)
for key, data in simulation.items():
    ax1.plot(data.t,((data.mp[0]/data.mp[0][0])), label=key, linestyle=data.style, color=data.color, linewidth=data.width+0.2)
    print('mp/mp0 final',key,(data.mp[0][-1]/data.mp[0][0])*100. -100.,"%")

#ax1.axhline(y=1., xmin=0., xmax=1000, linewidth=0.35, color='black', linestyle='solid')
plt.grid(True)
plt.minorticks_on()
plt.grid(which='major', linestyle='dashed', linewidth='0.15', color='gray')
plt.legend(loc='lower right', shadow=False, fontsize=5.)
plt.savefig('mass' + '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close('all')




ansorge = np.fromfile('./ansorge.txt', dtype=float, count=-1, sep=' ').reshape((193, 6))

U = ansorge[:, 0]
V = -ansorge[:, 1]
yp = ansorge[:, 2]
ym = ansorge[:, 3]
tke = ansorge[:, 4]
y = ansorge[:, 5]

fig0, ax0 = plt.subplots()
fig0.set_size_inches(4., 3.)
plt.ylabel(r'$y$')
plt.plot(U,y, linewidth=1., linestyle='solid', color='black', label='$u$')
plt.plot(V,y, linewidth=1., linestyle='dotted', color='gray', label='$v$')
plt.plot(tke/np.max(tke),y, linewidth=1., linestyle='dashed', color='gray', label='$TKE/TKEmax$')
plt.legend(loc='upper center', shadow=False, fontsize='small')
plt.savefig('ansorge.pdf', format = 'pdf')
plt.close('All')

#x = np.linspace(0, Lx, nx, endpoint = True)
y = np.linspace(0., Ly, ny, endpoint = True)
#z = np.linspace(0, Lz, nz, endpoint = True)


#pfront = 1.
#ug=1.

#y = y/57.
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.ylabel(r'$y$')
# plt.yscale("log")
# plt.xscale("log")
#plt.xlim(0., 1.1)
plt.ylim(0, 57.)
plt.plot(U,y,               linewidth=0.8, linestyle='solid',  color='black', label='$u$')
#plt.plot(tke/np.max(tke),y, linewidth=0.8, linestyle='dashed', color='black',  label='$TKE$')
plt.plot(tke,y, linewidth=0.8, linestyle='dashed', color='black',  label='$TKE$')
plt.plot(V,y,               linewidth=0.8, linestyle='solid', color='gray',  label='$v$')

print(np.max(V))

#d_Ekman = 4.
#a_star = 4.#0.75*d_Ekman
#ri_b=0.4
#phi = 0.5 * a_star * ri_b * ((-pi/log(0.01))**0.5) * erf( (y/a_star) / (-log(0.01))**(-0.5) )




#ri_b=0.2
#phi = tanh(y * ri_b)
#plt.plot(phi,y, linewidth=0.8, linestyle='solid', color='blue', label='$Ri_b=0.2$')

phi = erf ( y*0.2)
plt.plot(phi,y, linewidth=0.8, linestyle='dashed', color='red', label='$\phi$')

plt.grid(True)
plt.grid(which='major', linestyle='dashed', linewidth='0.4', color='gray')
plt.grid(which='minor', linestyle='solid', linewidth='0.2', color='gray')
plt.legend(loc='upper center', shadow=False, fontsize='5')
plt.savefig('ansorge_phi'+'.'+formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close('All')



'''
x = np.linspace(0, xlx, nx, endpoint=True)
y = np.linspace(0, yly, ny, endpoint=True)
z = np.linspace(0, zlz, nz, endpoint=True)
X, Y = np.meshgrid(x, y)

u1_3d = (np.fromfile('../data/ux0032', dtype=np.float32)).reshape((nz, ny, nx))
u1_2d = scipy.integrate.simps(u1_3d, dx=dz, axis=0)/zlz
u1_1d = scipy.integrate.simps(u1_2d, dx=dx, axis=1)/xlx

utau = (np.fromfile('../data/utmap0032', dtype=np.float32)).reshape((nz, nx))
print(utau.shape,np.max(utau))
utau_2d = scipy.integrate.simps(utau, dx=dz, axis=0)/zlz
print(utau_2d.shape)

utau_1d = scipy.integrate.simps(utau_2d, dx=dx, axis=0)/xlx
print(utau_1d.shape,utau_1d)


fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.scatter(y*utau_1d*500.,u1_1d/utau_1d , c='black', marker='x', s=0.8)
yy = np.linspace(0, 2.5e3, 2e5, endpoint=True)
dyy = 2.5e3 / len(yy)
plt.plot(yy[600:int(200. / dyy)], (1. / 0.4) * np.log(yy[600:int(200. / dyy)]) + 5., linewidth=1.,
         linestyle='dashed', color='black')#, label=r'$log\ law$')
plt.plot(yy[0:775], yy[0:775], linewidth=1., linestyle='dotted', color='black')#, label=r'$U^{+}=y^{+}$')
plt.ylim(0., 20.)
plt.xlim(5e-1, 1000.)
plt.ylabel(r'$u_1^{+}$')
plt.xlabel(r'$x_2^{+}$')
plt.xscale('log')
#plt.grid(True)
#plt.minorticks_on()
plt.legend(loc='upper left')
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
plt.savefig('u1_log'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close()

u2_3d = (np.fromfile('../data/uy0032', dtype=np.float32)).reshape((nz, ny, nx))
u2_2d = scipy.integrate.simps(u2_3d, dx=dz, axis=0)/zlz
u2_1d = scipy.integrate.simps(u2_2d, dx=dx, axis=1)/xlx

u3_3d = (np.fromfile('../data/uz0032', dtype=np.float32)).reshape((nz, ny, nx))
u3_2d = scipy.integrate.simps(u3_3d, dx=dz, axis=0)/zlz
u3_1d = scipy.integrate.simps(u3_2d, dx=dx, axis=1)/xlx

phi_3d = (np.fromfile('../data/phi10032', dtype=np.float32)).reshape((nz, ny, nx))
phi_2d = scipy.integrate.simps(phi_3d, dx=dz, axis=0)/zlz
phi_1d = scipy.integrate.simps(phi_2d, dx=dx, axis=1)/xlx

diss_3d = (np.fromfile('../data/diss0032', dtype=np.float32)).reshape((nz, ny, nx))
diss_2d = scipy.integrate.simps(diss_3d, dx=dz, axis=0)/zlz
diss_1d = scipy.integrate.simps(diss_2d, dx=dx, axis=1)/xlx


plt.figure(1, figsize=(7, 7))
plt.contourf(X, Y, u1_3d[int(nz/2.),:,:], 256, extend='both', cmap=plt.cm.hot, corner_mask=False)
plt.axes().set_aspect('equal')
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.colorbar(orientation="horizontal",fraction=0.05)
plt.grid(which='major', linestyle='dashed', linewidth='0.5', color='black')
# plt.grid(which='minor', linestyle='solid', linewidth='0.25', color='gray')
#plt.xlim(0., 40.)
#plt.ylim(0., 8.)
plt.savefig('ux_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close(1)

plt.figure(2, figsize=(7, 7))
plt.contourf(X, Y, u1_2d, 256, extend='both', cmap=plt.cm.hot, corner_mask=False)
plt.axes().set_aspect('equal')
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.colorbar(orientation="horizontal",fraction=0.05)
plt.grid(which='major', linestyle='dashed', linewidth='0.5', color='black')
# plt.grid(which='minor', linestyle='solid', linewidth='0.25', color='gray')
#plt.xlim(0., 40.)
#plt.ylim(0., 8.)
plt.savefig('uxm_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close(2)


fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(u1_3d[int(nz/2.),:,int(nx/2.)],y,linewidth=1., color='black', linestyle='dotted',label='3D')
plt.plot(u1_2d[:,int(nx/2.)],y           ,linewidth=1., color='black', linestyle='dashed',label='2D')
plt.plot(u1_1d,y                         ,linewidth=1., color='black', linestyle='solid' ,label='1D')
plt.xlabel(r'$u_1$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.legend(loc='upper left')
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
#plt.xlim(0., 1.2)
plt.ylim(0.,yly)
plt.savefig('u1_profiles_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close()

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(u2_3d[int(nz/2.),:,int(nx/2.)],y,linewidth=1., color='black', linestyle='dotted',label='3D')
plt.plot(u2_2d[:,int(nx/2.)],y           ,linewidth=1., color='black', linestyle='dashed',label='2D')
plt.plot(u2_1d,y                         ,linewidth=1., color='black', linestyle='solid' ,label='1D')
plt.xlabel(r'$u_2$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.legend(loc='upper left')
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
#plt.xlim(0., 1.2)
plt.ylim(0.,yly)
plt.savefig('u2_profiles_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close()

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(u3_3d[int(nz/2.),:,int(nx/2.)],y,linewidth=1., color='black', linestyle='dotted',label='3D')
plt.plot(u3_2d[:,int(nx/2.)],y           ,linewidth=1., color='black', linestyle='dashed',label='2D')
plt.plot(u3_1d,y                         ,linewidth=1., color='black', linestyle='solid' ,label='1D')
plt.xlabel(r'$u_3$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.legend(loc='upper right')
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
#plt.xlim(0., 1.2)
plt.ylim(0.,yly)
plt.savefig('u3_profiles_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close()

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(phi_3d[int(nz/2.),:,int(nx/2.)],y,linewidth=1., color='black', linestyle='dotted',label='3D')
plt.plot(phi_2d[:,int(nx/2.)],y           ,linewidth=1., color='black', linestyle='dashed',label='2D')
plt.plot(phi_1d,y                         ,linewidth=1., color='black', linestyle='solid' ,label='1D')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.legend(loc='lower right')
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
#plt.xlim(0., 1.2)
plt.ylim(0.,yly)
plt.savefig('phi_profiles_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close()

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(diss_3d[int(nz/2.),:,int(nx/2.)],y,linewidth=1., color='black', linestyle='dotted',label='3D')
plt.plot(diss_2d[:,int(nx/2.)],y           ,linewidth=1., color='black', linestyle='dashed',label='2D')
plt.plot(diss_1d,y                         ,linewidth=1., color='black', linestyle='solid' ,label='1D')
plt.xlabel(r'$\epsilon$')
plt.ylabel(r'$x_2$')
plt.grid(True)
plt.minorticks_on()
plt.legend(loc='upper right')
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
#plt.xlim(0., 1.2)
plt.ylim(0.,3.)
plt.savefig('diss_profiles_t250'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close()

#############
## Reading
#############
s0 = Statistics(nphi, '../out/statistics')
s0e = Energy('../out/budget')
s0.key = "DNS"
s0.cor = 'black'
s0.tam = 1.
s0.sty = 'solid'

#############
## Suspended
#############
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.xlabel(r'$t$')
plt.ylabel(r'$m_p$')
#plt.xlim(0., 1000.)

D = s0
plt.plot(D.t, D.mp[0], color=D.cor, linewidth=D.tam, linestyle=D.sty)#, label=D.key)
#plt.plot(D.t, D.mp[0] / (D.mp[0][0]), color=D.cor, linewidth=D.tam, linestyle=D.sty, label=D.key)
#plt.axhline(y=1., xmin=0., xmax=1000, linewidth=0.5, color='black', linestyle='dotted')

# plt.axvline(x=2.4, ymin=0, ymax=1000, linewidth=0.5, color='black',linestyle='dotted')
# plt.axvline(x=4.7, ymin=0, ymax=1000, linewidth=0.5, color='black',linestyle='dotted')

plt.grid(True)
plt.minorticks_on()
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
plt.legend(loc='lower right')
plt.savefig('suspended' + '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close('all')

############
# Energy
############
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.xlabel(r'$t$')
#plt.ylabel(r'$m_p$')
#plt.xlim(0., 1000.)

sim1 = Energy('../out/budget')
#sim1.dep[0]=0.; sim1.dek[0]=0.
sim1_deki = [0.]; sim1_depi = [0.]
for t in range(1, sim1.t.shape[0]):
    sim1_deki.append((sim1.dek[t-1]+sim1.dek[t])*dt*float(imodule)/2.+sim1_deki[t-1])
    sim1_depi.append((sim1.dep[t-1]+sim1.dep[t])*dt*float(imodule)/2.+sim1_depi[t-1])


sz = 1.
ls = 'solid'
e0 = sim1.ek[0] + sim1.ep[0]
plt.plot(sim1.t, sim1.ek  , linestyle=ls, linewidth=sz, color='blue')
plt.plot(sim1.t, sim1_deki, linestyle=ls, linewidth=sz, color='orange')
plt.plot(sim1.t, sim1.ep  , linestyle=ls, linewidth=sz, color='green')
plt.plot(sim1.t, sim1_depi, linestyle=ls, linewidth=sz, color='red')
plt.plot(sim1.t, (sim1.ek+sim1_deki+sim1.ep+sim1_depi), color='black', linestyle=ls, linewidth=sz)
#plt.axhline(y=1., xmin=0., xmax=1000, linewidth=0.5, color='black', linestyle='dotted')
plt.grid(True)
plt.minorticks_on()
plt.grid(which='major', linestyle='dashed', linewidth='0.35', color='gray')
plt.legend(loc='lower right')
plt.savefig('energy' + '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
plt.close('all')
'''


############
# Profiles
############

class Probe(object):
    def __init__(self, nphi, filename):
        self.nphi = nphi
        with open(filename, 'r') as f:
            data = list(csv.reader(f, delimiter=' ', skipinitialspace=True,
                                   quoting=csv.QUOTE_NONNUMERIC))
        data = np.transpose(data)
        self.t = data[0]
        counter = 1
        self.uxm = []
        for n in range(final):
            self.phi.append(data[counter])
            counter += 1


for probe in range(1,191):

    data = Probe(1,'./out/probe'+str(probe).zfill(4) )
    fig = plt.figure()
    fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$t$')

    plt.plot(data.t,data.ux    ,label='$u_1$', linewidth=0.5)

    plt.legend(loc='lower right')
    print('Saving probe',str(probe).zfill(4))
    plt.savefig('probe'+str(probe).zfill(4)+'.'+formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
    plt.close()



############
# Probes
############

class Probe(object):
    def __init__(self, nphi, filename):
        self.nphi = nphi
        with open(filename, 'r') as f:
            data = list(csv.reader(f, delimiter=' ', skipinitialspace=True,
                                   quoting=csv.QUOTE_NONNUMERIC))
        data = np.transpose(data)
        self.t = data[0]
        self.ux = data[1]
        self.uy = data[2]
        self.uz = data[3]
        self.pre = data[4]
        self.diss = data[5]
        counter = 6
        self.phi = []
        for n in range(nphi):
            self.phi.append(data[counter])
            counter += 1


for probe in range(1,191):

    data = Probe(1,'./out/probe'+str(probe).zfill(4) )
    fig = plt.figure()
    fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$t$')
    #plt.xlim(0.,1000.)

    plt.plot(data.t,data.ux    ,label='$u_1$', linewidth=0.5)
    plt.plot(data.t,data.uy    ,label='$u_2$', linewidth=0.5)
    plt.plot(data.t,data.uz    ,label='$u_3$', linewidth=0.5)
    plt.plot(data.t,data.pre   ,label='$pre$', linewidth=0.5)
    plt.plot(data.t,data.diss  ,label='$diss$', linewidth=0.5)
    plt.plot(data.t,data.phi[0],label='$\phi$', linewidth=0.5)

    plt.legend(loc='lower right')
    print('Saving probe',str(probe).zfill(4))
    plt.savefig('probe'+str(probe).zfill(4)+'.'+formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp)
    plt.close()

