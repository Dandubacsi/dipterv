# importing necessary libraries
from sys import argv
import csv
from matplotlib import pyplot as plt
import numpy as np
from math import sin, pi, sqrt
from auxilary import TableLookup1D, PrefixCorrection

# model class for holding information about a model configuration
class Model:
    def __init__(self, width, resistance, amplitude, frequency, separation, parasitic_cap):
        # electrode parameters
        self.w = width                  # electrode width
        self.c10 = np.array([])         # measuring electrode capacitance
        self.c12 = np.array([])         # mutual capactiance
        self.c20 = np.array([])         # shield capacitance
        self.c = np.array([])           # charges on the sample's surface
        self.distances = np.array([])   # simulated distances

        # vibration parameters
        self.f = frequency              # signal frequency
        self.A = amplitude              # vibration amplitude
        self.d0 = separation            # rest separation

        # readout circuit parameters
        self.R = resistance             # measuring resistance
        self.Cp = parasitic_cap         # parasitic capacitance to ground

    def addData(self, distance, c10, c20, c12, c):
        # add data points to the given configuration
        self.distances = np.append(self.distances, distance)
        self.c10 = np.append(self.c10,c10)
        self.c12 = np.append(self.c12,c12)
        self.c20 = np.append(self.c20,c20)
        if np.size(self.c) == 0:
            self.c = c
        else:
            self.c = np.dstack((self.c, c))

    def simulate(self, potential_map, time, U0):
        # separation as a function of time
        separation = [self.d0 + self.A*sin(2*pi*self.f*tau) for tau in time]

        # auxilary variables
        size = np.size(self.c, 0)
        dt = time[1]-time[0]

        # capactiance-time functions
        c10_t = np.sum(self.c,0)
        c10_t = np.sum(c10_t,0)
        C10_t = TableLookup1D(self.distances,c10_t,separation)
        C12_t = TableLookup1D(self.distances,self.c12,separation)
        C10_i_t = np.zeros([size, size , np.size(time)])
        for i in range(size):
            for j in range(size):
                C10_i_t[i, j, :] = TableLookup1D(self.distances, self.c[i, j, :], separation)

        # and their numerical derivatives
        dC10_t = np.diff(C10_t)/dt
        dC12_t = np.diff(C12_t)/dt
        dC10_i_t = np.zeros([size, size, np.size(time)-1])
        for i in range(size):
            for j in range(size):
                dC10_i_t[i, j, :] = np.diff(C10_i_t[i, j, :])/dt

        # solutions functions
        u_t = np.zeros(np.size(separation))

        # solving the ODE with the backwards Euler method
        for ind, t in enumerate(time[:-1]):
            # surface potential contribution
            contribution = 0
            for i in range(size):
                for j in range(size):
                    contribution = contribution + dC10_i_t[i, j, ind]*(U0-potential_map[i, j])

            # update solution
            u_t[ind + 1] = (u_t[ind]*self.R*(C12_t[ind]+C10_t[ind]+self.Cp)+contribution*self.R*dt)/(self.R*(C12_t[ind]+C10_t[ind]+self.Cp)+(1+self.R*(dC12_t[ind]+dC10_t[ind]))*dt)

        return u_t

    def simulateFFT(self, potential_map, time, U0):
        # auxilary variables
        size = np.size(self.c, 0)
        dt = time[1]-time[0]
        N = np.size(time, 0)
        w0 = 2*pi*self.f

        # separation as a function of time
        separation = [self.d0 + self.A*sin(w0*tau) for tau in time]

        # capactiance-time functions
        C10_t = np.sum(self.c,0)
        C10_t = np.sum(C10_t,0)
        C10_t = TableLookup1D(self.distances,C10_t,separation)
        C12_t = TableLookup1D(self.distances,self.c12,separation)
        C10_i_t = np.zeros([size, size , np.size(time)])
        for i in range(size):
            for j in range(size):
                C10_i_t[i, j, :] = TableLookup1D(self.distances, self.c[i, j, :], separation)

        # FFT of capacitnace-time functions
        C10_fft = np.fft.fft(C10_t)
        C12_fft = np.fft.fft(C12_t)
        C10_i_fft = np.fft.fft(C10_i_t)

        # spectral derivatives
        if N%2 == 0:
            D = np.array([0])
            D = np.concatenate((D, np.arange(-N/2+1, N/2), D))
            D = np.diag(1j*w0*D)

            C10_fft = np.fft.fftshift(C10_fft)
            C10_fft = np.append(C10_fft, C10_fft[0]/2)
            C10_fft[0] = C10_fft[0]/2
            C12_fft = np.fft.fftshift(C12_fft)
            C12_fft = np.append(C12_fft, C12_fft[0]/2)
            C12_fft[0] = C12_fft[0]/2
            C10_i_fft = np.fft.fftshift(C10_i_fft, 2)
            C10_i_fft = np.dstack((C10_i_fft, C10_i_fft[:,:,0]/2))
            C10_i_fft[:,:,0] = C10_i_fft[:,:,0]/2

            dC10_fft = D.dot(C10_fft)
            dC10_fft = dC10_fft[0:-1]
            dC10_fft[0] = 2*dC10_fft[0]
            dC12_fft = D.dot(C12_fft)
            dC12_fft = dC12_fft[0:-1]
            dC12_fft[0] = 2*dC12_fft[0]
            dC10_i_fft = np.zeros([size, size, N], dtype=np.cdouble)
            for i in range(size):
                for j in range(size):
                    tmp = D.dot(C10_i_fft[i,j,:])
                    tmp = tmp[0:-1]
                    tmp[0] = 2*tmp[0]
                    dC10_i_fft[i,j,:] = tmp

            dC10_t = np.real(np.fft.ifft(np.fft.ifftshift(dC10_fft)))
            dC12_t = np.real(np.fft.ifft(np.fft.ifftshift(dC12_fft)))
            dC10_i_t = np.real(np.fft.ifft(np.fft.ifftshift(dC10_i_fft, 2)))

        else:
            D = np.diag(1j*w0*np.arange(-(N-1)/2, (N+1)/2))

            dC10_t = np.fft.ifft(np.fft.ifftshift(D.dot(np.fft.fftshift(C10_fft))))
            dC12_t = np.fft.ifft(np.fft.ifftshift(D.dot(np.fft.fftshift(C12_fft))))
            for i in range(size):
                for j in range(size):
                    dC10_i_t[i,j,:] = np.fft.ifft(np.fft.ifftshift(D.dot(np.fftshift(C10_i_fft[i,j,:]))))


        # driver voltage time function
        U = np.ones(N)*U0                       # constant driving voltage
        # U = np.cos(2*pi*self.f*time)*U0         # sinusoidally changing driving voltage

        # surface potential contribution
        surface_pot_cont = np.zeros(N)
        for ind in range(N):
            for i in range(size):
                for j in range(size):
                    surface_pot_cont[ind] = surface_pot_cont[ind] + dC10_i_t[i, j, ind]*(U[ind]-potential_map[i, j])

        # coefficient functions
        f = np.divide((1*np.ones(N)+self.R*(dC12_t+dC10_t)),(self.R*(C12_t+self.Cp*np.ones(N)+C10_t)))
        g = np.divide(surface_pot_cont,C12_t+self.Cp*np.ones(N)+C10_t)

        # their fft's
        F = np.fft.fftshift(np.fft.fft(f))
        G = np.fft.fftshift(np.fft.fft(g))

        if N%2 == 0:
            Fmod = [F[0]/2]
            Fmod = np.append(Fmod, F[1:])
            Fmod = np.append(Fmod, F[0]/2)
            Gmod = [G[0]/2]
            Gmod = np.append(Gmod, G[1:])
            Gmod = np.append(Gmod, G[0]/2)

            Fmat = np.zeros([N+1, N+1],dtype=np.cdouble)
            for k in range(N+1):
                Fmat = Fmat + np.diag(Fmod[k]*np.ones(int(1+N-abs(k-N/2))),int(N/2-k))
            Fmat = Fmat/N

            Y = np.linalg.solve(D+Fmat,Gmod)
            Y = Y[:-1]
            Y[0] = 2*Y[0]
            return np.real(np.fft.ifft(np.fft.ifftshift(Y)))
        else:
            Fmat = np.zeros([N, N],dtype=np.cdouble)
            for k in range(N):
                Fmat = Fmat + np.diag(F[k]*np.ones(int(N-abs(k-(N-1)/2))),int((N-1)/2-k))
            Fmat = Fmat/N

            Y = np.linalg.solve(D+Fmat,G)
            return np.real(np.fft.ifft(np.fft.ifftshift(Y)))

    def simulateNoise(self, time, power):
        # auxilary variables
        size = np.size(self.c, 0)
        dt = time[1]-time[0]
        N = np.size(time, 0)
        w0 = 2*pi*self.f

        # separation as a function of time
        separation = [self.d0 + self.A*sin(w0*tau) for tau in time]

        # capactiance-time functions
        C10_t = np.sum(self.c,0)
        C10_t = np.sum(C10_t,0)
        C10_t = TableLookup1D(self.distances,C10_t,separation)
        C12_t = TableLookup1D(self.distances,self.c12,separation)
        C10_i_t = np.zeros([size, size , np.size(time)])
        for i in range(size):
            for j in range(size):
                C10_i_t[i, j, :] = TableLookup1D(self.distances, self.c[i, j, :], separation)

        # FFT of capacitnace-time functions
        C10_fft = np.fft.fft(C10_t)
        C12_fft = np.fft.fft(C12_t)
        C10_i_fft = np.fft.fft(C10_i_t)

        # spectral derivatives
        if N%2 == 0:
            D = np.array([0])
            D = np.concatenate((D, np.arange(-N/2+1, N/2), D))
            D = np.diag(1j*w0*D)

            C10_fft = np.fft.fftshift(C10_fft)
            C10_fft = np.append(C10_fft, C10_fft[0]/2)
            C10_fft[0] = C10_fft[0]/2
            C12_fft = np.fft.fftshift(C12_fft)
            C12_fft = np.append(C12_fft, C12_fft[0]/2)
            C12_fft[0] = C12_fft[0]/2
            C10_i_fft = np.fft.fftshift(C10_i_fft, 2)
            C10_i_fft = np.dstack((C10_i_fft, C10_i_fft[:,:,0]/2))
            C10_i_fft[:,:,0] = C10_i_fft[:,:,0]/2

            dC10_fft = D.dot(C10_fft)
            dC10_fft = dC10_fft[0:-1]
            dC10_fft[0] = 2*dC10_fft[0]
            dC12_fft = D.dot(C12_fft)
            dC12_fft = dC12_fft[0:-1]
            dC12_fft[0] = 2*dC12_fft[0]
            dC10_i_fft = np.zeros([size, size, N], dtype=np.cdouble)
            for i in range(size):
                for j in range(size):
                    tmp = D.dot(C10_i_fft[i,j,:])
                    tmp = tmp[0:-1]
                    tmp[0] = 2*tmp[0]
                    dC10_i_fft[i,j,:] = tmp

            dC10_t = np.real(np.fft.ifft(np.fft.ifftshift(dC10_fft)))
            dC12_t = np.real(np.fft.ifft(np.fft.ifftshift(dC12_fft)))
            dC10_i_t = np.real(np.fft.ifft(np.fft.ifftshift(dC10_i_fft, 2)))

        else:
            D = np.diag(1j*w0*np.arange(-(N-1)/2, (N+1)/2))

            dC10_t = np.fft.ifft(np.fft.ifftshift(D.dot(np.fft.fftshift(C10_fft))))
            dC12_t = np.fft.ifft(np.fft.ifftshift(D.dot(np.fft.fftshift(C12_fft))))
            for i in range(size):
                for j in range(size):
                    dC10_i_t[i,j,:] = np.fft.ifft(np.fft.ifftshift(D.dot(np.fftshift(C10_i_fft[i,j,:]))))

        # coefficient functions
        a = np.divide((-1*np.ones(N)-self.R*(dC12_t+dC10_t)),(self.R*(C12_t+self.Cp*np.ones(N)+C10_t)))
        c = np.square(np.divide(1,(self.R*(C12_t+self.Cp*np.ones(N)+C10_t))))

        # their fft's
        A = np.fft.fftshift(np.fft.fft(a))
        C = np.fft.fftshift(np.fft.fft(c))

        if N%2 == 0:
            Amod = [A[0]/2]
            Amod = np.append(Amod, A[1:])
            Amod = np.append(Amod, A[0]/2)

            Cmod = [C[0]/2]
            Cmod = np.append(Cmod, C[1:])
            Cmod = np.append(Cmod, C[0]/2)

            Amat = np.zeros([N+1, N+1],dtype=np.cdouble)
            for k in range(N+1):
                Amat = Amat + np.diag(Amod[k]*np.ones(int(1+N-abs(k-N/2))),int(N/2-k))
            Amat = Amat/N

            V = np.linalg.solve(D-2*Amat, power*Cmod)
            V = V[:-1]
            V[0] = 2*V[0]
            return np.real(np.fft.ifft(np.fft.ifftshift(V)))
        else:
            Amat = np.zeros([N, N],dtype=np.cdouble)
            for k in range(N):
                Amat = Amat + np.diag(A[k]*np.ones(int(N-abs(k-(N-1)/2))),int((N-1)/2-k))
            Amat = Amat/N

            V = np.linalg.solve(D-2*Amat,power*C)
            return np.real(np.fft.ifft(np.fft.ifftshift(V)))

    def plotCapacitances(self):
        # plot capacitance as a function of distance

        # font sizes
        majorFS = 24
        minorFS = 20

        # create figure
        fig = plt.figure()
        ax1 = plt.subplot()
        plt.subplots_adjust(right=0.85)
        ax2 = ax1.twinx()
        ax3 = ax1.twinx()
        ax3.spines.right.set_position(("axes", 1.1))

        # plot
        p1, = ax1.plot(self.distances, self.c10, color='r', linewidth=2, label='Elektróda-minta kapacitás')
        p2, = ax2.plot(self.distances, self.c12, color='b', linewidth=2, label='Elektróda-fókuszálás kapacitás')
        p3, = ax3.plot(self.distances, self.c20, color='g', linewidth=2, label='Fókuszálás-minta kapacitás')

        # set title
        plt.title('Kapacitás együtthatók '+str(self.w) +' \u03BCm-es elektróda szélességnél\n',fontsize=majorFS)

        # set axis labels
        ax1.set_xlabel('Elektróda-minta távolság \u03BCm-ben',fontsize=minorFS)
        ax1.set_ylabel('fF',fontsize=minorFS)
        ax2.set_ylabel('fF',fontsize=minorFS)
        ax3.set_ylabel('fF',fontsize=minorFS)

        # get traces' colors
        tkw = dict(size=4, width=1.5)
        ax1.tick_params(axis='y', colors=p1.get_color(), **tkw, labelsize=minorFS)
        ax1.tick_params(axis='x', labelsize=minorFS)
        ax2.tick_params(axis='y', colors=p2.get_color(), **tkw, labelsize=minorFS)
        ax3.tick_params(axis='y', colors=p3.get_color(), **tkw, labelsize=minorFS)

        # set axis' colors
        ax1.yaxis.label.set_color(p1.get_color())
        ax2.yaxis.label.set_color(p2.get_color())
        ax3.yaxis.label.set_color(p3.get_color())

        # set legend
        ax1.legend(handles=[p1, p2, p3], loc='center right', prop={"size":minorFS})

    def plotStability(self):
        # plot the stability criterion
        # [mV, pA, Gohm, fF, MHz, us, um, m/s]

        # calculate auxilary variable
        c0 = self.c10+self.c12

        # create figure
        fig = plt.figure(figsize=(12.8, 8))
        ax1 = plt.subplot()

        # plot
        p1, = ax1.plot(self.distances[:-1], 2*pi*self.f*self.A*self.R*abs(np.diff(c0)), color='b', linewidth=2, label='Stabilitási kritérium')

        # set title
        plt.title('Stabilitás vizsgálat '+ str(self.w) +' \u03BCm-es elektróda szélességnél')

        # set axis' labels
        ax1.set_xlabel('Elektróda-minta távolság \u03BCm-ben')
        ax1.set_ylabel('1')

        # set legend
        ax1.legend(handles=[p1], loc='center right')

    def plotSignal(self):
        # [mV, pA, Gohm, fF, MHz, us, um, m/s]

        # auxilary variable
        T = 1/self.f
        periods = 2
        t = np.linspace(0,periods*T,periods*100)
        N = np.shape(self.c)[0]
        runs = 6

        distances = [self.d0 + self.A*sin(2*pi*self.f*tau) for tau in t]

        # font sizes
        majorFS = 22
        middleFS = 18
        minorFS = 16

        # create figure
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(3,4)
        plt.suptitle('Felületi potenciáleloszlások és a szimulált jelek\n ' + str(self.w) + '\u03BCm-es elektróda szélességnél', fontsize=majorFS)

        # potential distributions
        x = np.arange(N)*5
        y = np.arange(N)*5
        z0 = np.ones([N,N])*500
        z1 = np.array([[int(v > 5*N/3-1) for v in x] for v in x])*1000
        z2 = np.random.normal(333, 100,(N,N))

        # constant potential
        ax = fig.add_subplot(gs[0,0])
        plot = ax.pcolormesh(x, y, z0, cmap=plt.colormaps['gray'], vmin=0, vmax=1000)
        ax.set_title('Homogén potenciál', fontsize=middleFS, y=1.05)
        ax.set_aspect('equal')
        ax.set_xlim([-0.5,N-0.5])
        ax.set_ylim([-0.5,N-0.5])
        ax.set_xticks(np.arange(min(x), max(x)+1, (max(x)-min(x))/4))
        ax.set_yticks(np.arange(min(y), max(y)+1, (max(y)-min(y))/4))
        ax.tick_params(axis='x',labelsize=minorFS)
        ax.tick_params(axis='y',labelsize=minorFS)
        cb = fig.colorbar(plot, ax=ax, location='right')
        cb.ax.tick_params(labelsize=minorFS)
        cb.set_label(label='mV', size=minorFS)
        ax.axes.xaxis.set_visible(True)
        ax.axes.yaxis.set_visible(True)
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        ax = fig.add_subplot(gs[0,1:])
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax_t.set_ylabel('mV', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        for run in range(runs):
            u_t = self.simulate(z0, t, float(run)/runs*1000)
            lns = lns + ax_t.plot(t,u_t,label=f'{(float(run)/runs*1000):3.1f} mV')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)
        # ax.legend(lns, lbs)

        # potential step
        ax = fig.add_subplot(gs[1,0])
        plot = ax.pcolormesh(x, y, z1, cmap=plt.colormaps['gray'], vmin=0, vmax=1000)
        ax.set_title('Potenciál ugrás', fontsize=middleFS, y=1.05)
        ax.set_aspect('equal')
        ax.set_xlim([-0.5,N-0.5])
        ax.set_ylim([-0.5,N-0.5])
        ax.set_xticks(np.arange(min(x), max(x)+1, (max(x)-min(x))/4))
        ax.set_yticks(np.arange(min(y), max(y)+1, (max(y)-min(y))/4))
        ax.tick_params(axis='x',labelsize=minorFS)
        ax.tick_params(axis='y',labelsize=minorFS)
        cb = fig.colorbar(plot, ax=ax, location='right')
        cb.ax.tick_params(labelsize=minorFS)
        cb.set_label(label='mV', size=minorFS)
        ax.axes.xaxis.set_visible(True)
        ax.axes.yaxis.set_visible(True)
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        ax = fig.add_subplot(gs[1,1:])
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax_t.set_ylabel('mV', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        for run in range(runs):
            u_t = self.simulate(z1, t, float(run)/runs*1000)
            lns = lns + ax_t.plot(t,u_t,label=f'{(float(run)/runs*1000):3.1f} mV')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)
        ax.legend(lns, lbs, bbox_to_anchor=(1.1,0.5), loc='center left', prop={"size":minorFS})

        # noisy potential
        ax = fig.add_subplot(gs[2,0])
        plot = ax.pcolormesh(x, y, z2, cmap=plt.colormaps['gray'], vmin=0, vmax=1000)
        ax.set_title('Gauss potenciál', fontsize=middleFS, y=1.05)
        ax.set_aspect('equal')
        ax.set_xlim([-0.5,N-0.5])
        ax.set_ylim([-0.5,N-0.5])
        ax.set_xticks(np.arange(min(x), max(x)+1, (max(x)-min(x))/4))
        ax.set_yticks(np.arange(min(y), max(y)+1, (max(y)-min(y))/4))
        ax.tick_params(axis='x',labelsize=minorFS)
        ax.tick_params(axis='y',labelsize=minorFS)
        cb = fig.colorbar(plot, ax=ax, location='right')
        cb.ax.tick_params(labelsize=minorFS)
        cb.set_label(label='mV', size=minorFS)
        ax.axes.xaxis.set_visible(True)
        ax.axes.yaxis.set_visible(True)
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        ax = fig.add_subplot(gs[2,1:])
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax_t.set_ylabel('mV', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        for run in range(runs):
            u_t = self.simulate(z2, t, float(run)/runs*1000)
            lns = lns + ax_t.plot(t,u_t,label=f'{(float(run)/runs*1000):3.1f} mV')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)


        # ax.legend(lns, lbs)
        ax.set_xlabel('idő $\mu s$', fontsize=minorFS)

    def plotSensitivity(self):
        # font sizes
        majorFS = 24
        minorFS = 16

        # create figure
        fig, axs = plt.subplots(1,3,constrained_layout=True)
        fig.suptitle('Normalizált felületi érzékenységek '+ str(self.w) + '\u03BCm-es elektróda szélességnél', fontsize=majorFS)

        # auxilary variable
        d_min = self.d0-self.A
        d_min_index = np.absolute(self.distances-d_min).argmin()
        d_max = self.d0+self.A
        d_max_index = np.absolute(self.distances-d_max).argmin()
        d_mean = (d_min+d_max)/2
        d_mean_index = int((d_min_index+d_max_index)/2)
        N = np.shape(self.c)[0]

        # sensitivity distributions
        x = np.arange(N)*5
        y = np.arange(N)*5
        z_min = -np.diff(self.c)[:,:,d_min_index-1]
        z_min = z_min/np.sum(z_min)
        z_mean = -np.diff(self.c)[:,:,d_mean_index-1]
        z_mean = z_mean/np.sum(z_mean)
        z_max = -np.diff(self.c)[:,:,d_max_index-1]
        z_max = z_max/np.sum(z_max)

        # plot sensitivities
        plot = axs[0].pcolormesh(x,y,z_min,cmap=plt.colormaps['jet'])
        axs[0].set_title('Minimális elektróda távolság', fontsize=majorFS)
        cb = fig.colorbar(plot, ax = axs[0], location='right', fraction=0.05)
        cb.ax.tick_params(labelsize=minorFS)
        axs[0].axes.xaxis.set_visible(True)
        axs[0].axes.yaxis.set_visible(True)
        axs[0].set_xticks(np.arange(min(x), max(x)+2, (max(x)-min(x)+1)/4))
        axs[0].set_yticks(np.arange(min(y), max(y)+2, (max(y)-min(y)+1)/4))
        axs[0].tick_params(axis='x', labelsize=minorFS)
        axs[0].tick_params(axis='y', labelsize=minorFS)
        axs[0].set_aspect('equal')

        plot = axs[1].pcolormesh(x,y,z_mean,cmap=plt.colormaps['jet'])
        axs[1].set_title('Közepes elektróda távolság', fontsize=majorFS)
        cb = fig.colorbar(plot, ax = axs[1], location='right', fraction=0.05)
        cb.ax.tick_params(labelsize=minorFS)
        axs[1].axes.xaxis.set_visible(True)
        axs[1].axes.yaxis.set_visible(True)
        axs[1].set_xticks(np.arange(min(x), max(x)+2, (max(x)-min(x)+1)/4))
        axs[1].set_yticks(np.arange(min(y), max(y)+2, (max(y)-min(y)+1)/4))
        axs[1].tick_params(axis='x', labelsize=minorFS)
        axs[1].tick_params(axis='y', labelsize=minorFS)
        axs[1].set_aspect('equal')

        plot = axs[2].pcolormesh(x,y,z_max,cmap=plt.colormaps['jet'])
        axs[2].set_title('Maximális elektróda távolság', fontsize=majorFS)
        cb = fig.colorbar(plot, ax = axs[2], location='right', fraction=0.05)
        cb.ax.tick_params(labelsize=minorFS)
        axs[2].axes.xaxis.set_visible(True)
        axs[2].axes.yaxis.set_visible(True)
        axs[2].set_xticks(np.arange(min(x), max(x)+2, (max(x)-min(x)+1)/4))
        axs[2].set_yticks(np.arange(min(y), max(y)+2, (max(y)-min(y)+1)/4))
        axs[2].tick_params(axis='x', labelsize=minorFS)
        axs[2].tick_params(axis='y', labelsize=minorFS)
        axs[2].set_aspect('equal')

    def plotStationary(self):
        # [mV, pA, Gohm, fF, MHz, us, um, m/s]

        # auxilary variable
        T = 1/self.f
        t = np.linspace(0, T, 200)
        N = np.shape(self.c)[0]
        runs = 6

        distances = [self.d0 + self.A*sin(2*pi*self.f*tau) for tau in t]

        # font sizes
        majorFS = 22
        middleFS = 18
        minorFS = 16

        # create figure
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(3,4)
        plt.suptitle('Felületi potenciáleloszlások és az állandósult jelek\n ' + str(self.w) + '\u03BCm-es elektróda szélességnél', fontsize=majorFS)

        # potential distributions
        x = np.arange(N)*5
        y = np.arange(N)*5
        z0 = np.ones([N,N])*500
        z1 = np.array([[int(v > 5*N/3-1) for v in x] for v in x])*1000
        z2 = np.random.normal(333, 100,(N,N))

        # constant potential
        ax = fig.add_subplot(gs[0,0])
        plot = ax.pcolormesh(x, y, z0, cmap=plt.colormaps['gray'], vmin=0, vmax=1000)
        ax.set_title('Homogén potenciál', fontsize=middleFS, y=1.05)
        ax.set_aspect('equal')
        ax.set_xlim([-0.5,N-0.5])
        ax.set_ylim([-0.5,N-0.5])
        ax.set_xticks(np.arange(min(x), max(x)+1, (max(x)-min(x))/4))
        ax.set_yticks(np.arange(min(y), max(y)+1, (max(y)-min(y))/4))
        ax.tick_params(axis='x',labelsize=minorFS)
        ax.tick_params(axis='y',labelsize=minorFS)
        cb = fig.colorbar(plot, ax=ax, location='right')
        cb.ax.tick_params(labelsize=minorFS)
        cb.set_label(label='mV', size=minorFS)
        ax.axes.xaxis.set_visible(True)
        ax.axes.yaxis.set_visible(True)
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        ax = fig.add_subplot(gs[0,1:])
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax_t.set_ylabel('mV', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        for run in range(runs):
            u_t = self.simulateFFT(z0, t, float(run)/runs*1000)
            lns = lns + ax_t.plot(t,u_t,label=f'{(float(run)/runs*1000):3.1f} mV')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)


        # potential step
        ax = fig.add_subplot(gs[1,0])
        plot = ax.pcolormesh(x, y, z1, cmap=plt.colormaps['gray'], vmin=0, vmax=1000)
        ax.set_title('Potenciál ugrás', fontsize=middleFS, y=1.05)
        ax.set_aspect('equal')
        ax.set_xlim([-0.5,N-0.5])
        ax.set_ylim([-0.5,N-0.5])
        ax.set_xticks(np.arange(min(x), max(x)+1, (max(x)-min(x))/4))
        ax.set_yticks(np.arange(min(y), max(y)+1, (max(y)-min(y))/4))
        ax.tick_params(axis='x',labelsize=minorFS)
        ax.tick_params(axis='y',labelsize=minorFS)
        cb = fig.colorbar(plot, ax=ax, location='right')
        cb.ax.tick_params(labelsize=minorFS)
        cb.set_label(label='mV', size=minorFS)
        ax.axes.xaxis.set_visible(True)
        ax.axes.yaxis.set_visible(True)
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        ax = fig.add_subplot(gs[1,1:])
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax_t.set_ylabel('mV', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        for run in range(runs):
            u_t = self.simulateFFT(z1, t, float(run)/runs*1000)
            lns = lns + ax_t.plot(t,u_t,label=f'{(float(run)/runs*1000):3.1f} mV')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)
        ax.legend(lns, lbs, bbox_to_anchor=(1.1,0.5), loc='center left', prop={"size":minorFS})

        # noisy potential
        ax = fig.add_subplot(gs[2,0])
        plot = ax.pcolormesh(x, y, z2, cmap=plt.colormaps['gray'], vmin=0, vmax=1000)
        ax.set_title('Gauss potenciál', fontsize=middleFS, y=1.05)
        ax.set_aspect('equal')
        ax.set_xlim([-0.5,N-0.5])
        ax.set_ylim([-0.5,N-0.5])
        ax.set_xticks(np.arange(min(x), max(x)+1, (max(x)-min(x))/4))
        ax.set_yticks(np.arange(min(y), max(y)+1, (max(y)-min(y))/4))
        ax.tick_params(axis='x',labelsize=minorFS)
        ax.tick_params(axis='y',labelsize=minorFS)
        cb = fig.colorbar(plot, ax=ax, location='right')
        cb.ax.tick_params(labelsize=minorFS)
        cb.set_label(label='mV', size=minorFS)
        ax.axes.xaxis.set_visible(True)
        ax.axes.yaxis.set_visible(True)
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')

        ax = fig.add_subplot(gs[2,1:])
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax_t.set_ylabel('mV', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        for run in range(runs):
            u_t = self.simulateFFT(z2, t, float(run)/runs*1000)
            lns = lns + ax_t.plot(t,u_t,label=f'{(float(run)/runs*1000):3.1f} mV')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)

        # ax.legend(lns, lbs)
        ax.set_xlabel('idő $\mu s$', fontsize=minorFS)

    def plotNoise(self):
        # [mV, pA, Gohm, fF, MHz, us, um, m/s]

        # auxilary variable
        T = 1/self.f
        t = np.linspace(0, T, 200)
        N = np.shape(self.c)[0]
        k_B = 1.38e-2
        Temp = 300

        distances = [self.d0 + self.A*sin(2*pi*self.f*tau) for tau in t]

        # font sizes
        majorFS = 24
        minorFS = 18

        # create figure
        fig, ax = plt.subplots(1,1,constrained_layout=True)
        ax.set_title('Termikus zaj hatása ' + str(self.w) + '\u03BCm-es elektróda szélességnél', fontsize=majorFS)
        ax_t = ax.twinx()
        ax.set_ylabel('$\mu m$', fontsize=minorFS)
        ax.set_xlabel('idő $\mu s$', fontsize=minorFS)
        ax_t.set_ylabel('$\mu V$', fontsize=minorFS)
        lns = ax.plot(t,distances,label='Rezgés',color='b')
        sigma_t = np.sqrt(self.simulateNoise(t, 4*k_B*self.R*Temp*self.f*2))*1000
        lns = lns + ax_t.plot(t,sigma_t,label='Szórás',color='r')
        lbs = [l.get_label() for l in lns]
        ax.yaxis.label.set_color('b')
        ax.tick_params(axis='y', colors='b', labelsize=minorFS)
        ax.tick_params(axis='x', labelsize=minorFS)
        ax_t.yaxis.label.set_color('r')
        ax_t.tick_params(axis='y', colors='r', labelsize=minorFS)
        ax.legend(lns, lbs, prop={"size":minorFS})

def main(args):
    # default simulation parameters
    frequency = 0.05 # MHz
    resistance = 0.004 # Gohm
    rest_distance = 30 # um
    amplitude = 20 # um
    parasitic_cap = 100 # fF

    source_file_name = 'Simulation results.csv'

    # parse cli input to a coherent unit system
    # [mV, pA, Gohm, fF, MHz, us, um, m/s]
    if '-f' in args:
        frequency = PrefixCorrection(args[args.index('-f')+1])*1e-6
    if '-R' in args:
        resistance = PrefixCorrection(args[args.index('-R')+1])*1e-9
    if '-A' in args:
        amplitude = PrefixCorrection(args[args.index('-A')+1])*1e+6
    if '-d' in args:
        rest_distance = PrefixCorrection(args[args.index('-d')+1])*1e+6
    if '-s' in args:
        source_file_name = args[args.index('-s')+1]
    if '-Cp' in args:
        parasitic_cap = PrefixCorrection(args[args.index('-Cp')+1])*1e+15

    # open simulation results
    with open(source_file_name, newline='') as file:
        reader = csv.reader(file)
        data = np.array(list(reader))

    # simulation results should be in the following format:
    # sim. time; electrode width; electrode separation; shield separation; ...
    # electrode capacitance; shield capacitance; mutual capacitance; charges

    # create array of models
    models = [Model(data[1,1], resistance, amplitude, frequency, rest_distance, parasitic_cap)]
    for row in data[1:]:
        if row[1] not in [m.w for m in models]:
            models.append(Model(row[1], resistance, amplitude, frequency, rest_distance, parasitic_cap))
        m = next((x for x in models if x.w == row[1]), None)

        # extract capacitances bellow the electrode
        # 7 columns are for other parameters and data
        boxes = np.count_nonzero(~np.isnan(row[7:].astype('float')))
        N = int((sqrt(1+8*boxes)-1)/2)
        c = np.zeros([N, N])
        index = 0
        for i in range(N):
            for j in range(i,N):
                c[j,i] = float(row[7+index])
                c[i,j] = c[j,i]
                index = index + 1

        # add data to model
        m.addData(float(row[2]), float(row[4]), float(row[5]), float(row[6]), c)

    # plot capacitances
    for m in models:
        m.plotCapacitances()

    # plot stability
    # for m in models:
    #     m.plotStability(frequency, amplitude, resistance)

    # plot stationary signal
    # for m in models:
    #    m.plotStationary()

    # plot noise
    # for m in models:
    #    m.plotNoise()

    # plot signal
    # for m in models:
    #    m.plotSignal()

    # plot sensitivity
    # for m in models:
    #     m.plotSensitivity()

    plt.show()

if __name__ == '__main__':
    args = argv[1:]
    main(args)
