import numpy as np # numerical python--handles arrays
import scipy.optimize # optimization package within scipy -- use for curve fitting
import matplotlib.pyplot as plt # plotting package with all the same functions as matlab

def _f(Vgs, Vds, beta, Vt ,lamb):
    return beta*Vds+lamb
    #if (Vgs-Vt)<=0:                                      # off
    #    return 0
    #elif 0<=Vds<(Vgs-Vt):                                # linear
    #    return beta*Vds*(2*(Vgs-Vt)-Vds)*(1+lamb*Vds)
    #elif 0<(Vgs-Vt)<=Vds:                                # saturation
    #    return beta*(Vgs-Vt)**2*(1+lamb*Vds)
    #else:
    #    print 'error 01'
    #    return

f = np.vectorize(_f, excluded=["beta", "Vt", "lamb"])

def f2(xdata, beta, Vt, lamb):
    return f(xdata[0], xdata[1], beta, Vt, lamb)

def main(fname):
    with open(fname, 'r') as f:
        first_line = f.readline()
        second_line = f.readline()
    Nmax = second_line.split()
    ke2220Nmax = int(Nmax[0])
    ke2401Nmax = int(Nmax[1])
    NperPoint = int(Nmax[2])
    NumPerHEMT = ke2220Nmax*ke2401Nmax*NperPoint
    data = np.loadtxt(open(fname),skiprows=2)

    HEMTID = ['1A','1B','1C','1D','2A','2B','2C','2D','3A','3B','3C','3D','4A','4B','4C','4D','5A','5B','5C','5D']

    # run once for each HEMT
    for x in range(20):
        # read in HEMT x's data
        Vgs = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),0]*(-1)
        Vds = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),1]
        Ids = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),2]

        # Keithley 2220G-30-1 has 1 mV rms noise (in all ranges)
        #Vgserr = np.ones(NumPerHEMT)

        # Keithley 2401 has 5 microV peak to peak noise (in 200 mV range) when sourcing voltage. Round up to 2 microV rms
        #Vdserr = np.ones(NumPerHEMT)*0.002

        # Keithley 2401 has 0.035% uncertainty when measuring current (in 10 mA range)
        #Idserr = -Ids+4

        # do the minimization
        #xdata=np.vstack((Vgs,Vds))
        #p=[1,0,0]
        #v, pcov = scipy.optimize.curve_fit(f2, xdata, Ids, p0=p, sigma=Idserr, absolute_sigma=True)
        #verr = np.sqrt(np.diag(pcov))    # one standard deviation errors

        # get best-fit parameters and errors on those parameters
        #beta_opt = float(v[0])
        #Vt_opt = float(v[1])
        #lamb_opt = float(v[2])
        #beta_opterr = float(verr[0])
        #Vt_opterr = float(verr[1])
        #lamb_opterr = float(verr[2])

        plt.figure()
        plt.title('$I_{ds}$ vs $V_{ds}$ for all $V_{gs}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$[mA]')
        plt.xlabel('$V_{ds}$[mV]')

        #paramStringBeta = r'$\beta$ = ' + str(round(beta_opt,8)) + ' mA/mV'
        #paramStringVt = 'R = ' + str(round((1/beta_opt),3)) + ' Ohms'
        #paramStringLambda = '$\lambda$ = ' + str(round(lamb_opt,6)) + ' mA'

        #plt.figtext(.2,.13,paramStringBeta)
        #plt.figtext(.6,.18,paramStringVt)
        #plt.figtext(.6,.13,paramStringLambda)

        plt.plot(Vds, Ids, 'b.')
        #for i in range(ke2220Nmax):
        #    plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], f2(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], beta_opt, Vt_opt, lamb_opt), 'r')
        plt.show()

    return
