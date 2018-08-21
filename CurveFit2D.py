import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def JFET_piecewise(Vgs, Vds, beta, Vt, lamb):
    if (Vgs-Vt) <= 0:
        # off
        return 0.0
    elif 0 <= Vds < (Vgs-Vt):
        # linear
        return beta*Vds*(2*(Vgs-Vt)-Vds)*(1+lamb*Vds)
    elif 0 < (Vgs-Vt) <= Vds:
        # saturation
        return beta*((Vgs-Vt)**2)*(1+lamb*Vds)
    else:
        print 'error 01'

def dJFET_piecewise(Vgs, Vds, beta, Vt, lamb):
    if (Vgs-Vt) <= 0:
        # off
        return 0.0
    elif 0 <= Vds < (Vgs-Vt):
        # linear
        return beta*Vds*(2*(Vgs-Vt)-Vds)*(lamb) + beta*(2*(Vgs-Vt)-Vds)*(1+lamb*Vds) + beta*Vds*(-1)*(1+lamb*Vds)
    elif 0 < (Vgs-Vt) <= Vds:
        # saturation
        return beta*((Vgs-Vt)**2)*(lamb)
    else:
        print 'error 01'

JFET_vector = np.vectorize(JFET_piecewise, excluded=["beta", "Vt", "lamb"])

dJFET_vector = np.vectorize(dJFET_piecewise, excluded=["beta", "Vt", "lamb"])

def JFET_fit(x, beta, Vt, lamb):
    return JFET_vector(x[0], x[1], beta, Vt, lamb)

def dJFET_fit(x, beta, Vt, lamb):
    return dJFET_vector(x[0], x[1], beta, Vt, lamb)

def gm_fit(x, beta, Vt):
    return beta*(x-Vt)**2

def gm_func(x, beta, Vt):
    return 2*beta*(x-Vt)

def HEMTfit(fname):
    with open(fname, 'r') as fyle:
        first_line = fyle.readline()
        second_line = fyle.readline()
    Nmax = second_line.split()
    ke2220Nmax = int(Nmax[0])
    ke2401Nmax = int(Nmax[1])
    NperPoint = int(Nmax[2])
    NumPerHEMT = ke2220Nmax*ke2401Nmax*NperPoint
    data = np.loadtxt(open(fname),skiprows=2)

    HEMTID = ['1A','1B','1C','1D','2A','2B','2C','2D','3A','3B','3C','3D','4A','4B','4C','4D','5A','5B','5C','5D']

    # run once for each HEMT
    #for x in range(len(HEMTID)):
    for x in range(4):
        # read in HEMT x's data
        Vgs = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),0]*(-1)
        Vds = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),1]
        Ids = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),2]
        xdata = np.vstack((Vgs,Vds))

        # Keithley 2220G-30-1 has 1 mV rms noise (in all ranges)
        Vgserr = np.ones(NumPerHEMT)

        # Keithley 2401 has 5 microV peak to peak noise (in 200 mV range) when sourcing voltage. Round up to 2 microV rms
        Vdserr = np.ones(NumPerHEMT)*0.002

        # Keithley 2401 has (0.035% + 600 nA) uncertainty when measuring current (in 10 mA range)
        Idserr = (abs(Ids)*0.01*0.035) + 6E-4

        valu1s = np.array([])
        for k in range(ke2220Nmax):
            valu1s = np.append(valu1s, k*ke2401Nmax*NperPoint + np.arange(NperPoint))

        plt.figure(1, figsize=(8,7.5))
        plt.figure(2, figsize=(8,5))
        for i in range(ke2401Nmax):
            valu2s = i*NperPoint + valu1s
            valu3s = map(int,valu2s)
            gm_result, pcov3 = scipy.optimize.curve_fit(gm_fit, Vgs[valu3s], Ids[valu3s], p0=[2.56E-4,-120])
            plt.figure(1)
            if i == 0:
                plt.subplot(2,1,1) # Figure 1
                plt.plot(Vgs[valu3s], Ids[valu3s], 'g.')
                plt.plot(Vgs[valu3s], gm_fit(Vgs[valu3s], gm_result[0], gm_result[1]), 'g', linewidth=3, label='$V_{ds}$ = '+str(round(Vds[valu3s[0]],3))+' mV')
                plt.subplot(2,1,2) # Figure 1
                plt.plot(Ids[valu3s], gm_func(Vgs[valu3s], gm_result[0], gm_result[1]), 'g', linewidth=3)
            elif i == (ke2401Nmax-1):
                plt.subplot(2,1,1) # Figure 1
                plt.plot(Vgs[valu3s], Ids[valu3s], 'c.')
                plt.plot(Vgs[valu3s], gm_fit(Vgs[valu3s], gm_result[0], gm_result[1]), 'c', linewidth=3, label='$V_{ds}$ = '+str(round(Vds[valu3s[0]],3))+' mV')
                plt.subplot(2,1,2) # Figure 1
                plt.plot(Ids[valu3s], gm_func(Vgs[valu3s], gm_result[0], gm_result[1]), 'c', linewidth=3)
            else:
                plt.subplot(2,1,1) # Figure 1
                plt.plot(Vgs[valu3s], Ids[valu3s], 'b.')
                plt.plot(Vgs[valu3s], gm_fit(Vgs[valu3s], gm_result[0], gm_result[1]),'r')
                plt.subplot(2,1,2) # Figure 1
                plt.plot(Ids[valu3s], gm_func(Vgs[valu3s], gm_result[0], gm_result[1]),'r')
            plt.figure(2)
            plt.plot(Vds[valu3s[0]],2*np.sqrt(gm_result[0])*1E3,'.')

        plt.figure(1)
        plt.subplot(2,1,1)
        plt.title('$I_{ds}$($V_{gs}$) and $g_{m}$($I_{ds}$) for each $V_{ds}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$ [mA]')
        plt.xlabel('$V_{gs}$ [mV]')
        plt.legend(loc='best')
        plt.subplot(2,1,2)
        plt.ylabel('$g_{m}$ [mA/mV]')
        plt.xlabel('$I_{ds}$ [mA]')
        plt.plot(np.arange(0,max(Ids),0.01),2*np.sqrt(gm_result[0])*np.sqrt(np.arange(0,max(Ids),0.01)))

        plt.figure(2)
        plt.title(r'2*$\sqrt{\beta}$($V_{ds}$): HEMT ' + HEMTID[x])
        plt.ylabel(r'2*$\sqrt{\beta}$  [$\sqrt{mA}$/V]')
        plt.xlabel('$V_{ds}$ [mV]')
        plt.figtext(.2,.15,'$g_{m}$ = '+str(round(2*np.sqrt(gm_result[0])*1E3,3))+'$\sqrt{I_{ds}}$ mS')
        plt.figtext(.55,.15,r'$\beta$ = '+str(round(gm_result[0],8)) + ' mA/m$V^2$')
        plt.figtext(.55,.25,r'$V_{t}$ = '+str(round(gm_result[1],8)) + ' mV')

        plt.figure(4, figsize=(8,7.5))
        guess0 = [gm_result[0], gm_result[1], 0] # Beta, Vt, lambda
        for i in range(ke2220Nmax):
            #valu2s = i*NperPoint + valu1s
            #valu3s = map(int,valu2s)
            #gm_result, pcov3 = scipy.optimize.curve_fit(gm_fit, Vgs[valu3s], Ids[valu3s], p0=[2.56E-4,-120])
            lamb_result, pcov4 = scipy.optimize.curve_fit(JFET_fit, xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], p0=guess0)
            if i == 0:
                plt.subplot(2,1,1) # Figure 1
                plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], 'g.')
                plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'g', linewidth=3, label='$V_{gs}$ = '+str(round(Vgs[(ke2401Nmax*NperPoint*i)],3))+' mV')
                plt.subplot(2,1,2) # Figure 1
                plt.plot(Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], dJFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'g', linewidth=3)
            elif i == (ke2220Nmax-1):
                plt.subplot(2,1,1) # Figure 1
                plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], 'c.')
                plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'c', linewidth=3, label='$V_{gs}$ = '+str(round(Vgs[(ke2401Nmax*NperPoint*i)],3))+' mV')
                plt.subplot(2,1,2) # Figure 1
                plt.plot(Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], dJFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'c', linewidth=3)
            else:
                plt.subplot(2,1,1) # Figure 1
                plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], 'b.')
                plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'r')
                plt.subplot(2,1,2) # Figure 1
                plt.plot(Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], dJFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'r')

        plt.subplot(2,1,1)
        plt.title('$I_{ds}$($V_{ds}$) and $g_{d}$($I_{ds}$) for each $V_{gs}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$ [mA]')
        plt.xlabel('$V_{ds}$ [mV]')
        plt.legend(loc='best')
        plt.subplot(2,1,2)
        plt.ylabel('$g_{ds}$ [mA/mV]')
        plt.xlabel('$I_{ds}$ [mA]')

        # do the minimization
        #guess0 = [2.56E-4,-200,8.55E-4] # Beta, Vt, lambda
        guess0 = [gm_result[0],gm_result[1],0] # Beta, Vt, lambda
        guess1, pcov1 = scipy.optimize.curve_fit(JFET_fit, xdata[:,:(ke2401Nmax*NperPoint)], Ids[:(ke2401Nmax*NperPoint)], p0=guess0, sigma=Idserr[:(ke2401Nmax*NperPoint)], absolute_sigma=True)
        result, pcov2 = scipy.optimize.curve_fit(JFET_fit, xdata, Ids, p0=guess1, sigma=Idserr, absolute_sigma=True)
        result_error = np.sqrt(np.diag(pcov2))    # one standard deviation errors

        # get best-fit parameters and errors on those parameters
        beta_opt = float(result[0])
        Vt_opt = float(result[1])
        lamb_opt = float(result[2])
        beta_opterr = float(result_error[0])
        Vt_opterr = float(result_error[1])
        lamb_opterr = float(result_error[2])

        plt.figure(3, figsize=(10,6))
        plt.plot(Vds, Ids, 'b.', label='data')
        plt.plot(100, 1, 's', color='black')
        for i in range(ke2220Nmax):
            plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], beta_opt, Vt_opt, lamb_opt), linewidth=2, label='$V_{gs}$ = '+str(round(Vgs[(ke2401Nmax*NperPoint*i)+1],3))+' mV')

        paramStringBeta = r'$\beta$ = ' + str(round(beta_opt,8)) + ' mA/m$V^2$'
        paramStringVt = '$V_{threshold}$ = ' + str(round(Vt_opt,3)) + ' m$V$'
        paramStringLambda = '$\lambda$ = ' + str(round(lamb_opt,6)) + ' m$V^{-1}$'

        plt.title('$I_{ds}$ vs $V_{ds}$ for each $V_{gs}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$ [mA]')
        plt.xlabel('$V_{ds}$ [mV]')
        plt.xlim(0,250)
        plt.ylim(-0.5,)

        plt.figtext(.14,.11,paramStringBeta)
        plt.figtext(.4,.15,paramStringVt)
        plt.figtext(.4,.11,paramStringLambda)
        plt.legend()

        #plt.figure(4)
        #for k in range(NumPerHEMT):
        #    if (Vgs[k] - Vt_opt) <= 0:
        #        # off
        #        zed=k
        #    elif Vds[k] <= (Vgs[k] - Vt_opt):
        #        # linear
        #        zed=k
        #    elif (Vgs[k] - Vt_opt) <= Vds[k]:
        #        # saturation
        #        plt.plot(Vgs[k],Ids[k],'.')

        plt.show()

if __name__ == "__main__":
    HEMTfit(sys.argv[1])
