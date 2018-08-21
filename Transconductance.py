import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def target(x):
    for k in range(len(x)):
        if x[k] < 0:
            x[k] = 1000
        else:
            x[k] = x[k]
    return x

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

Ids_Goal = 1

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
    dataAvg = np.zeros((ke2220Nmax*ke2401Nmax,3))

    HEMTID = ['1A','1B','1C','1D','2A','2B','2C','2D','3A','3B','3C','3D','4A','4B','4C','4D','5A','5B','5C','5D']

    # run once for each HEMT
    #for x in range(len(HEMTID)):

    for k in range(ke2220Nmax*ke2401Nmax):
        kAvg = np.array((sum(data[k*NperPoint:(k+1)*NperPoint,0])/NperPoint, sum(data[k*NperPoint:(k+1)*NperPoint,1])/NperPoint, sum(data[k*NperPoint:(k+1)*NperPoint,2])/NperPoint))
        dataAvg[k,:] = kAvg

    for x in range(1):
        # read in HEMT x's data
        Vgs = dataAvg[((ke2220Nmax*ke2401Nmax)*x):((ke2220Nmax*ke2401Nmax)*(x+1)),0]*(-1)
        Vds = dataAvg[((ke2220Nmax*ke2401Nmax)*x):((ke2220Nmax*ke2401Nmax)*(x+1)),1]
        Ids = dataAvg[((ke2220Nmax*ke2401Nmax)*x):((ke2220Nmax*ke2401Nmax)*(x+1)),2]
        xdata = np.vstack((Vgs,Vds))

        closest = np.argmin(((Ids-Ids_Goal)**2)+((Vds-100)**2))
        upper = np.argmin((target(Ids-Ids_Goal)**2)+((Vds-100)**2))
        lower = np.argmin((target(Ids_Goal-Ids)**2)+((Vds-100)**2))
        #left = closest - 5
        #right = closest + 5

        # Keithley 2220G-30-1 has 1 mV rms noise (in all ranges)
        Vgserr = np.ones(NumPerHEMT)

        # Keithley 2401 has 5 microV peak to peak noise (in 200 mV range) when sourcing voltage. Round up to 2 microV rms
        Vdserr = np.ones(NumPerHEMT)*0.002

        # Keithley 2401 has (0.035% + 600 nA) uncertainty when measuring current (in 10 mA range)
        Idserr = (abs(Ids)*0.01*0.035) + 6E-4

        #valu1s = np.array([])
        #for k in range(ke2220Nmax):
        #    valu1s = np.append(valu1s, k*ke2401Nmax*NperPoint + np.arange(NperPoint))

        #plt.figure(1, figsize=(8,7.5))
        #plt.figure(2, figsize=(8,5))
        #for i in range(ke2401Nmax):
        #    valu2s = i*NperPoint + valu1s
        #    valu3s = map(int,valu2s)
        #    gm_result, pcov3 = scipy.optimize.curve_fit(gm_fit, Vgs[valu3s], Ids[valu3s], p0=[2.56E-4,-120])
        #    plt.figure(1)
        #    if i == 0:
        #        plt.subplot(2,1,1) # Figure 1
        #        plt.plot(Vgs[valu3s], Ids[valu3s], 'g.')
        #        plt.plot(Vgs[valu3s], gm_fit(Vgs[valu3s], gm_result[0], gm_result[1]), 'g', linewidth=3, label='$V_{ds}$ = '+str(round(Vds[valu3s[0]],3))+' mV')
        #        plt.subplot(2,1,2) # Figure 1
        #        plt.plot(Ids[valu3s], gm_func(Vgs[valu3s], gm_result[0], gm_result[1]), 'g', linewidth=3)
        #    elif i == (ke2401Nmax-1):
        #        plt.subplot(2,1,1) # Figure 1
        #        plt.plot(Vgs[valu3s], Ids[valu3s], 'c.')
        #        plt.plot(Vgs[valu3s], gm_fit(Vgs[valu3s], gm_result[0], gm_result[1]), 'c', linewidth=3, label='$V_{ds}$ = '+str(round(Vds[valu3s[0]],3))+' mV')
        #        plt.subplot(2,1,2) # Figure 1
        #        plt.plot(Ids[valu3s], gm_func(Vgs[valu3s], gm_result[0], gm_result[1]), 'c', linewidth=3)
        #    else:
        #        plt.subplot(2,1,1) # Figure 1
        #        plt.plot(Vgs[valu3s], Ids[valu3s], 'b.')
        #        plt.plot(Vgs[valu3s], gm_fit(Vgs[valu3s], gm_result[0], gm_result[1]),'r')
        #        plt.subplot(2,1,2) # Figure 1
        #        plt.plot(Ids[valu3s], gm_func(Vgs[valu3s], gm_result[0], gm_result[1]),'r')
        #    plt.figure(2)
        #    plt.plot(Vds[valu3s[0]],2*np.sqrt(gm_result[0])*1E3,'.')

        #plt.figure(1)
        #plt.subplot(2,1,1)
        #plt.title('$I_{ds}$($V_{gs}$) and $g_{m}$($I_{ds}$) for each $V_{ds}$: HEMT ' + HEMTID[x])
        #plt.ylabel('$I_{ds}$ [mA]')
        #plt.xlabel('$V_{gs}$ [mV]')
        #plt.legend(loc='best')
        #plt.subplot(2,1,2)
        #plt.ylabel('$g_{m}$ [mA/mV]')
        #plt.xlabel('$I_{ds}$ [mA]')
        #plt.plot(np.arange(0,max(Ids),0.01),2*np.sqrt(gm_result[0])*np.sqrt(np.arange(0,max(Ids),0.01)))

        #plt.figure(2)
        #plt.title(r'2*$\sqrt{\beta}$($V_{ds}$): HEMT ' + HEMTID[x])
        #plt.ylabel(r'2*$\sqrt{\beta}$  [$\sqrt{mA}$/V]')
        #plt.xlabel('$V_{ds}$ [mV]')

        #plt.figure(4, figsize=(8,7.5))
        #guess0 = [gm_result[0], gm_result[1], 0] # Beta, Vt, lambda
        #for i in range(ke2220Nmax):
            #valu2s = i*NperPoint + valu1s
            #valu3s = map(int,valu2s)
            #gm_result, pcov3 = scipy.optimize.curve_fit(gm_fit, Vgs[valu3s], Ids[valu3s], p0=[2.56E-4,-120])
        #    lamb_result, pcov4 = scipy.optimize.curve_fit(JFET_fit, xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], p0=guess0)
        #    if i == 0:
        #        plt.subplot(2,1,1) # Figure 1
        #        plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], 'g.')
        #        plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'g', linewidth=3, label='$V_{gs}$ = '+str(round(Vgs[(ke2401Nmax*NperPoint*i)],3))+' mV')
        #        plt.subplot(2,1,2) # Figure 1
        #        plt.plot(Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], dJFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'g', linewidth=3)
        #    elif i == (ke2220Nmax-1):
        #        plt.subplot(2,1,1) # Figure 1
        #        plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], 'c.')
        #        plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'c', linewidth=3, label='$V_{gs}$ = '+str(round(Vgs[(ke2401Nmax*NperPoint*i)],3))+' mV')
        #        plt.subplot(2,1,2) # Figure 1
        #        plt.plot(Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], dJFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'c', linewidth=3)
        #    else:
        #        plt.subplot(2,1,1) # Figure 1
        #        plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], 'b.')
        #        plt.plot(Vds[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], JFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'r')
        #        plt.subplot(2,1,2) # Figure 1
        #        plt.plot(Ids[(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], dJFET_fit(xdata[:,(ke2401Nmax*NperPoint*i):(ke2401Nmax*NperPoint*(i+1))], lamb_result[0], lamb_result[1], lamb_result[2]), 'r')

        #plt.subplot(2,1,1)
        #plt.title('$I_{ds}$($V_{ds}$) and $g_{d}$($I_{ds}$) for each $V_{gs}$: HEMT ' + HEMTID[x])
        #plt.ylabel('$I_{ds}$ [mA]')
        #plt.xlabel('$V_{ds}$ [mV]')
        #plt.legend(loc='best')
        #plt.subplot(2,1,2)
        #plt.ylabel('$g_{ds}$ [mA/mV]')
        #plt.xlabel('$I_{ds}$ [mA]')

        # do the minimization
        #guess0 = [2.56E-4,-200,8.55E-4] # Beta, Vt, lambda
        #guess0 = [gm_result[0],gm_result[1],0] # Beta, Vt, lambda
        #guess1, pcov1 = scipy.optimize.curve_fit(JFET_fit, xdata[:,:(ke2401Nmax*NperPoint)], Ids[:(ke2401Nmax*NperPoint)], p0=guess0, sigma=Idserr[:(ke2401Nmax*NperPoint)], absolute_sigma=True)
        #result, pcov2 = scipy.optimize.curve_fit(JFET_fit, xdata, Ids, p0=guess1, sigma=Idserr, absolute_sigma=True)
        #result_error = np.sqrt(np.diag(pcov2))    # one standard deviation errors

        # get best-fit parameters and errors on those parameters
        #beta_opt = float(result[0])
        #Vt_opt = float(result[1])
        #lamb_opt = float(result[2])
        #beta_opterr = float(result_error[0])
        #Vt_opterr = float(result_error[1])
        #lamb_opterr = float(result_error[2])

        gd_fit = np.polyfit(Vds[(closest/ke2401Nmax)*ke2401Nmax:(1+(closest/ke2401Nmax))*ke2401Nmax],Ids[(closest/ke2401Nmax)*ke2401Nmax:(1+(closest/ke2401Nmax))*ke2401Nmax],1)

        stringVgs = '$V_{gs}$: ' + str((round(Vgs[upper],2), round(Vgs[lower],2))) + ' mV'
        stringgm = '$g_{m}$ = ' + str(round(1000*(Ids[upper]-Ids[lower])/(Vgs[upper]-Vgs[lower]), 2)) + ' mS'
        stringgd = '$g_{d}$ = ' + str(round(1000*gd_fit[0], 4)) + ' mS'

        plt.figure()
        plt.plot(Vds, Ids, 'b.', label='data')
        plt.plot(100, Ids_Goal, 's', color='black', label='100 mV, 1 mA')
        plt.plot(100, Ids_Goal, color='white', label=stringVgs)
        plt.plot(100, Ids_Goal, color='white', label=stringgm)
        plt.plot(100, Ids_Goal, color='white', label=stringgd)
        plt.plot(Vds[lower], Ids[lower], 's', color='red')
        plt.plot(Vds[upper], Ids[upper], 's', color='green')
        #plt.plot(Vds[left], Ids[left], 'D', color='red')
        #plt.plot(Vds[right], Ids[right], 'D', color='green')

        plt.plot(Vds[(closest/ke2401Nmax)*ke2401Nmax:(1+(closest/ke2401Nmax))*ke2401Nmax], gd_fit[0]*Vds[(closest/ke2401Nmax)*ke2401Nmax:(1+(closest/ke2401Nmax))*ke2401Nmax]+gd_fit[1])

        plt.title('$I_{ds}$ vs $V_{ds}$ for each $V_{gs}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$ [mA]')
        plt.xlabel('$V_{ds}$ [mV]')
        plt.legend(loc='best')
        #plt.xlim(0,250)
        #plt.ylim(-0.5,)

        #plt.figtext(.14,.11,paramStringBeta)
        #plt.figtext(.4,.15,paramStringVt)
        #plt.figtext(.4,.11,paramStringLambda)

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

        plt.figure()
        plt.plot([min(Vgs),max(Vgs)], [Ids_Goal,Ids_Goal], 'k--', label='1 mA')
        plt.plot(Vgs, Ids, 'b.', label='data')
        plt.plot([Vgs[upper],Vgs[lower]], [Ids[upper],Ids[lower]])
        plt.plot(Vgs[lower], Ids[lower], 's', color='red')
        plt.plot(Vgs[upper], Ids[upper], 's', color='green')
        plt.title('$I_{ds}$ vs $V_{gs}$ for each $V_{ds}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$ [mA]')
        plt.xlabel('$V_{gs}$ [mV]')
        plt.legend(loc='best')

        plt.show()

if __name__ == "__main__":
    HEMTfit(sys.argv[1])
