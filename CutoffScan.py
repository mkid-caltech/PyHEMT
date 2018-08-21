from __future__ import division
import sys
import numpy as np
import datetime
import matplotlib.pyplot as plt
import scipy.optimize
from matplotlib.backends.backend_pdf import PdfPages

def saturation(Vgs,Vt,D,F,G,H,J,K,L):
    Ids = np.ones(len(Vgs))*0
    Ids[Vgs>Vt] = D*(Vgs[Vgs>Vt]-Vt)**2 + F*(Vgs[Vgs>Vt]-Vt)**4 + G*(Vgs[Vgs>Vt]-Vt)**6 + H*(Vgs[Vgs>Vt]-Vt)**8 + J*(Vgs[Vgs>Vt]-Vt)**10 + K*(Vgs[Vgs>Vt]-Vt)**12 + L*(Vgs[Vgs>Vt]-Vt)**14
    return Ids

#def saturation2(V,Vt,C,D,F,G,H,J,K,L,M):
#    return saturation(V[0],Vt,C,D,F,G,H,J,K,L)*(1+M*V[1])

def main(fname):
    with open(fname, 'r') as file:
        first_line = file.readline()
        second_line = file.readline()

    HEMTID = first_line.split()
    colorchoice = ['blue', 'green', 'red', 'cyan']
    Settings = second_line.split()
    NumPerHEMT = 105*int(Settings[2])*int(Settings[4])
    SearchIds = int(Settings[3])
    SearchVds = int(Settings[4])
    SearchVgs = int(Settings[5])
    CurveLength = int(Settings[6])
    data = np.loadtxt(open(fname),skiprows=2)
    Ids_Goal = [0.5,1.0,2.0,3.0]
    Vds_Goal = [100,175,250]
    device = fname.split("/")[-1].split("_")[0:2]
    location = fname.split("/")[:-1]

    chi2 = np.ones(SearchVds)
    SearchVds_vals = np.ones(SearchVds)
    Vcutoff_vals = np.ones(SearchVds)
    results = np.ones((SearchVds,8))

    plt.figure(1, figsize=(15,8.5))
    plt_back = plt.subplot(111)

    for x in range(len(HEMTID)):   # read in HEMT x's data

        Vgs = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),3]*(-1)
        Vds = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),1]
        Ids = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),2] - 0.00014503

        VgsHEMT = np.array([])
        VdsHEMT = np.array([])
        IdsHEMT = np.array([])
        IdsErrHEMT = np.array([])

        for y in range(SearchVds):
            VgsCut = Vgs[(105*y):(105*(y+1))]
            VdsCut = Vds[(105*y):(105*(y+1))]
            IdsCut = Ids[(105*y):(105*(y+1))]
            Cut = np.argmin(VdsCut)
            VgsCut = VgsCut[:Cut]
            VdsCut = VdsCut[:Cut]
            IdsCut = IdsCut[:Cut]
            IdsErr = IdsCut*(0.035/100) + 6E-4
            fit_result, pcov = scipy.optimize.curve_fit(saturation, VgsCut, IdsCut, p0=[-165, 3E-6, 0, 0, 0, 0.0001, 0, 0], sigma=IdsErr)
            results[y,:] = fit_result
            chi2[y] = sum((IdsCut-saturation(VgsCut, fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7]))**2)/len(IdsCut)

        nexttry = results[np.argmin(chi2)]

        figH = plt.figure(x+2, figsize=(15,8.5))

        figH.add_subplot(443)
        plt_zoom = plt.subplot2grid((4, 4), (1, 2), colspan=2, rowspan=3)

        figH.add_subplot(442)
        plt_main = plt.subplot2grid((4, 4), (0, 0), colspan=2, rowspan=4)

        for y in range(SearchVds):
            VgsCut = Vgs[(105*y):(105*(y+1))]
            VdsCut = Vds[(105*y):(105*(y+1))]
            IdsCut = Ids[(105*y):(105*(y+1))]
            Cut = np.argmin(VdsCut)
            VgsCut = VgsCut[:Cut]
            VdsCut = VdsCut[:Cut]
            IdsCut = IdsCut[:Cut]
            IdsErr = IdsCut*(0.035/100) + 6E-4
            fit_result, pcov = scipy.optimize.curve_fit(saturation, VgsCut, IdsCut, p0=nexttry, sigma=IdsErr)
            results[y,:] = fit_result
            Vcutoff_vals[y] = '%.1f' % fit_result[0]
            SearchVds_vals[y] = '%.1f' % np.mean(VdsCut)

            #chi2[y] = sum((IdsCut-saturation(VgsCut, fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7], fit_result[8]))**2)/len(IdsCut)
            fit_Vgs = np.arange(min(VgsCut), max(VgsCut)+0.25, 0.25)

            #plt.figure(1)

            if y == 0:
                plt_back.plot(VgsCut, IdsCut,'.', c=colorchoice[x], label='HEMT '+HEMTID[x])
                plt_back.plot(fit_Vgs, saturation(fit_Vgs, fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7]), c=colorchoice[x], label='Vt = %.1f mV ' % fit_result[0])
            else:
                plt_back.plot(VgsCut, IdsCut,'.', c=colorchoice[x])
                plt_back.plot(fit_Vgs, saturation(fit_Vgs, fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7]), c=colorchoice[x], label='Vt = %.1f mV ' % fit_result[0])

            #plt.figure(x+2)
            #plt.subplot(442)
            plt_main.semilogy(VgsCut, IdsCut,'.', c=colorchoice[y], label='Vds = %.1f mV' % np.mean(VdsCut))
            plt_main.semilogy(fit_Vgs, saturation(fit_Vgs, fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7]), c=colorchoice[y], label='Vt = %.1f mV ' % fit_result[0])
            plt_zoom.semilogy(VgsCut[(VgsCut>fit_result[0]-20) & (VgsCut<fit_result[0]+20)], IdsCut[(VgsCut>fit_result[0]-20) & (VgsCut<fit_result[0]+20)],'.', c=colorchoice[y], label='Vds = %.1f mV' % np.mean(VdsCut))
            plt_zoom.semilogy(fit_Vgs[(fit_Vgs>fit_result[0]-20) & (fit_Vgs<fit_result[0]+20)], saturation(fit_Vgs[(fit_Vgs>fit_result[0]-20) & (fit_Vgs<fit_result[0]+20)], fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7]), c=colorchoice[y], label='Vt = %.1f mV ' % fit_result[0])
            #plt.subplot(443)
            #plt.semilogy(VgsCut, IdsCut,'.', c=colorchoice[y], label='Vds = %.2f mV' % np.mean(VdsCut))
            #plt.semilogy(fit_Vgs, saturation(fit_Vgs, fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7], fit_result[8]), c=colorchoice[y], label='Vt = %.2f mV ' % fit_result[0])

        plt_zoom.set_title('$I_{ds}$ vs $V_{gs}$ (zoomed)')
        plt_zoom.set_xlabel('$V_{gs}$[mV]')
        plt_zoom.autoscale(tight=True)
        plt_main.set_title('$I_{ds}$ vs $V_{gs}$: HEMT ' + HEMTID[x])
        plt_main.set_ylabel('$I_{ds}[mA]$')
        plt_main.set_xlabel('$V_{gs}$[mV]')
        plt_main.legend(loc='best', fancybox=True)

        plt_table = figH.add_subplot(444)
        plt.subplot2grid((4, 4), (0, 2), colspan=2, rowspan=1)
        ytable = plt.table(cellText=np.transpose([SearchVds_vals,Vcutoff_vals]), loc='center', cellLoc='center', colLabels=['$V_{ds}$ [mV]', '$V_{t}$ [mV]'], cellColours=SearchVds*[['cyan', 'white']])
        ytable.scale(0.9, 1.5)
        plt.axis('off')
        plt.text(0,1, 'Data taken from HEMT ' + HEMTID[x] + ' in ' + fname.split("/")[-1] + '\n Plot made on ' + str(datetime.datetime.now()))


        if location:
            OutNamepdf = '/'.join(location)+'/'+device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_Vt.pdf'
        else:
            OutNamepdf = device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_Vt.pdf'
        Res_pdf = PdfPages(OutNamepdf)

        Res_pdf.savefig(x+2)
        Res_pdf.close()

        #nexttry = results[np.argmin(chi2)]
        #plt.legend(loc='best', fancybox=True)

        #plt.figure(x+2)
        ##plt.semilogy(VgsHEMT,IdsHEMT,'.')
        #fit_result, pcov = scipy.optimize.curve_fit(saturation2, [VgsHEMT,VdsHEMT], IdsHEMT, p0=[nexttry[0], nexttry[1], nexttry[2], nexttry[3], nexttry[4], nexttry[5], nexttry[6], nexttry[7], nexttry[8], 0.001], sigma=IdsErrHEMT)
        #print sum((IdsHEMT-saturation2([VgsHEMT,VdsHEMT], fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7], fit_result[8], fit_result[9]))**2)/len(IdsHEMT)
        #plt.semilogy(VgsHEMT,saturation2([VgsHEMT,VdsHEMT], fit_result[0], fit_result[1], fit_result[2], fit_result[3], fit_result[4], fit_result[5], fit_result[6], fit_result[7], fit_result[8], fit_result[9]), label='Vt = %.2f mV ' % fit_result[0])

    #plt.figure(1)
    plt_back.set_title('$I_{ds}$ vs $V_{gs}$: All HEMTs')
    plt_back.set_ylabel('$I_{ds}$[mA]')
    plt_back.set_xlabel('$V_{gs}$[mV]')
    plt_back.legend(loc='best', fancybox=True)
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1])
