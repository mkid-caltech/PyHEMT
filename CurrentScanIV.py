import sys
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages

def main(fname, fname2=None):

    with open(fname, 'r') as file:
        first_line = file.readline()
        second_line = file.readline()
    HEMTID = first_line.split()
    Settings = second_line.split()
    NumPerHEMT = int(Settings[0])*int(Settings[1])*int(Settings[2])
    SearchIds = int(Settings[3])
    SearchVds = int(Settings[4])
    SearchVgs = int(Settings[5])
    CurveLength = int(Settings[6])
    data = np.loadtxt(open(fname),skiprows=2)
    timing = str(datetime.now())
    Ids_Goal = [0.5,1.0,2.0,3.0]
    if SearchVds == 3:
        Vds_Goal = [100,175,250]
    elif SearchVds == 4:
        Vds_Goal = [75,100,175,250]
    device = fname.split("/")[-1].split("_")[0:2]
    location = fname.split("/")[:-1]

    cmap = plt.cm.jet
    norm = Normalize(vmin=-152, vmax=-62)
    FakePlot = np.arange(0,1)
    scalarmappaple = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalarmappaple.set_array(FakePlot)

    for x in range(len(HEMTID)):   # read in HEMT x's data
    #for x in range(1):
        timing_note = 'Data taken from HEMT ' + HEMTID[x] + ' in ' + fname.split("/")[-1] + '\n Plot made on ' + timing
        Vgs = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),3]*(-1)
        Vds = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),1]
        Ids = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),2] - 0.00014503

        gd = np.zeros((SearchIds, SearchVds+1))
        gm = np.zeros((SearchIds, SearchVds+1))
        bias = np.zeros((SearchVds*2*SearchIds, SearchVds+1))

        plt.figure(2, figsize=(15,8.5))
        plt.subplot(412)
        plt.subplot2grid((4, 4), (1, 0), colspan=4, rowspan=3)

        plt.title('$I_{ds}$ vs $V_{gs}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$[mA]')
        plt.xlabel('$V_{gs}$[mV]')
        plt.figure(1, figsize=(15,8.5))
        plt.subplot(412)
        plt.subplot2grid((4, 4), (1, 0), colspan=4, rowspan=3)
        plt.title('$I_{ds}$ vs $V_{ds}$ for all $V_{gs}$: HEMT ' + HEMTID[x])
        plt.ylabel('$I_{ds}$[mA]')
        plt.xlabel('$V_{ds}$[mV]')
        for y in range(SearchIds):
            gd[y,0] = Ids_Goal[y]
            gm[y,0] = Ids_Goal[y]
            plt.gca().set_prop_cycle(None)
            plt.figure(1)
            plt.gca().set_prop_cycle(None)
            for z in range(SearchVds):
                distance = np.sqrt(((Ids[(int(Settings[1])*(y+z*SearchIds))+SearchVgs:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+2*CurveLength]-Ids_Goal[y])**2)+((Vds[(int(Settings[1])*(y+z*SearchIds))+SearchVgs:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+2*CurveLength]-Vds_Goal[z])**2))
                closest = ((int(Settings[1])*(y+z*SearchIds))+SearchVgs) + np.argmin(distance)
                distance = np.sqrt(((Ids[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs]-Ids_Goal[y])**2)+((Vds[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs]-Vds_Goal[z])**2))
                upper = (int(Settings[1])*(y+z*SearchIds)) + np.argmin(abs(distance-min(distance[np.where(Ids[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs]-Ids_Goal[y] >= 0)])))
                lower = (int(Settings[1])*(y+z*SearchIds)) + np.argmin(abs(distance-min(distance[np.where(Ids[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs]-Ids_Goal[y] <= 0)])))
                gd_fit = np.polyfit(Vds[closest-2:closest+3],Ids[closest-2:closest+3],1)
                gd[y,z+1] = np.around(1000*gd_fit[0], decimals=3) # in mS
                gm[y,z+1] = np.around(1000*((Ids[upper]-Ids[lower])/(Vgs[upper]-Vgs[lower])), decimals=3) # in mS
                plt.figure(1)
                plt.plot(Vds[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], Ids[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs],'k.')
                Curve_Vgs = Vgs[(int(Settings[1])*(y+z*SearchIds))+SearchVgs:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+CurveLength]
                Curve_Vds = Vds[(int(Settings[1])*(y+z*SearchIds))+SearchVgs:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+CurveLength]
                Curve_Ids = Ids[(int(Settings[1])*(y+z*SearchIds))+SearchVgs:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+CurveLength]
                bias[SearchVds*2*y+2*z,0] = Curve_Ids[np.argmin(abs(Curve_Vds-Vds_Goal[z]))]
                bias[SearchVds*2*y+2*z,z+1] = np.mean(Curve_Vgs)
                plt.plot(Curve_Vds, Curve_Ids, linewidth=2.0, color=cmap(norm(np.mean(Curve_Vgs))), label="%.1f mV" % np.mean(Curve_Vgs))
                Curve_Vgs = Vgs[(int(Settings[1])*(y+z*SearchIds))+SearchVgs+CurveLength:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+2*CurveLength]
                Curve_Vds = Vds[(int(Settings[1])*(y+z*SearchIds))+SearchVgs+CurveLength:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+2*CurveLength]
                Curve_Ids = Ids[(int(Settings[1])*(y+z*SearchIds))+SearchVgs+CurveLength:(int(Settings[1])*(y+z*SearchIds))+SearchVgs+2*CurveLength]
                bias[SearchVds*2*y+2*z+1,0] = Curve_Ids[np.argmin(abs(Curve_Vds-Vds_Goal[z]))]
                bias[SearchVds*2*y+2*z+1,z+1] = np.mean(Curve_Vgs)
                plt.plot(Curve_Vds, Curve_Ids, linewidth=2.0, color=cmap(norm(np.mean(Curve_Vgs))), label="%.1f mV" % np.mean(Curve_Vgs))
                plt.plot(Vds[closest-2:closest+3], gd_fit[0]*Vds[closest-2:closest+3]+gd_fit[1], 'k', linewidth=3.0, zorder=100)
                plt.figure(2)
                if y == 0:
                    plt.plot(Vgs[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], Ids[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], '.-', label="Ids @ Vds = %.1f mV" % np.mean(Vds[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs]))
                else:
                    plt.plot(Vgs[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], Ids[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], '.-')
                #plt.plot(Vgs[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], 0.01*Vds[(int(Settings[1])*(y+z*SearchIds)):(int(Settings[1])*(y+z*SearchIds))+SearchVgs], ls='-', marker='.')
                plt.plot(Vgs[lower:upper+1], Ids[lower:upper+1], 'k', linewidth=3.0, zorder=100)
        plt.figure(1)
        for goal in Ids_Goal:
            plt.axhline(y=goal, c='k', ls='--', zorder=0)
        plt.xlim(0, 1.01*max(Vds))
        Vgsbar = plt.colorbar(scalarmappaple)
        Vgsbar.set_label('$V_{gs}$ [mV]')
        plt.figure(2)
        for goal in Ids_Goal:
            plt.axhline(y=goal, c='k', ls='--', zorder=0)
        plt.legend(loc='best', fancybox=True)


        #print "\n"
        #print "   HEMT %s transconductance (gm)" % HEMTID[x]
        #print "  ------------------------------------"
        #print "   Ids\Vds |    100 mV   |    175 mV   |    250 mV   "
        #print "  ---------+-------------+-------------+-------------"
        #for row in range(SearchIds):
        #    print "   %3.1f mA  |  %07.3f mS |  %07.3f mS |  %07.3f mS " % (gm[row,0],gm[row,1],gm[row,2],gm[row,3])

        #print "\n"
        #print "   HEMT %s drain conductance (gd)" % HEMTID[x]
        #print "  ------------------------------------"
        #print "   Ids\Vds |    100 mV   |    175 mV   |    250 mV   "
        #print "  ---------+-------------+-------------+-------------"
        #for row in range(SearchIds):
        #    print "   %3.1f mA |  %07.3f mS |  %07.3f mS |  %07.3f mS " % (gd[row,0],gd[row,1],gd[row,2],gd[row,3])

        #print "\n"
        #print "    HEMT %s bias voltages (Vgs)" % HEMTID[x]
        #print "  ------------------------------------"
        #print "    Ids\Vds  |    100 mV   |    175 mV   |    250 mV   "
        #print "  -----------+-------------+-------------+-------------"
        #for row in range(len(bias[:,0])):
        #    if bias[row,2] == 0 and bias[row,3] ==0:
        #        print "   %6.4f mA |  %07.2f mV |             |             " % (bias[row,0],bias[row,1])
        #    if bias[row,1] == 0 and bias[row,3] ==0:
        #        print "   %6.4f mA |             |  %07.2f mV |             " % (bias[row,0],bias[row,2])
        #    if bias[row,1] == 0 and bias[row,2] ==0:
        #        print "   %6.4f mA |             |             |  %07.2f mV " % (bias[row,0],bias[row,3])

        bias[bias == 0] = np.nan
        words = [map(str, np.around(bias[:,0], decimals=3))]
        for v in range(SearchVds):
            words = words + [map(str, np.around(bias[:,v+1], decimals=1))]

        for f in range(len(words)):
            for g in range(len(words[f])):
                if words[f][g] == 'nan':
                    words[f][g] = ''
        words = np.array(words).transpose()

        if location:
            OutNamegm = '/'.join(location)+'/'+device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_gm.csv'
            OutNamegd = '/'.join(location)+'/'+device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_gd.csv'
            OutNamebias = '/'.join(location)+'/'+device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_bias.csv'
            OutNamepdf = '/'.join(location)+'/'+device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'.pdf'
        else:
            OutNamegm = device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_gm.csv'
            OutNamegd = device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_gd.csv'
            OutNamebias = device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'_bias.csv'
            OutNamepdf = device[0]+'_'+device[1]+'_HEMT'+HEMTID[x]+'.pdf'

        np.savetxt(OutNamegm, gm, delimiter=', ', header=timing_note + '\nIds [mA],  gm [mS] @ Vds = 100 mV, gm [mS] @ Vds = 175 mV, gm [mS] @ Vds = 250 mV')
        np.savetxt(OutNamegd, gd, delimiter=', ', header=timing_note + '\nIds [mA],  gd [mS] @ Vds = 100 mV, gd [mS] @ Vds = 175 mV, gd [mS] @ Vds = 250 mV')
        np.savetxt(OutNamebias, bias, delimiter=', ', header=timing_note + '\nIds [mA],  bias [mV] @ Vds = 100 mV, bias [mV] @ Vds = 175 mV, bias [mV] @ Vds = 250 mV')

        pattern = ['cyan']+SearchVds*['white']

        labelsgd = ['$I_{ds}$ [mA]']
        labelsgm = ['$I_{ds}$ [mA]']
        labelsbias = ['$I_{ds}$ [mA]']
        for value in Vds_Goal:
            labelsgd = labelsgd + ['$g_{d}$ [mS] @ $V_{ds}$ = %.0f mV' % value]
            labelsgm = labelsgm + ['$g_{m}$ [mS] @ $V_{ds}$ = %.0f mV' % value]
            labelsbias = labelsbias + ['$V_{gs}$ [mS] @ $V_{ds}$ = %.0f mV' % value]

        plt.figure(1)
        plt.subplot(411)
        ytable = plt.table(cellText = gd, loc='center', cellLoc='center', colLabels=labelsgd, cellColours=SearchIds*[pattern])
        plt.text(0,1, timing_note)
        ytable.scale(0.9, 1.5)
        plt.axis('off')

        plt.figure(2)
        plt.subplot(411)
        ytable = plt.table(cellText = gm, loc='center', cellLoc='center', colLabels=labelsgm, cellColours=SearchIds*[pattern])
        ytable.scale(0.9, 1.5)
        plt.axis('off')

        plt.figure(3, figsize=(15,8.5*((1+2.*SearchIds*SearchVds)/25)))
        ytable = plt.table(cellText = words, loc='center', cellLoc='center', colLabels=labelsbias, cellColours=2*SearchIds*SearchVds*[pattern])
        ytable.scale(1, 2)
        plt.axis('off')

        Res_pdf = PdfPages(OutNamepdf)
        Res_pdf.savefig(1)
        Res_pdf.savefig(2)
        Res_pdf.savefig(3)
        Res_pdf.close()
        plt.close('all')

if __name__ == "__main__":
    if len(sys.argv)>2:
        main(sys.argv[1],fname2=sys.argv[2])
    else:
        main(sys.argv[1])
