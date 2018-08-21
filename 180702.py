import sys
import numpy as np
import matplotlib.pyplot as plt

def main(fname, fname2=None):
    if fname2:
        with open(fname2, 'r') as f2:
            first_line2 = f2.readline()
            second_line2 = f2.readline()
        Nmax2 = second_line2.split()
        ke2220Nmax2 = int(Nmax2[0])
        ke2401Nmax2 = int(Nmax2[1])
        NperPoint2 = int(Nmax2[2])
        NumPerHEMT2 = ke2220Nmax2*ke2401Nmax2*NperPoint2
        data2 = np.loadtxt(open(fname2),skiprows=2)

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

    for x in range(1):
        # read in HEMT x's data
        Vgs = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),0]*(+1)
        Vds = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),1]
        Ids = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),2]
        Vds_34401 = data[(NumPerHEMT*x):(NumPerHEMT*(x+1)),3]

        if fname2:
            Vgs2 = data2[(NumPerHEMT2*x):(NumPerHEMT2*(x+1)),0]*(+1)
            Vds2 = data2[(NumPerHEMT2*x):(NumPerHEMT2*(x+1)),1]
            Ids2 = data2[(NumPerHEMT2*x):(NumPerHEMT2*(x+1)),2]
            Vgs2_34401 = data2[(NumPerHEMT2*x):(NumPerHEMT2*(x+1)),3]*(+1)

        plt.figure(1)
        plt.title('$V_{ds}$ comparison')
        plt.ylabel('ke 2401 $V_{ds} - $ HP 34401 $V_{ds}$[mV]')
        plt.xlabel('HP 34401 $V_{ds}$[mV]')
        #plt.plot(Vds, Vds_34401, '.')
        #plt.plot(Vds-4*range(261))
        #plt.plot(Vds_34401-4*range(261))
        #plt.plot(4*range(261), label='input Vds')
        #plt.plot(Vds_34401, Vds, label='Ke 2401 Vds')
        plt.plot(Vds_34401, Vds-Vds_34401, '.', label='(Ke 2401 Vds)-(HP 34401 Vds)')
        #plt.plot(Vds_34401, label='HP 34401 Vds')
        plt.ylim(-0.2,0.9)
        #plt.legend()
        #plt.figure(3)
        #plt.plot(np.array(4*range(261))-np.array(4*range(261)), label='input Vds')
        #plt.plot(np.array(4*range(261)),Vds-4*range(261), '.', label='Ke 2401 Vds')
        #plt.plot(np.array(4*range(261)),Vds_34401-4*range(261), '.', label='HP 34401 Vds')

        if fname2:
            plt.figure(2)
            plt.title('$V_{gs}$ comparison')
            plt.ylabel('Ke 2220 $V_{gs} - $ HP 34401 $V_{gs}$[mV]')
            plt.xlabel('HP 34401 $V_{gs}$[mV]')
            #plt.plot(Vgs2_34401, Vgs2, label='Ke 2220 Vgs')
            plt.plot(Vgs2_34401, Vgs2-Vgs2_34401, '.',label='(Ke 2220 Vgs)-(HP 34401 Vgs)')
            #plt.plot((np.arange(200,402,1)/2), label='input Vgs')
            #plt.plot(Vgs2, label='Ke 2220 Vgs')
            #plt.plot(Vgs2_34401, label='HP 34401 Vgs')
            plt.ylim(-0.2,2.2)
            #plt.legend()

        plt.show()

if __name__ == "__main__":
    if len(sys.argv)>2:
        main(sys.argv[1],fname2=sys.argv[2])
    else:
        main(sys.argv[1])
