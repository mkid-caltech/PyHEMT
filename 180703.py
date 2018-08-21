import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    fname1 = "../../Data/Run1_delete_1_02Jul2018_16h15.txt"
    fname2 = "../../Data/Run1_delete_1_03Jul2018_11h58.txt"
    fname3 = "../../Data/Run1_delete_1_03Jul2018_13h26.txt"
    fname4 = "../../Data/Run1_delete_1_03Jul2018_13h55.txt"
    fname5 = "../../Data/Run1_delete_1_03Jul2018_14h07.txt"
    fname6 = "../../Data/Run1_delete_1_03Jul2018_14h50.txt"
    fname7 = "../../Data/Run1_delete_1_03Jul2018_15h04.txt"

    fname = [fname1, fname2, fname3, fname4, fname5, fname6, fname7]

    for x in range(len(fname)):

        with open(fname[x], 'r') as f:
            first_line = f.readline()
            second_line = f.readline()
        Nmax = second_line.split()
        ke2220Nmax = int(Nmax[0])
        ke2401Nmax = int(Nmax[1])
        NperPoint = int(Nmax[2])
        NumPerHEMT = ke2220Nmax*ke2401Nmax*NperPoint
        data = np.loadtxt(open(fname[x]),skiprows=2)

        Vgs = data[:,0]*(+1)
        Vds = data[:,1]
        Ids = data[:,2]
        Vgs_34401 = data[:,3]


        plt.figure(1)
        plt.title('$V_{gs}$ comparison')
        plt.ylabel('Ke 2220 $V_{gs} - $ HP 34401 $V_{gs}$[mV]')
        plt.xlabel('HP 34401 $V_{gs}$[mV]')
        plt.plot(Vgs_34401, Vgs-Vgs_34401, '.',label=fname[x].split('_')[-1])

    plt.ylim(-0.2,2.2)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
