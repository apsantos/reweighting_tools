#!/usr/bin/env python
"""
Author : Andrew P. Santos
Date   : July, 2015
Institution: Princeton University

Purpose: Compare the energy and N from the partition function in "pvt.dat" file output from entropy or entropy2 with that used in simulations

"""
import os, sys, argparse
import subprocess
import numpy as np
import math as ma
from scipy import interpolate
from scipy.optimize import fsolve
import matplotlib as mpl
import matplotlib.pyplot as plt

def runEntropy(T, mu, N):
    feed_file = open('./tmp.dat', 'w')
    for i in range( len(N) ):
        feed_file.write('%f %f %f\n' % (T[i], mu[i], N[i]))
        if i == len(N)-1:
            feed_file.write('stop\n')


    proc = subprocess.Popen("entropy.x < tmp.dat", shell=True, stdout=subprocess.PIPE)
    feed_file.close()
    proc.wait()
    subprocess.Popen("ls", shell=True, stdout=subprocess.PIPE)
    os.remove('tmp.dat')

def runEntropy2(T, mu, N):
    feed_file = open('./tmp.dat', 'w')
    for i in range( len(N) ):
        feed_file.write('%f %f %f\n' % (T[i], mu[i], N[i]))
        if i == len(N)-1:
            feed_file.write('stop\n')


    proc = subprocess.Popen("entropy2.x < tmp.dat", shell=True, stdout=subprocess.PIPE)
    feed_file.close()
    proc.wait()
    subprocess.Popen("ls", shell=True, stdout=subprocess.PIPE)
    os.remove('tmp.dat')

def runFSPVT(T, mu1, mu2):
    feed_file = open('inp_pvt.dat', 'w')
    for i in range( len(mu1) ):
        feed_file.write('%f %f %f\n' % (T[i], mu1[i], mu2[i]))

    proc = subprocess.Popen("fspvt.x", shell=True, stdout=subprocess.PIPE)
    feed_file.close()
    proc.wait()
    subprocess.Popen("ls", shell=True, stdout=subprocess.PIPE)

def readRunsFile(filename, version, skiplines=0):
    if version == "entropy":
        return readOneComponentRunsFile(filename, skiplines)
    elif version == "entropy2":
        return readOneComponentRunsFile(filename, skiplines)
    elif version == "fspvt":
        return readTwoComponentRunsFile(filename, skiplines)

def readOneComponentRunsFile(filename, skiplines=0):
    # open the file
    try:
        ifile = open('./' + filename, 'r')

    except IOError:
        raise IOError('cannot find: %s' % filename)

    n_endings = int(ifile.readline().strip().split()[0])
    endings = []
    for i in range(n_endings):
        end = ifile.readline().strip().split()[0]
        endings.append(end)

    runs = []
    for line in ifile:
        data = line.strip().split()
        for dat in data:
            for end in endings:
                runs.append(dat + end)
        
    return runs
        
def readTwoComponentRunsFile(filename, skiplines=0):
    # open the file
    try:
        ifile = open('./' + filename, 'r')

    except IOError:
        raise IOError('cannot find: %s' % filename)

    # read the header
    ifile.readline()
    runs = []
    for line in ifile:
        runs.append( line.strip().split()[0] )
        
    return runs
        
class hisFile(object):
    """
    Read and plot histogram(s) as a COO sparse matrix
    """
    def __init__(self, file_roots, temperature):
        # can only take 7
        self.runs = file_roots
        self.temp_list = []
        if (temperature):
            self.temp_calc = temperature
        else:
            self.temp_calc = []
        self.Tconv = 1.0
        self.Muconv = 1.0
        self.Econv = 1.0
        return

    def getmu1(self):
        return self.mu1

    def getmu2(self):
        return self.mu2

    def getT(self):
        return self.temp

    def getEave(self, n_components):
        if n_components == 1:
            ave = 0
            for i in range(len(self.histogram[0,:])):
                ave += self.E[i] * sum(self.histogram[:,i]) 
            return ave / float(self.his_sum)
        elif n_components == 2:
            return np.mean( self.E )

    def getN1ave(self, n_components):
        if n_components == 1:
            ave = 0
            for iN in range(len(self.histogram[:,0])):
                ave += iN * sum(self.histogram[iN,:]) / float(self.his_sum)
            return ave
        elif n_components == 2:
            return np.mean( self.N1 )

    def getN2ave(self, n_components):
        if n_components == 2:
            return np.mean( self.N2 )
        else:
            print 'Cannot getN2ave, if there is only 1 component'
            return None

    def getNmaxEmin(self, filename):
        # open the file
        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            raise IOError('cannot find: %s' % filename)

        header_text = ifile.readline()
        header = ifile.readline()
        data = header.strip().split()
        if (len(data) < 6):
            header2 = ifile.readline()

        N_line = True
        N_max = 0 
        E_min = 10000
        E_max = -5000 
        for line in ifile:
            data = line.strip().split()
            if (N_line):
                N = float(data[0])
                N_max = max(N, N_max)

                n_E_bin = int(data[1])

                E = float(data[2])
                E_min = min(E, E_min)
                E_max = max(E + (self.his_width * (n_E_bin-1)), E_max)

                n_E_bin_read = 0
                N_line = False
            else:
                n_E_bin_read += len(data)
                if (n_E_bin_read == n_E_bin):
                    N_line = True
                    
        ifile.close()

        return N_max, E_min, E_max

    def getN12maxEmin(self, filename):
        # open the file
        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            raise IOError('cannot find: %s' % filename)

        header_text = ifile.readline()
        header = ifile.readline()

        N_line = True
        N1_max = 0 
        N2_max = 0 
        E_min = 10000
        E_max = -5000 
        for line in ifile:
            data = line.strip().split()
            N1 = float(data[0])
            N1_max = max(N1, N1_max)

            N2 = float(data[1])
            N2_max = max(N2, N2_max)

            E = float(data[2])
            E_max = max(E, E_max)
            E_min = min(E, E_min)
                    
        ifile.close()

        return N1_max, N2_max, E_min, E_max

    def setJ(self, J):
        self.J = J

    def setBox(self, x, y, z):
        self.box = [x, y, z]

    def setTemp(self, temperature):
        self.temp = temperature

    def setMu(self, mu1):
        self.mu1 = mu1[0]

    def setWidth(self, width):
        self.width = width

    def setHistogram(self, histogram):
        self.his = histogram

    def setNbins(self, n_bins):
        self.n_bins  = n_bins

    def setNeMin(self, energy_start):
        self.e_start = energy_start 

    def setNe(self, n_e_bins):
        self.n_e = n_e_bins

    def setTconv(self, unit):
        if unit == 'K':
            self.Tconv = 1.0 / 0.8314472

    def setEconv(self, unit):
        if unit == 'kJ/mol':
            self.Econv = 0.01
        elif unit == 'K':
            self.Econv = 1.2027221933

    def setMuconv(self, unit):
        if unit == 'kJ/mol':
            self.Muconv = 0.01
        elif unit == 'K':
            self.Muconv = 1.2027221933

    def write(self):
        """ There are 2 types of data lines, the 1st has the:
            N (number of particles)
            n_E_bin (number of energy bins)
            E_min (the min value of energy for non-zero observations)
            
            The next line(s) have the tally of occurances in each bin starting from
        """
        for file_root in self.runs:
            # open the file
            filename = 'his' + file_root + 'a.dat'
            ofile = open('./' + filename, 'w')
    
            # write header information
            ofile.write('T       mu          width     x- y- zdim  \n')
            ofile.write("%12.6f %14.6f %22.12f %10f %10f %10f\n" % ( self.temp, self.mu1, self.width, self.box[0], self.box[1], self.box[2] ))
    
            for i_bin in range( len(self.his[:,0]) ):
                ofile.write("%7i %7i %7.5f\n" % (self.n_bins[i_bin], self.n_e[i_bin], self.e_start[i_bin]))
                l_bin = 1
                writing = False
                for e_bin in range(len(self.his[i_bin,:])):
                    if (self.his[i_bin,e_bin] != 0):
                        writing = True
                    if (writing):
                        ofile.write("%7.f" % self.his[i_bin,e_bin])
                        l_bin += 1
                        if (l_bin > self.n_e[i_bin]):
                            break
                        elif ( l_bin % 9 == 0):
                            ofile.write("\n")
    

    
                ofile.write("\n")
    
            ofile.close()

    def readOneComponent(self, file_root):
        """ There are 2 types of data lines, the 1st has the:
            N (number of particles)
            n_E_bin (number of energy bins)
            E_min (the min value of energy for non-zero observations)
            
            The next line(s) have the tally of occurances in each bin starting from
        """
        # open the file
        filename = 'his' + file_root + '.dat'

        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            raise IOError('cannot find: %s' % filename)
            return 1
        
        # read header information
        header_text = ifile.readline()
        header = ifile.readline()
        data = header.strip().split()
        self.temp = float(data[0])
        if (self.temp not in self.temp_list):
            self.temp_list.append(self.temp)
        self.mu1 = float(data[1])
        self.his_width = float(data[2])

        # if the header was cutoff to a new line
        if (len(data) < 6):
            header2 = ifile.readline()
            data2 = header2.strip().split() 
            for dat in data2:
                data.append(dat)
            
        self.L = []
        for dim in range(3):
            self.L.append(float(data[dim + 3]) )

        self.N_max, self.E_min, self.E_max = self.getNmaxEmin(filename)
        E = np.arange(self.E_min, self.E_max + self.his_width, self.his_width)
        N = np.arange(0, self.N_max+1, dtype=np.int)
        #self.phi = N * self.J / float(self.L[0] * self.L[1] * self.L[2])

        E_int = []
        N_int = []
        his = np.zeros( (len(N), len(E)) )

        i_E_bin = 0
        N_line = True
        # read in the data
        for line in ifile:
            data = line.strip().split()
            # if it is a line with the N and starting E line
            if (N_line):
                i_N = int(data[0])
                E_bin_min = float(data[2])
                n_E_bin = int(data[1])
                # if it is a line with the N and starting E line
                for i_E_bin in range(n_E_bin):
                    i_E = E_bin_min + (i_E_bin * self.his_width)
                    #i_E = E_min + ((n_E_bin - i_E_bin - 1) * self.his_width)
                    #i_E_int = np.where(E == i_E)[0][0]
                    i_E_int = i_E_bin + int( np.floor((E_bin_min - self.E_min)/self.his_width+.1) )
                    E_int.append(i_E_int)
                    N_int.append(i_N)
                i_E_bin = 0 
                n_E_bin_read = 0
                N_line = False

            else:
                for counts in data:
                    i_E = E_bin_min + (i_E_bin * self.his_width)
                    i_E_int = i_E_bin + int( np.floor((E_bin_min - self.E_min)/self.his_width+.1) )
                    #i_E_int = np.where(E == i_E)[0][0]
                    his[i_N, i_E_int] = float(counts)
                    i_E_bin += 1

                n_E_bin_read += len(data)
                if (n_E_bin_read == n_E_bin):
                    i_E_bin = 0
                    N_line = True

        self.histogram = his[:i_N,:]
        self.his_sum = sum(sum(self.histogram))
        N_int = np.array(N_int)
        E_int = np.array(E_int)
        #self.histogram = coo_matrix( (his, (N_int, E_int)), shape=(max(N_int)+1, max(E_int)+1) ).toarray()
        self.N = N
        self.N_int = N_int
        self.E = E
        self.E_int = E_int
        if min(N) != 0:
            print 'min(N) != 0; you should probably get a histogram that does ;)'
        return 0
                    
    def readTwoComponent(self, file_root):
        """ two-component histogram with a line for each instance of 
            N1 (number of particles type 1)
            N2 (number of particles type 2)
            Energy
        """
        # open the file
        filename = 'his' + file_root + '.dat'

        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            raise IOError('cannot find: %s' % filename)
            return 1
        
        header_text = ifile.readline()
        header = ifile.readline()
        data = header.strip().split()
        self.temp = float(data[0])
        if (self.temp not in self.temp_list):
            self.temp_list.append(self.temp)
        self.mu1 = float(data[1])
        self.mu2 = float(data[2])
        self.L = []
        for dim in range(3):
            self.L.append(float(data[dim + 3]) )

        N1 = []
        N2 = []
        E = []
        #N1_max = 0 
        #N2_max = 0 
        #E_min = 10000
        #E_max = -5000 
        for line in ifile:
            data = line.strip().split()
            N1.append( int(data[0]) )
            #N1_max = max(N1, N1_max)

            N2.append( int(data[1]) )
            #N2_max = max(N2, N2_max)

            E.append( float(data[2]) )
            #E_max = max(E, E_max)
            #E_min = min(E, E_min)

        self.N1 = np.array(N1)
        self.N2 = np.array(N2)
        self.E = np.array(E)

        return 0

    def plot(self, show=False, legend=False, plot_phi=False):
        legend_temp = False
        colors = ['black', 'red', 'green', 'blue', 'magenta', 'cyan', 'yellow', 'gold',
                  'lawngreen', 'mediumaquamarine', 'ivory', 'fuchsia']
        labels = [0] * len(self.runs)
        plt.figure()
        plt.xlabel("E", fontsize=15)
        if (plot_phi):
            plt.ylabel("$\phi$", fontsize=18)
        else:
            plt.ylabel("N", fontsize=18)

        if (legend):
            if (len(self.runs) > 6):
                print("Legend cannot show more than 7 entries, turning off")
                legend = False
                if (len(self.temp_list) <= 13):
                    legend_temp = True
                    temp_labels = [0] * len(self.temp_list)

        levels = [100]
        plot_color = 'k'
        plotted_temps = []
        for i in range( len(self.runs) ):
            read_err = self.readOneComponent(self.runs[i])
            if (read_err == 0):
                if (plot_phi):
                    X, Y = np.meshgrid(self.E, self.phi)
                else:
                    X, Y = np.meshgrid(self.E, self.N)

                if (legend):
                    labels[i] = '%s: T=%s, $\mu$=%s' % (self.runs[i], self.temp, self.mu1)
                    plot_color = colors[i]

                elif (legend_temp):
                    if (self.temp not in plotted_temps):
                        plot_color = colors[len(plotted_temps)]

                    else:
                        plot_color = colors[plotted_temps.index(self.temp)]

                CS = plt.contour(X, Y, self.histogram, levels, colors=plot_color)
                
                if (legend):
                    CS.collections[0].set_label(labels[i])

                if (legend_temp):
                    if (self.temp not in plotted_temps):
                        plotted_temps.append(self.temp)
                        CS.collections[0].set_label('T=%s' % self.temp)

                #plt.clabel(CS, inline=1, fontsize=10)

        if (legend or legend_temp):
            plt.legend()

        if (show):
            plt.show()

        return

    def generateMu(self, mu_max, T):
        """
        Generate a list of mus
        Procedure:
            assume mu[i] = mu_step * i**exponent
                (increasing exponent will get larger chem spacing 
                as you go to lower chemical potentials)
        """
        mu_step_cut = [0.0002,5]
        temp = []
        mu = []
        N = []
        mu.append(mu_max)
        temp.append(T)
        N.append(5)
        exponent = 1.01
        for i_mu in range(1, self.mu_len):
            temp.append(T)
            N.append(5)
            mu_step_temp = i_mu**exponent * self.mu_step
            if (mu_step_temp <= mu_step_cut[0]):
                mu_step_temp = 0.0005 * i_mu
            elif (mu_step_temp > mu_step_cut[1]):
                mu_step_temp = 0.05 * i_mu

            mu_test = mu[i_mu-1] - mu_step_temp
            if (mu_test > self.mu_low):
                mu.append(mu_test)
            else:
                mu.append(mu_test)

        return temp[::-1], mu[::-1], N[::-1]

    def getMuMaxT(self):
        """
        Get the runs' temperature list and number of mu
        """
        T = []
        mu_max = []
        for i in range( len(self.runs) ):
            read_err = self.readOneComponent(self.runs[i])
            T_temp = self.getT()
            mu_temp = self.getmu1()
                
            if (T_temp not in T):
                if (self.temp_calc):
                    for i_T in self.temp_calc:
                        if (abs(T_temp - i_T) / float(i_T) < 0.05):
                            T.append(T_temp)
                            mu_max.append(self.mu_low)
                else:
                    T.append(T_temp)
                    mu_max.append(self.mu_low)

            for i_T in range( len(T) ):
                if (T_temp == T[i_T]):
                    if (mu_temp > mu_max[i_T]):
                        mu_max[i_T] = mu_temp
                    break
    
        for T_calc in self.temp_calc:
            not_in_array = True
            for i_T in T:
                if (abs(i_T - T_calc) / float(T_calc) < 0.05):
                    not_in_array = False
                    break

            if not_in_array:
                print "The temperature in the input (", str(T_calc), ") is not used in any simulation run given to entropy. It may not be as accurate."
                T.append(T_calc)
                mu_max.append(self.mu_low)

        return mu_max, T

    def getMuMinMaxToneComponent(self):
        """
        Get the runs' temperature list and number of mu
        """
        T = []
        mu_min = []
        mu_max = []
        mu_low = -100000
        mu_high = 100000
        for i in range( len(self.runs) ):
            read_err = self.readOneComponent(self.runs[i])
            T_temp = self.getT()
            mu_temp = self.getmu1()
                
            if (T_temp not in T):
                if (self.temp_calc):
                    for i_T in self.temp_calc:
                        if (abs(T_temp - i_T) / float(i_T) < 0.05):
                            T.append(T_temp)
                            mu_min.append(mu_high)
                            mu_max.append(mu_low)
                else:
                    T.append(T_temp)
                    mu_min.append(mu_high)
                    mu_max.append(mu_low)

            for i_T in range( len(T) ):
                if (T_temp == T[i_T]):
                    if (mu_temp < mu_min[i_T]):
                        mu_min[i_T] = mu_temp
                    if (mu_temp > mu_max[i_T]):
                        mu_max[i_T] = mu_temp
    
        for T_calc in self.temp_calc:
            not_in_array = True
            for i_T in T:
                if (abs(i_T - T_calc) / float(T_calc) < 0.05):
                    not_in_array = False
                    break

            if not_in_array:
                print "The temperature in the input (", str(T_calc), ") is not used in any simulation run given to entropy. It may not be as accurate."
                T.append(T_calc)
                mu_min.append(mu_high)
                mu_max.append(mu_low)

        return mu_min, mu_max, T

    def getMuMinMaxTtwoComponent(self):
        """
        Get the runs' temperature list and number of mu
        """
        T = []
        mu1_min = []
        mu1_max = []
        mu2_min = []
        mu2_max = []
        mu_low = [-100000, -100000]
        mu_high = [100000, 100000]
        for i in range( len(self.runs) ):
            read_err = self.readTwoComponent(self.runs[i])
            T_temp = self.getT()
            mu1_temp = self.getmu1()
            mu2_temp = self.getmu2()
                
            if (T_temp not in T):
                if (self.temp_calc):
                    for i_T in self.temp_calc:
                        if (abs(T_temp - i_T) / float(i_T) < 0.05):
                            T.append(T_temp)
                            mu1_min.append(mu_high[0])
                            mu1_max.append(mu_low[0])
                            mu2_min.append(mu_high[1])
                            mu2_max.append(mu_low[1])
                else:
                    T.append(T_temp)
                    mu1_min.append(mu_high[0])
                    mu1_max.append(mu_low[0])
                    mu2_min.append(mu_high[1])
                    mu2_max.append(mu_low[1])

            for i_T in range( len(T) ):
                if (T_temp == T[i_T]):
                    if (mu1_temp < mu1_min[i_T]):
                        mu1_min[i_T] = mu1_temp
                    if (mu1_temp > mu1_max[i_T]):
                        mu1_max[i_T] = mu1_temp
                    if (mu2_temp < mu2_min[i_T]):
                        mu2_min[i_T] = mu2_temp
                    if (mu2_temp > mu2_max[i_T]):
                        mu2_max[i_T] = mu2_temp
    
        for T_calc in self.temp_calc:
            not_in_array = True
            for i_T in T:
                if (abs(i_T - T_calc) / float(T_calc) < 0.05):
                    not_in_array = False
                    break

            if not_in_array:
                print "The temperature in the input (", str(T_calc), \
                       ") is not used in any simulation run given to entropy. ", \
                       "It may not be as accurate.", \
                       "Assuming this for prediction."
                T_near = T_calc*1000
                for i_T in range(len(T)):
                    if (abs(T[i_T] - T_calc) < abs(T_near - T_calc)):
                        T_near = T[i_T]
                        mu1_min.append(mu1_min[i_T])
                        mu1_max.append(mu1_max[i_T])
                        mu2_min.append(mu2_min[i_T])
                        mu2_max.append(mu2_max[i_T])
                
                T.append(T_calc)
                # since there are no chemical potential data at this temperature
                # find the the nearest tempearture and use those
                if T_near == T_calc*1000:
                    mu1_min.append(mu_high[0])
                    mu1_max.append(mu_low[0])
                    mu2_min.append(mu_high[1])
                    mu2_max.append(mu_low[1])

                print "Assuming this for prediction.", \
                      " using chemical potential initial guess from T = ", str(T_near)

        return mu1_min, mu1_max, mu2_min, mu2_max, T

    def generate(self, entropy_version="entropy", save_file=False, mu_len=30, N_lo=[0.01,0.1],N_hi=[60, 130], max_iter=200, n1_eq_n2=False):
        if (entropy_version == "entropy" or entropy_version == "entropy2"):
            self.generateCMCpvtOne(entropy_version, save_file, mu_len, N_lo, N_hi, max_iter)

        elif (entropy_version == "fspvt"):
            #self.generateCMCpvtTwo(entropy_version, save_file, mu_len, mu_len, [0, 280], [320, 330], [30,35], [], True, max_iter)
            self.generateCMCpvtTwo(entropy_version, save_file, mu_len, mu_len, N_lo, N_hi, N_lo, N_hi, n1_eq_n2, 30)

    def generateCMCpvtTwo(self, entropy_version="fspvt", save_file=False, mu1_len=3, mu2_len=3, 
                                N1_lo=[0.01,0.1], N1_hi=[60, 130], N2_lo=[0.01,0.1], N2_hi=[60, 130],
                                N1_eq_N2=False, max_iter=50):
        """
        Generate a P-N1-N2 that covers ideal gas and high pressure
        for micellization
        self-optimizing process:
            for each temperature
            generate N1(mu1,mu2) & N2(mu1,mu2) surfaces 
            if N1_eq_N2==True, then identify the line which intersects the surfaces
        """
        import shutil 
        griddata = False

        if (entropy_version == "entropy" or entropy_version == "entropy2"):
            self.n_components = 1    
        elif (entropy_version == "fspvt"):
            self.n_components = 2

        # The range that we want the N values to be in
        # N_lim = [<N min/max>, <upper/lower range bound>,<species>]
        N_lim = np.array([ [ [N1_lo[0], N2_lo[0]],
                             [N1_lo[1], N2_lo[1]] ],
                           [ [N1_hi[0], N2_hi[0]],
                             [N1_hi[1], N2_hi[1]] ],
                           [ [N1_lo[0], N2_lo[0]],
                             [N1_lo[1], N2_lo[1]] ],
                           [ [N1_hi[0], N2_hi[0]],
                             [N1_hi[1], N2_hi[1]] ] ])
        mu1_min, mu1_max, mu2_min, mu2_max, T = self.getMuMinMaxTtwoComponent()
        print T

        # for modifying the range of chemical potentials, so that we get the
        # range of N values we want
        o_mu_min = [0, 0]
        o_mu_max = [0, 0]
        mu_min = [0, 0]
        mu_max = [0, 0]
        N_min = [1000, 1000]
        N_max = [0, 0]
        o_N_min = [1000, 1000]
        o_N_max = [0, 0]
        for i_T in range(len(T)):
            mu_min[0] = mu1_min[i_T]
            mu_max[0] = mu1_max[i_T]
            mu_min[1] = mu2_min[i_T]
            mu_max[1] = mu2_max[i_T]
            if griddata:
                # generate a grid of values calculated from histogram reweighting
                print '# grid_iteration N1_low N2_low N1_hi N2_hi'
                for mu_iter in range(max_iter):
                    
                    a1 = -np.logspace(ma.log(-mu_min[0], 10.), 
                                     ma.log(-mu_max[0], 10.), mu1_len, base=10.)
                    a2 = -np.logspace(ma.log(-mu_min[1], 10.), 
                                     ma.log(-mu_max[1], 10.), mu2_len, base=10.)

                    matrix_mu1, matrix_mu2 = np.meshgrid(a1, a2)

                    array_mu1 = matrix_mu1.ravel()
                    array_mu2 = matrix_mu2.ravel()
                    array_T = np.ones(len(array_mu1)) * T[i_T]
                    part = partition()
                    part.setNcomponents(2)
                    runFSPVT(array_T, array_mu1, array_mu2)
                    part.readPVTsimple(self.n_components)
                    N_min[0] = np.amin(part.N1)
                    N_max[0] = np.amax(part.N1)
                    N_min[1] = np.amin(part.N2)
                    N_max[1] = np.amax(part.N2)

                    # Refine the chemical potentials to generate the surfaces
                    for i in range(2):
                        tm_min, tm_max =self.modifyMuGuess(o_N_min[i], N_min[i],
                                                           o_N_max[i], N_max[i],
                                                           o_mu_min[i], mu_min[i],
                                                           o_mu_max[i], mu_max[i],
                                                           N_lim[0,0,i], N_lim[1,0,i],
                                                           N_lim[0,1,i], N_lim[1,1,i])
                        o_mu_max[i] = mu_max[i]
                        o_mu_min[i] = mu_min[i]
                        mu_max[i] = max(tm_max, mu_min[i]*0.99)
                        mu_min[i] = min(tm_min, mu_max[i]*1.01)
                        o_N_min[i] = N_min[i]
                        o_N_max[i] = N_max[i]
                        
                    converge_fail, in_N_range = self.checkFailureTwo(o_N_min, N_min,
                                                                 o_N_max, N_max,
                                                                 N_lim[0,0,:], N_lim[1,0,:],
                                                                 N_lim[0,1,:], N_lim[1,1,:])

                    print mu_iter, N_min[0], N_min[1], N_max[0], N_max[1]

                    if in_N_range: break
                    if converge_fail: break

            else:
                # first loop over the corners of the matrix, i.e.
                # set the max and min mu1 and mu2, 
                # results in 8 chemical potentials
                # mu_extrema[<mu1_minmax>, <mu2_minmax>, <species_mu>]
                mu_extrema = np.array( [[ [mu_min[0], mu_min[1]],
                                         [mu_min[0], mu_max[1]] ],
                                        [ [mu_max[0], mu_min[1]],
                                         [mu_max[0], mu_max[1]] ]] )
                # Set the mu1 extrema to be in range at the current mu2 extrema
                o_mu_min = mu_extrema[0,0]
                o_mu_max = mu_extrema[0,1]
                mu_min = mu_extrema[0,0]
                mu_max = mu_extrema[0,1]
                for a1i in range(2):
                  if a1i == 1:
                    # use the mu1 that you used previously to get N1 hi or low
                    mu_extrema[1, 0, 1] = mu_extrema[0, 0, 1]
                    mu_extrema[1, 1, 1] = mu_extrema[0, 1, 1]

                  for a2i in range(2):
                    if a2i == 0:
                      # use the mu1 and mu2 from the preset
                      o_mu = [mu_extrema[a1i,a2i, 0], mu_extrema[a1i,a2i, 1]]
                      mu = [mu_extrema[a1i,a2i, 0], mu_extrema[a1i,a2i, 1]]
                    elif a2i == 1:
                      # use the mu1 that you used previously to get N1 hi or low
                      o_mu = [mu_extrema[a1i, 0, 0], mu_extrema[a1i, a2i, 1]]
                      mu = [mu_extrema[a1i, 0, 0], mu_extrema[a1i, a2i, 1]]

                    o_N = [0, 0]
                    N = [0, 0]
                    print '# iter 1_hilo 2_hilo N1 N2 mu1 mu2'
                    for mu1_iter in range(max_iter):
                      part = partition()
                      part.setNcomponents(2)
                      runFSPVT([T[i_T]], [mu[0]], [mu[1]])
                      part.readPVTsimple(self.n_components)
                      N[0] = part.N1[0]
                      N[1] = part.N2[0]

                      print '  ', mu1_iter, a1i, a2i, N[0], N[1], mu[0], mu[1]

                      ii = 0
                      jj = 0
                      for i in range(2):
                        if i == 0: aii = a1i
                        if i == 1: aii = a2i
                        converge_fail, in_N_range = self.checkFailureOneOne(
                                                       o_N[i], N[i], 
                                                       N_lim[aii,0,i], 
                                                       N_lim[aii,1,i])
                        if in_N_range: ii += 1
                        if converge_fail: jj += 1
  
                      if ii == 2: break
                      if jj == 2: break
  
                      for i in range(2):
                        # Refine the chemical potentials to generate the surfaces
                        if i == 0: aii = a1i
                        if i == 1: aii = a2i
                        t_mu = self.modifyMuGuessOne(N[i], mu[i],
                                                     N_lim[aii,0,i], 
                                                     N_lim[aii,1,i])
  
                        o_mu[i] = mu[i]
                        mu[i] = t_mu
                        o_N[i] = N[i]
  
                    mu_extrema[a1i,a2i,0] = mu[0]
                    mu_extrema[a1i,a2i,1] = mu[1]

                if True:
                    a1 = -np.logspace(ma.log(-mu_extrema[0,0,0], 10.), 
                                     ma.log(-mu_extrema[1,1,0], 10.), mu1_len, base=10.)
                    a2 = -np.logspace(ma.log(-mu_extrema[0,0,1], 10.), 
                                     ma.log(-mu_extrema[1,1,1], 10.), mu2_len, base=10.)

                    matrix_mu1, matrix_mu2 = np.meshgrid(a1, a2)
    
                    array_mu1 = matrix_mu1.ravel()
                    array_mu2 = matrix_mu2.ravel()
                    array_T = np.ones(len(array_mu1)) * T[i_T]
                    part = partition()
                    part.setNcomponents(2)
                    runFSPVT(array_T, array_mu1, array_mu2)
                    part.readPVTsimple(self.n_components)
                shutil.copyfile('./pvt.dat','./pvt_'+str(T[i_T])+'.dat')

            # Find the curve where these two intersect
            if N1_eq_N2:
                # interpolate to get a full function
                style = 'CloughTocher2DInterpolator'
                N1_interpolate = self.getInterpolation(array_mu1, array_mu2, part.N1, style)
                N2_interpolate = self.getInterpolation(array_mu1, array_mu2, part.N2, style)

                n_test = 15
                mu1_test_min = o_mu_min[0]
                mu1_test_max = o_mu_max[0]
                o_max_N_intersection = -100
                o_min_N_intersection = 100000
                print '# interpolate_iter min_N_intersection max_N_intersection'
                for i_iter in range(max_iter):
                    mu1_test = np.linspace(mu1_test_min, mu1_test_max, n_test)
                    N_intersection = np.zeros(n_test)
                    mu2_intersection = np.zeros(n_test)
                    for im1 in range(n_test):
                        if im1 < 5:
                            mu2_guess = mu1_test[im1]
                        mu2_intersection[im1] = self.findIntersection( N1_interpolate, 
                                                                       N2_interpolate, 
                                                                       mu1_test[im1], 
                                                                       mu2_guess )
                        N_intersection[im1] = N1_interpolate( mu1_test[im1], 
                                                              mu2_intersection[im1] )
                        #print mu1_test[im1], mu2_guess , N_intersection[im1]
                        if ma.isnan( N_intersection[im1] ):
                            break
                        mu2_guess = mu2_intersection[im1]

                    if ma.isnan( N_intersection[im1] ):
                        break
        
                    # check if the N is in the range
                    min_N_intersection = min(N_intersection)
                    max_N_intersection = max(N_intersection)
                    print i_iter, min_N_intersection, max_N_intersection

                    if (( N_low_max[0] > min_N_intersection > N_low_min[0]) and
                        ( N_hi_max[0] > max_N_intersection > N_hi_min[0])):
                        break
                    elif ((abs(o_min_N_intersection - min_N_intersection) / 
                                        min_N_intersection < 0.001) and
                          (abs(o_max_N_intersection - max_N_intersection) / 
                                        max_N_intersection < 0.001)):
                        break
        
                    # Otherwise, let's refine the guess
                    if min_N_intersection > N_low_max[0]:
                        mu1_test_min = mu1_test_min - abs(mu1_test_min*0.1)
                    elif min_N_intersection < N_low_min[0]:
                        mu1_test_min = mu1_test_min + abs(mu1_test_min*0.1)
        
                    if min_N_intersection > N_hi_max[0]:
                        mu1_test_max = mu1_test_max - abs(mu1_test_max*0.1)
                    elif min_N_intersection < N_hi_min[0]:
                        mu1_test_max = mu1_test_max + abs(mu1_test_max*0.1)
    
                    o_min_N_intersection = min_N_intersection
                    o_max_N_intersection = max_N_intersection
        
                shutil.copyfile('./pvt.dat','./pvt_grid.dat')
                # tweak prediction
                on1 = -1
                on2 = -1
                for im1 in range(n_test):
                    for j_iter in range(max_iter):
                        part = partition()
                        part.setNcomponents(2)
                        runFSPVT([T[i_T]], [mu1_test[im1]], [mu2_intersection[im1]])
                        part.readPVTsimple(self.n_components)

                        n1 = part.N1[0]
                        n2 = part.N2[0]

                        dn = (n1 - n2) / float(n1)
                        if abs(dn) > 0.01:
                            mu2_intersection[im1] = (mu2_intersection[im1] 
                                                   + dn * 0.1 * abs(mu2_intersection[im1]))
                        else:
                            break

                        if (abs((on1 - n1) / on1) < 0.0001 and
                            abs((on2 - n2) / on2) < 0.0001):
                            break
                        on1 = n1
                        on2 = n2

                    print mu1_test[im1], mu2_intersection[im1], n1, n2, dn, abs(on1 - n1) / on1, abs(on2 - n2) / on2
     
                array_T = np.ones(n_test) * T[i_T]
                part = partition()
                part.setNcomponents(2)
                runFSPVT(array_T, mu1_test, mu2_intersection)
                # now starting from intersection of planes, find the actual intersection
                shutil.copyfile('./pvt.dat','./pvt_n1_eq_n2.dat')

    def findIntersection(self, fun1, fun2, x, y0):
        return fsolve(lambda y : fun1(x, y) - fun2(x, y), y0, xtol=1.e-02, maxfev=10000)
            
    def getInterpolation(self, x, y, z, style='CloughTocher2DInterpolator'):
        """
        x, y and z are 1-dimensional data, that are preferably on a grid
        To use RectBivariateSpline the data must be on a grid
        """
        if style == 'RectBivariateSpline': #2
            # reshape x list to a grid 
            xx, counts = np.unique(x, return_counts=True)
            yy, counts = np.unique(y, return_counts=True)

            zz = np.reshape(z, (len(xx), len(yy)))
            return interpolate.RectBivariateSpline(xx, yy, zz).ev
        elif style == 'interp2d': #4
            # reshape x list to a grid 
            unique, counts = np.unique(x, return_counts=True)
            n_x = counts[0]
            unique, counts = np.unique(y, return_counts=True)
            n_y = counts[0]
            
            xx = np.reshape(x, (n_x, n_y))
            yy = np.reshape(y, (n_x, n_y))
            zz = np.reshape(z, (n_x, n_y))
            return interpolate.interp2d(xx, yy, zz, kind='cubic')
        elif style == 'LinearNDInterpolator': #3
            crd = list(zip(x, y))
            return interpolate.LinearNDInterpolator(crd, z)
        elif style == 'NearestNDInterpolator': #5
            crd = list(zip(x, y))
            return interpolate.NearestNDInterpolator(crd, z)
        elif style == 'CloughTocher2DInterpolator': #1
            crd = list(zip(x, y))
            return interpolate.CloughTocher2DInterpolator(crd, z)
        else:
            print 'scipy.interpolate method not available'
            return None

    def checkFailureTwo(self, o_N_min, N_min, o_N_max, N_max,
                              N_low_min, N_low_max, N_hi_min, N_hi_max):
        converge_fail = False
        in_N_range = False
        for i in range(2):
            if (abs(o_N_min[i] - N_min[i])/N_min[i] < 0.0001):
                if N_min[i] < N_low_min[i] or N_min[i] > N_low_max[i]:
                    converge_fail = True
            if (abs(o_N_max[i] - N_max[i])/N_max[i] < 0.0001):
                if N_max[i] < N_hi_min[i] or N_max[i] > N_hi_max[i]:
                    converge_fail = True

            if( (N_low_max[i] > N_min[i] > N_low_min[i]) and 
                  (N_hi_max[i] > N_max[i] > N_hi_min[i]) ):
                if i == 0:
                    in_N_range = True
            else:
                in_N_range = False

        return converge_fail, in_N_range

    def checkFailureOne(self, o_N_min, N_min, o_N_max, N_max,
                              N_low_min, N_low_max, N_hi_min, N_hi_max):

        converge_fail = False
        in_N_range = False
        if (abs(o_N_min - N_min)/N_min < 0.0001):
            if N_min < N_low_min or N_min > N_low_max:
                converge_fail = True
        if (abs(o_N_max - N_max)/N_max < 0.0001):
            if N_max < N_hi_min or N_max > N_hi_max:
                converge_fail = True

        if( (N_low_max > N_min > N_low_min) and 
              (N_hi_max > N_max > N_hi_min) ):
            in_N_range = True

        return converge_fail, in_N_range

    def checkFailureOneOne(self, o_N, N, N_min, N_max):

        converge_fail = False
        in_N_range = False
        if (abs(o_N - N)/N < 0.0001):
            if N < N_min or N > N_max:
                converge_fail = True

        if N_max > N > N_min:
            in_N_range = True

        return converge_fail, in_N_range

    def modifyMuGuess(self, o_N_min, N_min, o_N_max, N_max,
                            o_mu_min, mu_min, o_mu_max, mu_max,
                            N_low_min, N_low_max, N_hi_min, N_hi_max):

        mu_max_temp = mu_max
        mu_min_temp = mu_min

        if N_min > N_low_max:
            if True:
                mu_diff = mu_max - mu_min
                mu_min_temp = ( mu_min
                                - abs(mu_diff * 0.7 *
                                 (N_min-N_low_min) / N_min) )
            else:
                mu_min_temp = (( (o_mu_min - mu_min)
                                  / (o_N_min - N_min) 
                                  * (N_low_min - o_N_min) ) 
                                  + o_mu_min)
        elif N_min < N_low_min:
            if True:
                mu_diff = mu_max - mu_min
                mu_min_temp = ( mu_min
                                + abs(mu_diff * 0.3 *
                                  (N_low_max-N_min) / N_low_max))
            else:
                mu_min_temp = (( (o_mu_min - mu_min)
                                  / (o_N_min - N_min) 
                                  * (N_low_max - o_N_min) ) 
                                  + o_mu_min)
   
        if N_max > N_hi_max:
            if True:
                mu_diff = mu_max - mu_min
                mu_max_temp = ( mu_max
                                - abs(mu_diff * 0.3 *
                                 (N_max-N_hi_min) / N_max) )
            else:
                mu_max_temp = (( (o_mu_max - mu_max)
                                  / (o_N_max - N_max) 
                                  * (N_hi_min - o_N_max) ) 
                                  + o_mu_max)
        elif N_max < N_hi_min:
            if True:
                mu_diff = mu_max - mu_min
                mu_max_temp = ( mu_max
                                + abs(mu_diff * 0.7 * 
                                  (N_hi_max-N_max) / N_hi_min) )
            else:
                mu_max_temp = (( (o_mu_max - mu_max)
                                  / (o_N_max - N_max) 
                                  * (N_hi_max - o_N_max) ) 
                                  + o_mu_max)

        return mu_min_temp, mu_max_temp
    
    def modifyMuGuessOne(self, N, mu, N_min, N_max):
        if N > N_max:
            return ( mu - abs(mu * 0.1 * (N_min-N) / N_min) )

        elif N < N_min:
            return ( mu + abs(mu * 0.1 * (N_max-N) / N_max) )
   
        else:
            return mu
    
    def anotherone(self):
        temp_tot = []
        mu_tot = []
        N_tot = []
        print 'generating for T =', T[i_T]
        # Find the max chemical potential
        mu1_m = [mu1_max[i_T]]
        mu1_range = mu1_max[i_T] - mu1_min[i_T]
        mu2_m = [mu2_max[i_T]]
        mu2_range = mu2_max[i_T] - mu2_min[i_T]
        temp = []
        mu1 = []
        N1 = []
        mu2 = []
        N2 = []
        for i in range(max_iter):
            done = True
            part = partition()

            if (entropy_version == "fspvt"):
                runFSPVT([T[i_T]], mu1_m, mu2_m)
                part.setNcomponents(2)

            part.readPVTsimple(self.n_components)
            if (part.N1[0] < N_hi_min[0]):
                mu1_m[0] = max(mu1_m[0]+(mu1_range*0.01), mu1_m[0]+0.5)
                done = False
            elif (part.N1[0] > N_hi_max[0]):
                mu1_m[0] = min(mu1_m[0]-(mu1_range*0.05), mu1_m[0]-0.01)
                done = False
            else:
                N1_mu_m[0] = part.N1[0]

            if (part.N2[0] < N_hi_min[1]):
                mu2_m[0] = max(mu2_m[0]+(mu2_range*0.01), mu2_m[0]+0.5)
                done = False
            elif (part.N2[0] > N_hi_max[1]):
                mu2_m[0] = min(mu2_m[0]-(mu2_range*0.05), mu2_m[0]-0.01)
                done = False
            else:
                N2_mu_m[0] = part.N2[0]

            if N1_eq_N2:
                if done == True:
                    if abs(part.N2[0] - part.N1[0]) > 1.0:
                        done = False
                        if (part.N2[0] > part.N1[0]):
                            mu2_m[0] = max(mu2_m[0]+(mu2_range*0.01), mu2_m[0]+0.5)
                            mu1_m[0] = min(mu1_m[0]-(mu1_range*0.05), mu1_m[0]-0.01)
                        elif (part.N2[0] < part.N1[0]):
                            mu2_m[0] = min(mu2_m[0]-(mu2_range*0.05), mu2_m[0]-0.01)
                            mu1_m[0] = max(mu1_m[0]+(mu1_range*0.01), mu1_m[0]+0.5)

            if done == True: 
                break
            print i, mu1_m[0], part.N1[0], mu2_m, part.N2[0]

        if (i == max_iter):
            print 'change the maximum N range!'
            return
        print mu1_m[0], mu2_m[0]

        temp.append(T[i_T])
        mu1.append(mu1_m[0])
        mu2.append(mu2_m[0])

        # Find the min chemical potential
        mu1_n = [mu1_min[i_T]]
        mu2_n = [mu2_min[i_T]]

        part = partition()
        if (entropy_version == "fspvt"):
            runFSPVT([T[i_T]], mu1_n, mu2_n)
            part.setNcomponents(2)

        part.readPVTsimple(self.n_components)
        pmN1 = part.N1[0]
        pi_mu1 = mu1_n[0]
        if (part.N1[0] < N_low_min[0]):
            mu1_n[0] = max(mu1_n[0]*1.3, mu1_n[0]+0.1)
        elif (part.N1[0] > N_low_max[0]):
            mu1_n[0] = min(mu1_n[0]/3.0, mu1_n[0]-2.0)

        pmN2 = part.N2[0]
        pi_mu2 = mu2_n[0]
        if (part.N2[0] < N_low_min[1]):
            mu2_n[0] = max(mu2_n[0]*1.3, mu2_n[0]+0.1)
        elif (part.N2[0] > N_low_max[1]):
            mu2_n[0] = min(mu2_n[0]/3.0, mu2_n[0]-2.0)

        for i in range(max_iter):
            if (N_low_min < part.N1[0] < N_low_max):
                N_mu_n = part.N1[0]
                break

            part = partition()
            if (entropy_version == "fspvt"):
                runFSPVT([T[i_T]], mu1_n, mu2_n)
                part.setNcomponents(2)

            part.readPVTsimple(self.n_components)

            m = (part.N1[0] - pmN1) / (mu1_n[0] - pi_mu1)
            b = part.N1[0] - (m * mu1_n[0])
            pmN1 = part.N1[0]
            pi_mu1 = mu1_n[0]
            mu1_n[0] = (N_low_min[0] - b) / m 

            m = (part.N2[0] - pmN2) / (mu2_n[0] - pi_mu2)
            b = part.N2[0] - (m * mu2_n[0])
            pmN2 = part.N2[0]
            pi_mu2 = mu2_n[0]
            mu2_n[0] = (N_low_min[1] - b) / m 

        if (i == max_iter):
            print 'change the minimum N range!'
            return

        # Find the chemical potentials in between
        mu1_step = abs(mu1_m[0] - mu1_n[0]) / float(mu1_len)**1.0
        exponent1 = 1.0 + float(1.0 / mu1_len)
        mu1_step_temp = mu1_step
        i_mu1 = mu1_m
        for im in range(mu1_len):
            i_mu1[0] -= mu1_step_temp
            mu1_step_temp = im**exponent1 * mu1_step

            temp.append(T[i_T])
            mu1.append(i_mu1[0])

        mu2_step = abs(mu2_m[0] - mu2_n[0]) / float(mu2_len)**2.0
        exponent2 = 1.0 + float(1.0 / mu2_len)
        mu2_step_temp = mu2_step
        i_mu2 = mu2_m
        for im in range(mu2_len):
            i_mu2[0] -= mu2_step_temp
            mu2_step_temp = im**exponent2 * mu2_step

            temp.append(T[i_T])
            mu2.append(i_mu2[0])

        temp.append(T[i_T])
        mu1.append(mu1_n[0])
        mu2.append(mu2_n[0])

        temp_tot.extend(temp)
        mu1_tot.extend(mu1)
        mu2_tot.extend(mu2)

        temp_tot.reverse()
        mu1_tot.reverse()
        mu2_tot.reverse()

        if (entropy_version == "fspvt"):
            runEntropy(temp_tot, mu1_tot, mu2_tot)
            part.setNcomponents(2)
        
        if (save_file):
            s_file = open( 'TmuN_sVp.dat', 'w')
            for i in range( len(mu1_tot) ):
                s_file.write('%f %f %f\n' % (temp_tot[i], mu1_tot[i], mu2_tot[i]))
                if i == len(mu1_tot)-1:
                    s_file.write('stop\n')

    def generateCMCpvtOne(self, entropy_version="entropy", save_file=False, mu_len=30, N_lo=[0.01,0.1],N_hi=[60, 130], max_iter=200):
        """
        Generate a PVT that covers ideal gas and high pressure
        for micellization
        self-optimizing process:
            for each temperature
            find a chemical potential that lands in the range of upper N
            find a chemical potential that lands in the range of lower N
            starting from the top step down ramping up the step size to scale somewhat with the N(\mu) function
        """
        if (entropy_version == "entropy" or entropy_version == "entropy2"):
            self.n_components = 1    
        elif (entropy_version == "fspvt"):
            self.n_components = 2

        N_low_min = N_lo[0]
        N_low_max = N_lo[1]
        N_hi_min = N_hi[0]
        N_hi_max = N_hi[1]
        mu_min, mu_max, T = self.getMuMinMaxToneComponent()
        temp_tot = []
        mu_tot = []
        N_tot = []
        for i_T in range(len(T)):
            print 'generating for T =', T[i_T]
            # Find the max chemical potential
            mu_m = [mu_max[i_T]]
            mu_range = mu_max[i_T] - mu_min[i_T]
            temp = []
            mu = []
            N = []
            for i in range(max_iter):
                part = partition()

                if (entropy_version == "entropy"):
                    runEntropy([T[i_T]], mu_m, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "entropy2"):
                    runEntropy2([T[i_T]], mu_m, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "fspvt"):
                    part.setNcomponents(2)

                part.readPVTsimple(self.n_components)
                if (part.N1[0] < N_hi_min):
                    mu_m[0] = max(mu_m[0]+(mu_range*0.01), mu_m[0]+0.5)
                elif (part.N1[0] > N_hi_max):
                    mu_m[0] = min(mu_m[0]-(mu_range*0.05), mu_m[0]-0.01)
                else:
                    N_mu_m = part.N1[0]
                    break

            if (i == max_iter):
                print 'change the maximum N range!'
                return

            temp.append(T[i_T])
            N.append(5)
            mu.append(mu_m[0])

            # Find the min chemical potential
            mu_n = [mu_min[i_T]]

            part = partition()
            if (entropy_version == "entropy"):
                runEntropy([T[i_T]], mu_n, [5])
                part.setNcomponents(1)
            elif (entropy_version == "entropy2"):
                runEntropy2([T[i_T]], mu_n, [5])
                part.setNcomponents(1)
            elif (entropy_version == "fspvt"):
                part.setNcomponents(2)

            part.readPVTsimple(self.n_components)
            pmN = part.N1[0]
            pi_mu = mu_n[0]
            if (part.N1[0] < N_low_min):
                mu_n[0] = max(mu_n[0]*1.3, mu_n[0]+0.1)
            elif (part.N1[0] > N_low_max):
                mu_n[0] = min(mu_n[0]/3.0, mu_n[0]-2.0)

            for i in range(max_iter):
                if (N_low_min < part.N1[0] < N_low_max):
                    N_mu_n = part.N1[0]
                    break

                part = partition()
                if (entropy_version == "entropy"):
                    runEntropy([T[i_T]], mu_n, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "entropy2"):
                    runEntropy2([T[i_T]], mu_n, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "fspvt"):
                    part.setNcomponents(2)

                part.readPVTsimple(self.n_components)

                m = (part.N1[0] - pmN) / (mu_n[0] - pi_mu)
                b = part.N1[0] - (m * mu_n[0])
                pmN = part.N1[0]
                pi_mu = mu_n[0]
                mu_n[0] = (N_low_min - b) / m 

            if (i == max_iter):
                print 'change the minimum N range!'
                return

            # Find the chemical potentials in between
            mu_step = abs(mu_m[0] - mu_n[0]) / float(mu_len)**2.0
            exponent = 1.0 + float(1.0 / mu_len)
            mu_step_temp = mu_step
            i_mu = mu_m
            for im in range(mu_len):
                i_mu[0] -= mu_step_temp
                mu_step_temp = im**exponent * mu_step

                temp.append(T[i_T])
                N.append(5)
                mu.append(i_mu[0])

            temp.append(T[i_T])
            N.append(5)
            mu.append(mu_n[0])

            temp_tot.extend(temp)
            mu_tot.extend(mu)
            N_tot.extend(N)

        temp_tot.reverse()
        mu_tot.reverse()
        N_tot.reverse()

        if (entropy_version == "entropy"):
            runEntropy(temp_tot, mu_tot, N_tot)
            part.setNcomponents(1)
        elif (entropy_version == "entropy2"):
            runEntropy2(temp_tot, mu_tot, N_tot)
            part.setNcomponents(1)
        elif (entropy_version == "fspvt"):
            part.setNcomponents(2)

        
        if (save_file):
            s_file = open( 'TmuN_sVp.dat', 'w')
            for i in range( len(N_tot) ):
                s_file.write('%f %f %f\n' % (temp_tot[i], mu_tot[i], N_tot[i]))
                if i == len(N_tot)-1:
                    s_file.write('stop\n')

    def generateCurve(self, entropy_version=1, save_file=False, mu_len=30, N_lo=[0.01,0.1],N_hi=[60, 130], max_iter=100):
        """
        Generate a PVT that covers ideal gas and high pressure
        for micellization
        self-optimizing process:
            for each temperature
            find a chemical potential that lands in the range of upper N
            find a chemical potential that lands in the range of lower N
            starting from the top step down in equal increments of N
        """
        if (entropy_version == "entropy" or entropy_version == "entropy2"):
            self.n_components = 1    
        elif (entropy_version == "fspvt"):
            self.n_components = 2

        N_low_min = N_lo[0]
        N_low_max = N_lo[1]
        N_hi_min = N_hi[0]
        N_hi_max = N_hi[1]
        mu_min, mu_max, T = self.getMuMinMaxToneComponent()
        temp_tot = []
        mu_tot = []
        N_tot = []
        for i_T in range(len(T)):
            print 'generating for T =', T[i_T]
            # Find the max chemical potential
            mu_m = [mu_max[i_T]]
            temp = []
            mu = []
            N = []
            for i in range(max_iter):
                part = partition()
                if (entropy_version == "entropy"):
                    runEntropy([T[i_T]], mu_n, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "entropy2"):
                    runEntropy2([T[i_T]], mu_n, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "fspvt"):
                    part.setNcomponents(2)

                if (part.N1[0] < N_hi_min):
                    mu_m[0] = max(mu_m[0]*1.1, mu_m[0]+0.5)
                elif (part.N1[0] > N_hi_max):
                    mu_m[0] = min(mu_m[0]/1.1, mu_m[0]-0.01)
                else:
                    N_mu_m = part.N1[0]
                    break
            if (i == max_iter):
                print 'change the maximum N range!'
                return

            temp.append(T[i_T])
            N.append(5)
            mu.append(mu_m[0])

            # Find the min chemical potential
            mu_n = [mu_min[i_T]]
            part = partition()
            if (entropy_version == "entropy"):
                runEntropy([T[i_T]], mu_n, [5])
                part.setNcomponents(1)
            elif (entropy_version == "entropy2"):
                runEntropy2([T[i_T]], mu_n, [5])
                part.setNcomponents(1)
            elif (entropy_version == "fspvt"):
                part.setNcomponents(2)

            part.readPVTsimple(self.n_components)
            pmN = part.N1[0]
            pi_mu = mu_n[0]
            if (part.N1[0] < N_low_min):
                mu_n[0] = max(mu_n[0]*1.3, mu_n[0]+0.1)
            elif (part.N1[0] > N_low_max):
                mu_n[0] = min(mu_n[0]/3.0, mu_n[0]-2.0)

            for i in range(max_iter):
                if (N_low_min < part.N1[0] < N_low_max):
                    N_mu_n = part.N1[0]
                    break

                part = partition()
                if (entropy_version == "entropy"):
                    runEntropy([T[i_T]], mu_n, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "entropy2"):
                    runEntropy2([T[i_T]], mu_n, [5])
                    part.setNcomponents(1)
                elif (entropy_version == "fspvt"):
                    part.setNcomponents(2)

                part.readPVTsimple(self.n_components)

                m = (part.N1[0] - pmN) / (mu_n[0] - pi_mu)
                b = part.N1[0] - (m * mu_n[0])
                pmN = part.N1[0]
                pi_mu = mu_n[0]
                mu_n[0] = (N_low_min - b) / m 
                #if (part.N1[0] < N_low_min):
                #    mu_n[0] = max(mu_n[0]*1.3, mu_n[0]+0.1)
                #elif (part.N1[0] > N_low_max):
                #    mu_n[0] = min(mu_n[0]/3.0, mu_n[0]-2.0)
                #else:

            if (i == max_iter):
                print 'change the minimum N range!'
                return

            # Find the chemical potentials in between
            mu_step = abs(mu_m[0] - mu_n[0]) / float(mu_len)
            i_mu = mu_m
            N_step = abs(N_mu_m - N_mu_n) / float(mu_len)
            Nim = N_mu_m
            pmN = Nim
            pi_mu = i_mu[0]
            for im in range(mu_len):
                Nim -= N_step
                i_mu[0] -= mu_step
                Nim_min = 0
                Nim_max = 10000
                for i in range(max_iter):
                    part = partition()
                    if (entropy_version == "entropy"):
                        runEntropy([T[i_T]], mu_n, [5])
                        part.setNcomponents(1)
                    elif (entropy_version == "entropy2"):
                        runEntropy2([T[i_T]], mu_n, [5])
                        part.setNcomponents(1)
                    elif (entropy_version == "fspvt"):
                        part.setNcomponents(2)

                    part.readPVTsimple(self.n_components)


                    if (abs((part.N1[0] - Nim) / Nim) < 0.01):
                        mu_step = abs(mu[im] - i_mu[0])
                        #pmN = part.N1[0]
                        #pi_mu = i_mu[0]
                        break

                    if i == 0:
                        if part.N1[0] > Nim:
                            top = 1
                            bottom = 0
                            above = -0.2
                            below = 0.01
                            factor = above
                        else:
                            top = 0
                            bottom = 1
                            above = -0.01
                            below = 0.2
                            factor = below
                        i_mu[0] += mu_step * factor
                        #m = (i_mu[0] - pi_mu) / (part.N1[0] - pmN)
                        #i_mu[0] = ( (Nim - pmN) * m ) + pi_mu
                    else:
                        if part.N1[0] > Nim:
                            factor = above
                            if top: 
                                # was top before bottom?
                                if bottom and top > bottom:
                                    mu_step = abs(mu[im] - i_mu[0])
                                    break
                            else:
                                top = 0.5
                            #i_mu[0] -= mu_step/above
                            #i_mu[0] -= mu_step * 0.01 * factor
                        else: 
                            factor = below
                            if bottom:
                                # was top before bottom?
                                if top and bottom > top:
                                    mu_step = abs(mu[im] - i_mu[0])
                                    break
                            else:
                                bottom = 0.5
                                
                        i_mu[0] += mu_step * factor
                        #i_mu[0] -= mu_step * (part.N1[0] - Nim) / Nim * factor
                            #i_mu[0] += mu_step/below

                temp.append(T[i_T])
                N.append(5)
                mu.append(i_mu[0])

            temp.append(T[i_T])
            N.append(5)
            mu.append(mu_n[0])

            temp_tot.extend(temp)
            mu_tot.extend(mu)
            N_tot.extend(N)

        temp_tot.reverse()
        mu_tot.reverse()
        N_tot.reverse()

        if (entropy_version == "entropy"):
            runEntropy(temp_tot, mu_tot, N_tot)
        elif (entropy_version == "entropy2"):
            runEntropy2(temp_tot, mu_tot, N_tot)
        
        if (save_file):
            s_file = open( 'TmuN_sVp.dat', 'w')
            for i in range( len(N_tot) ):
                s_file.write('%f %f %f\n' % (temp_tot[i], mu_tot[i], N_tot[i]))
                if i == len(N_tot)-1:
                    s_file.write('stop\n')


    def generateCurveOld(self, entropy_version=1, save_file=False, mu_step=0.0002, mu_low=-0.05, mu_len=30, N_hi=[60, 130]):
        """
        Generate a PVT that covers ideal gas and high pressure
        for micellization
        self-optimizing process:
            for each temperature
            generate a range of chemical potentials, 
            starting from the largest chemical potential in the histograms at the temperature

            if the N from the maximum chemical potential is:
                Below the range increase the maximum chemical potential (mu_m)
                Above the range decrease the maximum chemical potential (mu_m)
            if the N from the minimum chemical potential is:
                Below the range, decrease the step size (self.mu_step), so its to the bottom slower
                Above the range, increase the step size (self.mu_step), so it gets down faster
        """
        if (entropy_version == "entropy" or entropy_version == "entropy2"):
            self.n_components = 1    
        elif (entropy_version == "fspvt"):
            self.n_components = 2

        max_iter = 100
        N_low_min = 0.01
        N_low_max = 0.1
        N_hi_min = N_hi[0]
        N_hi_max = N_hi[1]
        mu_step_start = mu_step
        self.mu_step = mu_step
        self.mu_low = mu_low
        self.mu_len = mu_len
        mu_max, T = self.getMuMaxT()
        temp_tot = []
        mu_tot = []
        N_tot = []
        for i_T in range(len(T)):
            mu_m = mu_max[i_T]
            # change mu until the upper and lower N are within both ranges
            for i in range(max_iter):
                temp, mu, N = self.generateMu(mu_m, T[i_T])
                part = partition()
                part.setNcomponents(self.n_components)
                if (entropy_version == "entropy"):
                    runEntropy(temp, mu, N)
                elif (entropy_version == "entropy2"):
                    runEntropy2(temp, mu, N)
                elif (entropy_version == "fspvt"):
                    runFSPVT(temp, mu1, mu2)

                part.readPVTsimple(self.n_components)
                N_mu_high = part.N1[mu_len-1]
                if (N_mu_high < N_hi_min):
                    mu_m += 0.01
                elif (N_mu_high > N_hi_max):
                    mu_m -= 0.01

                N_m_low = part.N1[0]
                if (N_m_low > N_low_max):
                    self.mu_step *= 2.0
                elif (N_m_low < N_low_min):
                    self.mu_step /= 1.25
                else:
                    if (N_hi_min < N_mu_high < N_hi_max):
                        break

            self.mu_step = mu_step_start
            temp_tot.extend(temp)
            mu_tot.extend(mu)
            N_tot.extend(N)
        if (entropy_version == 1):
            runEntropy(temp_tot, mu_tot, N_tot)
        else:
            runEntropy2(temp_tot, mu_tot, N_tot)

        if (save_file):
            s_file = open( 'TmuN_sVp.dat', 'w')
            for i in range( len(N_tot) ):
                s_file.write('%f %f %f\n' % (temp_tot[i], mu_tot[i], N_tot[i]))
                if i == len(N_tot)-1:
                    s_file.write('stop\n')
            
    def calcError(self, version="entropy2", write_latex=False):
        """
        Calculate the relative error between the calculated
        pratition function and the
        simlulation
        """
        import shutil 
        T = []
        mu1 = []
        N1 = []
        E = []

        if (version == "entropy" or version == "entropy2"): 
            self.n_components = 1
        elif (version == "fspvt"):
            self.n_components = 2
            N2 = []
            mu2 = []

        # copy files that may have been run from entropy2, just to save the user in case they forget
        cp_pvt = True
        try:
            shutil.copyfile('./pvt.dat','./pvt_tmp.dat')
        except IOError:
            cp_pvt = False

        cp_hs2 = True
        try:
            shutil.copyfile('./input_hs2.dat','./input_hs2_tmp.dat')
        except IOError:
            cp_hs2 = False

        N1min = 1000
        for i in range( len(self.runs) ):
            if (version == "entropy" or version == "entropy2"): 
                read_err = self.readOneComponent(self.runs[i])
            elif (version == "fspvt"):
                read_err = self.readTwoComponent(self.runs[i])
                mu2.append(self.getmu2())
                N2.append(self.getN2ave(2))
            # find the <N> and <E> from histograms
            T.append(self.getT())
            mu1.append(self.getmu1())
            N1.append(self.getN1ave(self.n_components))
            E.append(self.getEave(self.n_components))
            if min(N1) < N1min:
                N1min = min(N1)

        if N1min != 0:
            print 'WARNING none of the histograms go down to N=0'

        # run entropy with runs used to develop partition function
        part = partition()
        if (version == "entropy"):
            runEntropy(T, mu1, N1)
            part.setNcomponents(1)
        elif (version == "entropy2"):
            runEntropy2(T, mu1, N1)
            part.setNcomponents(1)
        elif (version == "fspvt"):
            runFSPVT(T, mu1, mu2)
            part.setNcomponents(2)
        else:
            print 'version needs to entropy2, entropy or fspatch'
            return

        part.readPVTsimple(self.n_components)
        if (False): #version == 1):
            for i in range( len(self.runs) ):
                part.E[i] *= part.N1[i]

        if self.n_components == 1:
            print '       run      |  kT      mu   |  N_sim    E_sim  |  N_part   E_part |   %e(N)   %e(E)'
            print '----------------+---------------+------------------+------------------+-----------------'
            if write_latex:
                latex_file = open( 'error.tex', 'w')
                latex_file.write('\\begin{table}\n')
                latex_file.write('  \\begin{center}\n')
                latex_file.write('    \\begin{tabular}{c c | c c | c c | c c}\n')
                latex_file.write('    \hline\n')
                latex_file.write('    \hline\n')
                latex_file.write('    $kT$ & $\mu$ & $\left<N\\right>_{\\text{sim}}$ & '
                                                   '$\left<E\\right>_{\\text{sim}}$ & '
                                                   '$\left<N\\right>_{\ln\\text{Z}}$ & '
                                                   '$\left<E\\right>_{\ln \\text{Z}}$ & '
                                      '\%e($\left<N\\right>$) & \%e($\left<E\\right>$)\\\\\n')
                latex_file.write('     \hline\n')
    
            for i in range( len(self.runs) ):
                N1_err = float(part.N1[i] - N1[i]) / N1[i] * 100.0 
                E_err = float(part.E[i] - E[i]) / (E[i]+1E-8) * 100.0
                print ('%15s | %6.3f %6.3f | %6.2f %9.3f | %6.2f %9.3f | %7.2f %7.2f' % 
                      (self.runs[i], T[i]*self.Tconv, mu1[i]*self.Muconv, N1[i], E[i]*self.Econv, 
                      part.N1[i], part.E[i]*self.Econv, N1_err, E_err))
                if write_latex:
                    latex_file.write('    %6.3f %s %6.3f %s %6.2f %s %9.3f %s %6.2f %s %9.3f %s %7.2f %s %7.2f \\\\\n' % 
                                    ( T[i]*self.Tconv, '&', mu1[i]*self.Muconv, '&', N1[i], '&', E[i]*self.Econv, '&',
                                      part.N1[i], '&', part.E*self.Econv[i], '&', N1_err, '&', E_err))
    
            if write_latex:
                latex_file.write('    \hline\n')
                latex_file.write('    \hline\n')
                latex_file.write('    \end{tabular}\n')
                latex_file.write('  \end{center}\n')
                latex_file.write('  \\vspace{-0pt}\n')
                latex_file.write('  \\caption{Comparison of the $\left<N\\right>$ and $\left<E\\right>$ '
                                 'values from simulations and predicted from the partition function $\ln$Z.}\n')
                latex_file.write('  \\label{tab:error_table}\n')
                latex_file.write('\end{table}\n')
                latex_file.close()

        elif self.n_components == 2:
            print '       run      |  kT       mu1     mu2   | N1_sim  N2_sim   E_sim  | N1_part N2_part E_part  |   %e(N1)  %e(N2)  %e(E)'
            print '----------------+-------------------------+-------------------------+-------------------------+-------------------------'
            if write_latex:
                latex_file = open( 'error.tex', 'w')
                latex_file.write('\\begin{table}\n')
                latex_file.write('  \\begin{center}\n')
                latex_file.write('    \\begin{tabular}{c c | c c c | c c c | c c c}\n')
                latex_file.write('    \hline\n')
                latex_file.write('    \hline\n')
                latex_file.write('    $kT$ & $\mu_1$ & $\left<N_1\\right>_{\\text{sim}}$ & '
                                                   '$\left<N_2\\right>_{\\text{sim}}$ & '
                                                   '$\left<E\\right>_{\\text{sim}}$ & '
                                                   '$\left<N_1\\right>_{\ln\\text{Z}}$ & '
                                                   '$\left<N_2\\right>_{\ln\\text{Z}}$ & '
                                                   '$\left<E\\right>_{\ln \\text{Z}}$ & '
                                                   '\%e($\left<N_1\\right>$) & '
                                                   '\%e($\left<N_2\\right>$) & '
                                                   '\%e($\left<E\\right>$)\\\\\n')
                latex_file.write('     \hline\n')
    
            for i in range( len(self.runs) ):
                N1_err = float(part.N1[i] - N1[i]) / N1[i] * 100.0 
                N2_err = float(part.N2[i] - N2[i]) / N2[i] * 100.0 
                E_err = float(part.E[i] - E[i]) / (E[i]+1E-8) * 100.0
                print ('%15s | %6.3f %6.3f %6.3f | %6.2f %6.2f %9.3f | %6.2f %6.2f %9.3f | %7.2f %7.2f %7.2f' % 
                      (self.runs[i], T[i]*self.Tconv, mu1[i]*self.Muconv, mu2[i]*self.Muconv, N1[i], N2[i], E[i]*self.Econv, 
                      part.N1[i], part.N2[i], part.E[i]*self.Econv, N1_err, N2_err, E_err))
                if write_latex:
                    latex_file.write('    %6.3f %s %6.3f %s %6.3f %s %6.2f %s %6.2f %s %9.3f %s %6.2f %s %6.2f %s %9.3f %s %7.2f %s %7.2f %s %7.2f \\\\\n' % 
                                    ( T[i]*self.Tconv, '&', mu1[i]*self.Muconv, '&', mu2[i]*self.Muconv, '&', N1[i], '&', N2[i], '&', E[i]*self.Econv, '&',
                                      part.N1[i], '&', part.N2[i], '&', part.E*self.Econv[i], '&', N1_err, '&', N2_err, '&', E_err))
    
            if write_latex:
                latex_file.write('    \hline\n')
                latex_file.write('    \hline\n')
                latex_file.write('    \end{tabular}\n')
                latex_file.write('  \end{center}\n')
                latex_file.write('  \\vspace{-0pt}\n')
                latex_file.write('  \\caption{Comparison of the $\left<N_1\\right>$, $\left<N_2\\right>$ and $\left<E\\right>$ '
                                 'values from simulations and predicted from the partition function $\ln$Z.}\n')
                latex_file.write('  \\label{tab:error_table}\n')
                latex_file.write('\end{table}\n')
                latex_file.close()

        if (cp_hs2):
            shutil.move('./input_hs2_tmp.dat','./input_hs2.dat')
        #if (cp_pvt):
            #shutil.move('./pvt_tmp.dat','./pvt.dat')

class partition(object):
    """
    Read the PVT file from the entropy program suite
    """

    def setNcomponents(self, n_components):
        self.n_components = n_components

    def getE(self):
        return self.E

    def getN1(self):
        return self.N1

    def getN2(self):
        if self.n_components >= 2:
            return self.N2
        else:
            return None

    def readPVTsimple(self, n_components):
        pwd = os.getcwd()
        pvt_name = 'pvt.dat'
        # Check that file exists
        if not os.path.isfile(pvt_name):
            sys.exit("Error! file "+pvt_name+" does not exist.")
    
        params = {}
    
        # get the ln of partition function
        pvt_file = open(pvt_name,'r')
        # read header
        data = pvt_file.readline().strip().split()
        while ('T' not in data):
            data = pvt_file.readline().strip().split()
    
        self.mu1 = []
        self.N1 = []
        self.E = []
        self.lnZ = []
        
        if self.n_components >= 2:
            self.mu2 = []
            self.N2 = []

        elif self.n_components == 1:
            self.lnZliq = []
            self.lnZgas = []
    
        for line in pvt_file:
            data = line.strip().split()
            if self.n_components == 1:
                self.mu1.append(float(data[1]) )
                self.N1.append( float(data[2]) )
                self.E.append( float(data[3]) )
                self.lnZ.append( float(data[4]) )
                self.lnZliq.append( float(data[5]) )
                self.lnZgas.append( float(data[6]) )

            elif self.n_components == 2:
                self.mu1.append(float(data[1]) )
                self.mu2.append(float(data[2]) )
                self.E.append( float(data[3]) )
                self.N1.append( float(data[4]) )
                self.N2.append( float(data[5]) )
                self.lnZ.append( float(data[6]) )
    
        pvt_file.close()

def main(argv=None):
    parser = argparse.ArgumentParser(description='Analyze simulation histograms '
                                    'and compare or generate a partition function. ')
    parser.add_argument("-i","--input_file", metavar='inFile', type=str,
                   help='File containing list of historgram file runnames '
                        '(e.g. input_hs.dat.')
    parser.add_argument("-e","--error", action="store_true",
                   help='Calculate the error in the N and E by the '
                        'partition function compared to the simulation')
    parser.add_argument("-v","--version",  metavar='histogram reweighting version',
                        type=str, choices=["entropy","entropy2","fspvt"],
                   help='Generate the pvt.dat file for CMC calculation')
    parser.add_argument("-g","--generate", action="store_true",
                   help='Use the partiton function generated a nice series mu values for pvt.dat')
    parser.add_argument("-t","--temperature", type=float, nargs="+",
                   help='temperature you want to pvt to be calculated at '
                        'assumed to be all those in the input_file.')
    parser.add_argument("--Npoints", type=int,
                   help='Number of points to calculate')
    parser.add_argument("--N_low", type=float, nargs="+",
                   help='Lower bound range of N particles in the generation. eg:'
                        '--N_low 0.01 0.1')
    parser.add_argument("--N_hi", type=float, nargs="+",
                   help='Upper bound range of N particles in the generation. eg:'
                        '--N_hi 50 80')
    parser.add_argument("--n1_eq_n2", action="store_true",
                   help='Set components 1 and 2 equal eachother. '
                        'Useful for charge neutrality')
    parser.add_argument("-s","--save", action="store_true",
                   help='Save the T, mu, N generated.  Use the output without simulationVpartition by:'
                        ' entropy2.x < TmuN_sVp.dat')
    parser.add_argument("--latex", action="store_true",
                   help='Output error table in latex format/')
    parser.add_argument("--Tunits", type=str, choices=['K','kT'],
                   help='Temperature units for output')
    parser.add_argument("--Eunits", type=str, choices=['kJ/mol', 'K'],
                   help='Energy units for output (should match Energy units)')
    parser.add_argument("--Muunits", type=str, choices=['kJ/mol', 'K'],
                   help='Chemical potential units for output (should match Energy units)')

    if (parser.parse_args().input_file):
        runs = readRunsFile(parser.parse_args().input_file, parser.parse_args().version)
    else:
        print 'Must give input file e.g.:'
        print '-i input_hs.dat'
        return

    # read in the input_hs.dat
    HIS = hisFile(runs, parser.parse_args().temperature)

    if (parser.parse_args().error):
        if(parser.parse_args().Tunits): HIS.setTconv( parser.parse_args().Tunits )
        if(parser.parse_args().Muunits): HIS.setMuconv( parser.parse_args().Muunits )
        if(parser.parse_args().Eunits): HIS.setEconv( parser.parse_args().Eunits )
        HIS.calcError(parser.parse_args().version, parser.parse_args().latex)

    elif (parser.parse_args().generate):
        mu_step_size = 0.0001
        mu_min = -0.1
        if parser.parse_args().Npoints:
            mu_length = parser.parse_args().Npoints
        else:
            mu_length = 50
        if parser.parse_args().N_low:
            N_min_range = parser.parse_args().N_low
        else:
            N_min_range = [0.01, 0.1]
        if parser.parse_args().N_hi:
            N_max_range = parser.parse_args().N_hi
        else:
            N_max_range = [100, 150]

        max_iterations = 100
        HIS.generate(parser.parse_args().version, parser.parse_args().save, mu_length, N_min_range, N_max_range, max_iterations, parser.parse_args().n1_eq_n2)
        #HIS.generateCurveOld(parser.parse_args().generate, parser.parse_args().save, mu_step_size, mu_min, mu_length, N_max_range)

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
