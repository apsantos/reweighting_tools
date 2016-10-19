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
from scipy.optimize import curve_fit
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

def readRunsFile(filename, skiplines=0):
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
        return

    def getmu(self):
        return self.mu

    def getT(self):
        return self.temp

    def getEave(self):
        ave = 0
        for i in range(len(self.histogram[0,:])):
            ave += self.E[i] * sum(self.histogram[:,i]) / float(self.his_sum)
        return ave

    def getNave(self):
        ave = 0
        for iN in range(len(self.histogram[:,0])):
            ave += iN * sum(self.histogram[iN,:]) / float(self.his_sum)
        return ave

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

    def setJ(self, J):
        self.J = J

    def setBox(self, x, y, z):
        self.box = [x, y, z]

    def setTemp(self, temperature):
        self.temp = temperature

    def setMu(self, mu):
        self.mu = mu[0]

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
            ofile.write("%12.6f %14.6f %22.12f %10f %10f %10f\n" % ( self.temp, self.mu, self.width, self.box[0], self.box[1], self.box[2] ))
    
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

    def read(self, file_root):
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
        self.mu = float(data[1])
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
            read_err = self.read(self.runs[i])
            if (read_err == 0):
                if (plot_phi):
                    X, Y = np.meshgrid(self.E, self.phi)
                else:
                    X, Y = np.meshgrid(self.E, self.N)

                if (legend):
                    labels[i] = '%s: T=%s, $\mu$=%s' % (self.runs[i], self.temp, self.mu)
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
            read_err = self.read(self.runs[i])
            T_temp = self.getT()
            mu_temp = self.getmu()
                
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

    def getMuMinMaxT(self):
        """
        Get the runs' temperature list and number of mu
        """
        T = []
        mu_min = []
        mu_max = []
        mu_low = -1000
        mu_high = 1000
        for i in range( len(self.runs) ):
            read_err = self.read(self.runs[i])
            T_temp = self.getT()
            mu_temp = self.getmu()
                
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

    def generateCurve2(self, entropy_version=1, save_file=False, mu_len=30, N_lo=[0.01,0.1],N_hi=[60, 130], max_iter=100):
        """
        Generate a PVT that covers ideal gas and high pressure
        for micellization
        self-optimizing process:
            for each temperature
            find a chemical potential that lands in the range of upper N
            find a chemical potential that lands in the range of lower N
            starting from the top step down ramping up the step size to scale somewhat with the N(\mu) function
        """
        N_low_min = N_lo[0]
        N_low_max = N_lo[1]
        N_hi_min = N_hi[0]
        N_hi_max = N_hi[1]
        mu_min, mu_max, T = self.getMuMinMaxT()
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
                if (entropy_version == 1):
                    runEntropy([T[i_T]], mu_m, [5])
                else:
                    runEntropy2([T[i_T]], mu_m, [5])
                part = partition()
                part.readPVTsimple()
                if (part.N[0] < N_hi_min):
                    mu_m[0] = max(mu_m[0]*1.1, mu_m[0]+0.5)
                elif (part.N[0] > N_hi_max):
                    mu_m[0] = min(mu_m[0]/1.1, mu_m[0]-0.01)
                else:
                    N_mu_m = part.N[0]
                    break
            if (i == max_iter):
                print 'change the maximum N range!'
                return

            temp.append(T[i_T])
            N.append(5)
            mu.append(mu_m[0])

            # Find the min chemical potential
            mu_n = [mu_min[i_T]]
            if (entropy_version == 1):
                runEntropy([T[i_T]], mu_n, [5])
            else:
                runEntropy2([T[i_T]], mu_n, [5])
            part = partition()
            part.readPVTsimple()
            pmN = part.N[0]
            pi_mu = mu_n[0]
            if (part.N[0] < N_low_min):
                mu_n[0] = max(mu_n[0]*1.3, mu_n[0]+0.1)
            elif (part.N[0] > N_low_max):
                mu_n[0] = min(mu_n[0]/3.0, mu_n[0]-2.0)

            for i in range(max_iter):
                if (N_low_min < part.N[0] < N_low_max):
                    N_mu_n = part.N[0]
                    break

                if (entropy_version == 1):
                    runEntropy([T[i_T]], mu_n, [5])
                else:
                    runEntropy2([T[i_T]], mu_n, [5])
                part = partition()
                part.readPVTsimple()

                m = (part.N[0] - pmN) / (mu_n[0] - pi_mu)
                b = part.N[0] - (m * mu_n[0])
                pmN = part.N[0]
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
        N_low_min = N_lo[0]
        N_low_max = N_lo[1]
        N_hi_min = N_hi[0]
        N_hi_max = N_hi[1]
        mu_min, mu_max, T = self.getMuMinMaxT()
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
                if (entropy_version == 1):
                    runEntropy([T[i_T]], mu_m, [5])
                else:
                    runEntropy2([T[i_T]], mu_m, [5])
                part = partition()
                part.readPVTsimple()
                if (part.N[0] < N_hi_min):
                    mu_m[0] = max(mu_m[0]*1.1, mu_m[0]+0.5)
                elif (part.N[0] > N_hi_max):
                    mu_m[0] = min(mu_m[0]/1.1, mu_m[0]-0.01)
                else:
                    N_mu_m = part.N[0]
                    break
            if (i == max_iter):
                print 'change the maximum N range!'
                return

            temp.append(T[i_T])
            N.append(5)
            mu.append(mu_m[0])

            # Find the min chemical potential
            mu_n = [mu_min[i_T]]
            if (entropy_version == 1):
                runEntropy([T[i_T]], mu_n, [5])
            else:
                runEntropy2([T[i_T]], mu_n, [5])
            part = partition()
            part.readPVTsimple()
            pmN = part.N[0]
            pi_mu = mu_n[0]
            if (part.N[0] < N_low_min):
                mu_n[0] = max(mu_n[0]*1.3, mu_n[0]+0.1)
            elif (part.N[0] > N_low_max):
                mu_n[0] = min(mu_n[0]/3.0, mu_n[0]-2.0)

            for i in range(max_iter):
                if (N_low_min < part.N[0] < N_low_max):
                    N_mu_n = part.N[0]
                    break

                if (entropy_version == 1):
                    runEntropy([T[i_T]], mu_n, [5])
                else:
                    runEntropy2([T[i_T]], mu_n, [5])
                part = partition()
                part.readPVTsimple()

                m = (part.N[0] - pmN) / (mu_n[0] - pi_mu)
                b = part.N[0] - (m * mu_n[0])
                pmN = part.N[0]
                pi_mu = mu_n[0]
                mu_n[0] = (N_low_min - b) / m 
                #if (part.N[0] < N_low_min):
                #    mu_n[0] = max(mu_n[0]*1.3, mu_n[0]+0.1)
                #elif (part.N[0] > N_low_max):
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
                    if (entropy_version == 1):
                        runEntropy([T[i_T]], i_mu, [5])
                    else:
                        runEntropy2([T[i_T]], i_mu, [5])
                    part = partition()
                    part.readPVTsimple()


                    #print part.N[0]
                    if (abs((part.N[0] - Nim) / Nim) < 0.01):
                        mu_step = abs(mu[im] - i_mu[0])
                        #pmN = part.N[0]
                        #pi_mu = i_mu[0]
                        break

                    if i == 0:
                        if part.N[0] > Nim:
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
                        #m = (i_mu[0] - pi_mu) / (part.N[0] - pmN)
                        #i_mu[0] = ( (Nim - pmN) * m ) + pi_mu
                    else:
                        if part.N[0] > Nim:
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
                        #i_mu[0] -= mu_step * (part.N[0] - Nim) / Nim * factor
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
                if (entropy_version == 1):
                    runEntropy(temp, mu, N)
                else:
                    runEntropy2(temp, mu, N)
                part = partition()
                part.readPVTsimple()
                N_mu_high = part.N[mu_len-1]
                if (N_mu_high < N_hi_min):
                    mu_m += 0.01
                elif (N_mu_high > N_hi_max):
                    mu_m -= 0.01

                N_m_low = part.N[0]
                if (N_m_low > N_low_max):
                    self.mu_step *= 2.0
                elif (N_m_low < N_low_min):
                    self.mu_step /= 1.25
                else:
                    if (N_hi_min < N_mu_high < N_hi_max):
                        break

                #print N_mu_high, N_m_low

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
            
    def calcError(self):
        """
        Calculate the relative error between the calculated
        pratition function and the
        simlulation
        """
        import shutil 
        T = []
        mu = []
        N = []
        E = []

        # copy files that mave have been run from entropy2, just to save the user in case they forget
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

        for i in range( len(self.runs) ):
            read_err = self.read(self.runs[i])
            # find the <N> and <E> from histograms
            T.append(self.getT())
            mu.append(self.getmu())
            N.append(self.getNave())
            E.append(self.getEave())

        # run entropy with runs used to develop partition function
        runEntropy(T, mu, N)
        # read pvt.dat 
        part = partition()
        part.readPVTsimple()
        # write %error
        print '       run      |   T      mu   |  N_sim    E_sim  |  N_part   E_part |   %e(N)   %e(E)'
        print '----------------+---------------+------------------+------------------+-----------------'
        for i in range( len(self.runs) ):
            N_err = float(part.N[i] - N[i]) / N[i] * 100.0 
            E_err = float(part.E[i]*part.N[i] - E[i]) / (E[i]+1E-8) * 100.0
            print ('%15s | %6.3f %6.2f | %6.2f %9.3f | %6.2f %9.3f | %7.2f %7.2f' % 
                  (self.runs[i], T[i], mu[i], N[i], E[i], 
                  part.N[i], part.E[i]*part.N[i], N_err, E_err))
        if (cp_hs2):
            shutil.move('./input_hs2_tmp.dat','./input_hs2.dat')
        if (cp_pvt):
            shutil.move('./pvt_tmp.dat','./pvt.dat')

class partition(object):
    """
    Read the PVT file from the entropy program suite
    """

    def getE(self):
        return self.E

    def getN(self):
        return self.N

    def readPVTsimple(self):
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
    
        self.mu = []
        self.N = []
        self.E = []
        self.lnZ = []
        self.lnZliq = []
        self.lnZgas = []
    
        for line in pvt_file:
            data = line.strip().split()
            self.mu.append(float(data[1]) )
            self.N.append( float(data[2]) )
            self.E.append( float(data[3]) )
            self.lnZ.append( float(data[4]) )
            self.lnZliq.append( float(data[5]) )
            self.lnZgas.append( float(data[6]) )
    
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
    parser.add_argument("-g","--generate", metavar='entropy version', 
                        type=int, choices=[1,2],
                   help='Generate the pvt.dat file for CMC calculation, '
                        'can either use the partiton function generated '
                        'by entropy[1,2].x.  Assumes 2 if none given.')
    parser.add_argument("-t","--temperature", type=float, nargs="+",
                   help='temperature you want to pvt to be calculated at '
                        'assumed to be all those in the input_file.')
    parser.add_argument("-s","--save", action="store_true",
                   help='Save the T, mu, N generated.  Use the output without simulationVpartition by:'
                        ' entropy2.x < TmuN_sVp.dat')

    if (parser.parse_args().input_file):
        runs = readRunsFile(parser.parse_args().input_file)
    else:
        print 'Must give input file e.g.:'
        print '-i input_hs.dat'
        return

    # read in the input_hs.dat
    HIS = hisFile(runs, parser.parse_args().temperature)

    if (parser.parse_args().error):
        HIS.calcError()
    elif (parser.parse_args().generate):
        mu_step_size = 0.0001
        mu_min = -0.1
        mu_length = 50
        N_min_range = [0.01, 0.1]
        N_max_range = [70, 90]
        #N_max_range = [400, 450]
        HIS.generateCurve2(parser.parse_args().generate, parser.parse_args().save, mu_length, N_min_range, N_max_range)
        #HIS.generateCurveOld(parser.parse_args().generate, parser.parse_args().save, mu_step_size, mu_min, mu_length, N_max_range)

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
