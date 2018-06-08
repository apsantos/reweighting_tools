#!/usr/bin/env python
"""
Author : Andrew P. Santos
Date   : May, 2014
Institution: Princeton University

Purpose: Calculate the pressure from the resulting "pvt.dat" file output from entropy or entropy2 

	The top two def sections are called from the main section which is below
"""
import os, sys, argparse
import numpy as np
import math as ma
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1.2
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1.2

def sortArray(x, y):
    order_indicies = []
    delete_indicies = []
    for i in range( 1, len(x) ):
        if (x[i] == x[i-1]):
            print 'duplicate at', i, x[i]
            delete_indicies.append(i)

        if (len(delete_indicies) > 4):
            return x, y, True

    for i in delete_indicies:
        x = np.delete(x, i)
        y = np.delete(y, i)

    for i in range( 1, len(x) ):
        if (x[i] < x[i-1]):
            x_temp = x[i]
            y_temp = y[i]
            x[i] = x[i-1]
            x[i-1] = x_temp
            y[i] = y[i-1]
            y[i-1] = y_temp
            print 'reorder data at', i, x[i]
            order_indicies.append(i)

        if (( len(order_indicies)+len(delete_indicies) ) > 4):
            return x, y, True

    return x, y, False

#def sigmoidIntegral(x, F0, A1, A2, x0, dx):
#    return F0 + A1*x + dx*(A2 - A1)*np.log((1.0 + np.exp((x - x0)/dx))/(1.0 + np.exp(-x0/dx)))
def sigmoidIntegral(x, F0, A2, x0, dx):
    return F0 + x + dx*(A2 - 1.0)*np.log((1.0 + np.exp((x - x0)/dx))/(1.0 + np.exp(-x0/dx)))

def getSpline(x, y, n_points=1000):
    x, y, order_error = sortArray(x, y)
    if (order_error):
        print 'too many entries are not in order'
        print 'The program assumes that the lowest N is first in the pvt.dat file for each temperature.'
        sys.exit(main())

    tck = interpolate.splrep(x, y, s=0)
    x_spl = np.linspace(0, max(x), n_points)
    return x_spl, interpolate.splev(x_spl,tck, der=0), interpolate.splev(x_spl,tck, der=1), interpolate.splev(x_spl,tck, der=2)

def maxCurvature(x, y, plot=False):
    kappa_max = 0.0
    x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y, 1000)
    dx = x_spl[2] - x_spl[0]
    dx2 = dx * dx
    kappa_max_e = dx2
    x_max_e = dx2
    kappa = np.zeros( (len(x_spl)) )
    for i in range(1,len(x_spl)-1):
        kappa[i] = -(y_spl_d2[i]) / (1.0 + (y_spl_d1[i])**2.0)**(1.5)
        if (kappa[i] > kappa_max):
            kappa_max = kappa[i]
            y_max = y_spl[i]
            x_max = x_spl[i]
            i_max = i

    if (plot):
        fig = plt.figure()
        plt.plot(x_spl, y_spl, 'r', label='Spline fit')
        plt.plot([0, x_max], [0, x_max], 'k')
        plt.plot([x_max, x_max], [0, x_max], 'k')
        plt.plot(x, y, 'og', label='Histogram rewieghting data')
        border = plt.legend( numpoints=1, prop={'size':12}, loc=4)
        border.draw_frame(False)
        plt.xlabel("$\phi_{tot}$",fontsize=20)
        plt.ylabel("$\Pi$",fontsize=20)
        
        plt.show()

    return x_max, x_max_e, x_spl, y_spl, kappa

def zero2ndDerivativeIntercept(x, y, method='spline', plot=False):
    d2ydx2_min = 10.0
    d2ydx2_min_e = 0.0
    x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y)
    dx = x_spl[2] - x_spl[0]
    dx2 = dx * dx
    d2ydx2_inter_e = dx2
    for i in range(1,len(x_spl)-1):
        d2ydx2 = y_spl_d2[i]
        if (d2ydx2 < d2ydx2_min):
            d2ydx2_min = d2ydx2
            y_min = y_spl[i]
            i_min = i

    x_zero = 0
    x_inter = 0
    for i in range(i_min,len(x_spl)-1):
        d2ydx2 = y_spl_d2[i]
        if (d2ydx2 > 0.0000):
            d2ydx2_zero = d2ydx2
            x_zero = x_spl[i]
            x_inter = y_spl[i]
            break
            
    if (plot):
        fig = plt.figure()
        plt.plot(x_spl, y_spl, 'r', label='Spline fit')
        plt.plot([0, x_inter], [0, x_inter], 'k')
        plt.plot([x_inter, x_zero], [x_inter, x_inter], 'k',linestyle='-')
        plt.plot([0.02738, x[len(x)-7]], [0.02738, y[len(x)-7]], 'k',linestyle='--', label='Method 1') #Slope-intercept method')
        plt.plot([0.02738, 0.02738], [0, 0.02738], 'k',linestyle='--')
        plt.plot([0.02804457, 0.02804457], [0, 0.02804457], 'k',linestyle=':', label='Method 2 w/ spline')#Min-in-2nd-derivative method w/ spline')
        plt.plot([0.02434607, 0.02434607], [0, 0.02434607], 'k',linestyle='-.', label='Method 2 w/out spline')#Min-in-2nd-derivative method w/out spline')
        plt.plot([x_inter, x_inter], [0, x_inter], 'k',linestyle='-', label='Method 3')#Zero-derivative-gas-intercept method')
        plt.plot(x, y, 'og', label='Histogram rewieghting data')
        border = plt.legend( numpoints=1, prop={'size':12}, loc=4)
        border.draw_frame(False)
        plt.xlabel("$\phi_{tot}$",fontsize=20)
        plt.ylabel("$\Pi$",fontsize=20)
        
        plt.show()

    x_inter_e = d2ydx2_inter_e

    return x_inter, x_inter_e

def max2ndDerivative(x, y, method='spline', plot=False):
    """
    Central finite difference numerical 2nd derivative
    """
    d2ydx2_min = 10.0
    d2ydx2_min_e = 0.0
    x, y, order_error = sortArray(x, y)
    if (order_error):
        sys.exit(main())

    if (method == 'spline'):
        tck = interpolate.splrep(x, y, s=0)
        x_spl = np.linspace(0, max(x), 1000)
        y_spl = interpolate.splev(x_spl,tck, der=0)
        y_spl_d2 = interpolate.splev(x_spl,tck, der=2)
        dx = x_spl[2] - x_spl[0]
        dx2 = dx * dx
        d2ydx2_min_e = dx2
        for i in range(1,len(x_spl)-1):
            d2ydx2 = y_spl_d2[i]
            if (d2ydx2 < d2ydx2_min):
                d2ydx2_min = d2ydx2
                x_min = x_spl[i]
                y_min = y_spl[i]
            
    if (plot):
        fig = plt.figure()
        plt.plot(x_spl, y_spl, 'r', label='spline fit')
        plt.plot([0, x_min], [0, x_min], 'k')
        plt.plot([x_min, x_min], [0, x_min], 'k',linestyle='dashed', label='minimum in 2nd derivative via spline')
        #plt.plot([0.01654329, 0.01654329], [0, 0.01654329], 'k',linestyle='dotted', label='minimum in 2nd derivative via finite difference in data')
        plt.plot(x, y, 'og', label='histogram rewieghting data')
        border = plt.legend( numpoints=1, prop={'size':12}, loc=2)
        border.draw_frame(False)
        plt.xlabel("$\phi_{tot}$",fontsize=18)
        plt.ylabel("$\Pi$",fontsize=18)
        
        plt.show()

    elif (method == 'finite'):
        for i in range(2, len(x)-1 ):
            dx = x[i+1] - x[i-1]
            if (abs(dx) > 0.0005):
                d2ydx2 = (y[i+1] - 2*y[i] + y[i-1]) / (dx * dx)
                if (d2ydx2 < d2ydx2_min):
                    d2ydx2_min = d2ydx2
                    d2ydx2_min_e = (dx * dx)
                    x_min = x[i]
                    y_min = y[i]

    x_min_e = d2ydx2_min_e
    y_min_e = d2ydx2_min_e

    return x_min, y_min, d2ydx2_min, x_min_e

def linFit(x, y):
    """
    Using Harris Quantitative Chemical Analysis chapter 4
    """
    x_sum = np.sum(x)
    xx_sum = np.sum( x * x )
    y_sum = np.sum(y)
    yy_sum = np.sum( y * y )
    xy_sum = np.sum( x * y )
    n = len(x)

    D = np.linalg.det([[xx_sum, x_sum], 
                       [x_sum,      n]])
    if D == 0:
        print ('lnZ(N) needs to be an independent function, i.e.' +
              ' change the chemical potential values used')
        return 0, 0, 0, 0, 0, -1
    slope = np.linalg.det([[xy_sum, x_sum], 
                           [y_sum,      n]])/D
    intercept = np.linalg.det([[xx_sum, xy_sum], 
                               [x_sum,   y_sum]])/D

    dev = y - ( slope * x + intercept )
    dev2_sum = np.sum( dev * dev )
    
    y_std2 = dev2_sum / float(n-2)
    slope_std = ( y_std2 * n / D )**0.5
    intercept_std = ( y_std2 * xx_sum / D )**0.5
    r = ( ( n * xy_sum - x_sum * y_sum ) / 
          ( ( n*xx_sum - x_sum**2 ) * ( n*yy_sum - y_sum**2 )  )**0.5)
    return slope, intercept, slope_std, intercept_std, r, 0

class partition2pressure(object):
    """
    Read and process PVT and cmc files from the entropy program suite
    """
    def __init__(self):
        self.kB = 1.38064852E-23 # J K-1
        self.convert_kPa = 1.E27 # J / A^3 to kPa
        self.A3_to_m3 = 1.E-30 # A^3 to m^3
        self.atomic_to_kJmol = 100.0
        self.bar_to_kPa = 100.0
        self.N1a = 6.022140857E23 # 1/ mol
        self.gas_const = 0.0083144621 # kJ / (K mol)
        self.calc_conc = False
        self.calc_rho = False
        self.calc_phi = False
        self.calc_num_dens = True
        self.calc_pressure = True
        self.calc_lnZ = False
        self.colors = ['red', 'blue', 'green', 'magenta', 'cyan', 
                       'yellow', 'black', 'darkgoldenrod','firebrick', 
                       'purple', 'burlywood', 'chartreuse', 'red', 
                       'blue', 'green', 'magenta', 'cyan']
        # including the points very near min or max N/density can confuse
        # the fitting process, so we as  a shorter set
        self.begin_fit = 10
        self.end_fit = 3

    def addParser(self, parser):
        """
        Get relevant values from an argparse parser
        """
        self.vol = 1.0
        n_sides = len(parser.parse_args().box_length)
        for i in range(3):
            if n_sides == 1:
                self.vol *= parser.parse_args().box_length[0]
            elif n_sides == 3:
                self.vol *= parser.parse_args().box_length[i]
        if (n_sides not in [1,3]):
            print 'You can only give 1 or 3 entries for the box size'
            
        self.calc_conc = parser.parse_args().concentration
        self.calc_num_dens = parser.parse_args().num_density

        if (parser.parse_args().n_components == 1):
            self.n_spec = 1
        elif (parser.parse_args().n_components == 2):
            self.n_spec = 2
        else:
            self.n_spec = 1

        if (parser.parse_args().density):
            self.calc_rho = True
            self.mass = parser.parse_args().density

        elif (parser.parse_args().phi):
            self.calc_phi = True
            self.mono_vol = parser.parse_args().phi

        self.calc_pressure = parser.parse_args().pressure
        self.calc_partition = parser.parse_args().partition
        self.calc_t_norm = parser.parse_args().t_norm

        self.show = parser.parse_args().show_plot
        self.save = parser.parse_args().save_plot
        self.legend = parser.parse_args().legend
        self.show_gas = parser.parse_args().show_gas
        self.cmc_method = parser.parse_args().cmc_method

    def readPVT(self):

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
        
        self.temp, self.t_max = self.getTempRange()

        self.lnZ = np.zeros( (self.t_max, len(self.temp)) )
        self.mu1 = np.zeros( (self.t_max, len(self.temp)) )
        self.N1 = np.zeros( (self.t_max, len(self.temp)) )
        params["avgU"] = np.zeros( (self.t_max, len(self.temp)) )

        if self.n_spec == 1:
            params["lnZliq"] = np.zeros( (self.t_max, len(self.temp)) )
            params["lnZgas"] = np.zeros( (self.t_max, len(self.temp)) )

        elif self.n_spec == 2:
            self.mu2 = np.zeros( (self.t_max, len(self.temp)) )
            self.N2 = np.zeros( (self.t_max, len(self.temp)) )

        t_count = 0
        i = 0
        T_cur = self.temp[t_count]
        for line in pvt_file:
            data = line.strip().split()
            T = float(data[0])
            if (T != T_cur):
                i = 0
                t_count += 1
                T_cur =  self.temp[t_count]

            if self.n_spec == 1:
                self.mu1[i, t_count] = float(data[1])
                self.N1[i, t_count] = float(data[2])
                params["avgU"][i, t_count] = float(data[3])
                self.lnZ[i, t_count] = float(data[4])
                params["lnZliq"][i, t_count] = float(data[5])
                params["lnZgas"][i, t_count] = float(data[6])

            elif self.n_spec == 2:
                self.mu1[i, t_count] = float(data[1])
                self.mu2[i, t_count] = float(data[2])
                params["avgU"][i, t_count] = float(data[3])
                self.N1[i, t_count] = float(data[4])
                self.N2[i, t_count] = float(data[5])
                self.lnZ[i, t_count] = float(data[6])

            i += 1

        # calculate conversions
        self.conversions(t_count)
    
        pvt_file.close()
        return params

    def conversions(self, t_count):
        if (self.calc_phi):
            self.x = self.N1 * self.mono_vol / self.vol
            if self.calc_pressure:
                self.pressure = self.lnZ * self.mono_vol / self.vol
        elif (self.calc_conc):
            self.x = self.N1 / (self.vol * self.N1a * self.A3_to_m3)
            if self.calc_pressure:
                self.pressure = self.lnZ * self.temp[t_count] / (self.vol * self.N1a * self.A3_to_m3 * self.atomic_to_kJmol * self.bar_to_kPa) #* self.gas_const # temp [=] kJ/mol
                #self.pressure = self.lnZ * self.kB * self.temp[t_count] * self.convert_kPa / self.vol
                if self.calc_t_norm:
                    self.pressure = self.pressure / self.temp[t_count]
        elif (self.calc_rho):
            self.x = self.N1 * self.mass / (self.vol * self.N1a * self.A3_to_cm3)
            if self.calc_pressure:
                self.pressure = self.lnZ * self.temp[t_count] * self.convert_kPa / self.vol
                if self.calc_t_norm:
                    self.pressure = self.pressure / self.temp[t_count]
        elif (self.calc_num_dens):
            self.x = self.N1 / self.vol
            if self.calc_pressure:
                self.pressure = self.lnZ * self.temp[t_count] / self.vol
                if self.calc_t_norm:
                    self.pressure = self.pressure / self.temp[t_count]
                #self.pressure = self.lnZ * self.kB * self.temp[t_count] * self.convert_kPa / self.vol
        else:
            self.x = self.N1

        if not self.calc_pressure:
            self.pressure = self.lnZ

    def calculateCMC(self):
        if (self.cmc_method == 'spline' or self.cmc_method == 'finite'):
            self.calcCMC2ndD()

        elif (self.cmc_method == 'zero'):
            self.calcCMC2ndDIntercept()
    
        elif (self.cmc_method == 'curvature'):
            self.calcCMCcurvature()
    
        elif (self.cmc_method == 'intercept'):
            self.calcCMCslope()

        elif (self.cmc_method == 'percent'):
            self.calcCMCslope()

        elif (self.cmc_method == 'sigmoid'):
            self.calcCMCsigmoidIntegral()

        else :
            print 'CMCs are not being calcualted use the -m flag.'
    
    def getTempRange(self):
        pvt_name = 'pvt.dat'
        # get the ln of partition function
        pvt_file = open(pvt_name,'r')
        # read header
        data = pvt_file.readline().strip().split()
        while ('T' not in data):
            data = pvt_file.readline().strip().split()
        
        line = pvt_file.readline()

        data = line.strip().split()
        try:
            temp = [float(data[0])]
        except IndexError:
            sys.exit('you forgot to generate pvt.dat')
        t_range = [1]
        t_count = 0

        for line in pvt_file:
            data = line.strip().split()
            T = float(data[0])
            if (T != temp[len(temp)-1]):
                temp.append(T)
                t_range.append(0)
                t_count += 1
            t_range[t_count] += 1

        pvt_file.close()

        return temp, max(t_range)
    
    # READ the cmc.dat file #
    def readCMC(self):
        pwd = os.getcwd()
        cmc_name='cmc.dat'
        # Check that file exists
        if not os.path.isfile(cmc_name):
            sys.exit("Error! file "+cmc_name+" does not exist.")
    
        params = {}
    
        # get info on the min slope from cmc
        cmc_file = open(cmc_name,'r')
        # read header
        data = cmc_file.readline().strip().split()
        while ('T' not in data):
            data = cmc_file.readline().strip().split()

        params["temp"] = []
        params["min_slope"] = []
        params["intercept"] = []
        params["cmc"] = []

        for line in cmc_file:
            data = line.strip().split()
            params["temp"].append(float(data[0]))
            params["min_slope"].append(float(data[1]))
            params["cmc"].append(float(data[3])*self.mono_vol/self.vol)
            params["intercept"].append(float(data[3])*self.mono_vol/self.vol * (1 - float(data[1])))
        cmc_file.close()
        return params
    
    def convert(self):
        # keep track of how many temepratures there are
        m_gas = 1.0
        self.itemp_skip = []
        
        slope_find_range = 4
        for itemp in range(len(self.temp)):
            if len(self.pressure) < slope_find_range:
                print "not enough data points for each temperature"
                self.itemp_skip.append( itemp )
                continue

            if (self.calc_conc):
                m_gas = self.temp[itemp] #* self.gas_const*1000
            elif (self.calc_rho):
                m_gas = self.temp[itemp] * self.kB
            elif (self.calc_num_dens):
                m_gas = self.temp[itemp] * self.kB

            if (self.calc_t_norm):
                m_gas = 1.0

            m_min = m_gas
            m_best = 0
            ave_x = sum(self.x[:,itemp]) / float(len(self.pressure))
            for i in range(int( len(self.pressure) - slope_find_range )):
                x = self.x[i:i+slope_find_range, itemp]
                y = self.pressure[i:i+slope_find_range, itemp]
                # y = mx + b
                m, b, m_s, b_s, r, err = linFit(x, y)
                if err: continue
                if (i == 0 ):
                    if (self.x[i,itemp] < ave_x):
                        b_shift = b
                if (i == (len(self.pressure) - slope_find_range - 1) ): 
                    if (self.x[i+slope_find_range,itemp] < ave_x):
                        b_shift = b

                if ( abs(m - m_gas) < abs(m_best - m_gas) ):
                    m_best = m
                    i_best = [i, i+slope_find_range]

            if (m_best == 0):
                self.itemp_skip.append( itemp )
                print "The maximum slope in your system is either "
                print "greater than "+str(m_gas)+" or less than 0."
                print "This does not make sense for calculating cmcs."
                b_shift = b
            else:
                # check how close the the pressure obeys the IGEOS
                print ('min IGEOS deviation(T =%6.2f) = %f' % (self.temp[itemp], m_best))
                print ('   over a range of points %f - %f' % (self.x[ i_best[0], itemp], self.x[ i_best[1], itemp]))
                
            for i in range( len(self.pressure[:,itemp]) ):
                self.pressure[i, itemp] -= b_shift

    def calcCMCslope(self):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu1_cmc = []
        percentage = 0.1
        m_gas = 1.0
        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            if (self.calc_conc):
                m_gas = self.temp[itemp] #* self.gas_const
            elif (self.calc_rho):
                m_gas = self.temp[itemp] * self.kB
            m_min = m_gas
            b_min = 0
            m_max = 0
            for i in range(int( len(self.pressure) - 5 )):
                x = self.x[i:i+5, itemp]
                y = self.pressure[i:i+5, itemp]
                # y = mx + b
                m, b, m_s, b_s, r, err = linFit(x, y)
                if err: continue
                if (m < m_min):
                    m_min = m
                    b_min = b
                    m_s_min = m_s
                    b_s_min = b_s
    
                elif (1.1*m_gas >= m > m_max):
                    m_max = m
                    b_max = b
                    m_s_max = m_s
                    b_s_max = b_s

            if (self.cmc_method == 'percent'):
                x = self.x[6:len(self.x)-3, itemp]
                y = self.pressure[6:len(self.x)-3, itemp]
                x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y, 1000)
                for i in range( len(x_spl) ):
                    if ( m_min >= (y_spl_d1[i] * percentage) ):
                        x_per = x_spl[i]
                        break

                self.cmc.append( x_per )
                self.cmc_s.append( x_spl[2] - x_spl[0] )
                
            elif (self.cmc_method == 'intercept'):
                if m_min == m_gas:
                    self.itemp_skip.append(itemp)
                    print 'could not find the intercept'
                    return

                self.cmc.append( abs( b_min / (m_gas - m_min)) )
                self.cmc_s.append( self.cmc[itemp] * 
                                  ( (b_s_min * 100)**2.0 + (m_s_min * 100)**2.0)**0.5 )
                self.cmc_slope.append( m_min )
                self.cmc_intercept.append( b_min )

            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu1_cmc.append( self.mu1[cmc_index,itemp] )

        return self.cmc, self.cmc_s

    def calcCMC2ndDIntercept(self):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu1_cmc = []
        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            n_points = len(self.x)
            x = self.x[self.begin_fit:n_points-self.end_fit, itemp]
            y = self.pressure[self.begin_fit:n_points-self.end_fit, itemp]
            i_cmc, i_cmc_e = zero2ndDerivativeIntercept(x, y, self.cmc_method, False)
            self.cmc.append(i_cmc)
            self.cmc_s.append(i_cmc_e)
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu1_cmc.append( self.mu1[cmc_index,itemp] )
        
        return self.cmc, self.cmc_s
    
    def calcCMCsigmoidIntegral(self, n_pass=2):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu1_cmc = []
        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            n_points = len(self.x)
            x = self.x[self.begin_fit:n_points-self.end_fit, itemp]
            y = self.pressure[self.begin_fit:n_points-self.end_fit, itemp]
            x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y)
            funky=True
            f_int = 6
            while (funky):
                try:
                    popt, pcov = curve_fit(sigmoidIntegral, x_spl, y_spl)
                    funky=False

                except RuntimeError:
                    print 'funky data, shortening range'
                    return
                    x = self.x[f_int:n_points-f_int, itemp]
                    y = self.pressure[f_int:n_points-f_int, itemp]
                    f_int -= 3

            #x0_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - popt[2]))
            #x = self.x[x0_index+10:x0_index-10, itemp]
            #y = self.pressure[x0_index+10:x0_index-10, itemp]
            #x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y)
            #popt, pcov = curve_fit(sigmoidIntegral, x_spl, y_spl)
                
            y_sig = sigmoidIntegral(x_spl, popt[0], popt[1], popt[2], popt[3])
            #for i in range(len(y_sig)):
            #    print x_spl[i], y_sig[i]
            
            self.cmc.append(popt[2])
            self.cmc_s.append(np.sqrt(np.diag(pcov))[2])
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu1_cmc.append( self.mu1[cmc_index,itemp] )
            if (False):
                fig = plt.figure()
                x = self.x[self.begin_fit:n_points-self.end_fit, itemp]
                y = self.pressure[self.begin_fit:n_points-self.end_fit, itemp]
                plt.plot(x_spl, y_sig, 'g', label='Sigmoid fit')
                plt.plot(x_spl, y_spl, 'r', label='Spline fit')
                plt.plot(x, y, 'og', label='Histogram rewieghting data')
                border = plt.legend( numpoints=1, prop={'size':12}, loc=4)
                border.draw_frame(False)
                if (self.calc_phi):
                    plt.xlabel("$\phi_{tot}$",fontsize=20)
                    if self.calc_t_norm:
                        plt.ylabel("$\Pi$",fontsize=20)
                    else:
                        plt.ylabel("$\P$",fontsize=20)
                elif (self.calc_conc):
                    plt.xlabel("$\\rho_{tot}$[mM]",fontsize=20)
                    if self.calc_t_norm:
                        plt.ylabel("$P / kT$",fontsize=20)
                    else:
                        plt.ylabel("$P$[kPa]",fontsize=20)
                elif (self.calc_rho):
                    plt.xlabel("$\\rho_{tot}$[g/cm$^3$]",fontsize=20)
                    if self.calc_t_norm:
                        plt.ylabel("$p / kT$",fontsize=20)
                    else:
                        plt.ylabel("$p$[kPa]",fontsize=20)
                elif self.calc_num_dens:
                    plt.xlabel("$\\rho_{tot}$",fontsize=20)
                    if self.calc_t_norm:
                        plt.ylabel("$p / kT$",fontsize=20)
                    else:
                        plt.ylabel("$p$[kPa]",fontsize=20)
                else:
                    plt.xlabel("N$_{tot}$",fontsize=20)
                    plt.ylabel("ln$\Omega$",fontsize=20)
        
                plt.show()

        return self.cmc, self.cmc_s

    def calcCMCcurvature(self):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu1_cmc = []
        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            n_points = len(self.x)
            i_cmc, i_cmc_e, x_spl, p_spl, kappa = maxCurvature(self.x[self.begin_fit:n_points-self.end_fit, itemp], self.pressure[self.begin_fit:n_points-self.end_fit, itemp], False)
            self.cmc.append(i_cmc)
            self.cmc_s.append(i_cmc_e)
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu1_cmc.append( self.mu1[cmc_index,itemp] )
            curavture_file = open("Pcurvature" + str(self.temp[itemp]) + ".dat", 'w')
            curavture_file.write("# rho pressure curvature\n")
            for i in range( len(kappa)-1):
                curavture_file.write("%f %f %f\n" % (x_spl[i], p_spl[i], kappa[i]) )
            curavture_file.close()
        
            
        return self.cmc, self.cmc_s
    
    def calcCMC2ndD(self):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu1_cmc = []
        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            n_points = len(self.x)
            i_cmc, i_p, i_d2, i_cmc_e = max2ndDerivative(self.x[self.begin_fit:n_points-self.end_fit, itemp], self.pressure[self.begin_fit:n_points-self.end_fit, itemp], self.cmc_method, False)
            self.cmc.append(i_cmc)
            self.cmc_s.append(i_cmc_e)
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu1_cmc.append( self.mu1[cmc_index,itemp] )
        
        return self.cmc, self.cmc_s
    
    def writeCMC(self):
        cmc_file = open("cmc_" + self.cmc_method + ".dat", 'w')
        if (self.calc_phi):
            cmc_file.write("T   phi_cmc d_phi_cmc    mu_cmc\n")
        elif (self.calc_rho or self.calc_conc or self.calc_num_dens):
            cmc_file.write("T   rho_cmc d_rho_cmc    mu_cmc\n")
        else:
            cmc_file.write("T   N_cmc d_N_cmc    mu_cmc\n")

        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            cmc_file.write("%f %10.8f %10.8f %10.8f\n" % ( self.temp[itemp], 
                                                    self.cmc[itemp], self.cmc_s[itemp],
                                                    self.mu1_cmc[itemp] ))

    def writePN(self):
        for itemp in range(len(self.temp)):
            if (itemp in self.itemp_skip): continue

            mu_file = open("PVphiVmu" + str(self.temp[itemp]) + ".dat", 'w')
            if (self.calc_phi):
                mu_file.write("#phi     P/kT       mu\n")
            elif self.calc_rho:
                if self.calc_t_norm:
                    mu_file.write("#rho[g/cm3] P/kT       mu\n")
                else:
                    mu_file.write("#rho[g/cm3] P[]        mu\n")
            elif self.calc_conc:
                if self.calc_t_norm:
                    mu_file.write("#rho[mM] P/kT       mu\n")
                else:
                    mu_file.write("#rho[mM] P[bar]        mu\n")
            elif self.calc_num_dens:
                if self.calc_t_norm:
                    mu_file.write("#rho     P/kT       mu\n")
                else:
                    mu_file.write("#rho     P[atomistic] mu\n")
            else:
                mu_file.write("#N     lnZ          mu\n")

            for i in range( len(self.pressure[:,itemp]) ):
                mu_file.write("%10.8f %10.8f %f\n" % (self.x[i, itemp], 
                                                  self.pressure[i, itemp], 
                                                  self.mu1[i, itemp]) )
            mu_file.close()
    
    def plot(self):
        fig = plt.figure()
        handles = []
        for itemp in range(len(self.temp)):
            x = self.x[:, itemp]
            y = self.pressure[:, itemp]
            plt.plot(x, y, '-', c=self.colors[itemp], label='T = '+str(self.temp[itemp]))
        
            if (itemp in self.itemp_skip): continue
           
            if self.calc_t_norm:
                convert = 1.0
            elif (self.calc_rho):
                convert = self.temp[itemp] * self.kB 
            elif (self.calc_conc):
                convert = self.temp[itemp] #* self.gas_const
            else:
                convert = 1.0

            if (self.cmc_method == 'intercept'): 
                x = np.array([0.0, self.cmc[itemp]])
                plt.plot([self.cmc[itemp], self.cmc[itemp]], x * convert, '--k')
                x_inter = np.array([self.cmc[itemp], max(self.x[:, itemp])])
                plt.plot(x_inter, self.cmc_slope[itemp]*x_inter + self.cmc_intercept[itemp], 'k')
                plt.plot(x, x*convert, '-k', lw=2.1)
            elif (self.cmc_method != None):
                if (self.show_gas):
                    x = np.array([0.0, self.cmc[itemp]])
                    plt.plot(x, convert * x, '-k', lw=2.1)
                    plt.plot([self.cmc[itemp], self.cmc[itemp]], [0, max(y)], '--k')
                    #plt.plot([self.cmc[itemp], self.cmc[itemp]], convert * x, '--k')
                else:
                    x = np.array([0.0, self.cmc[itemp]*1.5])
                    plt.plot([self.cmc[itemp], self.cmc[itemp]], [0, max(y)], '--k')
                    #plt.plot([self.cmc[itemp], self.cmc[itemp]], convert * x, '--k')
                #x = np.array([self.cmc[itemp], max(self.x[:, itemp])])
                #plt.plot(x, self.cmc*np.ones((len(x),1)), 'k')

    
        if (self.legend):
            border = plt.legend( numpoints=1, prop={'size':12}, loc=2)
            border.draw_frame(False)
        if (self.calc_phi):
            plt.xlabel("$\phi_{tot}$",fontsize=20)
            if self.calc_t_norm:
                plt.ylabel("$\Pi$",fontsize=20)
            else:
                plt.ylabel("$\P$",fontsize=20)
        elif (self.calc_conc):
            plt.xlabel("$\\rho_{tot}$[mM]",fontsize=20)
            if self.calc_t_norm:
                plt.ylabel("$P / kT$",fontsize=20)
            else:
                plt.ylabel("$P$[kPa]",fontsize=20)
        elif (self.calc_rho):
            plt.xlabel("$\\rho_{tot}$[g/cm$^3$]",fontsize=20)
            if self.calc_t_norm:
                plt.ylabel("$p / kT$",fontsize=20)
            else:
                plt.ylabel("$p$[kPa]",fontsize=20)
        elif self.calc_num_dens:
            plt.xlabel("$\\rho_{tot}$",fontsize=20)
            if self.calc_t_norm:
                plt.ylabel("$p / kT$",fontsize=20)
            else:
                plt.ylabel("$p$[kPa]",fontsize=20)
        else:
            plt.xlabel("N$_{tot}$",fontsize=20)
            plt.ylabel("ln$\Omega$",fontsize=20)
        
        plt.gcf().subplots_adjust(bottom=0.12)
        plt.gcf().subplots_adjust(left=0.11)
        plt.gcf().subplots_adjust(right=0.96)
        plt.gcf().subplots_adjust(top=0.96)
        if (self.show):
            plt.show()

        if (self.save):
            fig.savefig('pvt.png', format='png',dpi=100)
        return


def main(argv=None):
    parser = argparse.ArgumentParser(description='Calculate the pressure from'
                                    ' the resulting "pvt.dat" file output from'
                                    '  entropy or entropy2.')
    parser.add_argument("--n_components", type=int, choices=[1,2],
                   help='number of components, can either be 1 or 2')
    parser.add_argument("-p","--save_plot", action="store_true",
                   help='Save a figure ploting the pressure versus x')
    parser.add_argument("-s","--show_plot", action="store_true",
                   help='Show the plot.')
    parser.add_argument("-g","--show_gas", action="store_true",
                   help='Show the ideal-gas law limiting curve')
    parser.add_argument("-l","--legend", action="store_true",
                   help='Show the legend of temperature curves')
    parser.add_argument("--partition", action="store_true",
                   help='Plot the partition function, not the pressure')
    parser.add_argument("--pressure", action="store_true",
                   help='Plot the pressure in units corresponding to the x-axis')
    parser.add_argument("--t_norm", action="store_true",
                   help='Normalize y-axis by the temperature')
    parser.add_argument("-c","--concentration", action="store_true",
                   help='Plot versus concentration [mM]')
    parser.add_argument("--phi", type=float,
                   help='Plot versus volume fraction, give the volume of 1 moiety')
    parser.add_argument("-r","--density", type=float,
                   help='Plot versus density [g/cm3] must give the molar mass')
    parser.add_argument("--num_density", action="store_true",
                   help='Plot versus the number density [1/V]')
    parser.add_argument("-o","--output", action="store_true",
                   help='Output files for the cmc and P vs N vs mu')
    parser.add_argument("-m","--cmc_method", type=str, 
                        choices=['finite', 'zero', 'spline', 'intercept',
                                 'curvature', 'sigmoid', 'percent'],
                   help='Method for calculating the CMC')
    parser.add_argument("-b","--box_length", type=float, nargs="+",
                   help='Simulation box length, can be 1 or 3 numbers')

    pressure = partition2pressure()

    pressure.addParser(parser)
    pressure.readPVT()
    pressure.convert()
    if (parser.parse_args().output):
        pressure.writePN()

    if (parser.parse_args().cmc_method):
        pressure.calculateCMC()
        if (parser.parse_args().output):
            pressure.writeCMC()

    if (parser.parse_args().save_plot or parser.parse_args().show_plot):
        pressure.plot()

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
