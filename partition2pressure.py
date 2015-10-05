#!/usr/bin/env python
"""
Author : Andrew P. Santos
Date   : May, 2014
Institution: Princeton University

Purpose: Calculate the pressure from the resulting "pvt.dat" file output from entropy or entropy2 

	The top two def sections are called from the main section which is below
"""
import os, sys, getopt
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
        sys.exit(main())

    tck = interpolate.splrep(x, y, s=0)
    x_spl = np.linspace(0, max(x), n_points)
    return x_spl, interpolate.splev(x_spl,tck, der=0), interpolate.splev(x_spl,tck, der=1), interpolate.splev(x_spl,tck, der=2)

def maxCurvature(x, y, plot=False):
    kappa_max = 0.0
    x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y)
    dx = x_spl[2] - x_spl[0]
    dx2 = dx * dx
    kappa_max_e = dx2
    x_max_e = dx2
    for i in range(1,len(x_spl)-1):
        kappa = -(y_spl_d2[i]) / (1 + (y_spl_d1[i])**2.0)**(1.5)
        if (kappa > kappa_max):
            kappa_max = kappa
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

    return x_max, x_max_e

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
    return slope, intercept, slope_std, intercept_std, r

class partition2pressure(object):
    """
    Read and process PVT and cmc files from the entropy program suite
    """
    def __init__(self, parms):
        self.L = parms["L"]
        self.vol = float(self.L * self.L * self.L)
        self.calc_conc = False
        if ("conc" in parms):
            self.calc_conc = True
            self.kB = 1.38064852E-23 # J K-1
            self.convert_kPa = 1.E27 # J / A^3 to kPa
            self.A3_to_m3 = 1.E-30 # A^3 to m^3
            self.Na = 6.022140857E23 # 1/ mol
            self.gas_const = 0.0083144621 # kJ / (K mol)
        self.calc_rho = False
        if ("mass" in parms):
            self.calc_rho = True
            self.mass = parms["mass"]
            self.kB = 1.38064852E-23 # J K-1
            self.convert_kPa = 1.E27 # J / A^3 to kPa
            self.A3_to_cm3 = 1.E-24 # A^3 to m^3
            self.Na = 6.022140857E23 # 1/ mol
        self.calc_phi = False
        if ("J" in parms):
            self.calc_phi = True
            self.J = parms["J"]

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

        self.mu = np.zeros( (self.t_max, len(self.temp)) )
        self.N = np.zeros( (self.t_max, len(self.temp)) )
        params["avgU"] = np.zeros( (self.t_max, len(self.temp)) )
        lnZ = np.zeros( (self.t_max, len(self.temp)) )
        params["lnZliq"] = np.zeros( (self.t_max, len(self.temp)) )
        params["lnZgas"] = np.zeros( (self.t_max, len(self.temp)) )

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
            self.mu[i, t_count] = float(data[1])
            self.N[i, t_count] = float(data[2])
            params["avgU"][i, t_count] = float(data[3])
            lnZ[i, t_count] = float(data[4])
            params["lnZliq"][i, t_count] = float(data[5])
            params["lnZgas"][i, t_count] = float(data[6])
            i += 1

        if (self.calc_phi):
            self.x = self.N * self.J / self.vol
            self.pressure = lnZ * self.J / self.vol
        elif (self.calc_conc):
            self.x = self.N / (self.vol * self.Na * self.A3_to_m3)
            self.pressure = lnZ * self.temp[t_count] * self.convert_kPa / self.vol /self.Na * 1000.0
            #self.pressure = lnZ * self.kB * self.temp[t_count] * self.convert_kPa / self.vol
        elif (self.calc_rho):
            self.x = self.N * self.mass / (self.vol * self.Na * self.A3_to_cm3)
            self.pressure = lnZ * self.temp[t_count] * self.convert_kPa / self.vol
            #self.pressure = lnZ * self.kB * self.temp[t_count] * self.convert_kPa / self.vol
        else:
            self.x = self.N
            self.pressure = lnZ
    
        pvt_file.close()
        return params

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
        temp = [float(data[0])]
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
            params["cmc"].append(float(data[3])*self.J/self.vol)
            params["intercept"].append(float(data[3])*self.J/self.vol * (1 - float(data[1])))
        cmc_file.close()
        return params
    
    def convert(self):
        # keep track of how many temepratures there are
        m_gas = 1.0
        
        for itemp in range(len(self.temp)):
            if (self.calc_conc):
                m_gas = self.temp[itemp] #* self.gas_const
            elif (self.calc_rho):
                m_gas = self.temp[itemp] * self.kB
            m_min = m_gas
            m_max = 0
            for i in range(int( len(self.pressure) - 5 )):
                x = self.x[i:i+5, itemp]
                y = self.pressure[i:i+5, itemp]
                # y = mx + b
                m, b, m_s, b_s, r = linFit(x, y)
                if (m_gas*1.0 >= m > m_max):
                #if (m > m_max):
                    m_max = m
                    b_max = b
            print m_max
            if (m_max == 0):
                sys.exit("The maximum slope in your system is either "
                         "greater than "+str(m_gas)+" or less than 0.\n"
                         "This does not make sense for calculating cmcs.")
                
            # make sure that the pressure obeys the IGEOS at the limit
            for i in range( len(self.pressure[:,itemp]) ):
                self.pressure[i, itemp] -= b_max

    def calcCMCslope(self, cmcParams, method='intersept'):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu_cmc = []
        percentage = 0.1
        m_gas = 1.0
        for itemp in range(len(self.temp)):
            if (self.calc_conc):
                m_gas = self.temp[itemp] #* self.gas_const
            elif (self.calc_rho):
                m_gas = self.temp[itemp] * self.kB
            m_min = m_gas
            m_max = 0
            for i in range(int( len(self.pressure) - 5 )):
                x = self.x[i:i+5, itemp]
                y = self.pressure[i:i+5, itemp]
                # y = mx + b
                m, b, m_s, b_s, r = linFit(x, y)
                if (m < m_min):
                    m_min = m
                    b_min = b
                    m_s_min = m_s
                    b_s_min = b_s
    
                elif (1.0 >= m > m_max):
                    m_max = m
                    b_max = b
                    m_s_max = m_s
                    b_s_max = b_s

            if (method == 'percent'):
                x = self.x[6:len(self.x)-3, itemp]
                y = self.pressure[6:len(self.x)-3, itemp]
                x_spl, y_spl, y_spl_d1, y_spl_d2 = getSpline(x, y, 1000)
                for i in range( len(x_spl) ):
                    if ( m_min >= (y_spl_d1[i] * percentage) ):
                        x_per = x_spl[i]
                        break

                self.cmc.append( x_per )
                self.cmc_s.append( x_spl[2] - x_spl[0] )
                
            elif (method == 'intercept'):
                print b_min, m_min
                self.cmc.append( abs( b_min / (m_gas - m_min)) )
                self.cmc_s.append( self.cmc[itemp] * 
                                  ( (b_s_min * 100)**2.0 + (m_s_min * 100)**2.0)**0.5 )
                self.cmc_slope.append( m_min )
                self.cmc_intercept.append( b_min )

            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu_cmc.append( self.mu[cmc_index,itemp] )

        return self.cmc, self.cmc_s

    def calcCMC2ndDIntercept(self, cmcParams, method):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu_cmc = []
        for itemp in range(len(self.temp)):
            x = self.x[6:len(self.x)-3, itemp]
            y = self.pressure[6:len(self.x)-3, itemp]
            i_cmc, i_cmc_e = zero2ndDerivativeIntercept(x, y, method, False)
            self.cmc.append(i_cmc)
            self.cmc_s.append(i_cmc_e)
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu_cmc.append( self.mu[cmc_index,itemp] )
        
        return self.cmc, self.cmc_s
    
    def calcCMCsigmoidIntegral(self, cmcParams, n_pass=2):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu_cmc = []
        for itemp in range(len(self.temp)):
            print self.temp[itemp]
            n_points = len(self.x)
            x = self.x[6:n_points-3, itemp]
            y = self.pressure[6:n_points-3, itemp]
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
            for i in range(len(y_sig)):
                print x_spl[i], y_sig[i]
            
            self.cmc.append(popt[2])
            self.cmc_s.append(np.sqrt(np.diag(pcov))[2])
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu_cmc.append( self.mu[cmc_index,itemp] )
            if (False):
                fig = plt.figure()
                x = self.x[6:n_points-3, itemp]
                y = self.pressure[6:n_points-3, itemp]
                plt.plot(x_spl, y_sig, 'g', label='Sigmoid fit')
                plt.plot(x_spl, y_spl, 'r', label='Spline fit')
                plt.plot(x, y, 'og', label='Histogram rewieghting data')
                border = plt.legend( numpoints=1, prop={'size':12}, loc=4)
                border.draw_frame(False)
                if (self.calc_phi):
                    plt.xlabel("$\phi_{tot}$",fontsize=20)
                    plt.ylabel("$\Pi$",fontsize=20)
                elif (self.calc_conc):
                    plt.xlabel("$\\rho_{tot}$[mM]",fontsize=20)
                    plt.ylabel("$p$[kPa]",fontsize=20)
                elif (self.calc_rho):
                    plt.xlabel("$\\rho_{tot}$[g/cm$^3$]",fontsize=20)
                    plt.ylabel("$p$[kPa]",fontsize=20)
                else:
                    plt.xlabel("N$_{tot}$",fontsize=20)
                    plt.ylabel("ln$\Omega$",fontsize=20)
        
                plt.show()

        return self.cmc, self.cmc_s

    def calcCMCcurvature(self, cmcParams):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu_cmc = []
        for itemp in range(len(self.temp)):
            i_cmc, i_cmc_e = maxCurvature(self.x[6:len(self.x)-3, itemp], self.pressure[6:len(self.x)-3, itemp], False)
            self.cmc.append(i_cmc)
            self.cmc_s.append(i_cmc_e)
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu_cmc.append( self.mu[cmc_index,itemp] )
        
        return self.cmc, self.cmc_s
    
    def calcCMC2ndD(self, cmcParams, method):
        self.cmc = []
        self.cmc_s = []
        self.cmc_slope = []
        self.cmc_intercept = []
        self.mu_cmc = []
        for itemp in range(len(self.temp)):
            i_cmc, i_p, i_d2, i_cmc_e = max2ndDerivative(self.x[6:len(self.x)-3, itemp], self.pressure[6:len(self.x)-3, itemp], method, False)
            self.cmc.append(i_cmc)
            self.cmc_s.append(i_cmc_e)
            cmc_index = min(range(len(self.x[:,itemp])), key=lambda i: abs(self.x[i,itemp] - self.cmc[itemp]))
            self.mu_cmc.append( self.mu[cmc_index,itemp] )
        
        return self.cmc, self.cmc_s
    
    def write(self, parms):
        cmc_file = open("cmc_APS.dat", 'w')
        if (self.calc_phi):
            cmc_file.write("T   phi_cmc d_phi_cmc    mu_cmc\n")
        elif (self.calc_rho or self.calc_conc):
            cmc_file.write("T   rho_cmc d_rho_cmc    mu_cmc\n")
        else:
            cmc_file.write("T   N_cmc d_N_cmc    mu_cmc\n")
        for itemp in range(len(self.temp)):
            cmc_file.write("%f %10.8f %10.8f %10.8f\n" % ( self.temp[itemp], 
                                                    self.cmc[itemp], self.cmc_s[itemp],
                                                    self.mu_cmc[itemp] ))
            mu_file = open("PVphiVmu" + str(self.temp[itemp]) + ".dat", 'w')
            if (self.calc_phi):
                mu_file.write("phi     Pressure   mu\n")
            elif (self.calc_rho or self.calc_conc):
                mu_file.write("rho     Pressure   mu\n")
            else:
                mu_file.write("N     Pressure   mu\n")
            for i in range( len(self.pressure[:,itemp]) ):
                mu_file.write("%10.8f %10.8f %f\n" % (self.x[i, itemp], 
                                                  self.pressure[i, itemp], 
                                                  self.mu[i, itemp]) )
            mu_file.close()
        return
    
    def plot(self, pvtParams, show, method):
        fig = plt.figure()
        colors = ['red', 'blue', 'green', 'magenta', 'cyan', 'yellow', 'black', 'darkgoldenrod','firebrick', 'purple', 'burlywood', 'chartreuse', 'red', 'blue', 'green', 'magenta', 'cyan']
        handles = []
        for itemp in range(len(self.temp)):
            x = self.x[:, itemp]
            y = self.pressure[:, itemp]
            plt.plot(x, y, '-', c=colors[itemp], label='T = '+str(self.temp[itemp]))
        
           
            if (method == 'intercept'): 
                x = np.array([0.0, self.cmc[itemp]])
                plt.plot([self.cmc[itemp], self.cmc[itemp]], [0, self.cmc[itemp]], '--k')
                x_inter = np.array([self.cmc[itemp], max(self.x[:, itemp])])
                plt.plot(x_inter, self.cmc_slope[itemp]*x_inter + self.cmc_intercept[itemp], 'k')
                plt.plot(x, x, '-k', lw=2.1)
            elif (method != None):
                if (self.calc_rho):
                    convert = self.temp[itemp] * self.kB 
                elif (self.calc_conc):
                    convert = self.temp[itemp] #* self.gas_const
                else:
                    convert = 1.0
                x = np.array([0.0, self.cmc[itemp]])
                plt.plot(x, convert * x, '-k', lw=2.1)
                plt.plot([self.cmc[itemp], self.cmc[itemp]], convert * x, '--k')
                #x = np.array([self.cmc[itemp], max(self.x[:, itemp])])
                #plt.plot(x, self.cmc*np.ones((len(x),1)), 'k')

    
        border = plt.legend( numpoints=1, prop={'size':12}, loc=2)
        border.draw_frame(False)
        if (self.calc_phi):
            plt.xlabel("$\phi_{tot}$",fontsize=20)
            plt.ylabel("$\Pi$",fontsize=20)
        elif (self.calc_conc):
            plt.xlabel("$\\rho_{tot}$[mM]",fontsize=20)
            plt.ylabel("$p$[kPa]",fontsize=20)
        elif (self.calc_rho):
            plt.xlabel("$\\rho_{tot}$[g/cm$^3$]",fontsize=20)
            plt.ylabel("$p$[kPa]",fontsize=20)
        else:
            plt.xlabel("N$_{tot}$",fontsize=20)
            plt.ylabel("ln$\Omega$",fontsize=20)
        
        if (show):
            plt.show()

        fig.savefig('pvt.png', format='png',dpi=100)
        return


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.getopt(argv[1:], "hl:j:opsm:r:c",
                     ["help", "box_length=", "num_beads=", "output", "plot", "show plot",
                       "maxmethod=","mass=", "concentration"])

    except getopt.error, msg:
        print msg
        print "for help use --help"
        return 2

    inparams = {}
    output_on = False
    plot_on = False
    show_on = False
    cmc_method = None
    for opt, arg in opts:
        if opt == '-h':
            print "partition2pressure.py -l <box_length> -j <n_beads> [optional] -p -s -o -m <method> -r <molecular_mass> -c"
            print "methods: intercept, spline, finite, zero, curvature and sigmoid"
            return 1

        elif opt == '-l':
            inparams["L"] = int(arg)

        elif opt == '-c':
            inparams["conc"] = True

        elif opt == '-r':
            inparams["mass"] = float(arg)

        elif opt == '-j':
            inparams["J"] = int(arg)

        elif opt == '-o':
            output_on = True

        elif opt == '-p':
            plot_on = True

        elif opt == '-s':
            show_on = True

        elif opt == '-m':
            cmc_method = arg

    pressure = partition2pressure(inparams)

    PVTparams = pressure.readPVT()
    pressure.convert()

    #CMCparams = pressure.readCMC()
    CMCparams = {}
    if (cmc_method == 'spline' or cmc_method == 'finite'):
        CMCparams["cmc"], CMCparams["cmc_s"] = pressure.calcCMC2ndD(PVTparams, cmc_method)

    elif (cmc_method == 'zero'):
        CMCparams["cmc"], CMCparams["cmc_s"] = pressure.calcCMC2ndDIntercept(PVTparams, cmc_method)

    elif (cmc_method == 'curvature'):
        CMCparams["cmc"], CMCparams["cmc_s"] = pressure.calcCMCcurvature(PVTparams)

    elif (cmc_method == 'intercept'):
        CMCparams["cmc"], CMCparams["cmc_s"] = pressure.calcCMCslope(PVTparams, cmc_method)
    elif (cmc_method == 'percent'):
        CMCparams["cmc"], CMCparams["cmc_s"] = pressure.calcCMCslope(PVTparams, cmc_method)
    elif (cmc_method == 'sigmoid'):
        CMCparams["cmc"], CMCparams["cmc_s"] = pressure.calcCMCsigmoidIntegral(PVTparams)
    else :
        print 'CMCs are not being calcualted use the -m flag.'

    if (output_on):
        pressure.write(PVTparams)

    if (plot_on):
        pressure.plot(PVTparams, show_on, cmc_method)

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
