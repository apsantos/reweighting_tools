#!/usr/bin/env python
import sys, getopt
import math as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# from scipy.sparse import coo_matrix 
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1.2
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1.2

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

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

class cassandraFile(object):
    """
    Read Cassandra *box.prp file
    """
    def __init__(self, file_root):
        self.runname = file_root[0]
        return

    def setWidth(self, width):
        self.width = float(width)

    def read(self, start_line=0):
        self.readInp()
        self.readPrp(start_line)

    def getNspecies(self):
        return self.nspecies

    def getMu(self):
        return self.mu

    def getTemp(self):
        return self.temperature * 0.0083144621

    def getBox(self):
        return self.box

    def readInp(self):
        # open the file
        filename = self.runname + '.inp'

        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            print('cannot find: %s' % filename)
            return 1
        
        # read header information
        while 1:
            line = ifile.readline()
            if line.strip():
                data = line.strip().split()
                if (data[0] == '#'):
                    if (data[1] == 'Nbr_Species'):
                        self.nspecies = int( ifile.readline().strip().split()[0] )
                    elif (data[1] == 'Chemical_Potential_Info'):
                        self.mu = []
                        mu_str = ifile.readline().strip().split()
                        for i in range(self.nspecies):
                            self.mu.append( float( mu_str[i] ) )
                    elif (data[1] == 'Temperature_Info'):
                        self.temperature = float( ifile.readline().strip().split()[0] )
                    elif (data[1] == 'Box_Info'):
                        ifile.readline()
                        ifile.readline()
                        box_str = ifile.readline().strip().split()
                        self.box = []
                        self.box.append(float( box_str[0] ) )
                        self.box.append(float( box_str[1] ) )
                        self.box.append(float( box_str[2] ) )
            if not line: break
        ifile.close()
                
    def readPrp(self, start_line=0):
        # open the file
        filename = self.runname + '.box1.prp1'

        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            print('cannot find: %s' % filename)
            return 1
        
        # read header information
        ifile.readline()
        header = ifile.readline().strip().split()
        n_collumn = []
        e_collumn = 0
        i_collumn = -1
        for property in header:
            if (property == 'Nmols'):
                n_collumn.append(i_collumn)
                
            elif (property == 'Energy_Total'):
                e_collumn = i_collumn
            i_collumn += 1

        if (len(n_collumn) != self.nspecies): 
            print('Need to have a Nmols collumn in the %s for all (%d) species' %
                  (filename, self.nspecies) )
            return 1

        elif (e_collumn == 0):
            print('Need to have Energy_Total in the %s' % filename)
            return 1

        self.nmols = []
        self.energy = []
        i_line = 0
        for line in ifile:
            data = line.strip().split()
            if (i_line >= start_line):

                if (data[0] != "#"):
                    self.nmols.append( int(float(data[n_collumn[0]])) )
                    #for ispecies in n_collumn[1:]:
                    #    self.nmols[len(self.nmols) - 1] += int(float(data[ispecies]))

                    self.energy.append( float(data[e_collumn]) )

            i_line += 1

    def generateHis(self):
        """ There are 2 types of data lines, the 1st has the:
            N (number of particles)
            n_E_bin (number of energy bins)
            E_min (the min value of energy for non-zero observations)
            
            The next line(s) have the tally of occurances in each bin starting from
        """
        n_min = min(self.nmols)
        n_max = max(self.nmols)
        e_min = min(self.energy)
        e_max = max(self.energy)
        n_e_bins = int(abs(e_max - e_min)/self.width) + 1
        his = np.zeros(((n_max - n_min)+1, n_e_bins), np.int)
        n_his = np.zeros((n_max - n_min)+1, np.int)
        n_his[:] = -1
        for i_step in range(len(self.nmols)): 
            n_bin = self.nmols[i_step] - n_min
            e_bin = int( abs(self.energy[i_step] - e_min) / self.width)
            his[n_bin, e_bin] += 1
            if (self.nmols[i_step] not in n_his):
                n_his[n_bin] = self.nmols[i_step]

        n_e = []
        e_min_his = []
        n_pop = 0
        for n_bin in range(len(n_his)): 
            # compact array for empty spaces in N
            if (n_his[n_bin-n_pop] == -1):
                n_his = np.delete(n_his, n_bin-n_pop)
                his = np.delete(his, n_bin-n_pop, 0)
                n_pop += 1
            # calculate the minimum energy for each N and number of energy bins
            else:
                nz_index = np.nonzero( his[n_bin-n_pop,:] )
                e_min_his.append( (nz_index[0][0] * self.width) + e_min)
                n_e.append(nz_index[0][len(nz_index[0])-1] - nz_index[0][0] + 1)

        self.his = his
        self.e_min_his = e_min_his
        self.n_his = n_his
        self.n_e = n_e

    def getHistogram(self):
        return self.his

    def getNbins(self):
        return self.n_his

    def getNe(self):
        return self.n_e

    def getNeMin(self):
        return self.e_min_his

class hisFile(object):
    """
    Read and plot histogram(s) as a COO sparse matrix
    """
    def __init__(self, file_roots):
        # can only take 7
        self.runs = file_roots
        self.temp_list = []
        return

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
            if (N_line):
                i_N = int(data[0])
                E_bin_min = float(data[2])
                n_E_bin = int(data[1])
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

        self.histogram = his
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

def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.getopt(argv[1:], "hr:i:s:plj:o:w:", 
                     ["help", "ifile=", "start=", "plot", "legend", "chain_length=", "write_out=", "width"])

    except getopt.error, msg:
        print msg
        print "for help use --help"
        return 2

    show_on = False
    legend_on = False
    write_on = False
    J = 8
    start_line = 0
    for opt, arg in opts:
        if opt == '-h':
            print("hisFile.py -r ''<file_root_name>'' -i '<file_w/_root_names>' --plot --legend -j <chain_length> -o <file_type> -w <width> -s <start_line>\n"
                  "EG plot_his.py -r ''01a','02b'' -s\n"
                  "\tOR\n"
                  "hisFile.py -i input_hs.dat")
            return 1

        elif opt == '-i':
            runs_filename = arg
            runs = readRunsFile(runs_filename)

        elif opt == '-r':
            runs = [x.strip() for x in  arg.split(',')]

        elif opt == '-s':
            start_line = int(arg)

        elif opt == '-p':
            show_on = True

        elif opt == '-l':
            legend_on = True

        elif opt == '-j':
            J = int(arg)

        elif opt == '-w':
            width = float(arg)

        elif opt == '-o':
            write_on = True
            file_type = arg

    HIS = hisFile(runs)
    HIS.setJ(J)

    if (write_on):
        if (file_type == "Cassandra" or file_type == "cassandra"):
            inFile = cassandraFile(runs)
        else:
            print(file_type +' is not a file_type option')

        inFile.read(start_line)
        inFile.setWidth(width)
        inFile.generateHis()
        box = inFile.getBox()
        HIS.setBox(box[0], box[1], box[2])
        HIS.setMu(inFile.getMu() )
        HIS.setTemp(inFile.getTemp() )
        HIS.setHistogram(inFile.getHistogram() )
        HIS.setNe(inFile.getNe() )
        HIS.setNeMin(inFile.getNeMin() )
        HIS.setNbins(inFile.getNbins() )
        HIS.setWidth(width)
        HIS.write()
    else:
        HIS.plot(show_on, legend_on)

if __name__ == '__main__':
    sys.exit(main())

# "The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
