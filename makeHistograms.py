#!/usr/bin/env python
import sys, argparse
import math as ma
import numpy as np

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

class inputFile(object):
    """
    Read Cassandra *box.prp file
    Or a regular .dat file
    """
    def __init__(self, file_root):
        self.runname = file_root[0]
        return

    def setWidth(self, width):
        self.width = float(width)

    def read(self, start_line=0):
        if (self.data_type == "cassandra"):
            self.readInp()
            error = self.readPrp(start_line)
        else:
            error = self.readData(start_line)
        return error

    def getNspecies(self):
        return self.nspecies

    def getMu(self):
        return self.mu

    def getTemp(self):
        return self.temperature * 0.83144621

    def getBox(self):
        return self.box

    def readInp(self):
        # open the file
        filename = self.runname + '.log'

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
                            self.mu.append( float( mu_str[i] ) * 100.0 )
                    elif (data[1] == 'Temperature_Info'):
                        self.temperature = float( ifile.readline().strip().split()[0] )
                    elif (data[1] == 'Box_Info'):
                        ifile.readline()
                        ifile.readline()
                        box_str = ifile.readline().strip().split()
                        self.box = []
                        self.box.append(float( box_str[0] ) )
			if (len(box_str) == 1):
                            self.box.append(float( box_str[0] ) )
                            self.box.append(float( box_str[0] ) )
			else:
                            self.box.append(float( box_str[1] ) )
                            self.box.append(float( box_str[2] ) )
            if not line: break
        ifile.close()
                
    def readPrp(self, start_line=0, end_line=1000000):
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
            return 2

        self.nsteps = []
        self.nmols = []
        self.energy = []
        i_line = 0
        for line in ifile:
            data = line.strip().split()
            if i_line >= start_line and i_line < 500000:
                if (data[0] != "#"):
                    	self.nsteps.append( int(data[0]) )
                    	self.nmols.append( int(float(data[n_collumn[0]])) )
                    	#for ispecies in n_collumn[1:]:
                    	#    self.nmols[len(self.nmols) - 1] += int(float(data[ispecies]))

                    	self.energy.append( float(data[e_collumn]) * 100.0)

            if (i_line >= end_line): break

            i_line += 1

        if (i_line < (start_line - 10)):
            print ('starting line too high, there are only %s lines' % i_line)
            return 3
        return 0

    def readData(self, start_line=0, end_line=1000000):
        # open the file
        filename = self.runname + '.' + self.data_type

        try:
            ifile = open('./' + filename, 'r')

        except IOError:
            print('cannot find: %s' % filename)
            return 1
        
        # read header information
        ifile.readline()
        header = ifile.readline().strip().split()

        self.nsteps = []
        self.nmols = []
        self.energy = []
        i_line = 0
        for line in ifile:
            data = line.strip().split()
            if i_line >= start_line and i_line < 500000:
                if (data[0] != "#"):
                    	self.nmols.append( int(float(data[self.n_collumn])) )
                    	#for ispecies in n_collumn[1:]:
                    	#    self.nmols[len(self.nmols) - 1] += int(float(data[ispecies]))

                    	self.energy.append( float(data[self.e_collumn]) )

            if (i_line >= end_line): break

            i_line += 1

        if (i_line < (start_line - 10)):
            print ('starting line too high, there are only %s lines' % i_line)
            return 3
        return 0

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
    def __init__(self):
        # can only take 7
        self.temp_list = []
        self.colors = ['black', 'red', 'green', 'blue', 'magenta', 'cyan', 
                      'lawngreen', 'gold', 'mediumaquamarine', 'fuchsia',
                      'black', 'red', 'green', 'blue', 'magenta', 'cyan', 
                      'lawngreen', 'gold', 'mediumaquamarine', 'fuchsia',
                      'black', 'red', 'green', 'blue', 'magenta', 'cyan', 
                      'lawngreen', 'gold', 'mediumaquamarine', 'fuchsia']
        return

    def addParser(self, parser):
        self.J = 0
        if (parser.parse_args().num_beads):
            self.J = parser.parse_args().num_beads

        if (parser.parse_args().input_file):
            self.runs = readRunsFile(parser.parse_args().input_file)
        elif (parser.parse_args().runs_list):
            self.runs = parser.parse_args().runs_list
        else:
            print 'No input files given!'

        if (parser.parse_args().write_his):
            self.inFile = inputFile(self.runs)
            self.inFile.data_type = parser.parse_args().write_his
            if (parser.parse_args().write_his != 'cassandra'):
                if (not parser.parse_args().e_collumn):
                    print 'Need to specify:' 
                    print '  e_collumn'
                    return -1
                elif (not parser.parse_args().n_collumn):
                    print 'Need to specify:' 
                    print '  n_collumn'
                    return -1
                elif (not parser.parse_args().temperature):
                    print 'Need to specify:' 
                    print '  temperature'
                    return -1
                elif (not parser.parse_args().chem_pot):
                    print 'Need to specify:' 
                    print '  chem_pot'
                    return -1
                elif (not parser.parse_args().box):
                    print 'Need to specify:' 
                    print '  box'
                    return -1
                self.inFile.e_collumn = parser.parse_args().e_collumn
                self.inFile.n_collumn = parser.parse_args().n_collumn
                self.inFile.temperature = parser.parse_args().temperature /0.83144621
                self.inFile.mu = [parser.parse_args().chem_pot]
                self.inFile.box = [parser.parse_args().box, parser.parse_args().box, parser.parse_args().box]

            #self.runname = parser.parse_args().write_his
            self.inFile.setWidth(parser.parse_args().his_width)
            self.setWidth(parser.parse_args().his_width)

        if (parser.parse_args().start_line):
            self.start_line = parser.parse_args().start_line
        else:
            self.start_line = 0

        if (parser.parse_args().end_line):
            self.end_line = parser.parse_args().start_line
        else:
            self.end_line = 100000000000000

        if (parser.parse_args().show_legend):
            self.show_legend = parser.parse_args().show_legend
        else:
            self.show_legend = False

        if (parser.parse_args().show_color):
            self.show_color = parser.parse_args().show_color
        else:
            self.show_color = False

        return 0

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
        readErr = self.inFile.read(self.start_line)
        if (readErr > 0):
            return
        self.inFile.generateHis()
        box = self.inFile.getBox()
        self.setBox(box[0], box[1], box[2])
        self.setMu(self.inFile.getMu() )
        self.setTemp(self.inFile.getTemp() )
        self.setHistogram(self.inFile.getHistogram() )
        self.setNe(self.inFile.getNe() )
        self.setNeMin(self.inFile.getNeMin() )
        self.setNbins(self.inFile.getNbins() )

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
                    
    def reduceHis(self):
        for i in range( len(self.runs) ):
            read_err = self.read(self.runs[i])

def main(argv=None):
    parser = argparse.ArgumentParser(description='Program to create and do some '
                                    'processing of histogram files.')
    parser.add_argument("-i","--input_file", metavar='inFile', type=str,
                   help='File containing list of historgram file runnames '
                        '(e.g. input_hs.dat.')
    parser.add_argument("-r","--runs_list", type=str, nargs="+",
                   help='Instead of a file with input runs, '
                        'here give a list of the roots, e.g. - 1a 1b 1c')
    parser.add_argument("-o","--write_his", type=str,
                   help='Histogram file type (''cassandra'' for Cassandra output).')
    parser.add_argument("-w","--his_width", type=float,
                   help='Histogram energy width [kJ/mol]')
    parser.add_argument("-s","--start_line", type=int,
                   help='Starting line to operate on property file.')
    parser.add_argument("-f","--end_line", type=int,
                   help='Ending line to operate on property file.')
    parser.add_argument("--e_collumn", type=int,
                   help='The collumn number with the energy (starts with 0) (dimensionless).')
    parser.add_argument("--n_collumn", type=int,
                   help='The collumn number with the no. of molecules (starts with 0.')
    parser.add_argument("--chem_pot", type=float,
                   help='The chemical potential (dimensionless)')
    parser.add_argument("--temperature", type=float,
                   help='The temperature (dimensionless)')
    parser.add_argument("--box", type=float,
                   help='box length')

    HIS = hisFile()
    err = HIS.addParser(parser)
    if (err < 0): 
        print 'Exiting program'
        return

    if (parser.parse_args().write_his):
        HIS.write()

if __name__ == '__main__':
    sys.exit(main())

#The greatest happiness is to know the source of unhappiness." -Fyodor Dostoevsky
