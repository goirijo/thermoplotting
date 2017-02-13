import numpy as np


class EciInFiler(object):

    """Contains all the information from an eci.in file, with ways
    to fix a cluster's basis functions on or off as needed."""

    def __init__(self,filename):
        """Load all the information from a eci.in file"""
        
        self._source=filename

        with open(filename) as f:
            self._header=f.readline().split()

            self._lines=[]
            for line in f:
                self._lines.append(line.split())

    def fixon(self, indx):
        """Fix basis function off for the given index

        :indx: int
        :returns: EciInFiler

        """
        self._lines[indx][1]="FixOn"
        return self

    def fixoff(self, indx):
        """Fix basis function off for the given index

        :indx: int
        :returns: EciInFiler

        """
        self._lines[indx][1]="FixOff"

    def write(self, outfile):
        """Write the data back into a file, eci.in style.

        :outfile: string
        :returns: EciInFiler

        """
        ecistream=open(outfile,"w")
        template = "{:10}"*len(self._header)
        ecistream.write(template.format(*self._header))
        ecistream.write("\n")
        for line in self._lines:
            template = "{:10}"*len(line)
            ecistream.write(template.format(*line))
            ecistream.write("\n")
        return self

    def weight_array(self):
        """Return an ndarray of whether the eci were selected or not"
        :returns: ndarray

        """
        bits=[]
        for line in self._lines:
            bits.append(line[1])

        bits=np.array(bits,dtype=int)
        return bits.astype(bool)


class EciOutFiler(object):

    """Contains the values of an eci.out file"""

    def __init__(self,filename):
        """Just loads ECI, ECI/mult and cluster number into an ndarray

        :filename: path to eci.out

        """
        self._filename = filename
        self._values=np.loadtxt(filename,skiprows=7)

        return
        

class CorrFiler(object):

    """Contains all the data from a corr.in file.
    Doesn't really do much though. Meant to find
    relevant indexes for ECI."""

    def __init__(self,filename):
        """Fill all the corr.in info"""

        self._source=filename

        self._header=[]
        with open(filename) as f:
            self._header.append(f.readline().split())
            self._header.append(f.readline().split())
            self._header.append(f.readline().split())
        
        self._matrix=np.loadtxt(filename,skiprows=3)

    def extract_eci_inds(self, indx):
        """Get the indexes of the given configuration where the correlations are greater than 0.
        Only works if the correlations are made up of only 0 and 1 (pure endstate).

        :indx: int
        :returns: ndarray

        """
        row=self._matrix[indx,:]
        print row
        oneinds=row==1.00
        zeroinds=row==0.00

        if not np.alltrue(oneinds == np.logical_not(zeroinds)):
            print "A row with values other than 0.0 and 1.0 has been selected!"
            exit()

        return oneinds

    def extract_energy_inds(self, indx):
        """Get the indexes of all the configurations that are binary, given the index
        of a binary index.

        :indx: int
        :returns: ndarray

        """
        row=self._matrix[indx,:]
        print row
        oneinds=row==1.00
        zeroinds=row==0.00

        if not np.alltrue(oneinds == np.logical_not(zeroinds)):
            print "A row with values other than 0.0 and 1.0 has been selected!"
            exit()

        corrsum=np.sum(abs(self._matrix,axis=0))


class EnergyFiler(object):

    """Hold the data from an energy file"""

    def __init__(self,filename):
        """Read and store header separately from data

        :filename: path to energy file

        """
        self._filename = filename
        self._header=[]
        with open(filename) as f:
            self._header.append(f.readline().split())

        self._values=np.loadtxt(filename,dtype=str)


class Blocker(object):

    """Manipulate ECI, correlations and energies so you can create eci.in, energy
    and corr.in files for only ternary interactions by subtracting a binary fit.
    """

    def __init__(self, eciin3, corrin3, energy3):
        """Load all the files needed to perform a ternary
        fit.

        :eciin3: string
        :corrin3: string
        :energy3: string
        :returns: Blocker

        """
        self._eciF=EciInFiler(eciin3)
        self._corrF=CorrFiler(corrin3)
        self._eF=EnergyFiler(energy3)

        self._eci=np.arange(len(self._eciF._lines))
        self._corr=np.array(self._corrF._matrix,dtype=float)
        self._e=np.array(self._eF._values[:,0],dtype=float)

        #print np.shape(self._eci)
        #print np.shape(self._corr)
        #print np.shape(self._e)

        
    def subtract(self, eciin2, indx):
        """Given the index of a binary configuration, construct a new
        Blocker that has the binary energies subtracted and has had
        all the binary rows and columns removed from corr, eci and energy.

        :eciin2: string (eci.in file for binary fit)
        :indx: int (index corresponding to endstate)
        :returns: Blocker

        """
        eci2=EciInFiler(eciin2)
        corr3inds=self._corrF.extract_eci_inds(indx)

        corr3=self._corr[:,corr3inds]
        eci3=self._eci[corr3inds]
        pass
