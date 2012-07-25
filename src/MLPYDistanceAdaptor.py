"""
Author: Timothy Tickle
Description: Allows KMedoids on a custom metric space..
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
from scipy.spatial.distance import squareform

class MLPYDistanceAdaptor:
    """
    Allows one to use custom distance metrics with KMedoids in the MLPY package.
    """

    npaMatrix = None
    """
    Distance matrix to reference.
    """

    def __init__(self, npaDistanceMatrix, fIsCondensedMatrix):
        """
        Constructor requires a matrix of distances, could be condensed or square matrices

    	:param	npaDistanceMatrix:	The distance matrix to be used
	    :type	Numpy array
	    :param	fIsCondensedMatrix:	Indicator of the matrix being square (true = condensed; false = square)
	    :type	Boolean
        """

        if(fIsCondensedMatrix):
            self.npaMatrix = squareform(npaDistanceMatrix)
        else:
            self.npaMatrix = npaDistanceMatrix

    def compute(self,x,y):
        """
        This is the only method required in the interface to MLPY to be a distance metric.
        Does NOT want values but positions, the positions will be used for accessing the distance matrix already provided.

	    :param	x:	X position as a array of 1 number
	    :type	Numpy array
	    :param	y:	Y position as a array of 1 number
	    :type	Boolean
        """

        if(self.npaMatrix == None):
            raise Exception("".join(["MLPYDistanceAdaptor. Attempted to compute distance with out a distance matrix passed in during construction."]))
        return self.npaMatrix[x[0],y[0]]
