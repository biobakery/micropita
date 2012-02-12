#######################################################
# Author: Timothy Tickle
# Description: Class to Allow KMedoids on a metric space.
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
from scipy.spatial.distance import squareform

#This allows one to use custom distance metrics with Kmediods
class MLPYDistanceAdaptor:

    #Class variables
    #Distance matrix to reference
    matrix = None

    #Requires a matrix of distances, could be condensed or square matrices
    #@params tempDistanceMatrix The distance matrix to be used
    #@params tempIsCondensedMatrix Boolean indicator of the matrix being square (true = condensed; false = square)
    def __init__(self, tempDistanceMatrix, tempIsCondensedMatrix):
        if(tempIsCondensedMatrix):
            self.matrix = squareform(tempDistanceMatrix)
        else:
            self.matrix = tempDistanceMatrix

    #This is the only method required in the interface to MLPY to be a distance metric
    #Does NOT want values but positions, the positions will be used for accessing the distance matrix already provided.
    #@params x X position as a array of 1 number
    #@params y Y position as a array of 1 number
    def compute(self,x,y):
        if(self.matrix == None):
            raise Exception("".join(["MLPYDistanceAdaptor. Attempted to compute distance with out a distance matrix passed in during construction."]))
        return self.matrix[x[0],y[0]]
