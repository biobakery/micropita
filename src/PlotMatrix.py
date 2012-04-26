#######################################################
# Author: Timothy Tickle
# Description: Class to Create Confusion Matrix
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#External libraries
import logging
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

#Plots a matrix
class PlotMatrix:

  #Given a matrix and labels consistent to the matrix, plot a matrix
  @staticmethod
  def funcPlotMatrix(npMatrix, lsLabels, strOutputFigurePath, strXTitle="X Axis", strYTitle="Y Axis", fFlipYLabels=False, fInvert=False):

    plt.clf()
    figConfusionMatrix = plt.figure()
    objAxis = figConfusionMatrix.add_subplot(111)

    #Get y labels
    lNewYLabels = lsLabels
    if fFlipYLabels:
        lsLabels.reverse()

    #Set x axis and position
    objAxis.xaxis.set_ticklabels([""]+lsLabels)
    objAxis.xaxis.set_ticks_position('top')

    #Set y axis
    objAxis.yaxis.set_ticklabels([""]+lNewYLabels)

    #Set axis titles
    ylabel(strYTitle)
    plt.suptitle(strXTitle)

    #Plot matrix values
    objPlot = objAxis.imshow(np.array(npMatrix), cmap=get_cmap("Blues"), interpolation='nearest')

    #Plot text values
    for yIndex, ldRow in enumerate(npMatrix):
        for xIndex, dValue in enumerate(ldRow):
            plt.text(xIndex, yIndex, dValue, fontdict = {'size':18,'weight':'bold'} )

    #Add color bar
    figConfusionMatrix.colorbar(objPlot)

    #Save to a file
    savefig(strOutputFigurePath)

#npMatrix = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]
#lsLabels = ["Label 1","Label 2","Label 3","Label 4"]
#strOutputPath = "confusion.png"
#PlotMatrix.funcPlotMatrix(npMatrix, lsLabels, strOutputPath, fFlipYLabels=True)
