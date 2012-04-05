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

#Creates a confusion matrix
class PlotMatrix:

  #Given a matrix and labels consistent to the matrix, plot a confustion matrix
  @staticmethod
  def funcPlotMatrix(npMatrix, lsLabels, strOutputFigurePath, strPredictedTitle="Predicted", strActualTitle="Actual", fFlipYLabels=False):

    plt.clf()
    figConfusionMatrix = plt.figure()
    objAxis = figConfusionMatrix.add_subplot(111)

    #Format labels to be every other label
    lNewXLabels = []
    for strLabel in lsLabels:
        lNewXLabels.append("")
        lNewXLabels.append(strLabel)

    #Get y labels
    lNewYLabels = lNewXLabels
    if fFlipYLabels:
        lsLabels.reverse()
        lNewYLabels = []
        for strLabel in lsLabels:
            lNewYLabels.append("")
            lNewYLabels.append(strLabel)

    #Set x axis and position
    objAxis.xaxis.set_ticklabels(lNewXLabels)
    objAxis.xaxis.set_ticks_position('top')

    #Set y axis
    objAxis.yaxis.set_ticklabels(lNewYLabels)

    #Set axis titles
    ylabel(strActualTitle)
    plt.suptitle(strPredictedTitle)

    #Plot matrix values
    objPlot = objAxis.imshow(np.array(npMatrix), cmap=cm.jet, interpolation='nearest')

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
#strOutputPath = "confmat.png"
#PlotMatrix.funcPlotMatrix(npMatrix, lsLabels, strOutputPath, fFlipYLabels=True)
