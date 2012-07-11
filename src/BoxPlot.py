#######################################################
# Author: Timothy Tickle
# Description: Class to create box plots
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
from Constants_Figures import Constants_Figures
import matplotlib.pyplot as plt
import numpy as np
import random
from pylab import *

#Plots a matrix
class BoxPlot:

  @staticmethod
  def funcPlot(ly, lsLabels, strOutputFigurePath, strTitle = "Title", strXTitle="X Axis", strYTitle="Y Axis", strColor = "#83C8F9", fJitter=False, fInvert=False):
    """
    Plot a box plot with optional jittering.
    """

    #Start plot
    #Get plot object
    imgFigure = plt.figure()

    #Get plot colorsstrOutFigure
    objFigureControl = Constants_Figures()
    objFigureControl.invertColors(fInvert=fInvert)

    #Color/Invert figure
    imgFigure.set_facecolor(objFigureControl.c_strBackgroundColorWord)
    imgSubplot = imgFigure.add_subplot(111,axisbg=objFigureControl.c_strBackgroundColorLetter)
    imgSubplot.set_xlabel(strXTitle)
    imgSubplot.set_ylabel(strYTitle)
    imgSubplot.spines['top'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.spines['bottom'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.spines['left'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.spines['right'].set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.xaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)

    #Adds light grid for numbers and puts them in the background
    imgSubplot.yaxis.grid(True, linestyle='-', which='major', color=objFigureControl.c_strGridLineColor, alpha=objFigureControl.c_dAlpha)
    imgSubplot.set_axisbelow(True)
    imgSubplot.yaxis.label.set_color(objFigureControl.c_strDetailsColorLetter)
    imgSubplot.tick_params(axis='x', colors=objFigureControl.c_strDetailsColorLetter)
    imgSubplot.tick_params(axis='y', colors=objFigureControl.c_strDetailsColorLetter)
    charMarkerEdgeColor = objFigureControl.c_strDetailsColorLetter

    #Make box plot
    bp = plt.boxplot(x=ly, notch=1, patch_artist=True)
    for iindex, ldData in enumerate(ly):
      ldX = None
      if fJitter:
        ldX = [float(iindex+1)+ uniform(-.05,.05) for x in xrange(len(ldData))]
      else:
        ldX = [float(iindex+1) for x in xrange(len(ldData))]
      plt.scatter(x=ldX,y=ly[iindex],c=strColor,marker="o",alpha=objFigureControl.c_dAlpha)

    #Color boxes
    plt.setp(bp['boxes'], color=objFigureControl.c_strDetailsColorLetter, facecolor=strColor, alpha=objFigureControl.c_dAlpha)
    plt.setp(bp['whiskers'], color=objFigureControl.c_strDetailsColorLetter)

    #Set ticks and title
    xtickNames = plt.setp(imgSubplot, xticklabels=lsLabels)
    imgSubplot.set_title(strTitle)

    #End plot
    #Save to a file
    imgFigure.savefig(strOutputFigurePath, facecolor=imgFigure.get_facecolor())
