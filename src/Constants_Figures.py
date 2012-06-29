#######################################################
#
#	Title:		Constants_Figures
#	Author:		Timothy Tickle
#	Date:		02/05/2012
#	Purpose:	Class to hold Constants associated with figures
#
#######################################################

##
#Used to test the FileIO class
from MicroPITA import MicroPITA

class Constants_Figures():
    #lag to record if the colors are inverted
    c_fInverted = False

    #Controls the global alpha for global transparent plotting
    c_dAlpha = 0.5

    #General colors
    c_strBackgroundColorName = "Invisible"
    c_strBackgroundColor = "255,255,255"
    c_strBackgroundColorTuple = (255,255,255)
    c_strBackgroundColorWord = "white"
    c_strBackgroundColorLetter = "w"

    c_strDetailsColor = "0,0,0"
    c_strDetailsColorTuple = (0,0,0)
    c_strDetailsColorWord = "black"
    c_strDetailsColorLetter = "k"

    #Standardized coloring
    invSimpsonColor = "#FF6600"#Dark orange
    invSimpsonColorN = "255,102,0"#Dark orange
    invSimpsonColorLetter = 'o'#Dark orange
    invSimpsonColorTuple = (255,102,0)
    chao1Color = "#66FF00"#Lawn green
    chao1ColorN = "102,255,0"#Lawn green
    chao1ColorLetter = 'g'
    chao1ColorTuple = (102,255,0)
    unifracColor = "#FFCC33"#Gold
    unifracColorN = "255,204,51"#Gold
    inUnifracColor = "#FA8072"#Salmon
    inUnifracColorN = "250,128,114"#Salmon
    weightedUnifracColor = "#0000FF"#Blue
    weightedUnifracColorN = "0,0,255"#Blue
    inWeightedUnifracColor = "#0099FF"#Doger blue
    inWeightedUnifracColorN = "0,153,255"#Doger blue
    brayCurtisColor = "#6600FF"#Purple
    brayCurtisColorN = "102,0,255"#Purple
    brayCurtisColorLetter = ''
    brayCurtisColorTuple = (102,0,255)
    inBrayCurtisColor = "#FF33FF"#Fushia
    inBrayCurtisColorN = "255,51,255"#Fushia
    inBrayCurtisColorLetter = 'm'
    inBrayCurtisColorTuple = (255,51,255)
    randomColor = "#000000"#Black
    randomColorN = "0,0,0"#Black
    randomColorLetter = 'k'
    randomColorTuple = (0,0,0)
    userRanked = "#83C8F9"#Blue
    userRankedN = "131,200,249"#Blue
    userRankedLetter = 'c'
    userRankedTuple = (131,200,249)
#    userRanked = "#CCCCCC"#Grey
#    userRankedN = "204,204,204"#Grey
#Grey just was not working for the figure 2 got to have a more vibriant color    userRankedN = "204,204,204"#Grey
    svmClose = "#FF0000"#Red
    svmCloseN = "255,0,0"#Red
    svmCloseLetter = 'r'
    svmCloseTuple = (255,0,0)
    svmFar = "#009900"#Green
    svmFarLetter = 'g'
    svmFarN = "0,153,0"#Green
    svmFarTuple = (0,153,0)

    c_SELECTED_SAMPLE = "#FF0000" #Indicates a sample that is selected

    #General plotting
    c_strGridLineColor = "#CCCCCC"

    #PCOA related
    c_charNoSelect = "#FFFFFF" # White
    c_charPCOAShape = "o"
    c_charPCOAMultSelectColor = "#CCCCCC"
    c_charPCOAMultSelectionShape = "+"
    c_strPCOAMultSelectionName = "Multiple Selection"
    c_strPCOANotSelected = "Not_Selected"
    iMarkerSize = 100
    iPieMarkerSize = 30
    c_charPCOAPieChart = "o"
    c_charPCOASquarePieChart = "s"

    def invertColors(self,fInvert):
        if fInvert==True:
            #General colors
            self.c_strBackgroundColor = "0,0,0"
            self.c_strBackgroundColorTuple = (0,0,0)
            self.c_strBackgroundColorWord = "black"
            self.c_strBackgroundColorLetter = "k"

            self.c_strDetailsColor = "255,255,255"
            self.c_strDetailsColorTuple = (255,255,255)
            self.c_strDetailsColorWord = "white"
            self.c_strDetailsColorLetter = "w"

            #Invert random because it is the same color as the background toggling
            self.randomColor = "#FFFFFF"#White
            self.randomColorN = "255,255,255"#White

            #Invert no select color
            self.c_charNoSelect = "#000000" # black

            #Record that it is inverted
            self.c_fInverted = True

            #Alpha looks best at full in inversion
            self.c_dAlpha = 1.0

        else:
            #General colors
            self.c_strBackgroundColor = "255,255,255"
            self.c_strBackgroundColorTuple = (255,255,255)
            self.c_strBackgroundColorWord = "white"
            self.c_strBackgroundColorLetter = "w"

            self.c_strDetailsColor = "0,0,0"
            self.c_strDetailsColorTuple = (0,0,0)
            self.c_strDetailsColorWord = "black"
            self.c_strDetailsColorLetter = "k"

            #Standardized coloring
            self.randomColor = "#000000"#Black
            self.randomColorN = "0,0,0"#Black

            #No select color
            self.c_charNoSelect = "#FFFFFF" # White

            #Record that it is not inverted
            self.c_fInverted = False

            #Alpha looks best at full in inversion
            self.c_dAlpha = 0.5

    #Can be used to convert method names to colors
    dictConvertMethodToHEXColor={MicroPITA.c_strDiversity1:invSimpsonColor,
                              MicroPITA.c_strDiversity2:chao1Color,
                              MicroPITA.c_strExtremeDissimiarity1:inBrayCurtisColor,
                              MicroPITA.c_strSVMClose:svmClose,
                              MicroPITA.c_strSVMFar:svmFar,
                              MicroPITA.c_strRandom:randomColor,
                              MicroPITA.c_strRepresentativeDissimilarity1:brayCurtisColor,
                              MicroPITA.c_strUserRanked:userRanked}

    #Can be used to convert method names to rgb tuples
    dictConvertMethodToRGBLetter={MicroPITA.c_strDiversity1:invSimpsonColorLetter,
                              MicroPITA.c_strDiversity2:chao1ColorLetter,
                              MicroPITA.c_strExtremeDissimiarity1:inBrayCurtisColorLetter,
                              MicroPITA.c_strSVMClose:svmCloseLetter,
                              MicroPITA.c_strSVMFar:svmFarLetter,
                              MicroPITA.c_strRandom:randomColorLetter,
                              MicroPITA.c_strRepresentativeDissimilarity1:brayCurtisColorLetter,
                              MicroPITA.c_strUserRanked:userRankedLetter}

    #Can be used to convert method names to rgb strings
    dictConvertMethodToRGBString={MicroPITA.c_strDiversity1:invSimpsonColorN,
                              MicroPITA.c_strDiversity2:chao1ColorN,
                              MicroPITA.c_strExtremeDissimiarity1:inBrayCurtisColorN,
                              MicroPITA.c_strSVMClose:svmCloseN,
                              MicroPITA.c_strSVMFar:svmFarN,
                              MicroPITA.c_strRandom:randomColorN,
                              MicroPITA.c_strRepresentativeDissimilarity1:brayCurtisColorN,
                              MicroPITA.c_strUserRanked:userRankedN}
