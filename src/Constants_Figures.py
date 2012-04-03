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
    c_strBackgroundColorWord = "white"
    c_strBackgroundColorLetter = "w"

    c_strDetailsColor = "0,0,0"
    c_strDetailsColorWord = "black"
    c_strDetailsColorLetter = "k"

    #Standardized coloring
    invSimpsonColor = "#FF6600"#Dark orange
    invSimpsonColorN = "255,102,0"#Dark orange
    chao1Color = "#009900"#Green
    chao1ColorN = "0,153,0"#Green
    unifracColor = "#FFCC33"#Gold
    unifracColorN = "255,204,51"#Gold
    inUnifracColor = "#66FF00"#Lawn green
    inUnifracColorN = "102,255,0"#Lawn green
    weightedUnifracColor = "#0000FF"#Blue
    weightedUnifracColorN = "0,0,255"#Blue
    inWeightedUnifracColor = "#0099FF"#Doger blue
    inWeightedUnifracColorN = "0,153,255"#Doger blue
    brayCurtisColor = "#6600FF"#Purple
    brayCurtisColorN = "102,0,255"#Purple
    inBrayCurtisColor = "#FF33FF"#Fushia
    inBrayCurtisColorN = "255,51,255"#Fushia
    randomColor = "#000000"#Black
    randomColorN = "0,0,0"#Black
    userRanked = "#83C8F9"#Blue
    userRankedN = "131,200,249"#Blue
#    userRanked = "#CCCCCC"#Grey
#    userRankedN = "204,204,204"#Grey
#Grey just was not working for the figure 2 got to have a more vibriant color    userRankedN = "204,204,204"#Grey
    svmClose = "#FF0000"#Red
    svmCloseN = "255,0,0"#Red
    svmFar = "#FA8072"#Salmon
    svmFarN = "250,128,114"#Salmon

    c_SELECTED_SAMPLE = "#FF0000" #Indicates a sample that is selected

    #PCOA related
    c_charNoSelect = "#FFFFFF" # White
    c_charPCOAShape = "o"
    c_charPCOAMultSelectColor = "#CCCCCC"
    c_charPCOAMultSelectionShape = "+"
    c_strPCOAMultSelectionName = "Multiple Selection"

    def invertColors(self,fInvert):
        if fInvert==True:
            #General colors
            self.c_strBackgroundColor = "0,0,0"
            self.c_strBackgroundColorWord = "black"
            self.c_strBackgroundColorLetter = "k"

            self.c_strDetailsColor = "255,255,255"
            self.c_strDetailsColorWord = "white"
            self.c_strDetailsColorLetter = "w"

            #Invert random because it is the same color as the background toggling
            self.randomColor = "#FFFFFF"#White
            self.randomColorN = "255,255,255"#White

            #Invert no select color
            self.c_charNoSelect = "#000000" # black

            #Record that it is inverted
            c_fInverted = True

        else:
            #General colors
            self.c_strBackgroundColor = "255,255,255"
            self.c_strBackgroundColorWord = "white"
            self.c_strBackgroundColorLetter = "w"

            self.c_strDetailsColor = "0,0,0"
            self.c_strDetailsColorWord = "black"
            self.c_strDetailsColorLetter = "k"

            #Standardized coloring
            self.randomColor = "#000000"#Black
            self.randomColorN = "0,0,0"#Black

            #No select color
            self.c_charNoSelect = "#FFFFFF" # White

            #Record that it is not inverted
            c_fInverted = False

    #Can be used to convert method names to colors
    dictConvertMethodToHEXColor={MicroPITA.c_DIVERSITY_1:invSimpsonColor,
                              MicroPITA.c_DIVERSITY_2:chao1Color,
                              MicroPITA.c_EXTREME_DISSIMILARITY_1:inBrayCurtisColor,
                              MicroPITA.c_SVM_CLOSE:svmClose,
                              MicroPITA.c_SVM_FAR:svmFar,
                              MicroPITA.c_RANDOM:randomColor,
                              MicroPITA.c_REPRESENTATIVE_DISSIMILARITY_1:brayCurtisColor,
                              MicroPITA.c_USER_RANKED:userRanked}
