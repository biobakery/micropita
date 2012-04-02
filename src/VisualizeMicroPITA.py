#!/usr/bin/env python
"""
Author: Timothy Tickle
Description: Class to Visualize Analysis for the MicroPITA paper
"""

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"


class VisualizeMicroPITA:
    """
    Micropita class
    """


    c_progHC = "./external/hclust/hclust.py"

    #Figure 1b
    funcHierchicalClusterAll()
    {
        #Do not attempt to generate if you only have one selection criteria, there will be nothing to cluster/show as columns
        if(len(selectedSamples.keys())>1):
            #Figure 1b: Top ranked samples and how they are different
            microPITA.redefineMetricSelectionsAsAbsPres(tempSelectedSamplesDict=selectedSamples, tempOutputFile=absPresOutputFile)
            #Hierarchical cluster absence presence metrix matrix
            CommandLine.CommandLine().runCommandLine([c_progHC, "--in", absPresOutputFile, "--out", figure1bFile, "--label2cols", figure1ColorFilePath, "-l", figure1LabelFilePath, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1"])
    }

    #Figure 2
    funcCladogramSelection()
    {
        #Create files for circulator
        figureTreeFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_PhylotypicTaxaSelection_2.txt"
        figureColorFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_PhylotypicColor_2.txt"
        figureCircleFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_PhylotypicCircle_2.txt"
        figureTickFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_PhylotypicTick_2.txt"
        figureHighLightFilePath=Constants.Constants.INPUT_DATA_DIRECTORY+"microPITA/"+"PhylotypicHighLight_2.txt"
        figureStyleFilePath=Constants.Constants.INPUT_DATA_DIRECTORY+"microPITA/"+"PhylotypicStyle_2.txt"

        #Flags to make sure all files were generated or exist before calling script
        circleFileExist = False
        colorFileExist = False
        highLightFileExist = False
        styleFileExist = False
        tickFileExist = False
        treeFileExist = False

        #Do not remove the style file, it is static
        if(os.path.exists(figureStyleFilePath)):
            styleFileExist = True
        #Remove files if already exist
        #If the color file exists, delete
        if(os.path.exists(figureTreeFilePath)):
            os.remove(figureTreeFilePath)
        #If the label file exists, delete
        if(os.path.exists(figureColorFilePath)):
            os.remove(figureColorFilePath)
        #If the circle file exists, delete
        if(os.path.exists(figureCircleFilePath)):
            os.remove(figureCircleFilePath)
        #If the highlight file exists, delete
        if(os.path.exists(figureHighLightFilePath)):
            os.remove(figureHighLightFilePath)
        #If the tick file exists, delete
        if(os.path.exists(figureTickFilePath)):
            os.remove(figureTickFilePath)

        #Generate color file
        #Create file handle to write to files
        colorFileWriter = FileIO.FileIO(figureColorFilePath,False,True,False)
        colorList = list()
        #Color by technique
        colorList.append(microPITA.toLegendName(microPITA.c_CHAO1_A_DIVERSITY)+Constants.Constants.TAB+microPITA.chao1ColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_INV_SIMPSON_A_DIVERSITY)+Constants.Constants.TAB+microPITA.invSimpsonColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_UNIFRAC_B_DIVERSITY)+Constants.Constants.TAB+microPITA.unifracColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_WEIGHTED_UNIFRAC_B_DIVERSITY)+Constants.Constants.TAB+microPITA.weightedUnifracColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_BRAY_CURTIS_B_DIVERSITY)+Constants.Constants.TAB+microPITA.brayCurtisColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY)+Constants.Constants.TAB+microPITA.inBrayCurtisColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_INVERSE_UNIFRAC_B_DIVERSITY)+Constants.Constants.TAB+microPITA.inUnifracColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_INVERSE_WEIGHTED_UNIFRAC_B_DIVERSITY)+Constants.Constants.TAB+microPITA.inWeightedUnifracColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_RANDOM)+Constants.Constants.TAB+microPITA.randomColorN)
        colorList.append(microPITA.toLegendName(microPITA.c_USER_RANKED)+Constants.Constants.TAB+microPITA.userRankedN)
        colorList.append(microPITA.toLegendName(microPITA.c_SVM_CLOSE)+Constants.Constants.TAB+microPITA.svmCloseN)
        colorList.append(microPITA.toLegendName(microPITA.c_SVM_FAR)+Constants.Constants.TAB+microPITA.svmFarN)
        colorList.append(c_WHITEOUT+Constants.Constants.TAB+c_WHITEOUT_COLOR)
        colorFileWriter.writeToFile(Constants.Constants.ENDLINE.join(colorList))
        colorFileWriter.close()
        colorList = None
        colorFileExist = True

        #Start writing the tick file, we know the taxonmy levels will always be the same.
        #Then add the method names
        #Then close the writer.
        tickFileWriter = FileIO.FileIO(figureTickFilePath,False,True,True)
        tickFileWriter.writeToFile("1\tPhyla\n2\tClasses\n3\tOrders\n4\tFamilies\n5\tGenera")
        tickFileExist = True
        tickLevelCount = 6

        #Get the id names being used
        circleContent = [[id] for id in list(abundance[sampleID])]

        #Make sure dendrogram data is available before moving forward.
        if(len(circleContent)>0):
            #Indicate which taxa were in selected samples
            #Get selected samples of a method
            #Stores the taxa that are important to different methods
            taxaMembership = dict()
            for selectedSampleMethod in selectedSamples:
#TODO Waiting on circlader update                            #Write method name to tick file so that it will be displayed
#                            tickFileWriter.writeToFile("\n"+str(tickLevelCount)+"\t"+selectedSampleMethod)
#                            tickLevelCount=tickLevelCount+1
                #Get taxa interesting to a selection technique
                taxaMethodSet = set()
                #Holds the u-test score for each taxa comparing taxon abundances from selected samples and not selected samples
                taxaTScores = dict()
                selectedTaxaAbundancies = selectedSamples[selectedSampleMethod]
                selectedSamplesUTest = [] #Start with false to exclude sample ids
                notSelectedSamplesUTest = [] #Start with false to exclude sample ids
                #Set up boolean list to compress array
                for sample in sampleNames:
                    wasSelected = sample in selectedTaxaAbundancies
                    selectedSamplesUTest.append(wasSelected)
                    notSelectedSamplesUTest.append(not wasSelected)

                #Compress arrays to one or the other distribution
                #Conduct wilcoxon tests on all taxa
                maxScore = 0
                for taxonIndex in xrange(0,len(abundance[sampleID])):
                    taxaData = list(abundance[taxonIndex,])
                    taxaId = taxaData[0]
                    distribution = np.array(taxaData[1:])
                    selectedDistribution = np.compress(selectedSamplesUTest,distribution)
                    notSelectedDistribution = np.compress(notSelectedSamplesUTest,distribution)
                    score, pvalue = stats.ranksums(selectedDistribution,notSelectedDistribution)
#                    score, pvalue = stats.t_two_sample(selectedDistribution,notSelectedDistribution)
                    taxaTScores[taxaId]=[score,pvalue]
                    #Store the maxScore for my own use
                    if score > maxScore:
                        maxScore = score

                #This cuts down the taxa to only selected taxa of some abundance
                #For each selected sample, go through abundance data and select taxa which had more than 0 abundance
                for sample in selectedSamples[selectedSampleMethod]:
                    #Get taxa of each sample
                    taxaAbundanciesList = abundance[sample]
                    taxaOfInterest = list()
                    for abundanceIndex in xrange(0,len(taxaAbundanciesList)):
                        currentTaxon = abundance[sampleID][abundanceIndex]
                        #Automatically add if it is a user defined taxa
                        if(currentTaxon in userDefinedTaxa):
                            taxaOfInterest.append(currentTaxon)
                        #Otherwise add based on abundance
                        else: 
                            if(float(taxaAbundanciesList[abundanceIndex]) > 0.0):
                                taxaOfInterest.append(currentTaxon)
                    #Combine (union) taxa to a set
                    taxaMethodSet = taxaMethodSet | set(taxaOfInterest)
#                 #Check to see if the taxa was of interest
#                interest = microPITA.toLegendName(selectedSampleMethod)+":1.0!"
#                noInterest = microPITA.toLegendName(selectedSampleMethod)+":0.0!"
                #For every taxa in abundance data
                for circleIndex in xrange(0,len(circleContent)):
                    currentCircleElement = circleContent[circleIndex]
                    if(currentCircleElement[0] in taxaMethodSet):
                      #IF not showing occurence, showing enrichment
                      if(not c_FIGURE_2_ABSOLUTE_TAXA_OCCURENCE):
                          #Indicate direction of triangle based on score sign
                          taxaSignificance = taxaTScores[currentCircleElement[0]]
                          #If score is above 0 (Enriched)
                          if(taxaSignificance[microPITA.c_SCORE_INDEX]> 0.0):
                              #IF using score
                              if(c_FIGURE_2_USE_SCORE):
                                  #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                  if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                      currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(abs(taxaSignificance[microPITA.c_SCORE_INDEX])/5.2)+"!^#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                  else:
                                      currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(abs(taxaSignificance[microPITA.c_SCORE_INDEX])/5.2)+"!^#0.0")
                              #If Not using score, using pvalue
                              elif(c_FIGURE_2_USE_PVALUE):
                                  #IF only plotting significant p-values
                                  if(c_FIGURE_2_SIG_PVALUES_ONLY):
                                      #If pvalue is significant
                                      if(taxaSignificance[microPITA.c_PVALUE_INDEX]<0.05):
                                          #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                          if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                             currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!^#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                          else:
                                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!^#0.0")
                                     #If NOT pvalue is significant
                                     else:
                                          #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                          if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                              currentCircleElement.append(c_WHITEOUT+":1.0!^#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                          else:
                                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":0.0!^")
                                  #IF NOT only plotting significant p-values
                                  else:
                                      #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                      if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                          currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!^#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                      else:
                                          currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!^#0.0")

                          #IF score is below or equal to 0 (under enrichment)
                          else:
                              #If using score
                              if(c_FIGURE_2_USE_SCORE):
                                  #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                  if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                      currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(abs(taxaSignificance[microPITA.c_SCORE_INDEX])/5.2)+"!v#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                  else:
                                      currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(abs(taxaSignificance[microPITA.c_SCORE_INDEX])/5.2)+"!v#0.0")
                              #If NOT using score, using pvalue
                              elif(c_FIGURE_2_USE_PVALUE):
                                  #If using significant p-values only
                                  if(c_FIGURE_2_SIG_PVALUES_ONLY):
                                      #If the pvalue is significant
                                      if(taxaSignificance[microPITA.c_PVALUE_INDEX]<0.05):
                                          #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                          if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!v#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                          else:
                                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!v#0.0")
                                      #If the pvalue is not significant
                                      else:
                                          #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                          if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                              currentCircleElement.append(c_WHITEOUT+":1.0!v#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                          else:
                                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":0.0!v#0.0")
                                  #If using all pvalues
                                  else:
                                      #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                                      if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                                          currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!v#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                                      else:
                                          currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":"+str(1-taxaSignificance[microPITA.c_PVALUE_INDEX])+"!v#0.0")
                      #Performing absolute occurence so just indicate was selected
                      else:
                          #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                          if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":1.0!R#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                          else:
                              currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":1.0!R#0.0")
                #Not selected
                else:
                    #If the selection method is user ranked and the taxa is a user selected taxon, add a border
                    if((selectedSampleMethod == microPITA.c_USER_RANKED) and (currentCircleElement[0] in userDefinedTaxa)):
                        currentCircleElement.append(c_WHITEOUT+":1.0!R#"+c_FIGURE_2_TAXA_BORDER_WIDTH)
                    else:
                        currentCircleElement.append(microPITA.toLegendName(selectedSampleMethod)+":0.0!R#0.0")

            #Close the tick file writer
            tickFileWriter.close()

            #Add to treeContent
            #Build file
            for circleIndex in xrange(0,len(circleContent)):
                currentCircleData = circleContent[circleIndex]
                circleContent[circleIndex] = ".".join(re.split("\|",currentCircleData[0]))+Constants.Constants.TAB+Constants.Constants.TAB.join(currentCircleData[1:])
                    
            circleFileWriter = FileIO.FileIO(figureCircleFilePath,False,True,False)
            circleFileWriter.writeToFile(Constants.Constants.ENDLINE.join(circleContent))
            circleFileWriter.close()
            circleFileExist = True
            circleContent = None

        #Create tree file
        #Write to file and change | with .
        #Get the id names being used
        taxa = list(abundance[sampleID])
        treeContent = []
        treeContent = [".".join(re.split("\|",taxon)) for taxon in taxa]

        #Fill in dendrogram with ancester nodes
        fullTree = treeContent
        for node in treeContent:
            ancestryList = re.split("\.",node)
            currentNode = ancestryList[0]
            if(not currentNode in fullTree):
                fullTree.append(currentNode)
            if(len(ancestryList) > 1):
                for ancestryListNode in ancestryList[1:]:
                    currentNode = currentNode+"."+ancestryListNode
                    if(not currentNode in fullTree):
                        fullTree.append(currentNode)
        treeContent = None

        #Make sure dendrogram data is available before moving forward
        #If so write the file that will contain the tree structure
        if(len(fullTree)>0):
            treeFileWriter = FileIO.FileIO(figureTreeFilePath,False,True,False)
            treeFileWriter.writeToFile(Constants.Constants.ENDLINE.join(fullTree))
            treeFileWriter.close()
            treeFileExist = True

        #Generate Highlight file
        #Used to indicate the highlighting of the internal area or tree area
        #Attempting to highlight leaf nodes by the average abundance in the samples
        #Attempting to highlight internal biology at the class level
        #Holds the average abundance in the samples sum(abundace)/sample count 
        totalAbundanceInExperiment = dict()
        totalSampleCount = len(abundance[0])-1
        totalTaxaCount = len(abundance[sampleID])
        totalTaxaCountFloat = float(totalTaxaCount)
        #Get abundance data from each taxa
        #Holds the averageValues to get the median average abundance
#        averageAbundanceValues = []
        for taxaIndex in xrange(0,totalTaxaCount):
            taxaData = list(abundance[taxaIndex])
#            taxaPresence = 0.0
#This calculates the average occurence not abundance                        #In the abundace data, if the abundance is greater than one, count
#            for taxaCount in taxaData[1:]:
#                if float(taxaCount) > 0.0:
#                    taxaPresence = taxaPresence+1.0
#            totalAbundanceInExperiment[taxaData[0]] = taxaPresence/totalSampleCount
            averageTaxaAbundance = taxaData[1:]
            averageTaxaAbundance = sum(averageTaxaAbundance)/len(averageTaxaAbundance)
            #Update averages
#            averageAbundanceValues.append(averageTaxaAbundance)
            #Store average abundance by taxon id
            totalAbundanceInExperiment[taxaData[0]] = averageTaxaAbundance
           hLightFileWriter = FileIO.FileIO(figureHighLightFilePath,False,True,False)
           hLightWriteList = list()
        #Holds the taxa id information used to indicate what biology to highlight
        parentsToHighLightArea = dict()
        #For each taxon, write a entry indicating highlight color for leaf nodes
        #And store taxa id information for later processing for highlighting
        for totalTaxon in totalAbundanceInExperiment.keys():
            #Write leaf nodes for end node highlighting (or atleast append to the list that will be written)
#            currentHighLightColor = str(1-(totalAbundanceInExperiment[totalTaxon]/maxAverageTaxaAbundance))
            currentHighLightColor = str(1-totalAbundanceInExperiment[totalTaxon])
            highLightTaxonIdElements = re.split("\|",totalTaxon)
            highLightWriteList.append(Constants.Constants.TAB.join([".".join(highLightTaxonIdElements),"","","_c_["+currentHighLightColor+","+currentHighLightColor+","+currentHighLightColor+"]"]))

            #Store every level of the ancestor state/taxon id and an associated count
            previousAncestorLevel = highLightTaxonIdElements[0]
            if(not previousAncestorLevel in parentsToHighLightArea):
                parentsToHighLightArea[previousAncestorLevel] = 0
            parentsToHighLightArea[previousAncestorLevel] = parentsToHighLightArea[previousAncestorLevel]+1
            for ancestorLevel in highLightTaxonIdElements[1:]:
                previousAncestorLevel = ".".join([previousAncestorLevel,ancestorLevel])
                if(not previousAncestorLevel in parentsToHighLightArea):
                    parentsToHighLightArea[previousAncestorLevel] = 0
                parentsToHighLightArea[previousAncestorLevel] = parentsToHighLightArea[previousAncestorLevel]+1

        #Append to write list the different ancestor levels that have more leaves than a certain threshold
        parentCount = 1
        for level,count in parentsToHighLightArea.items():
            if(count > c_FIGURE_2_HIGHLIGHT_THRESHOLD):
                highLightWriteList.append(Constants.Constants.TAB.join([level,str(parentCount)+":"+re.split("\.",level)[-1],"","_c_[1.0,0.7,0.0]"]))#+currentHighLightColor+","+currentHighLightColor+","+currentHighLightColor+"]"]))
                parentCount = parentCount + 1
        highLightFileWriter.writeToFile(Constants.Constants.ENDLINE.join(highLightWriteList))
        highLightFileWriter.close()
        highLightFileExist = True

        #Call command
        print("treeFileExist:"+str(treeFileExist))
        print("colorFileExist:"+str(colorFileExist))
        print("circleFileExist:"+str(circleFileExist))
        print("styleFileExist:"+str(styleFileExist))
        print("tickFileExist:"+str(tickFileExist))
        print("highLightFileExist:"+str(highLightFileExist))
        if(treeFileExist and colorFileExist and circleFileExist and styleFileExist and highLightFileExist and tickFileExist):
            CommandLine.CommandLine().runCommandLine(["./external/circlader/circlader.py", figureTreeFilePath, figure2File, "--style_file", figureStyleFilePath, "--tree_format", "tabular", "--color_file", figureColorFilePath, "--circle_file", figureCircleFilePath, "--highlight_file", figureHighLightFilePath, "--tick_file", figureTickFilePath])
    }

    #Figure 3
    funcFigure3()
    {
        pass
    }

    #Figure 4
    funcFigure4()
    {
        pass
    }

    #Figure 5
    funcFigure5()
    {
        pass
    }

    #Figure 6
    funcFigure6()
    {
        pass
    }

    #Supplemental Figure 1
    funcFigureS1()
    {
        #Figure S1
        #Create file names
        figureColorFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarity_Color_S1.txt"
        figureLabelFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarity_Label_S1.txt"
        figureDataFile = inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarity_Data_S1.txt"
        outputFileName= inputFilePrefix+"_"+str(sampleSelectionCount)++"_ExtremeDissimilarityS1.pdf"

        #Sample selection amount
        selectCount = sampleSelectionCount

        #Transpose data
        tAbundance = rawData.transposeDataMatrix(tempMatrix=abundance, tempRemoveAdornments=True)

        #Get samples
        selectedSamplesBC = []
        selectedSamplesIBC = []
        selectedSamplesBoth = []
        if(microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY in selectedSamples):
            selectedSamplesIBC = selectedSamples[microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY]
        if(microPITA.c_BRAY_CURTIS_B_DIVERSITY in selectedSamples):
            selectedSamplesBC = selectedSamples[microPITA.c_BRAY_CURTIS_B_DIVERSITY]
        selectedSamplesBoth = list(set(selectedSamplesBC) & set(selectedSamplesIBC))

        #Generate HClust files
        #If the color file exists, delete
        if(os.path.exists(figureColorFilePath)):
            os.remove(figureColorFilePath)
        #If the label file exists, delete
        if(os.path.exists(figureLabelFilePath)):
            os.remove(figureLabelFilePath)
        #If the data file exists, delete
        if(os.path.exists(figureDataFile)):
            os.remove(figureDataFile)

        #Create color / label files
        #Create file handle to write to files
        colorFileWriter = FileIO.FileIO(figureColorFilePath,False,True,False)
        labelFileWriter = FileIO.FileIO(figureLabelFilePath,False,True,False)
        colorList = list()
        labelList = list()
        #Color representative dissimilarity samples
        for selectedSampleName in selectedSamplesBC:
            colorList.append(selectedSampleName+Constants.Constants.TAB+microPITA.brayCurtisColor)
            labelList.append(selectedSampleName+Constants.Constants.TAB+selectedSampleName)
        #Color extreme dissimilarity samples
        for selectedSampleName in selectedSamplesIBC:
             colorList.append(selectedSampleName+Constants.Constants.TAB+microPITA.inBrayCurtisColor)
             labelList.append(selectedSampleName+Constants.Constants.TAB+selectedSampleName)
        #Color dually selected samples
        for selectedSampleName in selectedSamplesIBC:
            colorList.append(selectedSampleName+Constants.Constants.TAB+microPITA.c_SELECTED_SAMPLE)
            labelList.append(selectedSampleName+Constants.Constants.TAB+selectedSampleName)
        colorFileWriter.writeToFile(Constants.Constants.ENDLINE.join(colorList))
        labelFileWriter.writeToFile(Constants.Constants.ENDLINE.join(labelList))
        colorFileWriter.close()
        labelFileWriter.close()

        #Create data file
        dataWriteContent = ["ID"+Constants.Constants.TAB+Constants.Constants.TAB.join(sampleNames)]
        for abundanceLine in abundance:
            abundanceLineList = []
            for abundanceLineElement in abundanceLine:
                abundanceLineList.append(str(abundanceLineElement))
            dataWriteContent.append(Constants.Constants.TAB.join(abundanceLineList))
        dataFileWriter = FileIO.FileIO(figureDataFile,False,True,False)
        dataFileWriter.writeToFile(Constants.Constants.ENDLINE.join(dataWriteContent))
        dataFileWriter.close()

        #Call command
        CommandLine.CommandLine().runCommandLine(["./external/hclust/hclust.py", "--in", figureDataFile, "--out", outputFileName, "--label2cols", figureColorFilePath, "-l", figureLabelFilePath, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "-d", "correlation", "-f", "braycurtis","-y","0.0075", "--grid","0"])
    }

    #Supplemental Figure 2
    funcFigureS2()
    {
        #Create file names
        figureColorFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarityColorS2.txt"
        figureLabelFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarityLabelS2.txt"
        figureDataFile = inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarityDataS2.txt"
        figureDistanceFile = inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarityDistancesS2.txt"
        outputFileName= inputFilePrefix+"_"+str(sampleSelectionCount)+"_ExtremeDissimilarityS2.pdf"

        #Sample selection amount
        selectCount = sampleSelectionCount

        #Transposed data
        tAbundance = rawData.transposeDataMatrix(tempMatrix=abundance, tempRemoveAdornments=True)

        #Generate distance matrix
        distanceMatrix = squareform(microPITA.getBetaMetric(tempAbundancies=tAbundance, tempMetric=microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY))
        distanceMatrix.transpose()

        #Get samples
        selectedSamples = selectedSamples[microPITA.c_INVERSE_BRAY_CURTIS_B_DIVERSITY]

        #Generate HClust files
        #If the color file exists, delete
        if(os.path.exists(figureColorFilePath)):
            os.remove(figureColorFilePath)
        #If the label file exists, delete
        if(os.path.exists(figureLabelFilePath)):
            os.remove(figureLabelFilePath)
        #If the distance file exists, delete
        if(os.path.exists(figureDistanceFile)):
            os.remove(figureDistanceFile)
        #If the data file exists, delete
        if(os.path.exists(figureDataFile)):
            os.remove(figureDataFile)

        #Create color / label files
        #Create file handle to write to files
        colorFileWriter = FileIO.FileIO(figureColorFilePath,False,True,True)
        labelFileWriter = FileIO.FileIO(figureLabelFilePath,False,True,True)
        colorList = list()
        labelList = list()
        for selectedSampleName in selectedSamples:
            colorList.append(selectedSampleName+Constants.Constants.TAB+microPITA.brayCurtisColor)
            labelList.append(selectedSampleName+Constants.Constants.TAB+selectedSampleName)
        colorFileWriter.writeToFile(Constants.Constants.ENDLINE.join(colorList))
        labelFileWriter.writeToFile(Constants.Constants.ENDLINE.join(labelList))
        colorFileWriter.close()
        labelFileWriter.close()

        #Create data file
        dataWriteContent = ["ID"+Constants.Constants.TAB+Constants.Constants.TAB.join(sampleNames)]
        for abundanceLine in abundance:
            abundanceLineList = []
            for abundanceLineElement in abundanceLine:
                abundanceLineList.append(str(abundanceLineElement))
            dataWriteContent.append(Constants.Constants.TAB.join(abundanceLineList))
        dataFileWriter = FileIO.FileIO(figureDataFile,False,True,False)
        dataFileWriter.writeToFile(Constants.Constants.ENDLINE.join(dataWriteContent))
        dataFileWriter.close()

        #Create distance matrix file
        dataWriteContent = []
        for distanceLine in distanceMatrix:
            distanceLineList = []
            for distanceLineElement in distanceLine:
                distanceLineList.append(str(distanceLineElement))
            dataWriteContent.append(Constants.Constants.TAB.join(distanceLineList))
        dataFileWriterDistances = FileIO.FileIO(figureDistanceFile,False,True,False)
        dataFileWriterDistances.writeToFile(Constants.Constants.ENDLINE.join(dataWriteContent))
        dataFileWriterDistances.close()

        #Call command
        CommandLine.CommandLine().runCommandLine(["./external/hclust/hclust.py", "--in", figureDataFile, "--out", outputFileName, "--label2cols", figureColorFilePath, "-l", figureLabelFilePath, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "-f", "correlation", "--dms", figureDistanceFile,"-y","0.0075", "--grid","0"])
    }

    #Supplemental Figure 3
    funcFigureS3()
    {
        #Generates a figure that shows the samples selected by which method.
        #Currently generates one file per selection method.
        #The title of the file has the selection method in it for reference.
        for selectionMethod in selectedSamples.keys():
            #Create file names
            figureColorFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_"+str(selectionMethod)+"_Color_S1.txt"
            figureLabelFilePath=inputFilePrefix+"_"+str(sampleSelectionCount)+"_"+str(selectionMethod)+"_Label_S1.txt"
            figureDataFile = inputFilePrefix+"_"+str(sampleSelectionCount)+"_"+str(selectionMethod)+"_Data_S1.txt"
            outputFileName= inputFilePrefix+"_"+str(sampleSelectionCount)+"_"+str(selectionMethod)+"S1.pdf"

            #Sample selection amount
            selectCount = sampleSelectionCount

            #Transpose data
            tAbundance = rawData.transposeDataMatrix(tempMatrix=abundance, tempRemoveAdornments=True)

            #Get samples
            currentSelectedSamples = selectedSamples[selectionMethod]

            #Generate HClust files
            #If the color file exists, delete
            if(os.path.exists(figureColorFilePath)):
                os.remove(figureColorFilePath)
            #If the label file exists, delete
            if(os.path.exists(figureLabelFilePath)):
                os.remove(figureLabelFilePath)
            #If the data file exists, delete
            if(os.path.exists(figureDataFile)):
                os.remove(figureDataFile)

            #Create color / label files
            #Create file handle to write to files
            colorFileWriter = FileIO.FileIO(figureColorFilePath,False,True,False)
            labelFileWriter = FileIO.FileIO(figureLabelFilePath,False,True,False)
            colorList = list()
            labelList = list()
            #Color representative dissimilarity samples
            for selectedSampleName in currentSelectedSamples:
                colorList.append(selectedSampleName+Constants.Constants.TAB+microPITA.c_SELECTED_SAMPLE)
                labelList.append(selectedSampleName+Constants.Constants.TAB+selectedSampleName)
            colorFileWriter.writeToFile(Constants.Constants.ENDLINE.join(colorList))
            labelFileWriter.writeToFile(Constants.Constants.ENDLINE.join(labelList))
            colorFileWriter.close()
            labelFileWriter.close()

            #Create data file
            dataWriteContent = ["ID"+Constants.Constants.TAB+Constants.Constants.TAB.join(sampleNames)]
            for abundanceLine in abundance:
                abundanceLineList = []
                for abundanceLineElement in abundanceLine:
                    abundanceLineList.append(str(abundanceLineElement))
                dataWriteContent.append(Constants.Constants.TAB.join(abundanceLineList))
            dataFileWriter = FileIO.FileIO(figureDataFile,False,True,False)
            dataFileWriter.writeToFile(Constants.Constants.ENDLINE.join(dataWriteContent))
            dataFileWriter.close()

            #Call command
            CommandLine.CommandLine().runCommandLine(["./external/hclust/hclust.py", "--in", figureDataFile, "--out", outputFileName, "--label2cols", figureColorFilePath, "-l", figureLabelFilePath, "--legend", "1", "--legend_ncol", "1", "--pad_inches", "1.5", "--fdend_w", "0", "--font_size", "12", "--cm_h", "0", "-c", "Blues", "--grid","1"])
    }

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MicroPITA.py", 
    description = """Selects samples from abundance tables based on various selection schemes.""" )

#Arguments
#Outputfile
argp.add_argument( "strOutFile", metavar = "output.txt", nargs = "?", help = "An optional output file" )
#Abundance file
argp.add_argument( "strFileAbund", metavar = "Abundance_file", help = "An abundance table." )
#Taxa file
argp.add_argument( "strFileTaxa", metavar = "Taxa_file",
    help = "A file containing taxa to be used in taxa-directed selection." )
#Temporary folder
argp.add_argument( "strTMPDir", metavar = "Temporary_Directory", help = "Directory to place temporary and intermediate files.")
#Select count
argp.add_argument( "icount", metavar = "number", type = int, help = "The number of samples to select (An integer greater than 0.)." )
#Selection parameter
argp.add_argument("strSelection", metavar = "Selection_Methods", help = "Select techniques listed one after another.", nargs="*")

__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
    args = argp.parse_args( )
    dictSelectedSamples = MicroPITA().run(strOutputFile=args.strOutFile, strInputAbundanceFile=args.strFileAbund, strUserDefinedTaxaFile=args.strFileTaxa, strTemporaryDirectory=args.strTMPDir, iSampleSelectionCount=args.icount, strSelectionTechnique=args.strSelection)
    strOutputContent = ""
    for sKey in dictSelectedSamples:
        strOutputContent = "".join([strOutputContent,sKey,Constants.COLON,Constants.ENDLINE,", ".join(dictSelectedSamples[sKey]),Constants.ENDLINE])

    #Write to file
    if(not strOutputContent == ""):
        fHndlOutput = open(args.strOutFile,'w')
        fHndlOutput.write(str(strOutputContent))
        fHndlOutput.close()

if __name__ == "__main__":
    _main( )

