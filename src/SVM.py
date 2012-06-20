#######################################################
# Author: Timothy Tickle
# Description: Class to Allow Support Vector Machine 
# analysis and to contain associated scripts.
#######################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2011"
__credits__ = ["Timothy Tickle"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

#Libraries
from AbundanceTable import AbundanceTable
from Constants import Constants
from CommandLine import CommandLine
import math
import operator
import os
from random import shuffle
from ValidateData import ValidateData

class SVM:

    #SVM programs to call
    c_SVM_SCALE = "svm-scale"
    c_SVM_SCALE_LOCATION = "./external/libsvm-3.11/python/"+c_SVM_SCALE+".py"
    c_SVM_TRAIN = "svm-train"
    c_SVM_TRAIN_LOCATION = "/usr/bin/"+c_SVM_TRAIN
    c_SVM_PREDICT = "svm-predict"
    c_SVM_PREDICT_LOCATION = "/usr/bin/"+c_SVM_PREDICT
    c_GRID_LOCATION = "./external/libsvm-3.11/tools/grid.py"
    c_GNUPLOT = "/usr/bin/gnuplot"

    #SVM Parameters
    c_LOG_C = "-log2c"
    c_LINEAR_KERNEL = "0"
    c_C_SVC = "0"
    c_SVC_PROBABILITY_ESTIMATES = "1"

    #Dictionary keywords for files
    c_KEYWORD_INPUT_FILE = "INPUT"
    c_KEYWORD_SCALED_FILE = "SCALED"
    c_KEYWORD_RANGE_FILE = "RANGE"
    c_KEYWORD_CV_OUT_FILE = "CV_OUT"
#    c_KEYWORD_CV_PNG_FILE = "CV_PNG"
    c_KEYWORD_MODEL_FILE = "MODEL"
    c_KEYWORD_PREDICTION_SCALED_FILE = "SCALED_FOR_PREDICTION"
    c_KEYWORD_PREDICTION_FILE = "PREDICTION"
    c_COST_VALUE = "C"
    c_ACCURACY = "ACCURACY"

#TODO Unfreeze two tests
    #Tested 2/4 test cases
    #Creates a model using a linear SVM
    #1. Scales data if indicated
    #2. Uses cross validation to define the gamma, cost, and range
    #3. Uses the optimal gamma to create a model
    #@params tempInputFileName String Must be a valid file  and must have data in rows in the format 
    #<label> <index1>:<value1> <index2>:<value2> ...
    #This method crates a system of outputfiles.
    #@params tempScaling -1 or 0 indicates that scaling should be performed on the data first.
    #The given value 0 or -1 will be used as the lower bound for scaling.
    #1 is always used for the upper bound.
    #This will perform X scaling not y.
    #A Scaled file will be saved with the extension ".scaled"
    #Given this is a linear kernel the process should be invariant to scale but this is given as an option.
    #@params tempPrediction Indicates prediction should occur.
    def createLinearModel(self, tempInputFileName=None, tempTMPDirectory=None, tempScaling=None, tempLogC="-1,2,1", tempCrossValidationFold = 10, tempProbabilistic=False):
        #ValidateData
        if(not ValidateData.isValidFileName(tempInputFileName)):
            print "".join(["SVM.createLinearModel: ","The input file ",tempInputFileName," is not a valid file name."])
            return False

        #TODO Need to validate
        costList = [filter(None,strPiece) for strPiece in tempLogC.split(",")]
        costList = sorted(costList)

        #Create output file names
        inputFile = Constants.QUOTE+tempInputFileName+Constants.QUOTE
        fileNamePrefix = os.path.splitext(tempInputFileName)[0]

        #Creat file names and delete of they exist
        if(tempTMPDirectory == None):
            scaledFile = "".join([Constants.QUOTE,fileNamePrefix,Constants.c_SCALED_FILE_EXT,Constants.QUOTE])
            rangeFile = "".join([Constants.QUOTE,fileNamePrefix,Constants.c_SCALING_PARAMETERS,Constants.QUOTE])
            cvOutFile = "".join([Constants.QUOTE,fileNamePrefix,Constants.c_CV_FILE_EXT,Constants.QUOTE])
            modelFile = "".join([Constants.QUOTE,fileNamePrefix,Constants.c_MODEL_FILE_EXT,Constants.QUOTE])
        else:
            scaledFile = "".join([Constants.QUOTE,tempTMPDirectory,Constants.c_SCALED_FILE_EXT,Constants.QUOTE])
            rangeFile = "".join([Constants.QUOTE,tempTMPDirectory,Constants.c_SCALING_PARAMETERS,Constants.QUOTE])
            cvOutFile = "".join([Constants.QUOTE,tempTMPDirectory,Constants.c_CV_FILE_EXT,Constants.QUOTE])
            modelFile = "".join([Constants.QUOTE,tempTMPDirectory,Constants.c_MODEL_FILE_EXT,Constants.QUOTE])
        if os.path.exists(scaledFile):
            os.remove(scaledFile)
        if os.path.exists(rangeFile):
            os.remove(rangeFile)
        if os.path.exists(cvOutFile):
            os.remove(cvOutFile)
        if os.path.exists(modelFile):
            os.remove(modelFile)

        #Dict contains generated files for return value and later reference
        generatedFiles = dict()
        generatedFiles[self.c_KEYWORD_INPUT_FILE] = tempInputFileName

        #Create Commandline interface
        shell = CommandLine()

        #Scale data
        if((tempScaling == -1)or(tempScaling == 0)):
            cmd = [self.c_SVM_SCALE,"-l",str(tempScaling),"-u","1","-s",rangeFile,inputFile,">",scaledFile]
            print cmd
            noError = shell.runPipedCommandLine(tempCommand = cmd)
            if(noError == False):
                print "".join(["Received error during scaling. Command=",cmd])
                return False
            generatedFiles[self.c_KEYWORD_SCALED_FILE] = scaledFile.strip(Constants.QUOTE)
            generatedFiles[self.c_KEYWORD_RANGE_FILE] = rangeFile.strip(Constants.QUOTE)
        
        #Cross validate for the c parameter
        #If scaling was used, use the scaled file, otherwise use the original file.
        if((tempScaling == -1)or(tempScaling == 0)):
            #Indicate file to perform CV on
            cv_file = scaledFile
        else:
            cv_file = inputFile

        #Perform cross validation
        #Holds the cost and accuracy
        lldCostAccuracy = []
        #Holds the highest cost
        bestCost = 0.0
        #Holds the associated accuracy
        highestAccuracy = -1
        with open(cvOutFile.strip(Constants.QUOTE),'a') as f:
            for cost in costList:
                cost = math.pow(2,int(cost))
                cmd = [Constants.QUOTE+self.c_SVM_TRAIN_LOCATION+Constants.QUOTE,"-s",self.c_C_SVC,"-t",self.c_LINEAR_KERNEL,"-c",str(cost),"-v",str(tempCrossValidationFold),"-q",cv_file]
                print "Cross Validation"
                print cmd
                result = shell.runPipedCommandLine(tempCommand = cmd)

                #TODO This is not very clean
                accuracy = result[0].split(" = ")
                if(len(accuracy) == 2):
                    accuracy = accuracy[1][:-2]
                else:
                    accuracy = False
                    break

                lldCostAccuracy.append([cost,float(accuracy)])
                f.write("".join(["Cost=",str(cost)," with ",str(accuracy),"% Cross Validation Accuracy\n"]))

                #Get best cost parameter
                if(not accuracy == False):
                    if(highestAccuracy < float(accuracy)):
                        bestCost = cost
                        highestAccuracy = float(accuracy)
        f.close()

        if(accuracy == False):
            print "".join(["Received error during cross validation on scaled file. Command=",cmd])
            return False

        #Sort [cost, accuracy] by accuracy
        lldCostAccuracy = sorted(lldCostAccuracy, key=operator.itemgetter(1))
        #If the accuracy does not change
        if len(lldCostAccuracy) > 1:
            #If the accuracy does not change, use the middle entry for the c parameter,
            #this assumes the c value are sorted in the function, of which they are.
            #Otherwise leave best cost and accuracy as is
            if lldCostAccuracy[0][1] == lldCostAccuracy[len(lldCostAccuracy)-1][1]:
                bestCost,highestAccuracy = lldCostAccuracy[int(len(lldCostAccuracy)/2)]

        print "LIBSVM::Using cost=", bestCost, " log cost=", math.log(bestCost,2)
        generatedFiles[self.c_COST_VALUE] = str(math.log(bestCost,2))
        generatedFiles[self.c_ACCURACY] = str(highestAccuracy)
        generatedFiles[self.c_KEYWORD_CV_OUT_FILE] = cvOutFile.strip(Constants.QUOTE)

        #Create model
        if(tempProbabilistic):
            cmd = [Constants.QUOTE+self.c_SVM_TRAIN_LOCATION+Constants.QUOTE,"-s",self.c_C_SVC,"-t",self.c_LINEAR_KERNEL,"-c",str(bestCost),"-b",self.c_SVC_PROBABILITY_ESTIMATES,scaledFile,modelFile]
        else:
            cmd = [Constants.QUOTE+self.c_SVM_TRAIN_LOCATION+Constants.QUOTE,"-s",self.c_C_SVC,"-t",self.c_LINEAR_KERNEL,"-c",str(bestCost),scaledFile,modelFile]
        print cmd
        noError = shell.runPipedCommandLine(tempCommand = cmd)
        if(noError == False):
            print "".join(["Received error during cross validation on scaled file. Command=",str(cmd)])
            return False
        generatedFiles[self.c_KEYWORD_MODEL_FILE] = modelFile.strip(Constants.QUOTE)

        #Return generated files
        return generatedFiles

#TODO Test
    def predictFromLinearModel(self, tempDataFileName=None, tempModelFileName=None, tempRangeFileName=None, tempProbabilistic=False):
        #ValidateData
        if(not ValidateData.isValidFileName(tempDataFileName)):
            print "".join(["The data file ",tempDataFileName," is not a valid file name."])
            return False
        if(not ValidateData.isValidString(tempModelFileName)):
            print "".join(["The model file ",tempModelFileName," is not a valid file name."])
            return False
        if(not ValidateData.isValidString(tempRangeFileName)):
            print "".join(["The range file ",tempRangeFileName," is not a valid file name."])
            return False

        #Get file prefix
        dataFileNamePrefix = os.path.splitext(tempDataFileName)[0]

        #Associated files
        generatedFiles = dict()

        #Create output file names
        modelFile = Constants.QUOTE+tempModelFileName+Constants.QUOTE
        dataFile = Constants.QUOTE+tempDataFileName+Constants.QUOTE
        scaledFile = Constants.QUOTE+dataFileNamePrefix+Constants.c_SCALED_FOR_PREDICTION_FILE_EXT+Constants.QUOTE
        rangeFile = Constants.QUOTE+tempRangeFileName+Constants.QUOTE
        predictFile = Constants.QUOTE+dataFileNamePrefix+Constants.c_PREDICT_FILE_EXT+Constants.QUOTE

        #Create Commandline interface
        shell = CommandLine()

        #Scale data
        cmd = [self.c_SVM_SCALE,"-r",rangeFile,dataFile,">",scaledFile]
        print "SCALE "+str(cmd)
        noError = shell.runPipedCommandLine(tempCommand = cmd)
        if(noError == False):
            print "".join(["Received error during scaling. Command=",cmd])
            return False
        generatedFiles[self.c_KEYWORD_PREDICTION_SCALED_FILE] = scaledFile.strip(Constants.QUOTE)
        
        #Predict data
        if(tempProbabilistic):
            cmd = [self.c_SVM_PREDICT_LOCATION,"-b",self.c_SVC_PROBABILITY_ESTIMATES,scaledFile,modelFile,predictFile]
        else:
            cmd = [self.c_SVM_PREDICT_LOCATION,scaledFile,modelFile,predictFile]
        print "PREDICT "+str(cmd)
        noError = shell.runPipedCommandLine(tempCommand = cmd)
        if(noError == False):
            print "".join(["Received error during cross validation on scaled file. Command=",cmd])
            return False
        generatedFiles[self.c_KEYWORD_PREDICTION_FILE] = predictFile.strip(Constants.QUOTE)
        return generatedFiles

    #1 Happy Path tested
    #Converts abundance files to input SVM files.
    #@tempInputFile Abundance file to read (should be a standard Qiime output abundance table)
    #@tempOutputSVMFile File to save SVM data to when converted from teh abundance table
    #@tempDelimiter Delimiter of the Abundance table
    #@tempLabels Ordered labels to use to classify the samples in the abundance table
    #@sLastMetadataName The name of the last row in the abundance table representing metadata
    #@tempSkipFirstColumn Boolean Indicates to skip the first column (true) (for instance if it contains taxonomy identifiers)
    #@tempNormalize Boolean to indicate if the abundance data should be normalized (true) before creating the file (normalized by total sample abundance)
    @staticmethod
    def convertAbundanceTableToSVMFile(abndAbundanceTable, tempOutputSVMFile, sMetadataLabel):
        #Validate parameters
        if abndAbundanceTable == None:
            print "Error, invalid Abundance table."
            return False
        if(not ValidateData.isValidString(tempOutputSVMFile)):
            print "Error, file not valid. File:"+str(tempOutputSVMFile)
            return False

        #If output file exists, delete
        if(os.path.exists(tempOutputSVMFile)):
            os.remove(tempOutputSVMFile)

        #Create data matrix
        dataMatrix = zip(*abndAbundanceTable.funcGetAbundanceCopy())

        #Add labels
        llData = []
        lsLabels = abndAbundanceTable.funcGetMetadata(sMetadataLabel)
        lsUniqueLabels = list(set(lsLabels))
        dictLabels = dict([[str(lenuLabels[1]),str(lenuLabels[0])] for lenuLabels in enumerate(lsUniqueLabels)])
        lsLabels = [dictLabels[sLabel] for sLabel in lsLabels]

        iRowIndex = 0
        for dataRow in dataMatrix[1:]:
            llData.append(" ".join([lsLabels[iRowIndex]]+[Constants.COLON.join([str(enuSamples[0]+1),str(enuSamples[1])])
                            for enuSamples in enumerate(dataRow)])+Constants.ENDLINE)
            iRowIndex = iRowIndex + 1

        #Output file
        with open(tempOutputSVMFile,'a') as f:
            (f.write("".join(llData)))
        f.close()
        return lsUniqueLabels

    #Tested
    @staticmethod
    def funcScaleFeature(npdData):
        if sum(npdData) == 0 or len(set(npdData))==1:
            return npdData
        dMin = min(npdData)
        return (npdData-dMin)/float(max(npdData-dMin))

    #Tested
    @staticmethod
    def funcWeightLabels(lLabels):
        #Convert to dict
        lUniqueLabels = list(set(lLabels))
        dictLabels = dict(zip(lUniqueLabels, range(len(lUniqueLabels))))

        #Build a dict of weights per label {label:weight, label:weight}
        #Get the occurence of each label
        dictWeights = dict()
        for sLabelKey in dictLabels:
            sCurLabel = dictLabels[sLabelKey]
            dictWeights[sCurLabel] = lLabels.count(sLabelKey)

        #Divide the highest occurence each occurence
        iMaxOccurence = max(dictWeights.values())
        for sWeightKey in dictWeights:
            dictWeights[sWeightKey]=iMaxOccurence/float(dictWeights[sWeightKey])

        return [dictWeights,lUniqueLabels]

    #Tested 3/4 cases TODO could add in test 12 with randomize True
    def func10FoldCrossvalidation(self, iTotalSampleCount, fRandomise = False):
        """
        Generates the indexes for a 10 fold crossvalidation given a sample count.
        If there are less than 10 samples, it uses the sample count as the K-fold crossvalidation
        as a leave one out method.

	:param	iTotalSampleCount:	Total Sample Count
	:type	int:	Sample Count
	:param	fRandomise:	Random sample indices
	:type	boolean:	True indicates randomise (Default False)
        """

        #Make indices and shuffle if needed
        liindices = range(iTotalSampleCount)
        if fRandomise:
            shuffle(liindices)

        #For 10 times
        iKFold = 10
        if iTotalSampleCount < iKFold:
            iKFold = iTotalSampleCount
        for iiteration in xrange(iKFold):
            lfTraining = [iindex % iKFold != iiteration for iindex in liindices]
            lfValidation = [not iindex for iindex in lfTraining]
            yield lfTraining, lfValidation
