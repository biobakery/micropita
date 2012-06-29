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
from SVM import SVM
from random import shuffle
from ValidateData import ValidateData

class LIBSVM(SVM):

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
    def createLinearModel(self, tempInputFileName, tempTMPDirectory, tempScaling, tempLogC, tempCrossValidationFold = 10, tempProbabilistic=False):
        #ValidateData
        if(not ValidateData.funcIsValidFileName(tempInputFileName)):
            print "".join(["SVM.createLinearModel: ","The input file ",tempInputFileName," is not a valid file name."])
            return False

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
        for f in [scaledFile,rangeFile,cvOutFile,modelFile]:
            if os.path.exists(f):
                os.remove(f)

        #Dict contains generated files for return value and later reference
        generatedFiles = dict()
        generatedFiles[Constants.c_strKeywordInputFile] = tempInputFileName

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
            generatedFiles[Constants.c_strKeywordScaledFile] = scaledFile.strip(Constants.QUOTE)
            generatedFiles[Constants.c_strKeywordRangeFile] = rangeFile.strip(Constants.QUOTE)
        
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
            for cost in tempLogC:
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
        generatedFiles[Constants.c_strKeywordCostValue] = str(math.log(bestCost,2))
        generatedFiles[Constants.c_strKeywordAccuracy] = str(highestAccuracy)
        generatedFiles[Constants.c_strKeywordCVOutFile] = cvOutFile.strip(Constants.QUOTE)

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
        generatedFiles[Constants.c_strKeywordModelFile] = modelFile.strip(Constants.QUOTE)

        #Return generated files
        return generatedFiles

#TODO Test
    def predictFromLinearModel(self, tempDataFileName=None, tempModelFileName=None, tempRangeFileName=None, tempProbabilistic=False):
        #ValidateData
        if(not ValidateData.funcIsValidFileName(tempDataFileName)):
            print "".join(["The data file ",tempDataFileName," is not a valid file name."])
            return False
        if(not ValidateData.funcIsValidString(tempModelFileName)):
            print "".join(["The model file ",tempModelFileName," is not a valid file name."])
            return False
        if(not ValidateData.funcIsValidString(tempRangeFileName)):
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
        generatedFiles[Constants.c_strKeywordScaledPredFile] = scaledFile.strip(Constants.QUOTE)
        
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
        generatedFiles[Constants.c_strKeywordPredFile] = predictFile.strip(Constants.QUOTE)
        return generatedFiles
