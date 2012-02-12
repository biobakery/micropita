#######################################################
#
#	Title:		Utility_File
#	Author:		Timothy Tickle
#	Date:		06/28/2010
#	Purpose:	Utility class for file operations not
#                       neccessarily tied to reading and writing.
#
#######################################################

#Import libaries
from FileIO import FileIO
from ValidateData import ValidateData
#import Log
import os
import shutil
import string

##
#Utility function for file manipulation
class Utility_File():

    ##
    #Contructor
    def __init__(): pass

    ##
    #Tested 6-29-10
    #Splits a file to to files of the given amount of lines.
    #If there is not an equal split, the last file is smaller
    #Files are stored given the fileName with _# to the end of the fileName
    #@param tempNumberOfLines Number of lines to put in each split (smaller) file. IF the amount if greater than the file the file will be generated as one file.
    #@param tempReadFileName File to read lines from
    #@param tempOutputFileName Base file name to use for the split files
    @staticmethod
    def splitFile(tempNumberOfLines,tempReadFileName,tempOutputFileName):

        #Validate parameters
        if(not ValidateData.isValidPositiveInteger(tempNumberOfLines)):return False
        if(not ValidateData.isValidString(tempOutputFileName)):return False
        if(not ValidateData.isValidString(tempReadFileName)):return False

        #If file exists
        if(False == os.path.isfile(tempReadFileName)):
            return False

        #Get the file name part
        outputFileName = Utility_File.incrementFileName(tempOutputFileName)
        if(outputFileName == False): return False
        #Create file writers
        fileWriter = None
        fileReader = None
        try:
            fileWriter = FileIO(outputFileName,False,True,False)
        except:
            return False

        try:
            fileReader = FileIO(tempReadFileName,True,False,False)
        except:
            fileWriter.close()
            return False

        #Read file and split into files of the given line number amount
        lineCount = 0
        currentLine = None

        while ((not currentLine == False) and (not currentLine == "")):
            while((lineCount<tempNumberOfLines) and (not currentLine == False) and (not currentLine == "")):
                currentLine = fileReader.readline()
                fileWriter.writeToFile(currentLine)
                lineCount = lineCount + 1

            lineCount = 0
            fileWriter.close()
            outputFileName = Utility_File.incrementFileName(outputFileName)
            fileWriter = FileIO(outputFileName,False,True,False)

        fileWriter.close()
        fileReader.close()

        return True

    #Tested 4/4/2011
    #Abridges a file to a given amount of lines starting at the beginning and removing the end if needed.
    #@param tempNumberOfLines Number of lines to abridge to
    #@param tempReadFileName File to read lines from
    #@param tempOutputFileName Base file name to use for the split files
    @staticmethod
    def abridgeFile(tempNumberOfLines,tempReadFileName,tempOutputFileName):

        #Validate parameters
        if(not ValidateData.isValidPositiveInteger(tempNumberOfLines)):return False
        if(not ValidateData.isValidString(tempOutputFileName)):return False
        if(not ValidateData.isValidString(tempReadFileName)):return False

        #If file exists
        if(False == os.path.isfile(tempReadFileName)):
            return False

        #Create file writers
        fileWriter = None
        fileReader = None
        try:
            fileWriter = FileIO(tempOutputFileName,False,True,False)
        except:
            return False

        try:
            fileReader = FileIO(tempReadFileName,True,False,False)
        except:
            fileWriter.close()
            return False

        #Read file and split into files of the given line number amount
        lineCount = 0
        currentLine = None

        while((lineCount<tempNumberOfLines) and (not currentLine == False) and (not currentLine == "")):
            currentLine = fileReader.readline()
            fileWriter.writeToFile(currentLine)
            lineCount = lineCount + 1
#        Log.Log().logMessage("Utility_File","AbridgeFile","Abridged file to "+str(lineCount)+" lines.")
        print(" ".join(["Utility_File","AbridgeFile","Abridged file to "+str(lineCount)+" lines."]))

        fileReader.close()
        fileWriter.close()
        return True

    #Take the base file name and add a _# to the name before the extension
    #If _# is already there then # is incremented
    #Tested 1/28/2010
    #@param tempFileName String file name to be updated
    @staticmethod
    def incrementFileName(tempFileName):
        #Validate parameters
        if(not ValidateData.isValidString(tempFileName)):return False
        
        #Split to get extension and body of file name
        fileNamePieces = tempFileName.split(".")
        fileNamePiecesCount = len(fileNamePieces)
        fileNameBody = tempFileName
        fileNameExtension = ""
        if(fileNamePiecesCount == 2):
            fileNameBody = fileNamePieces[0]
            fileNameExtension = "."+fileNamePieces[1]
        elif(fileNamePiecesCount>2):
            fileNameBody = ".".join(fileNamePieces[:(fileNamePiecesCount-1)])
            fileNameExtension = "."+fileNamePieces[fileNamePiecesCount-1]

        #Check to see if the body has already been updated in the past and update
        fileBodyPieces = fileNameBody.split("_")
        fileBodyPiecesCount = len(fileBodyPieces)

        if(fileBodyPieces>1):
            numberCheck = fileBodyPieces[fileBodyPiecesCount-1]
            if(numberCheck.isdigit()):
                number = int(numberCheck)
                fileNameBody = "_".join(fileBodyPieces[:len(fileBodyPieces)-1])
                fileNameBody = fileNameBody +"_"+str(number+1)
            else:
                fileNameBody = fileNameBody+"_1"
        else:
            fileNameBody = fileNameBody+"_1"

        return fileNameBody +fileNameExtension

    #Tested 4/4/2011
    #Given a file name return the prefix of the name (everything but the extension)
    #If a "." appears in the name, the last is considered to indicate the seperator between the prefix and extension.
    #@param tempFileName String file name to be updated
    @staticmethod
    def getFileNamePrefix(tempFileName):
        #Validate parameters
        if(not ValidateData.isValidString(tempFileName)):return False
        
        #Split to get extension and body of file name
        fileNamePieces = tempFileName.split(".")
        fileNamePiecesCount = len(fileNamePieces)
        fileNameBody = tempFileName
        if(fileNamePiecesCount == 2):
            fileNameBody = fileNamePieces[0]
        elif(fileNamePiecesCount>2):
            fileNameBody = ".".join(fileNamePieces[:(fileNamePiecesCount-1)])
        return fileNameBody

    #Tested 2-10-11
    #Combines the list of given files into one file
    #@params outputFileName Name of the file to output the combined results
    #@params listOfFiles List of strings which are names of files to be combined
    #params tempRemoveHeaderLine boolean True indicates removing the first line of each file to be combined.
    @staticmethod
    def combineFiles(tempOutputFileName = None, tempListOfFiles = None, tempRemoveHeaderLine = False):
        #check parameters
        if( not ValidateData.isValidString(tempOutputFileName)):
            return False
        if( not ValidateData.isValidList(tempListOfFiles)):
            return False
        if( not ValidateData.isValidBoolean(tempRemoveHeaderLine)):
            return False

        #Create file writers
        fileWriter = None
        try:
            fileWriter = open(tempOutputFileName, 'w')
        except:
            return False

        #Go through each file and combine
        for fileName in tempListOfFiles:
            if not fileName == None:
                if(os.path.isfile(fileName)):
                    fileReader = open(fileName, 'r')
                    if(tempRemoveHeaderLine):
                        fileReader.readline()
                    for line in fileReader:
                        fileWriter.write(line)
                    fileReader.close()

        #Close files
        fileWriter.close()
        return True

    ##
    #Tested 4-6-2011
    #Checks if two files are equal.
    #@params String tempFile1 First file
    #@params String tempFile2 Second file
    #@return boolean True (Equal), false (not equal), None (error)
    @staticmethod
    def areEqual(tempFile1 = None, tempFile2 = None):
        #Check parameters
        if(not ValidateData.isValidString(tempFile1)):
            return False
        if(not ValidateData.isValidString(tempFile2)):
            return False

        #Check that the files exist
        if(not os.path.exists(tempFile1)):
            return None
        if(not os.path.exists(tempFile2)):
            return None

        #Open files to read
        file1Reader = open(tempFile1, 'r')
        file2Reader = open(tempFile2, 'r')
        tempEOF = False
        tempEqual = True
        #Stop at an empty string return.
        #Return false and stop on inequality
        line = 0
        while(not tempEOF):
            line = line + 1
            line1 = file1Reader.readline()
            line2 = file2Reader.readline()
            if(line1 == line2):
                if(line1 == ""):
                    tempEOF = True
            else:
                tempEOF = False
                tempEqual = False
        #Close files
        file1Reader.close()
        file2Reader.close()
        #Return
        return tempEqual

    #Created 6-24-2012
    #Tested 6-24-2012
    #Clears the given directory, with no exception no matter the files
    #@params tempDirectory Directroy to be emptied, no take backs!
    #@return Returns boolean to indicate if it was cleared and None on an error.
    @staticmethod
    def clearDirectory(tempDirectory = None):
        #Check if is valid parameter
        if(not ValidateData.isValidString(tempDirectory)):
            return False
        #Check if is a directory path
        if(not os.path.exists(tempDirectory)):
            return False
        #Delete directory
        shutil.rmtree(tempDirectory, ignore_errors=True)
        #Make sure the directory is deleted
        if(os.path.exists(tempDirectory)):
            return False
        os.mkdir(tempDirectory)
        return True
    
    #Tested 6-24-2012
    #Gives the line count of a file, counts blank lines as lines
    #@params tempDirectory File to be counted
    #@return Returns Count of file lines or false on error.
    @staticmethod
    def fileLineCount(tempDirectory = None):
        #Check if is valid parameter
        if(not ValidateData.isValidString(tempDirectory)):
            return False
        #Check if is a directory path
        if(not os.path.exists(tempDirectory)):
            return False
        #Count lines
        i = -1
        f = open(tempDirectory)
        for i, l in enumerate(f):
            pass
        #Close file and return
        f.close()
        return i + 1
