#######################################################
#
#	Title:		CommandLine.py
#	Author:		Timothy Tickle
#	Date:		10/03/2011
#	Purpose:	Parent class for inheriting methods for command line functionality.
#
#######################################################

#Import libaries
from subprocess import call, Popen, PIPE
from ValidateData import ValidateData
import traceback

##
#Class containing commandline functionality.
class CommandLine():

    ##
    #Contructor
    def __init__(self): pass

    #Sends a command to command line interface
    #@parameter tempCommand Must be an list of command key word and string arguments, no whitespaces
    #Create new array elements instead of white spacing
    #@return boolean indicator of success (True = Success)
    def runCommandLine(self,tempCommand = None):
        print "Command="+str(tempCommand)
        #Makes sure the the input data is a list of strings
        if(not ValidateData.funcIsValidStringList(tempCommand)):
            print "Error:: tempCommand must be an array of strings. Received="+str(tempCommand)+"."
            return False

        #Run command
        try:
            returnCode = call(tempCommand)
            print "Return="+str(returnCode)
            if returnCode > 0:
                print "Error:: Error during command call. Script stopped."
                print "Error:: Error Code "+str(returnCode)+"."
                print "Error:: Command ="+str(tempCommand)+"."
                return False
        except (OSError,TypeError), e: 
                print "Error:: Error during command call. Script stopped."
                print "Error:: Command ="+str(tempCommand)+"."
                print "Error:: OS error: "+str(traceback.format_exc(e))+"."
                return False
        return True

    #Sends a command to command line interface
    #@parameter tempCommand Must be an list of command key word and string arguments, no whitespaces
    #Create new array of string elements instead of white spacing
    #Put file names in escaped quotation marks.
    #This uses shell == true so make sure the commandline is not malicious
    #This should wait for process completion
    #@return False = Failure or the return code from the subprocess
    def runPipedCommandLine(self,tempCommand = None):
        #Makes sure the the input data is a list of strings
        if(not ValidateData.funcIsValidStringList(tempCommand)):
            print "Error:: tempCommand must be an array of strings. Received="+str(tempCommand)+"."
            return False

        #Run command
        tempCommand = " ".join(tempCommand)
        try:
            returnCode = Popen(tempCommand, shell = True, stdout = PIPE).communicate()
            return returnCode
        except (OSError,TypeError), e: 
                print "Error:: Error during command call. Script stopped."
                print "Error:: Command ="+str(tempCommand)+"."
                print "Error:: OS error: "+str(traceback.format_exc(e))+"."
                return False

    #Sends a an array of commands to the commandline
    #@parameter tempCommand Must be an list of commands, parsing and removing whitespace is handled internally
    #Do not send mkdir and rm commands, use the appropriate os.* method call
    #@return boolean indicator of success (True = Success)
    def runBatchCommandline(self,tempArrayOfCommands = None):

        #Holds commands
        parsedCommmands = []

        #Indicates if success or error occured
        success = True

        #Makes sure the the input data is list of strings
        if(not ValidateData.funcIsValidStringList(tempArrayOfCommands)):
            print "Error:: tempCommand must be an array of strings. Received="+str(tempArrayOfCommands)+"."
            return False

        #Parse commands into an array and call
        #On an error break and return False
        for command in tempArrayOfCommands:
            commandElements = command.split(" ")
            if(not self.runCommandLine(commandElements)):
                success = False
                break
        return success
                
