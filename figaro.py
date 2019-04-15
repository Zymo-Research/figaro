import logging
import figaroSupport


def getApplicationParameters():
    import sys
    if len(sys.argv) > 1:
        return getApplicationParametersFromCommandLine()
    parameters = figaroSupport.environmentParameterParser.EnvParameters()
    parameters.addParameter("outputFileName", str, default=default.outputFileName, externalValidation=True)
    parameters.addParameter("ampliconLength", int, lowerBound=0, required=True)
    parameters.addParameter("forwardPrimerLength", int, required=True, lowerBound=0, upperBound=50)
    parameters.addParameter("reversePrimerLength", int, required=True, lowerBound=0, upperBound=50)
    parameters.addParameter("inputDirectory", str, default=default.inputFolder, expectedDirectory=True)
    parameters.addParameter("outputDirectory", str, default = default.outputFolder, expectedDirectory=True)
    parameters.addParameter("minimumOverlap", int, default=default.minOverlap, lowerBound=5, upperBound=30)
    parameters.addParameter("subsample", int, default=default.subsample, lowerBound=-1)
    parameters.addParameter("percentile", int, default = default.percentile, lowerBound=1, upperBound=100)
    parameters.checkCreatedFileStructures()
    combinedReadLengths = parameters.ampliconLength.value + parameters.minimumOverlap.value
    parameters.sideLoadParameter("minimumCombinedReadLength", combinedReadLengths)
    for character in parameters.outputFileName.value:
        if character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_.-":
            raise ValueError("Unusual character detected for output file name.  Contains %s" %character)
    if parameters.subsample.value == -1:
        totalFileSize = figaroSupport.fastqAnalysis.getEstimatedFastqSizeSumFromDirectory(parameters.inputDirectory.value)
        fastqGigabytes = totalFileSize / 1000000000
        parameters.subsample.value = round(fastqGigabytes * 10)
    return parameters


def getApplicationParametersFromCommandLine():
    import argparse
    import os
    # parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outputDirectory", help = "Directory for outputs", default = os.getcwd())
    parser.add_argument("-a", "--ampliconLength", help = "Length of amplicon (not including primers)", required=True, type=int)
    parser.add_argument("-f", "--forwardPrimerLength", help = "Length of forward primer", required=True, type=int)
    parser.add_argument("-r", "--reversePrimerLength", help = "Length of reverse primer", required=True, type=int)
    parser.add_argument("-i", "--inputDirectory", help = "Directory with Fastq files to analyze", default = os.getcwd())
    parser.add_argument("-n", "--outputFileName", help = "Output file for trim site JSON", default=default.outputFileName)
    parser.add_argument("-m", "--minimumOverlap", help = "Minimum overlap between the paired-end reads", default=default.minOverlap, type=int)
    parser.add_argument("-s", "--subsample", help = "Subsampling level (will analyze approximately 1/x reads", default=default.subsample, type=int)
    parser.add_argument("-p", "--percentile", help = "Percentile to use for expected error model", default=default.percentile, type=int)
    args = parser.parse_args()
    # validate args
    outputFileName = args.outputFileName
    ampliconLength = args.ampliconLength
    if not ampliconLength > 0:
        raise ValueError("Amplicon length must be a positive integer. %s was given" %ampliconLength)
    forwardPrimerLength = args.forwardPrimerLength
    if not forwardPrimerLength > 0:
        raise ValueError("Forward primer length must be a positive integer. %s was given" %forwardPrimerLength)
    reversePrimerLength = args.reversePrimerLength
    if not reversePrimerLength > 0:
        raise ValueError("Reverse primer length must be a positive integer. %s was given" %reversePrimerLength)
    inputDirectory = args.inputDirectory
    if not os.path.isdir(inputDirectory):
        raise NotADirectoryError("Unable to find input directory at %s" %inputDirectory)
    outputDirectory = args.outputDirectory
    if not os.path.isdir(outputDirectory):
        raise NotADirectoryError("Unable to find output directory at %s" %outputDirectory)
    minimumOverlap = args.minimumOverlap
    if not minimumOverlap > 0:
        raise ValueError("Minimum overlap must be a positive integer. %s was given." %minimumOverlap)
    subsample = args.subsample
    if subsample < 0:
        totalFileSize = figaroSupport.fastqAnalysis.getEstimatedFastqSizeSumFromDirectory(inputDirectory)
        fastqGigabytes = totalFileSize / 1000000000
        subsample = round(fastqGigabytes * 10)
    percentile = args.percentile
    if percentile < 0 or percentile > 100:
        raise ValueError("Percentile must be an integer value between 0 and 100. %s was given." %percentile)
    combinedReadLengths = ampliconLength + minimumOverlap
    # side-load args into parameter types
    parameters = figaroSupport.environmentParameterParser.EnvParameters()
    parameters.sideLoadParameter("outputFileName", outputFileName)
    parameters.sideLoadParameter("ampliconLength", ampliconLength)
    parameters.sideLoadParameter("forwardPrimerLength", forwardPrimerLength)
    parameters.sideLoadParameter("reversePrimerLength", reversePrimerLength)
    parameters.sideLoadParameter("inputDirectory", inputDirectory)
    parameters.sideLoadParameter("outputDirectory", outputDirectory)
    parameters.sideLoadParameter("minimumOverlap", minimumOverlap)
    parameters.sideLoadParameter("subsample", subsample)
    parameters.sideLoadParameter("percentile", percentile)
    parameters.sideLoadParameter("minimumCombinedReadLength", combinedReadLengths)
    return parameters


def getLoggingParameters():
    loggingParameters = figaroSupport.environmentParameterParser.EnvParameters()
    loggingParameters.addParameter("logFile", str, default=default.logFile, createdFile=True)
    loggingParameters.addParameter("logLevel", str, default=default.loggingLevel, logLevel=True)
    loggingParameters.addParameter("streamOff", bool, default=False)
    loggingParameters.addParameter("streamLoglevel", str, default=default.loggingLevel, logLevel=True)
    loggingParameters.addParameter("fileLogLevel", str, default=default.loggingLevel, logLevel=True)
    logFilePath = os.path.split(loggingParameters.logFile.value)[0]
    if not os.path.isdir(logFilePath):
        os.makedirs(logFilePath)
    loggingParameters.checkCreatedFileStructures()
    return loggingParameters


def loadDefaultPackage():
    defaultParameters = figaroSupport.environmentParameterParser.EnvParameters()
    defaultParameters.addParameter("defaultPackageName", str, default="standard", externalValidation=True)
    return figaroSupport.defaultParser.loadDefaultModule(defaultParameters.defaultPackageName.value)


def setLogging():
    loggingParameters = getLoggingParameters()
    formatter = logging.Formatter(loggingFormat)
    logStreamHandle = logging.StreamHandler()
    logStreamHandle.setFormatter(formatter)
    if not loggingParameters.streamLogLevel.usingDefaultValue:
        logStreamHandle.setLevel(loggingParameters.streamLogLevel.value)
    else:
        logStreamHandle.setLevel(loggingParameters.logLevel.value)
    logFileHandle = logging.FileHandler(loggingParameters.logFile.value)
    logFileHandle.setFormatter(formatter)
    if not loggingParameters.fileLogLevel.usingDefaultValue:
        logFileHandle.setLevel(loggingParameters.fileLogLevel.value)
    else:
        logFileHandle.setLevel(loggingParameters.logLevel.value)
    logger.addHandler(logFileHandle)
    if not loggingParameters.streamOff:
        logger.addHandler(logStreamHandle)


def makeResultJSON(resultTable:list, indent:int=0):
    import json
    resultDictList = []
    for result in resultTable:
        resultDictList.append(result.toDict())
    return json.dumps(resultDictList, indent=indent)


def saveResultOutput(outputDirectory:str, outputResultTableFileName:str, resultTable:list, forwardCurve, reverseCurve):
    import os
    outputResultTablePath = os.path.join(outputDirectory, outputResultTableFileName)
    outputResultTableFile = open(outputResultTablePath, 'w')
    outputResultTableFile.write(makeResultJSON(resultTable, indent=4))
    outputResultTableFile.close()
    outputForwardCurvePath = None
    outputReverseCurvePath = None
    if forwardCurve.curvePNG:
        import base64
        outputForwardCurvePath = os.path.join(outputDirectory, "forwardExpectedError.png")
        outputForwardCurveFile = open(outputForwardCurvePath, 'wb')
        outputForwardCurveFile.write(base64.b64decode(forwardCurve.curvePNG))
        outputForwardCurveFile.close()
    if reverseCurve.curvePNG:
        import base64
        outputReverseCurvePath = os.path.join(outputDirectory, "reverseExpectedError.png")
        outputReverseCurveFile = open(outputReverseCurvePath, 'wb')
        outputReverseCurveFile.write(base64.b64decode(reverseCurve.curvePNG))
        outputReverseCurveFile.close()
    return outputResultTablePath, outputForwardCurvePath, outputReverseCurvePath


if __name__ == "__main__":
    import datetime
    import os
    startTime = datetime.datetime.now()
    default = loadDefaultPackage()
    loggingFormat = "%(levelname)s:%(name)s:%(message)s"
    logger = logging.getLogger(__name__)
    logger.setLevel(
        logging.DEBUG)  # Do not change this line unless you know exactly what you are doing any why you are doing it. This will mess up logging in a way that can be hard to trace back.
    setLogging()
    parameters = getApplicationParameters()
    logger.debug("Starting analysis")
    resultTable, forwardCurve, reverseCurve = figaroSupport.trimParameterPrediction.performAnalysisLite(parameters.inputDirectory.value, parameters.minimumCombinedReadLength.value, subsample =  parameters.subsample.value, percentile = parameters.percentile.value, forwardPrimerLength=parameters.forwardPrimerLength.value, reversePrimerLength=parameters.reversePrimerLength.value)
    for result in resultTable:
        print(result)
    resultTableFileName = os.path.join(parameters.outputDirectory.value, parameters.outputFileName.value)
    saveResultOutput(parameters.outputDirectory.value, parameters.outputFileName.value, resultTable, forwardCurve, reverseCurve)
    print("Run time: %s" %(datetime.datetime.now() - startTime))
    exit(0)