#!/usr/bin/env python

import logging

try:
    from . import (
        environmentParameterParser,
        fileNamingStandards,
        fastqAnalysis,
        trimParameterPrediction,
    )
except ImportError:
    import environmentParameterParser, fileNamingStandards, fastqAnalysis, trimParameterPrediction

try:
    from figaro.defaults import standard as default
except ImportError:
    import defaults.standard as default


def getApplicationParameters():
    import sys

    if len(sys.argv) > 1:
        return getApplicationParametersFromCommandLine()
    parameters = environmentParameterParser.EnvParameters()
    parameters.addParameter(
        "outputFileName", str, default=default.outputFileName, externalValidation=True
    )
    parameters.addParameter("ampliconLength", int, lowerBound=0, required=True)
    parameters.addParameter(
        "forwardPrimerLength", int, required=True, lowerBound=0, upperBound=50
    )
    parameters.addParameter(
        "reversePrimerLength", int, required=True, lowerBound=0, upperBound=50
    )
    parameters.addParameter(
        "inputDirectory", str, default=default.inputFolder, expectedDirectory=True
    )
    parameters.addParameter(
        "outputDirectory", str, default=default.outputFolder, expectedDirectory=True
    )
    parameters.addParameter(
        "minimumOverlap", int, default=default.minOverlap, lowerBound=5, upperBound=30
    )
    parameters.addParameter("subsample", int, default=default.subsample, lowerBound=-1)
    parameters.addParameter(
        "percentile", int, default=default.percentile, lowerBound=1, upperBound=100
    )
    parameters.addParameter(
        "fileNamingStandard", str, default="nononsense", externalValidation=True
    )
    parameters.checkCreatedFileStructures()
    if (
        not parameters.fileNamingStandard.value.lower()
        in fileNamingStandards.aliasList.keys()
    ):
        raise ValueError(
            "%s is not a valid naming standard alias"
            % parameters.fileNamingStandard.value
        )
    combinedReadLengths = (
        parameters.ampliconLength.value + parameters.minimumOverlap.value
    )
    parameters.sideLoadParameter("minimumCombinedReadLength", combinedReadLengths)
    for character in parameters.outputFileName.value:
        if (
            character
            not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_.-"
        ):
            raise ValueError(
                "Unusual character detected for output file name.  Contains %s"
                % character
            )
    if parameters.subsample.value == -1:
        totalFileSize = fastqAnalysis.getEstimatedFastqSizeSumFromDirectory(
            parameters.inputDirectory.value, parameters.fileNamingStandard.value
        )
        fastqGigabytes = totalFileSize / 1000000000
        parameters.subsample.value = round(fastqGigabytes * 10)
    return parameters


def parseArgs():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--outputDirectory", help="Directory for outputs", default=os.getcwd()
    )
    parser.add_argument(
        "-a",
        "--ampliconLength",
        help="Length of amplicon (not including primers)",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-f",
        "--forwardPrimerLength",
        help="Length of forward primer",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-r",
        "--reversePrimerLength",
        help="Length of reverse primer",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-i",
        "--inputDirectory",
        help="Directory with Fastq files to analyze",
        default=os.getcwd(),
    )
    parser.add_argument(
        "-n",
        "--outputFileName",
        help="Output file for trim site JSON",
        default=default.outputFileName,
    )
    parser.add_argument(
        "-m",
        "--minimumOverlap",
        help="Minimum overlap between the paired-end reads",
        default=default.minOverlap,
        type=int,
    )
    parser.add_argument(
        "-s",
        "--subsample",
        help="Subsampling level (will analyze approximately 1/x reads",
        default=default.subsample,
        type=int,
    )
    parser.add_argument(
        "-p",
        "--percentile",
        help="Percentile to use for expected error model",
        default=default.percentile,
        type=int,
    )
    parser.add_argument(
        "-F",
        "--fileNamingStandard",
        help="File naming standard to use",
        default="nononsense",
    )
    parser.add_argument("-l", "--logFile", help="Log file path", default=None)
    return parser.parse_args()


def getApplicationParametersFromCommandLine():
    import os

    args = parseArgs()
    outputFileName = args.outputFileName
    ampliconLength = args.ampliconLength
    if not args.fileNamingStandard.lower() in fileNamingStandards.aliasList.keys():
        raise ValueError(
            "%s is not a valid naming standard alias" % args.fileNamingStandard
        )
    fileNamingStandard = args.fileNamingStandard
    if not ampliconLength > 0:
        raise ValueError(
            "Amplicon length must be a positive integer. %s was given" % ampliconLength
        )
    forwardPrimerLength = args.forwardPrimerLength
    if not forwardPrimerLength >= 0:
        raise ValueError(
            "Forward primer length must be a positive integer. %s was given"
            % forwardPrimerLength
        )
    reversePrimerLength = args.reversePrimerLength
    if not reversePrimerLength >= 0:
        raise ValueError(
            "Reverse primer length must be a positive integer. %s was given"
            % reversePrimerLength
        )
    inputDirectory = args.inputDirectory
    if not os.path.isdir(inputDirectory):
        raise NotADirectoryError(
            "Unable to find input directory at %s" % inputDirectory
        )
    outputDirectory = args.outputDirectory
    if not os.path.isdir(outputDirectory):
        raise NotADirectoryError(
            "Unable to find output directory at %s" % outputDirectory
        )
    minimumOverlap = args.minimumOverlap
    if not minimumOverlap > 0:
        raise ValueError(
            "Minimum overlap must be a positive integer. %s was given." % minimumOverlap
        )
    subsample = args.subsample
    if subsample < 0:
        totalFileSize = fastqAnalysis.getEstimatedFastqSizeSumFromDirectory(
            inputDirectory, fileNamingStandard
        )
        fastqGigabytes = totalFileSize / 1000000000
        subsample = round(fastqGigabytes * 10)
    percentile = args.percentile
    if percentile < 0 or percentile > 100:
        raise ValueError(
            "Percentile must be an integer value between 0 and 100. %s was given."
            % percentile
        )
    combinedReadLengths = ampliconLength + minimumOverlap
    # side-load args into parameter types
    parameters = environmentParameterParser.EnvParameters()
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
    parameters.sideLoadParameter("fileNamingStandard", fileNamingStandard)
    return parameters


def getLoggingParameters():
    import sys
    import os

    loggingParameters = environmentParameterParser.EnvParameters()
    if len(sys.argv) > 1:
        args = parseArgs()
        if args.logFile:
            loggingParameters.sideLoadParameter("logFile", args.logFile)
        else:
            import datetime

            timestamp = str(datetime.datetime.now().timestamp()).replace(".", "")
            loggingParameters.sideLoadParameter(
                "logFile",
                os.path.join(args.outputDirectory, "figaro.%s.log" % timestamp),
            )
    else:
        loggingParameters.addParameter(
            "logFile", str, default=default.logFile, createdFile=True
        )
    loggingParameters.addParameter(
        "logLevel", str, default=default.loggingLevel, logLevel=True
    )
    loggingParameters.addParameter("streamOff", bool, default=False)
    loggingParameters.addParameter(
        "streamLoglevel", str, default=default.loggingLevel, logLevel=True
    )
    loggingParameters.addParameter(
        "fileLogLevel", str, default=default.loggingLevel, logLevel=True
    )
    logFilePath = os.path.split(loggingParameters.logFile.value)[0]
    if not os.path.isdir(logFilePath):
        os.makedirs(logFilePath)
    loggingParameters.checkCreatedFileStructures()
    return loggingParameters


def setLogging():
    loggingFormat = "%(levelname)s:%(name)s:%(message)s"
    logger = logging.getLogger(__name__)
    logger.setLevel(
        logging.DEBUG
    )  # Do not change this line unless you know exactly what you are doing any why you are doing it. This will mess up logging in a way that can be hard to trace back.
    logger.debug("Starting analysis")
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


def makeResultJSON(resultTable: list, indent: int = 0):
    import json

    resultDictList = []
    for result in resultTable:
        resultDictList.append(result.toDict())
    return json.dumps(resultDictList, indent=indent)


def saveResultOutput(
    outputDirectory: str,
    outputResultTableFileName: str,
    resultTable: list,
    forwardCurve,
    reverseCurve,
):
    import os

    outputResultTablePath = os.path.join(outputDirectory, outputResultTableFileName)
    outputResultTableFile = open(outputResultTablePath, "w")
    outputResultTableFile.write(makeResultJSON(resultTable, indent=4))
    outputResultTableFile.close()
    outputForwardCurvePath = None
    outputReverseCurvePath = None
    if forwardCurve.curvePNG:
        import base64

        outputForwardCurvePath = os.path.join(
            outputDirectory, "forwardExpectedError.png"
        )
        outputForwardCurveFile = open(outputForwardCurvePath, "wb")
        outputForwardCurveFile.write(base64.b64decode(forwardCurve.curvePNG))
        outputForwardCurveFile.close()
    if reverseCurve.curvePNG:
        import base64

        outputReverseCurvePath = os.path.join(
            outputDirectory, "reverseExpectedError.png"
        )
        outputReverseCurveFile = open(outputReverseCurvePath, "wb")
        outputReverseCurveFile.write(base64.b64decode(reverseCurve.curvePNG))
        outputReverseCurveFile.close()
    return outputResultTablePath, outputForwardCurvePath, outputReverseCurvePath


def runAnalysis(
    inputDirectory: str,
    ampliconLength: int,
    forwardPrimerLength: int,
    reversePrimerLength: int,
    minimumOverlap: int = 20,
    fileNamingStandard: str = "nononsense",
    subsample: int = -1,
    percentile: int = 83,
):
    import os

    if not os.path.isdir(inputDirectory):
        raise NotADirectoryError("Unable to find directory at %s" % inputDirectory)
    if subsample == -1:
        totalFileSize = fastqAnalysis.getEstimatedFastqSizeSumFromDirectory(
            inputDirectory, fileNamingStandard
        )
        fastqGigabytes = totalFileSize / 1000000000
        subsample = round(fastqGigabytes * 10)
    resultTable, forwardCurve, reverseCurve = (
        trimParameterPrediction.performAnalysisLite(
            inputDirectory,
            ampliconLength + minimumOverlap,
            subsample=subsample,
            percentile=percentile,
            forwardPrimerLength=forwardPrimerLength,
            reversePrimerLength=reversePrimerLength,
            namingStandardAlias=fileNamingStandard,
        )
    )
    return resultTable, forwardCurve, reverseCurve


def main():
    import datetime
    import os

    startTime = datetime.datetime.now()
    setLogging()
    parameters = getApplicationParameters()
    fileNamingStandard = parameters.fileNamingStandard.value
    resultTable, forwardCurve, reverseCurve = (
        trimParameterPrediction.performAnalysisLite(
            parameters.inputDirectory.value,
            parameters.minimumCombinedReadLength.value,
            subsample=parameters.subsample.value,
            percentile=parameters.percentile.value,
            forwardPrimerLength=parameters.forwardPrimerLength.value,
            reversePrimerLength=parameters.reversePrimerLength.value,
            namingStandardAlias=fileNamingStandard,
        )
    )
    for result in resultTable:
        print(result)
    resultTableFileName = os.path.join(
        parameters.outputDirectory.value, parameters.outputFileName.value
    )
    saveResultOutput(
        parameters.outputDirectory.value,
        parameters.outputFileName.value,
        resultTable,
        forwardCurve,
        reverseCurve,
    )
    print("Run time: %s" % (datetime.datetime.now() - startTime))


if __name__ == "__main__":
    main()
