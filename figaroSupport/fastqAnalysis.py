import logging
logger = logging.getLogger(__name__)
from . import fileNamingStandards

def buildQualityMatrix(path:str):
    import numpy
    from .fastqHandler import FastqFile
    fastq = FastqFile(path, depth=1)
    qualityMatrix = []
    for read in fastq:
        qualityMatrix.append(read.quality.phredScores)
    fastq.close()
    return numpy.matrix(qualityMatrix, dtype='uint8') #Memory efficient, but if someone feeds in a phred score > 255, this will break. PacBio, I'm looking at you.


def buildQualityMatrixPaired(forward:str, reverse:str):
    return buildQualityMatrix(forward), buildQualityMatrix(reverse)


def buildExpectedErrorMatrix(path:str, superLean:bool = False, startPosition:int = 0, subsample:int=0, leftTrim:int=0, rightTrim:int=0):
    import numpy
    from . import qualityScoreHandler
    from .fastqHandler import FastqFile
    fastq = FastqFile(path, depth = 0, subsample = subsample, leftTrim=leftTrim, rightTrim=rightTrim)
    expectedErrorMatrix = []
    dataType = 'float16'
    if superLean:
        dataType = 'uint8'
    for line in fastq:
        expectedErrorLineList = qualityScoreHandler.cumulativeExpectedErrorArray(line.quality, fastq.qualityScoreScheme)[startPosition:]
        expectedErrorMatrix.append(expectedErrorLineList)  #low precision floating point. Usually users are looking for whole numbers anyway
    fastq.close()
    return numpy.array(expectedErrorMatrix, dataType, order='F')


def buildExpectedErrorMatrixPaired(forward:str, reverse:str, superLean:bool = False, startPositions:tuple = (0, 0), subsample:int=0):
    return buildExpectedErrorMatrix(forward, superLean, startPositions[0]), buildExpectedErrorMatrix(reverse, superLean, startPositions[1])


def findCutoffByPercentile(path:str, phredScore:int, percentile:int):
    '''
    This will analyze a fastq file to find where the given percentile of reads is at or below the given phred score (such as finding the read where the 10th percentile of reads is phred=10.
    Value returned is the position *INDEXED TO ZERO*
    :param path: path of the Fastq to analyze
    :param phredScore:  score to use in cutoff
    :param percentile:  percentile to use in cutoff
    :return:base position (integer)
    '''
    import numpy
    qualityMatrix = buildQualityMatrix(path).transpose() #faster calclation of percentiles if we have positions as rows and reads as columns
    for position, row in enumerate(qualityMatrix):
        nthPercentile = numpy.percentile(row, percentile)
        if nthPercentile < percentile:
            return position
    return numpy.size(qualityMatrix, 0)


def makeQualityMatrix(path:str):
    import numpy
    from . import fastqHandler
    readLength, variance = fastqHandler.estimateReadLength(path, getVariance=True)
    if variance != 0:
        readLength = fastqHandler.getLongestReadInFile(path)
    fastq = fastqHandler.FastqFile(path, depth=1)
    qualityRange = fastq.qualityScoreScheme.range
    readLengthMatrix = [0] * readLength
    qualityCountMatrix = []
    for i in range(qualityRange + 1):
        qualityCountMatrix.append(readLengthMatrix.copy())
    '''
    Building a matrix here where the correspond to all possibly quality scores and columns represent each base position of each read (indexed to zero)
    Calling a specific value is done by qualityMatrix[qualityScore][readPosition]
    '''
    for read in fastq:
        for position, phred in enumerate(read.quality.phredScores):
            qualityCountMatrix[phred][position] = qualityCountMatrix[phred][position] + 1
    fastq.close()
    qualityCountMatrix = numpy.matrix(qualityCountMatrix)
    return qualityCountMatrix
    # plt.imshow(qualityCountMatrix, origin='lower', aspect='auto')
    # plt.xlabel("Position")
    # plt.ylabel("Quality (Phred)")
    # plt.title("Read quality for %s" %path)
    # if not testingOnly:
    #     if outputFile:
    #         plt.savefig(outputFile)
    #     else:
    #         plt.show()
    # return qualityCountMatrix


def makeAverageExpectedErrorLine(path:str):
    import numpy
    expectedErrorMatrix = buildExpectedErrorMatrix(path)
    expectedErrorMatrix = expectedErrorMatrix.transpose()
    means = []
    for line in expectedErrorMatrix:
        means.append(numpy.mean(line))
    return means
    # plt.plot(means, 'k-')
    # plt.xlabel("Position")
    # plt.ylabel("Average Expected Error")
    # plt.show()


def getDataForFastqPlots(forwardFastq:fileNamingStandards.NamingStandard, reverseFastq:fileNamingStandards.NamingStandard = None):
    forwardQualityMatrix = makeQualityMatrix(forwardFastq.filePath)
    forwardExpectedErrorLine = makeAverageExpectedErrorLine(forwardFastq.filePath)
    if reverseFastq is None:
        reverseQualityMatrix = None
        reverseExpectedErrorLine = None
    else:
        reverseQualityMatrix = makeQualityMatrix(reverseFastq.filePath)
        reverseExpectedErrorLine = makeAverageExpectedErrorLine(reverseFastq.filePath)
    return forwardQualityMatrix, reverseQualityMatrix, forwardExpectedErrorLine, reverseExpectedErrorLine


def generateFastqPlotPaired(forwardFastq:fileNamingStandards.NamingStandard, reverseFastq:fileNamingStandards.NamingStandard, sampleTitle:str = None, outputFile:str = None, base64Format:str = None):
    import matplotlib.pyplot as plt
    if base64Format:
        import base64
    if outputFile and base64Format:
        outputFileFormat = outputFile.split(".")[-1]
        if not outputFileFormat == base64Format:
            logger.error(
                "Cannot save plot in one format and return base64 in a different format.  Returning file save format.  Save in %s. Return base64 %s" % (
                outputFileFormat, base64Format))
    if sampleTitle is None:
        sampleTitle = " ".join([str(item) for item in forwardFastq.sampleID])
    else:
        sampleTitle = str(sampleTitle)
    forwardQualityMatrix, reverseQualityMatrix, forwardExpectedErrorLine, reverseExpectedErrorLine = getDataForFastqPlots(forwardFastq, reverseFastq)
    plt.suptitle("Analysis of %s" % sampleTitle, horizontalalignment="center", fontsize=18, fontweight="bold")

    #make plots for forward reads
    plt.subplot(221)
    plt.imshow(forwardQualityMatrix, origin='lower', aspect='auto')
    plt.xlabel("Read 1 Position")
    plt.ylabel("Quality (Phred)")
    plt.title(" ", fontsize = 16) #making a whitespace buffer
    plt.subplot(222)
    plt.plot(forwardExpectedErrorLine, 'k-')
    plt.xlabel("Read 1 Position")
    plt.ylabel("Average Expected Error")
    plt.title(" ", fontsize = 16) #making a whitespace buffer

    #make plots for reverse reads
    plt.subplot(223)
    plt.imshow(reverseQualityMatrix, origin='lower', aspect='auto')
    plt.xlabel("Read 2 Position")
    plt.ylabel("Quality (Phred)")
    #plt.title("Read quality for %s" %reverseFastq.fileName)
    plt.subplot(224)
    plt.plot(reverseExpectedErrorLine, 'k-')
    plt.xlabel("Read 2 Position")
    plt.ylabel("Average Expected Error")
    #plt.title("Expected error for %s" % reverseFastq.fileName)

    plt.tight_layout()
    if outputFile:
        plt.savefig(outputFile)
        if base64Format:
            imageFile = open(outputFile)
            encodedFile = base64.b64encode(imageFile.read())
            imageFile.close()
            return encodedFile
    elif base64Format:
        import io
        byteStream = io.BytesIO()
        plt.savefig(byteStream, format=base64Format)
        byteStream.seek(0)
        encodedFile = base64.b64encode(byteStream.read())
        return encodedFile
    else:
        plt.show()


def generateFastqPlotSingle(forwardFastq: fileNamingStandards.NamingStandard, sampleTitle: str = None, outputFile: str = None, base64Format:str = None):
    import matplotlib.pyplot as plt
    if base64Format:
        import base64
    if outputFile and base64Format:
        outputFileFormat = outputFile.split(".")[-1]
        if not outputFileFormat == base64Format:
            logger.error("Cannot save plot in one format and return base64 in a different format.  Returning file save format.  Save in %s. Return base64 %s" %(outputFileFormat, base64Format))
    if sampleTitle is None:
        sampleTitle = " ".join([str(item) for item in forwardFastq.sampleID])
    else:
        sampleTitle = str(sampleTitle)
    forwardQualityMatrix, reverseQualityMatrix, forwardExpectedErrorLine, reverseExpectedErrorLine = getDataForFastqPlots(forwardFastq)
    plt.suptitle(sampleTitle, horizontalalignment="center", fontsize = 18, fontweight = "bold")

    # make plots for reads
    plt.subplot(211)
    plt.imshow(forwardQualityMatrix, origin='lower', aspect='auto')
    plt.xlabel("Position")
    plt.ylabel("Quality (Phred)")
    plt.title(" ", fontsize = 16) #making a whitespace buffer
    plt.subplot(212)
    plt.plot(forwardExpectedErrorLine, 'k-')
    plt.xlabel("Position")
    plt.ylabel("Average Expected Error")

    plt.tight_layout()
    if outputFile:
        plt.savefig(outputFile)
        if base64Format:
            imageFile = open(outputFile)
            encodedFile = base64.b64encode(imageFile.read())
            imageFile.close()
            return encodedFile
    elif base64Format:
        import io
        byteStream = io.BytesIO()
        plt.savefig(byteStream, format=base64Format)
        byteStream.seek(0)
        encodedFile = base64.b64encode(byteStream.read())
        return encodedFile
    else:
        plt.show()


class ParallelPlotAgent(object):

    def __init__(self, outputDirectory:str = None, base64Output:bool = False, outputFormat:str = None):
        self.outputDirectory = outputDirectory
        self.outputFormat = outputFormat
        self.base64Output = base64Output
        if outputDirectory or base64Output:
            if not outputFormat:
                raise ValueError("If output to file (directory) or base64 is set, an output format must be provided, but none was.")

    def parallelPlotter(self, fastq:[tuple, fileNamingStandards.NamingStandard]):
        import os
        if type(fastq) == tuple:
            sampleName =  "_".join([str(item) for item in fastq[0].sampleID])
            returnFastq = fastq[0]
        else:
            sampleName =  "_".join([str(item) for item in fastq.sampleID])
            returnFastq = fastq
        if self.outputDirectory:
            outputFileName = os.path.join(self.outputDirectory, sampleName + ".%s" %self.outputFormat)
        else:
            outputFileName = None
        if self.base64Output:
            base64Format = self.outputFormat
        else:
            base64Format = None
        if type(fastq) == tuple:
            base64EncodedPlot = generateFastqPlotPaired(fastq[0], fastq[1], outputFile=outputFileName, base64Format=base64Format)
        else:
            base64EncodedPlot = generateFastqPlotSingle(fastq, outputFile=outputFileName, base64Format=base64Format)
        return returnFastq, outputFileName, base64EncodedPlot #returnValue will be None unless a base64 encoded image was returned


def plotFastqFilesInFolder(directory:str, namingStandard:fileNamingStandards.NamingStandard, outputDirectory:str = None, base64Output:bool = False, outputFormat:str = None):
    import os
    from . import fastqHandler
    from ... import easyMultiprocessing
    if outputDirectory and not os.path.isdir(directory):
        raise NotADirectoryError("Unable to find a directory at %s" %directory)
    if outputDirectory and not os.path.isdir(outputDirectory):
        raise NotADirectoryError("Unable to find a directory at %s" %outputDirectory)
    if outputDirectory or base64Output:
        if not outputFormat:
            raise ValueError(
                "If output to file (directory) or base64 is set, an output format must be provided, but none was.")
    fastqTable = fastqHandler.getSamplePairTableFromFolder(directory, namingStandard)
    fastqSetList = []
    for key in fastqTable:
        if key == "unpaired":
            for fastq in fastqTable["unpaired"]:
                fastqSetList.append(fastq)
        else:
            fastqSetList.append(fastqTable[key])
    parallelPlotAgent = ParallelPlotAgent(outputDirectory=outputDirectory, base64Output=base64Output, outputFormat=outputFormat)
    if outputDirectory or base64Output:
        plotReturnValues = easyMultiprocessing.parallelProcessRunner(parallelPlotAgent.parallelPlotter, fastqSetList)
    else:
        plotReturnValues = [parallelPlotAgent.parallelPlotter(fastq) for fastq in fastqSetList]  #can't do parallel plotting if plotting to a display window
    returnTable = {}
    if outputDirectory and base64Output:
        for fastq, outputFile, base64EncodedPlot in plotReturnValues:
            returnTable[fastq] = (outputFile, base64EncodedPlot)
    elif base64Output:
        for fastq, outputFile, base64EncodedPlot in plotReturnValues:
            returnTable[fastq] = base64EncodedPlot
    elif outputDirectory:
        for fastq, outputFile, base64EncodedPlot in plotReturnValues:
            returnTable[fastq] = base64EncodedPlot
    return returnTable


def getEstimatedFastqFileSizeSumFromList(fastqList:list):
    import os
    from . import gzipIdentifier
    sum = 0
    for fastq in fastqList:
        fileSize = os.path.getsize(fastq.filePath)
        if gzipIdentifier.isGzipped(fastq.filePath):
            fileSize = round(fileSize * 3.5)  #best estimation without doing anything that will slow us down
        sum += fileSize
    return sum

def getEstimatedFastqSizeSumFromDirectory(path:str):
    import os
    from . import fastqHandler
    if not os.path.isdir(path):
        raise NotADirectoryError("Unable to find a directory at %s" %path)
    fastqList = fastqHandler.findSamplesInFolder(path)
    return getEstimatedFastqFileSizeSumFromList(fastqList)



