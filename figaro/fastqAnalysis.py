import logging

logger = logging.getLogger(__name__)


def buildQualityMatrix(path: str):
    import numpy
    from .fastqHandler import FastqFile

    fastq = FastqFile(path, depth=1)
    qualityMatrix = []
    for read in fastq:
        qualityMatrix.append(read.quality.phredScores)
    fastq.close()
    return numpy.matrix(
        qualityMatrix, dtype="uint8"
    )  # Memory efficient, but if someone feeds in a phred score > 255, this will break. PacBio, I'm looking at you.


def buildQualityMatrixPaired(forward: str, reverse: str):
    return buildQualityMatrix(forward), buildQualityMatrix(reverse)


def buildExpectedErrorMatrix(
    path: str,
    superLean: bool = False,
    startPosition: int = 0,
    subsample: int = 0,
    leftTrim: int = 0,
    rightTrim: int = 0,
):
    import numpy

    try:
        from . import qualityScoreHandler
        from .fastqHandler import FastqFile
    except ImportError:
        import qualityScoreHandler
        from fastqHandler import FastqFile
    fastq = FastqFile(
        path, depth=0, subsample=subsample, leftTrim=leftTrim, rightTrim=rightTrim
    )
    expectedErrorMatrix = []
    dataType = "float16"
    if superLean:
        dataType = "uint8"
    for line in fastq:
        expectedErrorLineList = qualityScoreHandler.cumulativeExpectedErrorArray(
            line.quality, fastq.qualityScoreScheme
        )[startPosition:]
        expectedErrorMatrix.append(
            expectedErrorLineList
        )  # low precision floating point. Usually users are looking for whole numbers anyway
    fastq.close()
    return numpy.array(expectedErrorMatrix, dataType, order="F")


def buildExpectedErrorMatrixPaired(
    forward: str,
    reverse: str,
    superLean: bool = False,
    startPositions: tuple = (0, 0),
    subsample: int = 0,
):
    return buildExpectedErrorMatrix(
        forward, superLean, startPositions[0]
    ), buildExpectedErrorMatrix(reverse, superLean, startPositions[1])


def findCutoffByPercentile(path: str, phredScore: int, percentile: int):
    """
    This will analyze a fastq file to find where the given percentile of reads is at or below the given phred score (such as finding the read where the 10th percentile of reads is phred=10.
    Value returned is the position *INDEXED TO ZERO*
    :param path: path of the Fastq to analyze
    :param phredScore:  score to use in cutoff
    :param percentile:  percentile to use in cutoff
    :return:base position (integer)
    """
    import numpy

    qualityMatrix = buildQualityMatrix(
        path
    ).transpose()  # faster calclation of percentiles if we have positions as rows and reads as columns
    for position, row in enumerate(qualityMatrix):
        nthPercentile = numpy.percentile(row, percentile)
        if nthPercentile < percentile:
            return position
    return numpy.size(qualityMatrix, 0)


def makeQualityMatrix(path: str):
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
    """
    Building a matrix here where the correspond to all possibly quality scores and columns represent each base position of each read (indexed to zero)
    Calling a specific value is done by qualityMatrix[qualityScore][readPosition]
    """
    for read in fastq:
        for position, phred in enumerate(read.quality.phredScores):
            qualityCountMatrix[phred][position] = (
                qualityCountMatrix[phred][position] + 1
            )
    fastq.close()
    qualityCountMatrix = numpy.matrix(qualityCountMatrix)
    return qualityCountMatrix


def makeAverageExpectedErrorLine(path: str):
    import numpy

    expectedErrorMatrix = buildExpectedErrorMatrix(path)
    expectedErrorMatrix = expectedErrorMatrix.transpose()
    means = []
    for line in expectedErrorMatrix:
        means.append(numpy.mean(line))
    return means


def getEstimatedFastqFileSizeSumFromList(fastqList: list):
    import os

    try:
        from . import gzipIdentifier
    except ImportError:
        import gzipIdentifier
    sum = 0
    for fastq in fastqList:
        fileSize = os.path.getsize(fastq.filePath)
        if gzipIdentifier.isGzipped(fastq.filePath):
            fileSize = round(
                fileSize * 3.5
            )  # best estimation without doing anything that will slow us down
        sum += fileSize
    return sum


def getEstimatedFastqSizeSumFromDirectory(path: str, fileNamingStandardAlias: str):
    import os

    try:
        from . import fileNamingStandards
        from . import fastqHandler
    except ImportError:
        import fileNamingStandards, fastqHandler
    fileNamingStandard = fileNamingStandards.loadNamingStandard(fileNamingStandardAlias)
    if not os.path.isdir(path):
        raise NotADirectoryError("Unable to find a directory at %s" % path)
    fastqList = fastqHandler.findSamplesInFolder(path, fileNamingStandard)
    return getEstimatedFastqFileSizeSumFromList(fastqList)
