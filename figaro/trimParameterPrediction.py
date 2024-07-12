import logging

logger = logging.getLogger(__name__)
try:
    from . import fileNamingStandards
    from . import fastqHandler
    from . import fastqAnalysis
    from . import expectedErrorCurve
except ImportError:
    import fileNamingStandards, fastqHandler, fastqAnalysis, expectedErrorCurve
import typing
import numpy


class TrimParameterSet(object):

    __slots__ = [
        "forwardTrimPosition",
        "reverseTrimPosition",
        "forwardMaxExpectedError",
        "reverseMaxExpectedError",
        "readRetention",
        "score",
    ]

    def __init__(
        self,
        forwardTrimPosition: int,
        reverseTrimPosition: int,
        forwardMaxExpectedError: int,
        reverseMaxExpectedError: int,
        readRetention: float,
    ):
        self.forwardTrimPosition = forwardTrimPosition
        self.reverseTrimPosition = reverseTrimPosition
        self.forwardMaxExpectedError = forwardMaxExpectedError
        self.reverseMaxExpectedError = reverseMaxExpectedError
        self.readRetention = readRetention
        self.score = self.calculateScore()

    def calculateScore(self):
        return (self.readRetention * 100) - (
            1
            * (
                ((self.forwardMaxExpectedError - 1) ** 2)
                + ((self.reverseMaxExpectedError - 1) ** 2)
            )
        )

    def toJson(self):
        import json

        valueDict = self.toDict()
        return json.dumps(valueDict)

    def toDict(self):
        valueDict = {
            "trimPosition": (self.forwardTrimPosition, self.reverseTrimPosition),
            "maxExpectedError": (
                self.forwardMaxExpectedError,
                self.reverseMaxExpectedError,
            ),
            "readRetentionPercent": round(100 * self.readRetention, 2),
            "score": self.score,
        }
        return valueDict

    def __str__(self):
        return self.toJson()


def calculateMaxExpectedErrorFromReadLength(readLength: int):
    dividedLength = readLength // 100
    maxExpectedError = 0
    for i in range(1, dividedLength + 1):
        maxExpectedError += i
    return maxExpectedError + 1


def calculateForwardExpectedErrorFromReadLength(readLength: int):
    import math

    calculatedValue = 0.0356 * (math.e ** (0.015 * readLength))
    roundedValue = round(calculatedValue)
    return roundedValue + 1


def calculateReverseExpectedErrorFromReadLength(readLength: int):
    import math

    calculatedValue = 0.0289 * (math.e ** (0.0203 * readLength))
    roundedValue = round(calculatedValue)
    return roundedValue + 1


def getFastqList(
    path: str, namingStandard: typing.Type[fileNamingStandards.NamingStandard]
):
    return fastqHandler.findSamplesInFolder(path, namingStandard)


def calculateLowestTrimBaseForPairedReads(
    forwardLength: int, reverseLength: int, minimumCombinedLength: int
):
    if forwardLength + reverseLength < minimumCombinedLength:
        logger.error(
            "Combined read lengths are less than the required combined length."
        )
        return forwardLength, reverseLength
    minimumForwardLength = minimumCombinedLength - reverseLength
    minimumReverseLength = minimumCombinedLength - forwardLength
    return minimumForwardLength, minimumReverseLength


def makeTrimLocations(
    forwardLength: int,
    reverseLength: int,
    minimumCombinedLength: int,
    numberOfIntermediateLocations: int = 10,
):
    minimumForwardLength, minimumReverseLength = calculateLowestTrimBaseForPairedReads(
        forwardLength, reverseLength, minimumCombinedLength
    )
    potentialTrimSpaceLength = forwardLength - minimumForwardLength
    if potentialTrimSpaceLength < numberOfIntermediateLocations:
        numberOfIntermediateLocations = potentialTrimSpaceLength - 2
    locationList = [(minimumForwardLength - 1, reverseLength - 1)]
    trimIncrement = potentialTrimSpaceLength // (numberOfIntermediateLocations + 1)
    for i in range(1, numberOfIntermediateLocations + 1):
        locationList.append(
            (
                minimumForwardLength - 1 + (i * trimIncrement),
                reverseLength - 1 - (i * trimIncrement),
            )
        )
    locationList.append((forwardLength - 1, minimumReverseLength - 1))
    return tuple(locationList)


def makeAllPossibleTrimLocations(
    forwardLength: int, reverseLength: int, minimumCombinedLength: int
):
    minimumForwardLength, minimumReverseLength = calculateLowestTrimBaseForPairedReads(
        forwardLength, reverseLength, minimumCombinedLength
    )
    forwardPosition = minimumForwardLength - 1
    reversePosition = reverseLength - 1
    trimPositions = []
    while forwardPosition < forwardLength:
        trimPositions.append((forwardPosition, reversePosition))
        forwardPosition += 1
        reversePosition -= 1
    return tuple(trimPositions)


class Q2ArrayParallelBuilderAgent(object):

    def __init__(self, subsample: int = 0, primerLength: int = 0):
        if subsample == 0:
            subsample = 1
        self.subsample = subsample
        self.primerLength = primerLength

    def makeQ2Array(self, fastqFileInfo: fileNamingStandards.NamingStandard):
        import numpy

        # print("Running %s" %fastq)
        fastq = fastqHandler.FastqFile(
            fastqFileInfo.filePath,
            depth=1,
            subsample=self.subsample,
            leftTrim=self.primerLength,
        )
        q2Locations = []
        for read in fastq:
            containedQ2 = False
            for position, qScore in enumerate(read.quality.phredScores):
                if qScore <= 2:
                    containedQ2 = True
                    q2Locations.append(position)
                    break
            if not containedQ2:
                q2Locations.append(len(read.sequence))
        firstQ2Array = numpy.array(q2Locations, "uint16")
        # print("%s Reads: %s. First Q2 Array: %s. First Q2 List: %s" %(fastqFileInfo.fileName, readCount, len(firstQ2Array), len(q2Locations)))
        return fastqFileInfo, firstQ2Array


def makeCombinedQ2ArrayForOneDirection(
    fastqList: list, sampleOrder: list, subsample: int = 0, primerLength: int = 0
):
    import numpy
    from . import easyMultiprocessing

    parallelBuildAgent = Q2ArrayParallelBuilderAgent(subsample, primerLength)
    firstQ2Arrays = easyMultiprocessing.parallelProcessRunner(
        parallelBuildAgent.makeQ2Array, fastqList
    )
    combinedArrayStarted = False
    combinedArray = None
    for array in firstQ2Arrays:
        if array[0].sameSample(sampleOrder[0]):
            combinedArray = array[1]
            combinedArrayStarted = True
            # print("Added %s" %array[0].fileName)
            # print(combinedArray.shape)
            break
    if not combinedArrayStarted:
        raise RuntimeError(
            "Did not find the initial combined matrix for first Q2. This requires debugging as it should not be possible."
        )
    for fastq in sampleOrder[1:]:
        for array in firstQ2Arrays:
            if fastq.sameSample(array[0]):
                combinedArray = numpy.concatenate((combinedArray, array[1]))
                # print("Added %s" % array[0].fileName)
                # print(combinedArray.shape)
                break
    return combinedArray


def makeCombinedQ2ArraysForBothEnds(
    fastqList: list,
    sampleOrder: list,
    subsample: int = 0,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    forwardFastqList = [fastq for fastq in fastqList if fastq.direction == 1]
    reverseFastqList = [fastq for fastq in fastqList if fastq.direction == 2]
    forwardQ2Array = makeCombinedQ2ArrayForOneDirection(
        forwardFastqList, sampleOrder, subsample, forwardPrimerLength
    )
    reverseQ2Array = makeCombinedQ2ArrayForOneDirection(
        reverseFastqList, sampleOrder, subsample, reversePrimerLength
    )
    # print("First Q2 array sizes:")
    # print("F: %s" %(len(forwardQ2Array)))
    # print("R: %s" %(len(reverseQ2Array)))
    return forwardQ2Array, reverseQ2Array


class NBaseArrayParallelBuilderAgent(object):

    def __init__(self, subsample: int = 0, primerLength: int = 0):
        if subsample == 0:
            subsample = 1
        self.subsample = subsample
        self.primerLength = primerLength

    def makeFirstNBaseArray(self, fastqFileInfo: fileNamingStandards.NamingStandard):
        import numpy

        # print("Running %s" %fastq)
        fastq = fastqHandler.FastqFile(
            fastqFileInfo.filePath, subsample=self.subsample, leftTrim=self.primerLength
        )
        nBaseLocations = []
        for read in fastq:
            containedN = False
            for position, base in enumerate(read.sequence):
                if base == "N":
                    containedN = True
                    nBaseLocations.append(position)
                    break
            if not containedN:
                nBaseLocations.append(len(read.sequence))
        firstNBaseArray = numpy.array(nBaseLocations, "uint16")
        # print("%s Reads: %s. First N Array: %s. First N List: %s" %(fastqFileInfo.fileName, readCount, len(firstNBaseArray), len(nBaseLocations)))
        return fastqFileInfo, firstNBaseArray


def makeCombinedFirstNBaseArrayForOneDirection(
    fastqList: list, sampleOrder: list, subsample: int = 0, primerLength: int = 0
):
    import numpy
    from . import easyMultiprocessing

    parallelBuildAgent = NBaseArrayParallelBuilderAgent(subsample, primerLength)
    firstNBaseArrays = easyMultiprocessing.parallelProcessRunner(
        parallelBuildAgent.makeFirstNBaseArray, fastqList
    )
    combinedArrayStarted = False
    combinedArray = None
    for array in firstNBaseArrays:
        if array[0].sameSample(sampleOrder[0]):
            combinedArray = array[1]
            combinedArrayStarted = True
            # print("Added %s" %array[0].fileName)
            # print(combinedArray.shape)
            break
    if not combinedArrayStarted:
        raise RuntimeError(
            "Did not find the initial combined matrix for first N base. This requires debugging as it should not be possible."
        )
    for fastq in sampleOrder[1:]:
        for array in firstNBaseArrays:
            if fastq.sameSample(array[0]):
                combinedArray = numpy.concatenate((combinedArray, array[1]))
                # print("Added %s" % array[0].fileName)
                # print(combinedArray.shape)
                break
    return combinedArray


def makeCombinedFirstNBaseArraysForBothEnds(
    fastqList: list,
    sampleOrder: list,
    subsample: int,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    forwardFastqList = [fastq for fastq in fastqList if fastq.direction == 1]
    reverseFastqList = [fastq for fastq in fastqList if fastq.direction == 2]
    forwardFirstNBaseArray = makeCombinedFirstNBaseArrayForOneDirection(
        forwardFastqList, sampleOrder, subsample, forwardPrimerLength
    )
    reverseFirstNBaseArray = makeCombinedFirstNBaseArrayForOneDirection(
        reverseFastqList, sampleOrder, subsample, reversePrimerLength
    )
    # print("First N base array sizes:")
    # print("F: %s" %(len(forwardFirstNBaseArray)))
    # print("R: %s" %(len(reverseFirstNBaseArray)))
    return forwardFirstNBaseArray, reverseFirstNBaseArray


class ExpectedErrorMatrixBuilderParallelAgent(object):

    def __init__(
        self, startPosition: int = 0, subsample: int = 0, primerLength: int = 0
    ):
        self.startPosition = startPosition
        self.subsample = subsample
        self.primerLength = primerLength

    def makeExpectedErrorMatrix(self, fastq: fileNamingStandards.NamingStandard):
        # print("Running %s" %fastq)
        expectedErrorMatrix = fastqAnalysis.buildExpectedErrorMatrix(
            fastq.filePath,
            superLean=True,
            startPosition=self.startPosition,
            subsample=self.subsample,
            leftTrim=self.primerLength,
        )
        return fastq, expectedErrorMatrix


def makeCombinedExpectedErrorMatrixForOneDirection(
    fastqList: list,
    sampleOrder: list,
    subsample: int,
    startPosition: int = 0,
    primerLength: int = 0,
):
    import numpy

    try:
        from . import easyMultiprocessing
    except ImportError:
        import easyMultiprocessing
    parallelBuildAgent = ExpectedErrorMatrixBuilderParallelAgent(
        startPosition, subsample, primerLength
    )
    expectedErrorMatrices = easyMultiprocessing.parallelProcessRunner(
        parallelBuildAgent.makeExpectedErrorMatrix, fastqList
    )
    combinedMatrixStarted = False
    combinedMatrix = None
    for matrix in expectedErrorMatrices:
        if matrix[0].sameSample(sampleOrder[0]):
            combinedMatrix = matrix[1]
            combinedMatrixStarted = True
            # print("Added %s" %matrix[0].fileName)
            # print(combinedMatrix.shape)
            break
    if not combinedMatrixStarted:
        raise RuntimeError(
            "Did not find the initial combined matrix. This requires debugging, as it should not be possible."
        )
    for fastq in sampleOrder[1:]:
        for matrix in expectedErrorMatrices:
            if fastq.sameSample(matrix[0]):
                combinedMatrix = numpy.concatenate((combinedMatrix, matrix[1]))
                # print("Added %s" % matrix[0].fileName)
                # print(combinedMatrix.shape)
                break
    # for matrix in expectedErrorMatrices:
    # print("%s, %s" %(matrix[0].fileName, matrix[1].size))
    return combinedMatrix.transpose()  # columns for reads, rows for positions


def makeCombinedErrorMatricesForBothEnds(
    fastqList: list,
    sampleOrder: list,
    subsample: int,
    minimumTrimPositions: tuple = (0, 0),
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    forwardMinimumTrimPosition, reverseMinimumTrimPosition = minimumTrimPositions
    forwardFastqList = [fastq for fastq in fastqList if fastq.direction == 1]
    reverseFastqList = [fastq for fastq in fastqList if fastq.direction == 2]
    forwardExpectedErrorMatrix = makeCombinedExpectedErrorMatrixForOneDirection(
        forwardFastqList,
        sampleOrder,
        subsample,
        forwardMinimumTrimPosition,
        forwardPrimerLength,
    )
    reverseExpectedErrorMatrix = makeCombinedExpectedErrorMatrixForOneDirection(
        reverseFastqList,
        sampleOrder,
        subsample,
        reverseMinimumTrimPosition,
        reversePrimerLength,
    )
    # print("Expected Error Matrix Sizes:")
    # print("F: %s" %(forwardExpectedErrorMatrix.size))
    # print("R: %s" %(reverseExpectedErrorMatrix.size))
    return forwardExpectedErrorMatrix, reverseExpectedErrorMatrix


def padMaxExpectedError(rawValue: float):
    roundedUpValue = -(int(-rawValue))
    return roundedUpValue + 1


def runTrimParameterTest(
    forwardExpectedErrorMatrix: numpy.ndarray,
    reverseExpectedErrorMatrix: numpy.ndarray,
    forwardFirstNBaseArray: numpy.ndarray,
    reverseFirstNBaseArray: numpy.ndarray,
    forwardQ2Array: numpy.ndarray,
    reverseQ2Array: numpy.ndarray,
    trimPositions: tuple,
    minimumTrimPositions: tuple = (0, 0),
    forwardCurve: expectedErrorCurve.ExponentialFit = None,
    reverseCurve: expectedErrorCurve.ExponentialFit = None,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    import operator

    forwardMinimumTrimPosition, reverseMinimumTrimPosition = minimumTrimPositions
    results = []
    for forwardTrimPosition, reverseTrimPosition in trimPositions:
        if not forwardCurve:
            forwardMaxExpectedError = calculateForwardExpectedErrorFromReadLength(
                forwardTrimPosition
            )
        else:
            forwardMaxExpectedError = padMaxExpectedError(
                forwardCurve.calculateValue(forwardTrimPosition)
            )
        if not reverseCurve:
            reverseMaxExpectedError = calculateReverseExpectedErrorFromReadLength(
                reverseTrimPosition
            )
        else:
            reverseMaxExpectedError = padMaxExpectedError(
                reverseCurve.calculateValue(reverseTrimPosition)
            )
        forwardExpectedErrors = forwardExpectedErrorMatrix[
            forwardTrimPosition - forwardMinimumTrimPosition
        ]
        reverseExpectedErrors = reverseExpectedErrorMatrix[
            reverseTrimPosition - reverseMinimumTrimPosition
        ]
        totalReads = 0
        keptReads = 0
        rejectedReads = 0
        for (
            forwardExpectedErrorValue,
            reverseExpectedErrorValue,
            forwardFirstNBasePosition,
            reverseFirstNBasePosition,
            forwardQ2Position,
            reverseQ2Position,
        ) in zip(
            forwardExpectedErrors,
            reverseExpectedErrors,
            forwardFirstNBaseArray,
            reverseFirstNBaseArray,
            forwardQ2Array,
            reverseQ2Array,
        ):
            totalReads += 1
            if (
                forwardExpectedErrorValue >= forwardMaxExpectedError
                or reverseExpectedErrorValue >= reverseMaxExpectedError
            ):  # Using this because I lose fractional values to save on memory. In theory, this would probably disagree with dada2's decision if the real value is exactly an integer with no fractional portion.  This is unlikely to happen often enough to be an issue.
                rejectedReads += 1
                continue
            elif (
                forwardTrimPosition >= forwardFirstNBasePosition
                or reverseTrimPosition >= reverseFirstNBasePosition
            ):
                rejectedReads += 1
                continue
            elif (
                forwardTrimPosition >= forwardQ2Position
                or reverseTrimPosition >= reverseQ2Position
            ):
                rejectedReads += 1
                continue
            else:
                keptReads += 1
        results.append(
            TrimParameterSet(
                forwardTrimPosition + 1 + forwardPrimerLength,
                reverseTrimPosition + 1 + reversePrimerLength,
                forwardMaxExpectedError,
                reverseMaxExpectedError,
                keptReads / totalReads,
            )
        )  # doing +1 to adjust for zero indexed matrices
        results.sort(key=operator.attrgetter("score"), reverse=True)
    return results


def runTrimParameterTestLite(
    forwardExpectedErrorMatrix: numpy.ndarray,
    reverseExpectedErrorMatrix: numpy.ndarray,
    trimPositions: tuple,
    minimumTrimPositions: tuple = (0, 0),
    forwardCurve: expectedErrorCurve.ExponentialFit = None,
    reverseCurve: expectedErrorCurve.ExponentialFit = None,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    import operator

    forwardMinimumTrimPosition, reverseMinimumTrimPosition = minimumTrimPositions
    results = []
    for forwardTrimPosition, reverseTrimPosition in trimPositions:
        if not forwardCurve:
            forwardMaxExpectedError = calculateForwardExpectedErrorFromReadLength(
                forwardTrimPosition
            )
        else:
            forwardMaxExpectedError = padMaxExpectedError(
                forwardCurve.calculateValue(forwardTrimPosition)
            )
        if not reverseCurve:
            reverseMaxExpectedError = calculateReverseExpectedErrorFromReadLength(
                reverseTrimPosition
            )
        else:
            reverseMaxExpectedError = padMaxExpectedError(
                reverseCurve.calculateValue(reverseTrimPosition)
            )
        forwardExpectedErrors = forwardExpectedErrorMatrix[
            forwardTrimPosition - forwardMinimumTrimPosition
        ]
        reverseExpectedErrors = reverseExpectedErrorMatrix[
            reverseTrimPosition - reverseMinimumTrimPosition
        ]
        totalReads = 0
        keptReads = 0
        rejectedReads = 0
        for forwardExpectedErrorValue, reverseExpectedErrorValue in zip(
            forwardExpectedErrors, reverseExpectedErrors
        ):
            totalReads += 1
            if (
                forwardExpectedErrorValue >= forwardMaxExpectedError
                or reverseExpectedErrorValue >= reverseMaxExpectedError
            ):  # Using this because I lose fractional values to save on memory. In theory, this would probably disagree with dada2's decision if the real value is exactly an integer with no fractional portion.  This is unlikely to happen often enough to be an issue.
                rejectedReads += 1
                continue
            else:
                keptReads += 1
        results.append(
            TrimParameterSet(
                forwardTrimPosition + 1 + forwardPrimerLength,
                reverseTrimPosition + 1 + reversePrimerLength,
                forwardMaxExpectedError,
                reverseMaxExpectedError,
                keptReads / totalReads,
            )
        )  # doing +1 to adjust for zero indexed matrices
        results.sort(key=operator.attrgetter("score"), reverse=True)
    return results


def getSampleOrder(fastqList: list):
    forwardFastqList = [fastq for fastq in fastqList if fastq.direction == 1]
    sampleOrder = forwardFastqList  # using this because I have to pick one
    return sampleOrder


def parallelReadLengthChecker(fastq: fileNamingStandards.NamingStandard):
    return fastq, fastqHandler.estimateReadLength(fastq.filePath, getVariance=True)


def checkReadLengths(fastqList: list):
    try:
        from . import easyMultiprocessing
    except ImportError:
        import easyMultiprocessing
    read1Data = []
    read2Data = []
    fastqReadLengthData = easyMultiprocessing.parallelProcessRunner(
        parallelReadLengthChecker, fastqList
    )
    for fastq, data in fastqReadLengthData:
        if fastq.direction == 1:
            read1Data.append(data)
        elif fastq.direction == 2:
            read2Data.append(data)
    read1DataSet = set(read1Data)
    read2DataSet = set(read2Data)
    filesPassCheck = True
    if not len(read1Data) == len(read2Data):
        logger.error(
            "There appears to be a different number of forward and reverse fastq files in the sequence folder. %s forward and %s reverse"
            % (len(read1Data), len(read2Data))
        )
        filesPassCheck = False
    if not len(read1DataSet) == 1:
        logger.error(
            "Forward read files appear to be of different lengths or of varied lengths. %s"
            % read1DataSet
        )
        filesPassCheck = False
    if not len(read2DataSet) == 1:
        logger.error(
            "Reverse read files appear to be of different lengths or of varied lengths. %s"
            % read2DataSet
        )
        filesPassCheck = False
    read1Length, read1Variance = list(read1DataSet)[0]
    read2Length, read2Variance = list(read2DataSet)[0]
    if read1Variance:
        logger.error(
            "Forward reads appear to not be of consistent length. %s" % read1DataSet
        )
        filesPassCheck = False
    if read2Variance:
        logger.error(
            "Reverse reads appear to not be of consistent length. %s" % read2DataSet
        )
        filesPassCheck = False
    if not filesPassCheck:
        raise fastqHandler.FastqValidationError(
            "Unable to validate fastq files enough to perform this operation. Please check log for specific error(s)."
        )
    return read1Length, read2Length


def performAnalysis(
    inputDirectory: str,
    minimumCombinedReadLength: int,
    subsample: int = 0,
    percentile: int = 83,
    fastqList: list = None,
    makeExpectedErrorPlots: bool = True,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    from . import expectedErrorCurve

    if not inputDirectory:
        if not fastqList:
            raise ValueError("No input directory and no fastq list were given.")
    if not fastqList:
        fastqList = getFastqList(inputDirectory)
        if not fastqList:
            raise ValueError("No fastq files found in input directory")
    sampleOrder = getSampleOrder(fastqList)
    forwardReadLength, reverseReadLength = checkReadLengths(fastqList)
    forwardReadLength = forwardReadLength - forwardPrimerLength
    reverseReadLength = reverseReadLength - reversePrimerLength
    forwardCurve, reverseCurve = (
        expectedErrorCurve.calculateExpectedErrorCurvesForFastqList(
            fastqList,
            subsample=subsample,
            percentile=percentile,
            makePNG=makeExpectedErrorPlots,
            forwardPrimerLength=forwardPrimerLength,
            reversePrimerLength=reversePrimerLength,
        )
    )
    minimumTrimmingPositions = calculateLowestTrimBaseForPairedReads(
        forwardReadLength, reverseReadLength, minimumCombinedReadLength
    )
    trimPositions = makeAllPossibleTrimLocations(
        forwardReadLength, reverseReadLength, minimumCombinedReadLength
    )
    forwardQ2Array, reverseQ2Array = makeCombinedQ2ArraysForBothEnds(
        fastqList, sampleOrder, subsample, forwardPrimerLength, reversePrimerLength
    )
    forwardFirstNBaseArray, reverseFirstNBaseArray = (
        makeCombinedFirstNBaseArraysForBothEnds(
            fastqList, sampleOrder, subsample, forwardPrimerLength, reversePrimerLength
        )
    )
    forwardExpectedErrorMatrix, reverseExpectedErrorMatrix = (
        makeCombinedErrorMatricesForBothEnds(
            fastqList,
            sampleOrder,
            subsample,
            minimumTrimmingPositions,
            forwardPrimerLength,
            reversePrimerLength,
        )
    )
    resultTable = runTrimParameterTest(
        forwardExpectedErrorMatrix,
        reverseExpectedErrorMatrix,
        forwardFirstNBaseArray,
        reverseFirstNBaseArray,
        forwardQ2Array,
        reverseQ2Array,
        trimPositions,
        minimumTrimmingPositions,
        forwardCurve,
        reverseCurve,
        forwardPrimerLength,
        reversePrimerLength,
    )
    return resultTable, forwardCurve, reverseCurve


def performAnalysisLite(
    inputDirectory: str,
    minimumCombinedReadLength: int,
    subsample: int = 0,
    percentile: int = 83,
    fastqList: list = None,
    makeExpectedErrorPlots: bool = True,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
    namingStandardAlias: str = "illumina",
):
    try:
        from . import expectedErrorCurve
    except:
        import expectedErrorCurve
    namingStandard = fileNamingStandards.loadNamingStandard(namingStandardAlias)
    if not inputDirectory:
        if not fastqList:
            raise ValueError("No input directory and no fastq list were given.")
    if not fastqList:
        fastqList = getFastqList(inputDirectory, namingStandard)
        if not fastqList:
            raise ValueError("No fastq files found in input directory")
    sampleOrder = getSampleOrder(fastqList)
    forwardReadLength, reverseReadLength = checkReadLengths(fastqList)
    print("Forward read length: %s" % forwardReadLength)
    print("Reverse read length: %s" % reverseReadLength)
    forwardReadLength = forwardReadLength - forwardPrimerLength
    reverseReadLength = reverseReadLength - reversePrimerLength
    forwardCurve, reverseCurve = (
        expectedErrorCurve.calculateExpectedErrorCurvesForFastqList(
            fastqList,
            subsample=subsample,
            percentile=percentile,
            makePNG=makeExpectedErrorPlots,
            forwardPrimerLength=forwardPrimerLength,
            reversePrimerLength=reversePrimerLength,
        )
    )
    minimumTrimmingPositions = calculateLowestTrimBaseForPairedReads(
        forwardReadLength, reverseReadLength, minimumCombinedReadLength
    )
    trimPositions = makeAllPossibleTrimLocations(
        forwardReadLength, reverseReadLength, minimumCombinedReadLength
    )
    forwardExpectedErrorMatrix, reverseExpectedErrorMatrix = (
        makeCombinedErrorMatricesForBothEnds(
            fastqList,
            sampleOrder,
            subsample,
            minimumTrimmingPositions,
            forwardPrimerLength,
            reversePrimerLength,
        )
    )
    resultTable = runTrimParameterTestLite(
        forwardExpectedErrorMatrix,
        reverseExpectedErrorMatrix,
        trimPositions,
        minimumTrimmingPositions,
        forwardCurve,
        reverseCurve,
        forwardPrimerLength,
        reversePrimerLength,
    )
    return resultTable, forwardCurve, reverseCurve
