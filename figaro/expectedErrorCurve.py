try:
    from . import fileNamingStandards
    from . import fastqHandler
    from . import fastqAnalysis
except ImportError:
    import fileNamingStandards, fastqHandler, fastqAnalysis
import numpy
import typing
import collections.abc as collections


class ExponentialFit(object):

    __slots__ = ["a", "b", "c", "covariance", "rSquared", "curvePNG"]

    def __init__(
        self,
        a: float,
        b: float,
        c: float,
        covariance: collections.Iterable = None,
        rSquared: float = None,
        curvePNG: str = None,
    ):
        self.a = a
        self.b = b
        self.c = c
        self.covariance = covariance
        self.rSquared = rSquared
        self.curvePNG = curvePNG

    def calculateValue(self, x: float):
        return self.a * numpy.exp(self.b * x) + self.c

    def __str__(self):
        sign = "+"
        if self.c < 0:
            sign = "-"
        return "%.4fe^(%.4fx) %s %.4f" % (self.a, self.b, sign, abs(self.c))

    def __bool__(self):
        return True


def exponentialPrototypeFunction(x, a, b, c):
    return a * numpy.exp(b * x) + c


def fitExponentialCurve(
    xValues: collections.Iterable,
    yValues: collections.Iterable,
    generateImage: bool = False,
    plotName: str = "Expected error by position",
):
    import scipy.optimize
    import scipy.stats

    coefficients, covariance = scipy.optimize.curve_fit(
        exponentialPrototypeFunction,
        xValues,
        yValues,
        p0=(0.03, 0.015, 0),
        bounds=((-2, -1, -8), (2, 1, 8)),
    )
    curve = ExponentialFit(*coefficients, covariance)
    modelPredictions = [curve.calculateValue(x) for x in xValues]
    pearson = scipy.stats.pearsonr(yValues, modelPredictions)
    rSquared = pearson[0] ** 2
    curve = ExponentialFit(*coefficients, covariance, rSquared)
    if generateImage:
        import matplotlib.pyplot as plt
        import io
        import base64

        plt.plot(xValues, yValues, "k-", label="Observed")
        plt.plot(xValues, modelPredictions, "b--", label="Predicted")
        plt.xlabel("Position in Read")
        plt.ylabel("Expected Error")
        plt.title(plotName)
        plt.legend(loc=2)
        text = "%s\nr^2=%.6f" % (curve, rSquared)
        textYPosition = max(modelPredictions) * 0.45
        plt.text(0, textYPosition, text)
        # plt.show()
        byteStream = io.BytesIO()
        plt.savefig(byteStream, format="png")
        byteStream.seek(0)
        encodedImage = base64.b64encode(byteStream.read())
        curve = ExponentialFit(*coefficients, covariance, rSquared, encodedImage)
        plt.clf()
    return curve


class ParallelExpectedErrorAverageAgent(object):

    def __init__(self, subsample: int = 0, primerLength: int = 0):
        if subsample == 0:
            subsample = 1
        self.subsample = subsample
        self.primerLength = primerLength

    def calculateAverageExpectedError(self, fastq: fileNamingStandards.NamingStandard):
        averageExpectedError = makeExpectedErrorAverageArrayForFastq(
            fastq.filePath, self.subsample, self.primerLength
        )
        return fastq, averageExpectedError


class ParallelExpectedErrorPercentileAgent(object):

    def __init__(self, subsample: int = 0, percentile: int = 83, primerLength: int = 0):
        if subsample == 0:
            subsample = 1
        self.subsample = subsample
        self.percentile = percentile
        self.primerLength = primerLength

    def calculateAverageExpectedError(self, fastq: fileNamingStandards.NamingStandard):
        percentileExpectedError = makeExpectedErrorPercentileArrayForFastq(
            fastq.filePath, self.subsample, self.percentile, self.primerLength
        )
        return fastq, percentileExpectedError


def makeExpectedErrorAverageArrayForFastq(
    path: str, subsample: int = 0, primerLength: int = 0
):
    expectedErrorMatrix = fastqAnalysis.buildExpectedErrorMatrix(
        path, subsample=subsample, leftTrim=primerLength
    )
    meanArray = numpy.mean(expectedErrorMatrix, axis=0)
    return meanArray


def makeExpectedErrorPercentileArrayForFastq(
    path: str, subsample: int = 0, percentile: int = 83, primerLength: int = 0
):
    expectedErrorMatrix = fastqAnalysis.buildExpectedErrorMatrix(
        path, subsample=subsample, leftTrim=primerLength
    )
    percentileList = []
    for positionArray in expectedErrorMatrix.transpose():
        percentileList.append(numpy.percentile(positionArray, percentile))
    return numpy.array(percentileList)


def makeExpectedErrorPercentileArrayForFastqList(
    fastqList: list, subsample: int = 0, percentile: int = 83, primerLength: int = 0
):
    try:
        from . import easyMultiprocessing
    except ImportError:
        import easyMultiprocessing
    parallelAgent = ParallelExpectedErrorPercentileAgent(
        subsample, percentile, primerLength
    )
    expectedErrorReturns = easyMultiprocessing.parallelProcessRunner(
        parallelAgent.calculateAverageExpectedError, fastqList
    )
    averageExpectedErrorMatrix = numpy.stack(
        [expectedErrorArray[1] for expectedErrorArray in expectedErrorReturns]
    )
    averageExpectedErrorArray = numpy.mean(averageExpectedErrorMatrix, axis=0)
    return averageExpectedErrorArray


def makeExpectedErrorPercentileArraysForDirectory(
    path: str,
    namingStandard: typing.Type[fileNamingStandards.NamingStandard],
    subsample: int = 0,
    percentile: int = 83,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    import os

    if not os.path.isdir(path):
        raise NotADirectoryError("Unable to find directory %s" % path)
    fastqList = fastqHandler.findSamplesInFolder(path, namingStandard)
    forwardFastqs = [fastq for fastq in fastqList if fastq.direction == 1]
    reverseFastqs = [fastq for fastq in fastqList if fastq.direction == 2]
    forwardExpectedErrorArray = makeExpectedErrorPercentileArrayForFastqList(
        forwardFastqs, subsample, percentile, forwardPrimerLength
    )
    reverseExpectedErrorArray = makeExpectedErrorPercentileArrayForFastqList(
        reverseFastqs, subsample, percentile, reversePrimerLength
    )
    return forwardExpectedErrorArray, reverseExpectedErrorArray


def makeXAndYValuesForPositionArray(positionArray: collections.Iterable):
    xValues = []
    yValues = []
    for position, value in enumerate(positionArray):
        xValues.append(position)
        yValues.append(value)
    xValues = numpy.asarray(xValues)
    yValues = numpy.asarray(yValues)
    return xValues, yValues


def getGroupName(
    path: str, namingStandard: typing.Type[fileNamingStandards.NamingStandard]
):
    import os

    if not os.path.isdir(path):
        raise NotADirectoryError("Unable to find directory %s" % path)
    fastqList = fastqHandler.findSamplesInFolder(path, namingStandard)
    return fastqList[0].group


def calculateExpectedErrorCurvesForFastqFolder(
    path: str,
    namingStandard: typing.Type[fileNamingStandards.NamingStandard],
    subsample: int = 0,
    percentile: int = 83,
    makePNG: bool = False,
    sampleGroupID: str = None,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    if not sampleGroupID:
        sampleGroupID = getGroupName(path, namingStandard)
    forwardExpectedErrorArray, reverseExpectedErrorArray = (
        makeExpectedErrorPercentileArraysForDirectory(
            path,
            namingStandard,
            subsample,
            percentile,
            forwardPrimerLength,
            reversePrimerLength,
        )
    )
    forwardPositions, forwardValues = makeXAndYValuesForPositionArray(
        forwardExpectedErrorArray
    )
    reversePositions, reverseValues = makeXAndYValuesForPositionArray(
        reverseExpectedErrorArray
    )
    forwardCurve = fitExponentialCurve(
        forwardPositions, forwardValues, makePNG, "%s forward reads" % sampleGroupID
    )
    reverseCurve = fitExponentialCurve(
        reversePositions, reverseValues, makePNG, "%s reverse reads" % sampleGroupID
    )
    return forwardCurve, reverseCurve


def calculateExpectedErrorCurvesForFastqList(
    fastqList,
    subsample: int = 0,
    percentile: int = 83,
    makePNG: bool = False,
    sampleGroupID: str = None,
    forwardPrimerLength: int = 0,
    reversePrimerLength: int = 0,
):
    if not sampleGroupID:
        sampleGroupID = fastqList[0].group
    forwardFastqs = [fastq for fastq in fastqList if fastq.direction == 1]
    reverseFastqs = [fastq for fastq in fastqList if fastq.direction == 2]
    forwardExpectedErrorArray = makeExpectedErrorPercentileArrayForFastqList(
        forwardFastqs, subsample, percentile, forwardPrimerLength
    )
    reverseExpectedErrorArray = makeExpectedErrorPercentileArrayForFastqList(
        reverseFastqs, subsample, percentile, reversePrimerLength
    )
    forwardPositions, forwardValues = makeXAndYValuesForPositionArray(
        forwardExpectedErrorArray
    )
    reversePositions, reverseValues = makeXAndYValuesForPositionArray(
        reverseExpectedErrorArray
    )
    forwardCurve = fitExponentialCurve(
        forwardPositions,
        forwardValues,
        makePNG,
        "%s forward reads. %s percentile" % (sampleGroupID, ordinal(percentile)),
    )
    reverseCurve = fitExponentialCurve(
        reversePositions,
        reverseValues,
        makePNG,
        "%s reverse reads. %s percentile" % (sampleGroupID, ordinal(percentile)),
    )
    return forwardCurve, reverseCurve


def ordinal(number: int):
    onesDigit = number % 10
    append = {1: "st", 2: "nd", 3: "rd"}
    if onesDigit in append:
        return "%s%s" % (number, append[onesDigit])
    else:
        return "%s%s" % (number, "th")
