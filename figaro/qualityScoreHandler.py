import logging
import math
import typing

logger = logging.getLogger(__name__)


class EncodingScheme(object):

    def __init__(
        self,
        name: str,
        base: int,
        startCharacter: str,
        endCharacter: str,
        pErrorToScore: typing.Callable,
        scoreToPError: typing.Callable,
    ):
        self.name = name
        self.base = base
        self.characterSet = self.makeCharacterSet(startCharacter, endCharacter)
        self.range = self.calculateRange(startCharacter, endCharacter)
        self.fromPErrorFormula = pErrorToScore
        self.toPErrorFormula = scoreToPError

    def makeCharacterSet(self, start: str, end: str):
        rangeStart = ord(start)
        rangeEnd = ord(end) + 1
        return [chr(asciiValue) for asciiValue in range(rangeStart, rangeEnd)]

    def calculateRange(self, start: str, end: str):
        rangeStart = ord(start)
        rangeEnd = ord(end)
        return rangeEnd - rangeStart

    def toPError(self, score: [int, str]):
        if isinstance(score, str):
            if len(score) == 1:
                score = convertCharacterToScore(score, self.base)
            else:
                logger.critical(
                    "Attempt to convert multiple characters to error probability. Function can only handle one conversion per call."
                )
                raise ValueError(
                    "Attempt to get pError for entire string. Need one value at a time. String: %s"
                    % score
                )
        return self.toPErrorFormula(score)

    def scoreFromPError(self, pError: float, round: bool = True):
        return self.fromPErrorFormula(pError, round)

    def encodedFromPError(self, pError: float):
        return chr(self.scoreFromPError(pError, round=True) + self.base)

    def qualifyWithQualityString(self, qualityString: str):
        try:
            throwaway = self.eliminated
        except AttributeError:
            self.eliminated = False
        if not self.eliminated:
            qualityString = str(qualityString)
            for character in qualityString:
                if character not in self.characterSet:
                    self.eliminated = True
                    break

    def __str__(self):
        return self.name

    def __eq__(self, other: [str]):
        if not isinstance(other, (str, EncodingScheme)):
            raise TypeError(
                "Unable to compare encoding scheme types with anything but string or other EncodingScheme objects"
            )
        return self.name == str(other)


def makeEncodingTable():
    encodingTable = [  # In order of likelihood
        EncodingScheme(
            "Sanger/Illumina 1.8+", 33, "!", "I", pErrorToPhred, phredToPError
        ),
        EncodingScheme("Illumina 1.8+", 33, "!", "J", pErrorToPhred, phredToPError),
        EncodingScheme("Illumina 1.5-7", 64, "B", "i", pErrorToPhred, phredToPError),
        EncodingScheme("Illumina 1.3-4", 64, "@", "h", pErrorToPhred, phredToPError),
        EncodingScheme("Solexa", 64, ";", "h", pErrorToSolexa, solexaToPError),
        EncodingScheme("Pacbio", 33, "!", "~", pErrorToPhred, phredToPError),
    ]
    return encodingTable


def convertCharacterToScore(character, base: int = 33):
    return ord(character) - base


def convertToNumericArray(qualityString, base: int = 33):
    phredScores = []
    for character in qualityString:
        phredScores.append(convertCharacterToScore(character, base))
    return tuple(phredScores)


def pErrorToPhred(pError: float, roundValue: bool = True):
    score = -10 * (math.log(pError, 10))
    if roundValue:
        score = round(score)
    return score


def phredToPError(phred: [int, float]):
    return 10 ** (-phred / 10)


def pErrorToSolexa(
    pError: float, roundValue: bool = True
):  # google the definition of "arcane"
    score = -10 * (math.log(pError / (1 - pError), 10))
    if roundValue:
        score = round(score)
    return score


def solexaToPError(
    solexa: [int, float]
):  # seriously, who uses this encoding anymore, and who realizes that it's a slightly different formula?
    return 1 / (
        (10 ** (solexa / 10)) + 1
    )  # Let's hope I don't have to derive that one again


class _Encodings(object):

    def __init__(self):
        self.sanger = EncodingScheme(
            "Sanger/Illumina 1.8+", 33, "!", "I", pErrorToPhred, phredToPError
        )
        self.illumina = EncodingScheme(
            "Illumina 1.8+", 33, "!", "J", pErrorToPhred, phredToPError
        )
        self.illumina1_8 = self.illumina
        self.illumina1_5 = EncodingScheme(
            "Illumina 1.5-7", 64, "B", "i", pErrorToPhred, phredToPError
        )
        self.illumina1_3 = EncodingScheme(
            "Illumina 1.3-4", 64, "@", "h", pErrorToPhred, phredToPError
        )
        self.solexa = EncodingScheme(
            "Solexa", 64, ";", "h", pErrorToSolexa, solexaToPError
        )
        self.pacbio = EncodingScheme(
            "Pacbio", 33, "!", "~", pErrorToPhred, phredToPError
        )


encodingSchemes = _Encodings()


def cumulativeExpectedErrorArray(
    qualityString: str, encoding: EncodingScheme = encodingSchemes.illumina
):
    cumulativeExpectedErrorArray = []
    cumulativeExpectedError = 0.0  # ask me no questions, I'll tell you no lies/errors
    qualityString = str(qualityString)
    for character in qualityString:
        cumulativeExpectedError += encoding.toPError(character)
        cumulativeExpectedErrorArray.append(cumulativeExpectedError)
    return cumulativeExpectedErrorArray


def cumulativeExpectedErrorArrayDada2Exact(
    qualityString: str, encoding: EncodingScheme = encodingSchemes.illumina
):
    cumulativeExpectedErrorArray = []
    cumulativeExpectedError = 0.0  # ask me no questions, I'll tell you no lies/errors
    qualityString = str(qualityString)
    for character in qualityString:
        score = ord(character) - encoding.base
        cumulativeExpectedError += 10 ** (-score / 10)
        cumulativeExpectedErrorArray.append(cumulativeExpectedError)
    return cumulativeExpectedErrorArray


def convertQualityString(
    qualityString: str, inputScheme: EncodingScheme, outputScheme: EncodingScheme
):
    qualityString = str(qualityString)
    if inputScheme.fromPErrorFormula == outputScheme.fromPErrorFormula:
        baseDifference = inputScheme.base - outputScheme.base
        outputString = ""
        for character in qualityString:
            outputString += chr(ord(character) - baseDifference)
        return outputString
    else:
        outputString = ""
        for character in qualityString:
            pError = inputScheme.toPError(character)
            outputString += outputScheme.encodedFromPError(pError)
        return outputString
