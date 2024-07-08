import os
import logging
import typing

logger = logging.getLogger(__name__)
try:
    from . import qualityScoreHandler
    from . import fileNamingStandards
except ImportError:
    import qualityScoreHandler, fileNamingStandards


class ReadMetadataLine(object):

    def __init__(self, rawMetadata):
        self.rawMetadata = rawMetadata
        if not rawMetadata.startswith("@"):
            logger.warning(
                "Got a metadata line that did not start with an @ symobol. This goes against the fastq standard and may suggest a corrupt file. Line: %s"
                % rawMetadata
            )
        metadataSplit = rawMetadata.strip().split(" ")
        if not len(metadataSplit) == 2:
            errorMessage = (
                "Got a metadata line that appears to have more than two elements divided by space. %s"
                % rawMetadata
            )
            logger.critical(errorMessage)
            raise FastqFormatError(errorMessage)
        equipmentInfo, readInfo = metadataSplit
        self.validEquipmentInfo = self.processEquipmentInfo(equipmentInfo, rawMetadata)
        self.validReadInfo = self.processReadInfo(readInfo, rawMetadata)
        self.allValidInfo = self.validEquipmentInfo and self.validReadInfo

    def processReadInfo(self, readInfo: str, rawMetadata: str = ""):
        validFields = True
        readInfo = readInfo.split(":")
        if not len(readInfo) == 4:
            errorMessage = (
                "Got a read info section of metadata that did not have 4 elements. Line: %s"
                % rawMetadata
            )
            logger.critical(errorMessage)
            raise FastqFormatError(errorMessage)
        self.direction, self.filtered, self.controlBits, self.index = readInfo
        try:
            self.direction = int(self.direction)
            if self.direction not in [1, 2]:
                validFields = False
                logger.error(
                    "Read direction found that was not 1 or 2. Line: %s" % rawMetadata
                )
        except ValueError:
            validFields = False
            logger.error(
                "Read direction could not be cast to integer. Line: %s" % rawMetadata
            )
        if self.filtered.upper() == "Y":
            self.filtered = True
            self.passedFilter = False
        elif self.filtered.upper() == "N":
            self.filtered = False
            self.passedFilter = True
        else:
            self.passedFilter = None
            validFields = False
            logger.error(
                "Got a value for filtered that was not Y or N. Line: %s" % rawMetadata
            )
        try:
            self.controlBits = int(self.controlBits)
            if not self.controlBits % 2 == 0:
                validFields = False
                logger.error(
                    "Got a control bits value of %s. Control bits should be an even number. Line: %s "
                    % (self.controlBits, rawMetadata)
                )
        except ValueError:
            validFields = False
            logger.error(
                "Unable to cast control bits to an integer. Line: %s " % rawMetadata
            )
        return validFields

    def processEquipmentInfo(self, equipmentInfo: str, rawMetadata: str = ""):
        validFields = True
        equipmentInfo = equipmentInfo.replace("@", "")
        equipmentInfo = equipmentInfo.split(":")
        if not len(equipmentInfo) == 7:
            logger.critical(
                "Equipment info section of metadata did not have 7 elements. Line: %s"
                % rawMetadata
            )
            raise FastqFormatError(
                "Equipment info section of metadata did not have 7 elements. Line: %s"
                % rawMetadata
            )
        (
            self.instrumentName,
            self.runID,
            self.flowcellID,
            self.tileNumber,
            self.laneNumber,
            self.xCoordinate,
            self.yCoordinate,
        ) = equipmentInfo
        try:
            self.runID = int(self.runID)
        except ValueError:
            validFields = False
            logger.error(
                "Run ID number could not be cast to integer. Metadata line: %s"
                % rawMetadata
            )
        try:
            self.laneNumber = int(self.laneNumber)
        except ValueError:
            validFields = False
            logger.error(
                "Lane number could not be cast to integer. Metadata line: %s"
                % rawMetadata
            )
        try:
            self.tileNumber = int(self.tileNumber)
        except ValueError:
            validFields = False
            logger.error(
                "Tile number could not be cast to integer. Metadata line: %s"
                % rawMetadata
            )
        try:
            self.xCoordinate = int(self.xCoordinate)
        except ValueError:
            validFields = False
            logger.error(
                "X-coordinate could not be cast to integer. Metadata line: %s"
                % rawMetadata
            )
        try:
            self.yCoordinate = int(self.yCoordinate)
        except ValueError:
            validFields = False
            logger.error(
                "Y-coordinate could not be cast to integer. Metadata line: %s"
                % rawMetadata
            )
        return validFields

    def __str__(self):
        return self.rawMetadata


class QualityScoreLine(object):

    def __init__(self, rawQualityLine: str, base: int = 33):
        self.qualityString = rawQualityLine
        self.phredScores = self.calculatePhredScores(base)

    def calculatePhredScores(self, base: int = 33):
        return qualityScoreHandler.convertToNumericArray(self.qualityString, base)

    def __str__(self):
        return self.qualityString

    def __getitem__(self, item):
        return self.phredScores[item]

    def __iter__(self):
        for value in self.phredScores:
            yield value


class SequenceLine(object):

    def __init__(self, rawSequence, runAnalysis: bool = False):
        self.sequence = rawSequence.strip().upper().replace(".", "N")
        self.length = len(self.sequence)
        if runAnalysis:
            self.baseFrequency = self.getBaseFrequencyTable()
            self.gcContent = self.calculateGCContent()

    def getBaseFrequencyTable(self):
        freq = {"A": 0, "G": 0, "C": 0, "T": 0, "N": 0}
        for base in self.sequence:
            try:
                freq[base] += 1
            except KeyError:
                logger.error(
                    "Found a sequence with an invalid character. Character: %s  Sequence: %s"
                    % (base, self.sequence)
                )
        return freq

    def calculateGCContent(self):
        totalReadBases = 0
        gcBases = 0
        for base in "ATGC":
            totalReadBases += self.baseFrequency[base]
            if base in "GC":
                gcBases += self.baseFrequency[base]
        if totalReadBases == 0:
            return 0
        return gcBases / totalReadBases

    def __len__(self):
        return self.length

    def __str__(self):
        return self.sequence

    def __eq__(self, other):
        if type(other) == SequenceLine:
            return self.sequence == other.sequence
        elif type(other) == str:
            return self.sequence == SequenceLine(other).sequence
        else:
            logger.critical(
                "Attempted to compare a sequence to something that is not a sequence line type or string. Value in question was type %s: %s"
                % (type(other), other)
            )


class FastqLineSet(object):

    def __init__(
        self,
        metadata: str,
        sequence: str,
        spacer: str,
        quality: str,
        depth: int = 0,
        analyzeMetadata: bool = False,
        analyzeSequence: bool = False,
        analyzeSequenceInDepth: bool = False,
        analyzeQuality: bool = False,
        qualityBase: int = 33,
    ):
        self.metadata = metadata.strip()
        self.sequence = sequence.strip()
        self.spacer = spacer.strip()
        self.quality = quality.strip()
        if depth >= 1 or analyzeQuality:
            self.quality = QualityScoreLine(quality, qualityBase)
        if depth >= 2 or analyzeSequence or analyzeSequenceInDepth:
            if depth >= 4 or analyzeSequenceInDepth:
                self.sequence = SequenceLine(self.sequence, runAnalysis=True)
            else:
                self.sequence = SequenceLine(self.sequence)
        if depth >= 3 or analyzeMetadata:
            self.metadata = ReadMetadataLine(self.metadata)

    def __str__(self):
        return "%s\n%s\n%s\n%s" % (
            self.metadata,
            self.sequence,
            self.spacer,
            self.quality,
        )


def reanalyzeFastqLineSet(
    fastqLineSet: FastqLineSet,
    depth: int = 0,
    analyzeMetadata: bool = False,
    analyzeSequence: bool = False,
    analyzeSequenceInDepth: bool = False,
    analyzeQuality: bool = False,
    qualityBase: int = 33,
):
    return FastqLineSet(
        str(fastqLineSet.metadata),
        str(fastqLineSet.sequence),
        str(fastqLineSet.spacer),
        str(fastqLineSet.quality),
        depth,
        analyzeMetadata,
        analyzeSequence,
        analyzeSequenceInDepth,
        analyzeQuality,
        qualityBase,
    )


class FastqFile(object):

    def __init__(
        self,
        path: str,
        depth: int = 0,
        analyzeMetadata: bool = False,
        analyzeSequence: bool = False,
        analyzeSequenceInDepth: bool = False,
        analyzeQuality: bool = False,
        fullValidation: bool = False,
        qualityScoreScheme: [qualityScoreHandler.EncodingScheme, None] = None,
        subsample: int = 0,
        leftTrim: int = 0,
        rightTrim: int = 0,
    ):
        self.path = path
        if not os.path.isfile(path):
            logger.critical("Unable to find fastq file at %s" % path)
            raise FileNotFoundError("Unable to find fastq file at %s" % path)
        if not qualityScoreScheme:
            qualityScoreScheme = findQualityScoreEncoding(path)
        if type(qualityScoreScheme) == qualityScoreHandler.EncodingScheme:
            self.qualityScoreScheme = qualityScoreScheme
        else:
            raise TypeError(
                "Quality score scheme must be of qualityScoreHandler.EncodingScheme type. Passed: %s of type %s."
                % (qualityScoreScheme, type(qualityScoreScheme))
            )
        self.depth = depth
        self.leftTrim = leftTrim
        if rightTrim == 0:
            self.rightTrim = None
        elif rightTrim < 0:
            self.rightTrim = -rightTrim
        else:
            raise ValueError("Right trim can only be zero or a positive integer.")
        self.analyzeMetadata = analyzeMetadata
        self.analyzeSequence = analyzeSequence
        self.analyzeSequenceInDepth = analyzeSequenceInDepth
        self.analyzeQuality = analyzeQuality
        self.fullValidation = fullValidation
        self.reachedEnd = False
        self.gzipped = self.checkGzip(path)
        if self.gzipped:
            import gzip

            self.filehandle = gzip.open(path, "rt")
        else:
            self.filehandle = open(path, "r")
        self.open = True
        subsample = int(subsample)
        if subsample == 0:
            subsample = 1
        self.subsample = subsample
        self.currentLine = 0

    def checkGzip(self, path):
        try:
            from . import gzipIdentifier
        except ImportError:
            import gzipIdentifier
        return gzipIdentifier.isGzipped(path)

    def getNextRead(self):

        def read4Lines():
            readBuffer = []
            for i in range(4):
                nextLine = self.filehandle.readline()
                if not nextLine:
                    self.reachedEnd = True
                    break
                nextLine = nextLine.strip()
                if nextLine:
                    readBuffer.append(nextLine)
            if self.reachedEnd:
                if readBuffer:
                    logger.error(
                        "Fastq file at %s appears to me missing lines (found something not a multiple of 4."
                        % self.path
                    )
                    for i in range(4 - len(readBuffer)):
                        readBuffer.append("")
            if readBuffer:
                readBuffer[1] = readBuffer[1][self.leftTrim : self.rightTrim]
                readBuffer[3] = readBuffer[3][self.leftTrim : self.rightTrim]
            return readBuffer

        if not self.open:
            logger.critical(
                "Attempting to read from a closed fastq file at %s" % self.path
            )
            raise ValueError("I/O operation on a closed file")
        readBuffer = None
        includedLine = False
        while not includedLine:
            readBuffer = read4Lines()
            self.currentLine += 1
            includedLine = (
                self.currentLine - 1
            ) % self.subsample == 0 or self.reachedEnd
        if not readBuffer:
            return readBuffer
        else:
            fastqLineSet = FastqLineSet(
                *readBuffer,
                depth=self.depth,
                analyzeMetadata=self.analyzeMetadata,
                analyzeSequence=self.analyzeSequence,
                analyzeSequenceInDepth=self.analyzeSequenceInDepth,
                analyzeQuality=self.analyzeQuality,
                qualityBase=self.qualityScoreScheme.base
            )
            if self.fullValidation:
                if not len(readBuffer[1]) == len(readBuffer[3]):
                    raise FastqValidationError(
                        "Got mismatched sequence and quality line lengths for line %s"
                        % readBuffer
                    )
                if type(fastqLineSet.metadata) == str:
                    metadata = ReadMetadataLine(str(fastqLineSet.metadata))
                else:
                    metadata = fastqLineSet.metadata
                if not metadata.allValidInfo:
                    raise FastqValidationError(
                        "Got some invalid metadata for line %s" % readBuffer
                    )
            return fastqLineSet

    def close(self):
        if not self.filehandle.closed:
            self.filehandle.close()

    def __iter__(self):
        return self

    def __next__(self):
        returnValue = self.getNextRead()
        if self.reachedEnd:
            self.close()
            raise StopIteration
        else:
            return returnValue

    def __str__(self):
        return "Fastq file object at %s" % self.path


class FastqFilePair(object):

    def __init__(
        self,
        pe1Path: str,
        pe2Path: str,
        depth: int = 0,
        analyzeMetadata: bool = False,
        analyzeSequence: bool = False,
        analyzeSequenceInDepth: bool = False,
        analyzeQuality: bool = False,
        fullValidation: bool = False,
        qualityScoreScheme: qualityScoreHandler = None,
        subsample: int = 0,
    ):
        self.pe1Path = pe1Path
        if not os.path.isfile(pe1Path):
            logger.critical("Unable to find fastq file at %s" % pe1Path)
            raise FileNotFoundError(
                "Unable to find paired-end 1 fastq file at %s" % pe1Path
            )
        self.pe2Path = pe2Path
        if not os.path.isfile(pe2Path):
            logger.critical("Unable to find fastq file at %s" % pe2Path)
            raise FileNotFoundError(
                "Unable to find paired-end 1 fastq file at %s" % pe2Path
            )
        self.depth = depth
        self.analyzeMetadata = analyzeMetadata
        self.analyzeSequence = analyzeSequence
        self.analyzeSequenceInDepth = analyzeSequenceInDepth
        self.analyzeQuality = analyzeQuality
        self.fullValidation = fullValidation
        self.reachedEnd = False
        if subsample == 0:
            subsample = 1
        self.subsample = subsample
        self.pe1FileHandle = FastqFile(
            pe1Path,
            depth=depth,
            analyzeMetadata=analyzeMetadata,
            analyzeSequence=analyzeSequence,
            analyzeSequenceInDepth=analyzeSequenceInDepth,
            analyzeQuality=analyzeQuality,
            fullValidation=fullValidation,
            qualityScoreScheme=qualityScoreScheme,
            subsample=subsample,
        )
        self.pe2FileHandle = FastqFile(
            pe2Path,
            depth=depth,
            analyzeMetadata=analyzeMetadata,
            analyzeSequence=analyzeSequence,
            analyzeSequenceInDepth=analyzeSequenceInDepth,
            analyzeQuality=analyzeQuality,
            fullValidation=fullValidation,
            qualityScoreScheme=qualityScoreScheme,
            subsample=subsample,
        )
        if (
            not self.pe1FileHandle.qualityScoreScheme
            == self.pe2FileHandle.qualityScoreScheme
        ):
            logger.warning(
                "Paired end files appear to have different quality score encodings. Pe1: %s:%s. Pe2: %s%s"
                % (
                    self.pe1FileHandle.qualityScoreScheme,
                    self.pe1FileHandle.path,
                    self.pe2FileHandle.qualityScoreScheme,
                    self.pe2FileHandle.path,
                )
            )
        self.open = True
        self.reportedReadMismatch = False

    def getNextReadPair(self):
        if not self.open:
            logger.critical(
                "Attempting to read from a closed fastq files at %s and %s"
                % (self.pe1Path, self.pe2Path)
            )
            raise ValueError("I/O operation on a closed file")
        nextPe1 = self.pe1FileHandle.getNextRead()
        nextPe2 = self.pe2FileHandle.getNextRead()
        if (nextPe1 and not nextPe2) or (not nextPe1 and nextPe2):
            if nextPe1:
                logger.error(
                    "Ran out of paired-end 2 reads with remaining paired-end 1 reads for file pair %s and %s"
                    % (self.pe1Path, self.pe2Path)
                )
            else:
                logger.error(
                    "Ran out of paired-end 1 reads with remaining paired-end 2 reads for file pair %s and %s"
                    % (self.pe1Path, self.pe2Path)
                )
            if self.fullValidation:
                raise FastqValidationError(
                    "Reached end of one paired-end file before the other. Files: %s and %s"
                    % (self.pe1Path, self.pe2Path)
                )
        if not nextPe1 and not nextPe2:
            self.reachedEnd = True
            return None
        if nextPe1 and nextPe2 and self.fullValidation:
            self.runValidation(nextPe1, nextPe2)
        return nextPe1, nextPe2

    def runValidation(self, pe1: FastqLineSet, pe2: FastqLineSet):
        if type(pe1.metadata) == str:
            pe1Metadata = ReadMetadataLine(str(pe1.metadata))
        elif type(pe1.metadata) == ReadMetadataLine:
            pe1Metadata = pe1.metadata
        else:
            raise TypeError(
                "Only able to compare metadata as string or metadata objects"
            )
        if type(pe2.metadata) == str:
            pe2Metadata = ReadMetadataLine(str(pe2.metadata))
        elif type(pe1.metadata) == ReadMetadataLine:
            pe2Metadata = pe2.metadata
        else:
            raise TypeError(
                "Only able to compare metadata as string or metadata objects"
            )
        if not pe1Metadata.allValidInfo or not pe2Metadata.allValidInfo:
            raise FastqValidationError(
                "Got invalid metadata field for at least one read in paired end mates:\n%s\n%s"
                % (pe1, pe2)
            )
        if not validPairedEndMetadata(pe1Metadata, pe2Metadata):
            raise FastqValidationError(
                "Got invalid metadata match for paired end mates:\n%s\n%s" % (pe1, pe2)
            )

    def close(self):
        self.pe1FileHandle.close()
        self.pe2FileHandle.close()
        self.open = False

    def __iter__(self):
        return self

    def __next__(self):
        returnValue = self.getNextReadPair()
        if self.reachedEnd:
            raise StopIteration
        else:
            return returnValue

    def __str__(self):
        return "Fastq file pair object at %s and %s" % (self.pe1Path, self.pe2Path)


class FastqValidationError(Exception):
    pass


class FastqFormatError(Exception):
    pass


def validPairedEndMetadata(pe1: ReadMetadataLine, pe2: ReadMetadataLine):
    matchFields = [
        "instrumentName",
        "runID",
        "flowcellID",
        "laneNumber",
        "tileNumber",
        "xCoordinate",
        "yCoordinate",
        "index",
    ]
    for field in matchFields:
        pe1Value = getattr(pe1, field)
        pe2Value = getattr(pe2, field)
        if not pe1Value == pe2Value:
            logger.error("Mismatch on %s" % matchFields)
            return False
    if not (
        (pe1.direction == 1 and pe2.direction == 2)
        or (pe2.direction == 1 and pe1.direction == 2)
    ):
        return False
    return True


def validFastqFile(path: str):
    readCount = 0
    fastq = FastqFile(path, fullValidation=True)
    read = fastq.getNextRead()
    while read:
        try:
            read = fastq.getNextRead()
            readCount += 1
        except Exception as error:
            logger.error(error)
            return False
    fastq.close()
    return readCount


def validFastqPair(pe1Path: str, pe2Path: str):
    readCount = 0
    fastqPair = FastqFilePair(pe1Path, pe2Path, fullValidation=True)
    read = fastqPair.getNextReadPair()
    while read:
        try:
            read = fastqPair.getNextReadPair()
            readCount += 1
        except Exception as error:
            logger.error(error)
            fastqPair.close()
            return False
    fastqPair.close()
    return readCount


def estimateReadLength(path: str, samplesize: int = 100, getVariance=False):
    lengths = []
    fastq = FastqFile(path)
    read = fastq.getNextRead()
    while read:
        lengths.append(len(read.sequence))
        if len(lengths) >= samplesize:
            break
        read = fastq.getNextRead()
    meanReadLength = sum(lengths) / len(lengths)
    if getVariance:
        import statistics

        if len(lengths) > 1:
            lengthVariance = statistics.variance(lengths)
        else:
            lengthVariance = 0
        return round(meanReadLength), lengthVariance
    return round(meanReadLength)


def getLongestReadInFile(path: str):
    longestReadLength = 0
    fastq = FastqFile(path)
    for read in fastq:
        if len(read.sequence) > longestReadLength:
            longestReadLength = len(read.sequence)
    fastq.close()
    return longestReadLength


def countReads(path: str):
    readCount = 0
    fastq = FastqFile(path)
    read = fastq.getNextRead()
    while read:
        readCount += 1
        read = fastq.getNextRead()
    fastq.close()
    return readCount


def findQualityScoreEncoding(path: str, lineLimit: int = 100):
    candidates = qualityScoreHandler.makeEncodingTable()
    for i in range(len(candidates)):
        candidates[i].eliminated = False
    fastq = FastqFile(
        path, qualityScoreScheme=qualityScoreHandler.encodingSchemes.sanger
    )
    line = fastq.getNextRead()
    lineCount = 0
    while line:
        for candidate in candidates:
            candidate.qualifyWithQualityString(line.quality)
        remaining = len([scheme for scheme in candidates if not scheme.eliminated])
        lineCount += 1
        if lineLimit > 0:
            if lineCount >= lineLimit:
                break
        if remaining == 0:
            logger.error(
                "No valid quality scoring scheme found for fastq file %s" % path
            )
            fastq.close()
            return None
        elif remaining == 1:
            break
    for candidate in candidates:
        if not candidate.eliminated:
            del candidate.eliminated
            fastq.close()
            return candidate


def findSamplesInFolder(
    directory: str,
    namingStandard: typing.Type[
        fileNamingStandards.NamingStandard
    ] = fileNamingStandards.IlluminaStandard,
):
    import os

    if not os.path.isdir(directory):
        raise NotADirectoryError("%s is not a directory or not found." % directory)
    fastqFileInfoList = []
    expectedEndings = fileNamingStandards.expectedEndings
    for item in os.listdir(directory):
        isFastqFile = False
        for expectedEnding in expectedEndings:
            if item.endswith(expectedEnding):
                isFastqFile = True
                break
        if not isFastqFile:
            continue
        filePath = os.path.join(directory, item)
        fastqFileInfoList.append(namingStandard(filePath))
    return fastqFileInfoList


def getSamplePairTableFromFolder(
    directory: str, namingStandard: typing.Type[fileNamingStandards.NamingStandard]
):
    def hasMate(fastq: fileNamingStandards.NamingStandard, potentialMates: list):
        for potentialMate in potentialMates:
            if fastq.sameSample(potentialMate):
                return potentialMate
        return False

    allFastqs = findSamplesInFolder(directory, namingStandard)
    pairedFastqs = {"unpaired": []}
    forwardFiles = [fastq for fastq in allFastqs if fastq.direction == 1]
    reverseFiles = [fastq for fastq in allFastqs if fastq.direction == 2]
    for fastq in forwardFiles:
        foundMate = hasMate(fastq, reverseFiles)
        if foundMate:
            reverseFiles.remove(foundMate)
            pairedFastqs[fastq.sampleID] = (fastq, foundMate)
        else:
            pairedFastqs["unpaired"].append(fastq)
    for fastq in reverseFiles:
        pairedFastqs["unpaired"].append(fastq)
    if not pairedFastqs["unpaired"]:
        del pairedFastqs["unpaired"]
    return pairedFastqs
