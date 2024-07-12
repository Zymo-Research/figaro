expectedEndings = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
aliasList = {
    "zymo": "zymo",
    "zymoservicesnamingstandard": "zymo",
    "zymoservices": "zymo",
    "illumina": "illumina",
    "keriksson": "keriksson",
    "nononsense": "nononsense",
    "fvieira": "fvieira",
    "yzhang": "yzhang",
}


class NamingStandard(object):

    __slots__ = [
        "fileName",
        "fileDirectory",
        "filePath",
        "sampleNumber",
        "group",
        "direction",
        "sampleID",
    ]

    def __init__(self, filePath: str):
        self.filePath = filePath
        self.fileDirectory, self.fileName = self.separateNameAndDirectory(filePath)
        self.group, self.sampleNumber, self.direction = self.getSampleInfo(
            self.fileName
        )
        self.sampleID = (self.group, self.sampleNumber)

    def separateNameAndDirectory(self, path: str):
        import os

        directory, name = os.path.split(path)
        return directory, name

    def getSampleInfo(self, fileName: str):
        raise RuntimeError(
            "This function should always be getting overridden. If you see this, someone called the base class by mistake."
        )

    def sameSample(self, other):
        if not isinstance(other, NamingStandard):
            raise TypeError(
                "Can only check for same sample in another naming standard type"
            )
        if self.group == other.group and self.sampleNumber == other.sampleNumber:
            return True
        return False

    def __str__(self):
        return self.filePath

    def __hash__(self):
        return hash(self.filePath)

    def __eq__(self, other):
        return (
            self.group == other.group
            and self.sampleNumber == other.sample
            and self.direction == other.direction
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __xor__(self, other):
        return self.sameSample(other)


class NoNonsenseNamingStandard(NamingStandard):

    def getSampleInfo(self, fileName: str):
        import re
        import os

        regex = r"_R?([12])(_\\d\\d\\d)?$"
        baseName = re.sub(r"\.(fq|fastq)(.gz)?$", "", os.path.basename(fileName))
        regexResult = re.search(regex, baseName)
        if not regexResult:
            raise ValueError(
                "Could not infer read orientation from filename: {}".format(fileName)
            )
        direction = int(regexResult[1])
        sample = group = re.sub(regex, "", baseName)
        return group, sample, direction


class ZymoServicesNamingStandard(NamingStandard):

    def getSampleInfo(self, fileName: str):
        baseName = fileName.split(".")[0]
        try:
            group, sample, direction = baseName.split("_")
        except ValueError:
            raise ValueError(
                "%s does not appear to be a valid Zymo Services file name. Please check file naming convention argument."
                % fileName
            )
        direction = int(direction.replace("R", ""))
        return group, sample, direction


class IlluminaStandard(NamingStandard):

    def getSampleInfo(self, fileName: str):
        try:
            baseName = fileName.split(".")[0]
            baseSplit = baseName.split("_")
            group = "_".join(baseSplit[:-4])
            sample = int(baseSplit[-4].replace("S", ""))
            direction = int(baseSplit[-2].replace("R", ""))
            return group, sample, direction
        except (ValueError, IndexError):
            raise ValueError(
                "%s does not appear to be a valid Illumina file name. Please check file naming convention argument."
                % fileName
            )


class KErickssonStandard(NamingStandard):

    def getSampleInfo(self, fileName: str):
        group, sampleAndDirection = fileName.split(".")[:2]
        try:
            sample, direction = sampleAndDirection.split("_")
            direction = direction.replace("R", "")
            direction = direction.replace("r", "")
            direction = int(direction)
        except ValueError:
            raise ValueError(
                "%s does not appear to be a valid file for this standard. Please check file naming convention argument."
                % fileName
            )

        return group, sample, direction


class FVieiraStandard(NamingStandard):

    def getSampleInfo(self, fileName: str):
        basename = fileName.split(".")[0]
        group = "default"
        try:
            sample, direction = basename.split("_")
            direction = direction.replace("R", "")
            direction = direction.replace("r", "")
            direction = int(direction)
        except ValueError:
            raise ValueError(
                "%s does not appear to be a valid file for this standard. Please check file naming convention argument."
                % fileName
            )

        return group, sample, direction


class YZhangStandard(NamingStandard):

    def getSampleInfo(self, fileName: str):
        basename = fileName.split(".")[0]
        group = "default"
        try:
            sample, seqType, direction = basename.split("_")
            direction = direction.replace("R", "")
            direction = direction.replace("r", "")
            direction = int(direction)
        except ValueError:
            raise ValueError(
                "%s does not appear to be a valid file for this standard. Please check file naming convention argument."
                % fileName
            )

        return group, sample, direction


class ManualNamingStandard(NamingStandard):
    __slots__ = [
        "fileName",
        "fileDirectory",
        "filePath",
        "sampleNumber",
        "group",
        "direction",
        "sampleID",
    ]

    def __init__(self, filePath: str, group: str, number: int, direction: int):
        self.filePath = filePath
        self.fileDirectory, self.fileName = self.separateNameAndDirectory(filePath)
        self.group = group
        self.sampleNumber = number
        self.direction = direction
        if direction not in [1, 2]:
            raise ValueError(
                "Read direction must be either 1 or 2. %s was given" % direction
            )
        self.sampleID = (self.group, self.sampleNumber)


def loadNamingStandard(name: str):
    aliasObjectKey = {
        "zymo": ZymoServicesNamingStandard,
        "illumina": IlluminaStandard,
        "keriksson": KErickssonStandard,
        "nononsense": NoNonsenseNamingStandard,
        "fvieira": FVieiraStandard,
        "yzhang": YZhangStandard,
    }
    nameLower = name.lower()
    if nameLower not in aliasList:
        raise ValueError("%s is not a valid naming standard identifier" % name)
    return aliasObjectKey[aliasList[nameLower]]
