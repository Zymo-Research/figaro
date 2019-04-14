expectedEndings = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]

class NamingStandard(object):

    __slots__ = ["fileName", "fileDirectory", "filePath", "sampleNumber", "group", "direction", "sampleID"]

    def __init__(self, filePath:str):
        self.filePath = filePath
        self.fileDirectory, self.fileName = self.separateNameAndDirectory(filePath)
        self.group, self.sampleNumber, self.direction = self.getSampleInfo(self.fileName)
        self.sampleID = (self.group, self.sampleNumber)

    def separateNameAndDirectory(self, path:str):
        import os
        directory, name = os.path.split(path)
        return directory, name

    def getSampleInfo(self, fileName:str):
        raise RuntimeError("This function should always be getting overridden. If you see this, someone called the base class by mistake.")

    def sameSample(self, other):
        if not isinstance(other, NamingStandard):
            raise TypeError("Can only check for same sample in another naming standard type")
        if self.group == other.group and self.sampleNumber == other.sampleNumber:
            return True
        return False

    def __str__(self):
        return self.filePath

    def __hash__(self):
        return hash(self.filePath)

    def __eq__(self, other):
        return self.group == other.group and self.sampleNumber == other.sample and self.direction == other.direction

    def __ne__(self, other):
        return not self.__eq__(other)

    def __xor__(self, other):
        return self.sameSample(other)


class ZymoServicesNamingStandard(NamingStandard):

    def getSampleInfo(self, fileName:str):
        baseName = fileName.split(".")[0]
        group, sample, direction = baseName.split("_")
        direction = int(direction.replace("R",""))
        return group, sample, direction


class IlluminaStandard(NamingStandard):

    def getSampleInfo(self, fileName:str):
        baseName = fileName.split(".")[0]
        baseSplit = baseName.split("_")
        group = baseSplit[0]
        sample = baseSplit[1]
        direction = int(baseSplit[2].replace("R",""))
        return group, sample, direction


class ManualNamingStandard(NamingStandard):
    __slots__ = ["fileName", "fileDirectory", "filePath", "sampleNumber", "group", "direction", "sampleID"]

    def __init__(self, filePath: str, group:str, number:int, direction:int):
        self.filePath = filePath
        self.fileDirectory, self.fileName = self.separateNameAndDirectory(filePath)
        self.group = group
        self.sampleNumber = number
        self.direction = direction
        if direction not in [1, 2]:
            raise ValueError("Read direction must be either 1 or 2. %s was given" %direction)
        self.sampleID = (self.group, self.sampleNumber)