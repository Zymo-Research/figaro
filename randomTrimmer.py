import figaroSupport
import os
import random


def getFileList(directory:str):
    if not os.path.isdir(directory):
        raise NotADirectoryError()
    fileList = os.listdir(directory)
    fullFileList = [os.path.join(directory, file) for file in fileList]
    return fullFileList


def randomTrimFile(inputFilePath:str, outputFilePath:str, maxTrim:int=10):
    fastq = figaroSupport.fastqHandler.FastqFile(inputFilePath)
    outputFile = open(outputFilePath, 'w')
    for read in fastq:
        trimLength = random.randint(0, maxTrim - 1)
        if trimLength:
            read.sequence = read.sequence[:-trimLength]
            read.quality = read.quality[:-trimLength]
        print(read, file=outputFile)
    outputFile.close()
    fastq.close()


def trimFiles(inputDir:str, outputDir:str, maxTrim:int=10):
    inputFiles = getFileList(inputDir)
    for filePath in inputFiles:
        outputFilePath = os.path.join(
            outputDir,
            os.path.split(filePath)[1]
        )
        randomTrimFile(filePath, outputFilePath, maxTrim)
        print("Random trimmed %s" %filePath)


def main():
    trimFiles("C:\\Users\\mweinstein\\dadatest\\input\\sequence", "randomTrim", 10)

if __name__ == "__main__":
    main()