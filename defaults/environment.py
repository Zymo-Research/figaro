import os
import datetime
timestamp = str(datetime.datetime.now().timestamp()).replace(".", "")
dataFolder = "/data"
inputFolder = os.path.join(dataFolder, "input")
outputFolder = os.path.join(dataFolder, "output")
outputFileName = "trimParameters.json"
logFile = os.path.join(outputFolder, "dada2.%s.log" %timestamp)
