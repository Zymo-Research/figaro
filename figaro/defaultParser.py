import os
import logging

logger = logging.getLogger(__name__)


def getDefaultsFolder():
    import sys

    mainModule = sys.modules["__main__"]
    if hasattr(mainModule, "__file__"):
        mainModulePath = os.path.abspath(mainModule.__file__)
        projectFolder = os.path.split(mainModulePath)[0]
        return os.path.join(projectFolder, "defaults")
    else:
        return os.path.join(os.getcwd(), "defaults")


def getDefaultPackageDict():
    defaultsFolder = getDefaultsFolder()
    dirContents = os.listdir(defaultsFolder)
    pythonPackages = [
        file
        for file in dirContents
        if file.endswith(".py") and os.path.isfile(os.path.join(defaultsFolder, file))
    ]
    pythonPackages = [file.replace(".py", "") for file in pythonPackages]
    defaultPackageDict = {}
    for package in pythonPackages:
        packageID = package.lower()
        if packageID == "environment" or packageID == "__init__":
            continue
        defaultPackageDict[packageID] = package
    return defaultPackageDict


def loadDefaultModule(name: str):
    import importlib

    packageID = name.lower()
    defaultPackages = getDefaultPackageDict()
    if packageID not in defaultPackages:
        logger.error(
            "Attempted to load default package %s. Only default packages found: %s"
            % (name, defaultPackages)
        )
        raise ValueError("Unable to find a default package called %s" % name)
    logger.info("Loading %s default package" % name)
    return importlib.import_module("defaults.%s" % defaultPackages[packageID])
