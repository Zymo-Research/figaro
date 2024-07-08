import os
import collections.abc as collections
import logging

logger = logging.getLogger(__name__)
typeHeirarchy = (int, float, str)


class EnvVariable(object):
    """
    Holder for individual environment variable parameters and associated tests and values.

    Accessing values:
    The value can always be directly accessed by calling EnvVariable.value.  Additionally, if the value is a boolean, trying to use EnvVariable in a logical statement will automatically call its value as the boolean for evaluation.  Casting it to a string will return the string of the value.
    """

    def __init__(
        self,
        name: str,
        typeRequirement: [type, list, tuple],
        default=None,
        flag: [int, str] = None,
        validationList: list = None,
        lowerBound: [int, float] = None,
        upperBound: [int, float] = None,
        expectedFile: bool = False,
        createdFile: bool = False,
        expectedDirectory: bool = False,
        createdDirectory: bool = False,
        logLevel: bool = False,
        required: bool = False,
        externalValidation: bool = False,
    ):
        """
        :param name: Name of the parameter and environment variable (environment variable should be all upper case when passed. Calling the name is functionally case insensitive
        :param typeRequirement: Can be either a single type or an iterable of types where multiple types may be acceptable.
        :param default: Default value for the parameter. Should be used for any non-required parameters. Must fit within type requirements
        :param flag: Flag value, used for automatically building argument strings
        :param validationList: Potential values for the parameter
        :param lowerBound: Minimum value for a numerical parameter
        :param upperBound: Maximum value for a numerical parameter
        :param expectedFile: File that is expected to exist at the time of checking. Will throw an exception if it is not there
        :param createdFile: File not expected to exist at the time of checking. This package will test to make sure the file can be created
        :param expectedDirectory: Folder that is expected to exxist at the time of checking. Will throw an exception if it is not there.
        :param createdDirectory: Folder that is being created as part of this run. Will be created at the time of checking to ensure files can be written there
        :param logLevel: Boolean value indicating that the parameter is to take in a logging level (sets some different checks for it)
        :param required: If true, this is a parameter that must be passed in and should not be using a default value.
        """
        self.name = name
        self.typeRequirement = typeRequirement
        self.default = default
        self.value = default
        self.flag = flag
        self.validationList = validationList
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.expectedFile = expectedFile
        self.createdFile = createdFile
        self.isFilePath = expectedFile or createdFile
        self.expectedDirectory = expectedDirectory
        self.createdDirectory = createdDirectory
        self.isDirectoryPath = expectedDirectory or createdDirectory
        self.required = required
        self.isArgument = not flag is None
        self.positionalArg = isinstance(flag, int)
        self.logLevel = logLevel
        self.externalValidation = externalValidation
        self.setValueValidations()
        self.runValidations()
        self.environmentVariableName = name.upper()
        self.usingDefaultValue = not self.setEnvironmentValue()

    def setValueValidations(self):
        if self.lowerBound is None and self.upperBound is None:
            if self.validationList:
                self.usingValidationBounds = False
                self.usingValidationList = True
            else:
                self.usingValidationBounds = False
                self.usingValidationList = False
                inherentlyBoundedCondition = (
                    self.typeRequirement == bool
                    or self.isFilePath
                    or self.isDirectoryPath
                    or self.logLevel
                )
                if not (inherentlyBoundedCondition or self.externalValidation):
                    logger.warning(
                        "Environment variable parameter %s has no validation list or bounds set."
                        % self.name
                    )
        else:
            self.usingValidationBounds = True
            infinity = float("inf")
            if self.lowerBound is None:
                self.lowerBound = -infinity
            if self.upperBound is None:
                self.upperBound = infinity
            if self.validationList:
                self.usingValidationList = True
                logger.warning(
                    "Environment variable parameter %s is using both bounds- and list-based validations. This is an unusual, but not impossible situation."
                    % self.name
                )
            else:
                self.usingValidationList = False

    def runValidations(self):
        self.validateTypeAndFlag()
        if not self.passedArgumentAssertions():
            raise ArgumentValueException(
                "Invalid values were given for one or more argument values when declaring environment variable paramater %s"
                % self.name
            )
        if self.expectedFile:
            self.validateExpectedFilePath()
        if self.expectedDirectory:
            self.validateExpectedDirectoryPath()

    def validLogLevelSettings(self):
        if not self.logLevel:
            logger.critical(
                "Hit the log level settings checker with self.logLevel as false for %s. This should never be able to happen and needs to be debugged."
                % self.name
            )
            return True
        if self.isDirectoryPath:
            logger.error(
                "Parameter %s is set as both directory and log level. This should not happen."
                % self.name
            )
            return False
        if self.isFilePath:
            logger.error(
                "Parameter %s is set as both file and log level.  This should never happen."
                % self.name
            )
            return False
        if self.usingValidationBounds:
            logger.error(
                "Parameter %s is set as a log level with validation bounds. This dev can't imagine a valid scenario for that. Let him know if you found one."
                % self.name
            )
            return False
        if not self.typeRequirement == str:
            logger.critical(
                "Logger level parameter %s is not set to expect a string value. This should be the only value type it takes."
                % self.name
            )
            return False
        logLevels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        if self.validationList:
            self.validationList = [item.upper() for item in self.validationList]
            testSet = set(self.validationList)
            levelSet = set(logLevels)
            if not testSet.issubset(levelSet):
                logger.error(
                    "Invalid logging levels given for log level validatio list on parameter %s"
                    % self.name
                )
                return False
        else:
            self.validationList = logLevels
        return True

    def validateExpectedFilePath(self):
        if self.expectedFile:
            if not os.path.isfile(self.value):
                logger.critical(
                    "Unable to find expected file for environment variable parameter %s at %s"
                    % (self.name, self.value)
                )
                raise FileNotFoundError("Unable to find expected file %s" % self.value)
        else:
            logger.error(
                "Validating expected file path when not an expected file for parameter %s"
                % self.name
            )

    def validateCreatedFile(self):
        if self.createdFile:
            if os.path.isfile(self.value):
                logger.critical(
                    "Environment variable parameter %s shows file %s as being created here, but it already exists.  If the file is expected to already exist and just be overwritten or modified, it should be flagged as an expected file only."
                    % (self.name, self.value)
                )
                raise FileExistsError(
                    "File %s created by parameter %s already exists."
                    % (self.value, self.name)
                )
            try:
                touchFile = open(self.value, "w")
                touchFile.close()
                os.remove(self.value)
            except (
                Exception
            ) as err:  # using a generic catch here because I don't care why it fails, just that it fails
                logger.critical(
                    "Unable to touch/create file set for environment variable parameter %s at %s. It generated an exception:\n%s"
                    % (self.name, self.value, err)
                )
                print("Walker started")
                for item in os.walk("/data"):
                    print(item)
                print("Walker done")
                raise ArgumentValueException(
                    "Unable to write to file %s as suggested by parameter %s"
                    % (self.value, self.name)
                )
        else:
            logger.error(
                "Testing file creation for what is not an expected file for parameter %s"
                % self.name
            )

    def validateExpectedDirectoryPath(self):
        if self.expectedDirectory:
            if not os.path.isdir(self.value):
                logger.critical(
                    "Unable to find expected directory for environment variable parameter %s at %s"
                    % (self.name, self.value)
                )
                raise NotADirectoryError("Unable to find expected file %s" % self.value)
        else:
            logger.error(
                "Validating directory existence for what is not an expected directory for parameter %s"
                % self.name
            )

    def createDirectory(self):
        if self.createdDirectory:
            try:
                os.makedirs(self.value)
            except FileExistsError:
                logger.info(
                    "Attempted to make directory %s for parameter %s, but it already appears to exist."
                    % (self.value, self.name)
                )
            except (
                Exception
            ) as err:  # using a generic catch here because I don't care why it fails, just that it fails
                logger.critical(
                    "Unable to touch/create directory set for environment variable parameter %s at %s. It generated an exception:\n%s"
                    % (self.name, self.value, err)
                )
                raise ArgumentValueException(
                    "Unable to create directory %s as suggested by parameter %s"
                    % (self.value, self.name)
                )

    def validateTypeAndFlag(self):
        if self.typeRequirement == bool:
            if self.positionalArg:
                logMessage = (
                    "Error on setting up environmental variable argument %s: boolean args cannot be positional."
                    % self.name
                )
                logger.exception(logMessage)
                raise ArgumentTypeValidationError()

    def passedArgumentAssertions(self):
        failures = {}
        assertion = isinstance(self.name, str)
        if assertionFails(assertion):
            failures["NameTypeCheck"] = (
                "Name value must be of string type. Given value was %s of type %s"
                % (self.name, type(self.name))
            )
        assertion = self.name != ""
        if assertionFails(assertion):
            failures["NameSetCheck"] = "Name value cannot be a blank"
        assertion = self.hasValidEnvironmentVariableName()
        if assertionFails(assertion):
            failures["EnvironmentVariableName"] = (
                "Name value must be directly translated to an environment variable. Valid characters include alphanumerics and underscores. The name cannot begin with a digit. %s was given as the name"
                % self.name
            )
        assertion = self.hasValidTypeRequirement()
        if assertionFails(assertion):
            failures["TypeRequirementCheck"] = (
                "Invalid type requirement given. Value passed: %s"
                % self.typeRequirement
            )
        if not self.default is None:
            assertion = self.fitsTypeRequirement(self.default)
            if assertionFails(assertion):
                failures["DefaultValueTypeCheck"] = (
                    "Invalid data type given for default. Valid types: %s. Default value: %s: %s"
                    % (self.typeRequirement, type(self.default), self.default)
                )
        if self.usingValidationBounds:
            assertion = self.lowerBound <= self.upperBound
            if assertionFails(assertion):
                failures["BoundsSanityCheck"] = (
                    "Invalid bounds given for environment variable parameter %s; lower value is greater than upper. Lower: %s.  Upper: %s"
                    % (self.name, self.lowerBound, self.upperBound)
                )
        if self.usingValidationList:
            assertion = self.validValidationListValues()
            if assertionFails(assertion):
                failures["ValidationValuesSanityCheck"] = (
                    "Got validation list elements that are not of an expected data type for environment variable parameter %s.  Validation list: %s.  Type requirement: %s."
                    % (self.name, self.validationList, self.typeRequirement)
                )
        if self.isArgument:
            assertion = type(self.flag) in [str, int]
            if assertionFails(assertion):
                failures["FlagDataType"] = (
                    "Flag data type for environment variable parameter %s was %s and value was %s. Acceptable types for this value are string and integer, or None for a non-argument parameter."
                    % (self.name, type(self.flag), self.flag)
                )
        assertion = not (self.expectedFile and self.createdFile)
        if assertionFails(assertion):
            failures["FilePathFlags"] = (
                "Environment variable parameter %s was set as both an expected file and a created file.  This should not be possible.  If the file is expected from the start and will be modified/overwritten, then mark it ONLY as expected."
                % self.name
            )
        assertion = not (self.expectedDirectory and self.createdDirectory)
        if assertionFails(assertion):
            failures["DirectoryPathFlags"] = (
                "Environment variable parameter %s was set as both an expected directory and a created directory.  This should not be possible.  If the directory is expected from the start and will be modified/overwritten, then mark it ONLY as expected."
                % self.name
            )
        assertion = not (self.isFilePath and self.isDirectoryPath)
        if assertionFails(assertion):
            failures["FileAndDirectoryPathFlags"] = (
                "Environment variable parameter %s has flags suggesting that it is both a file path and a directory path. This is not possible and needs to be corrected."
                % self.name
            )
        if self.logLevel:
            assertion = self.validLogLevelSettings()
            if assertionFails(assertion):
                failures["LogLevelParameterSetting"] = (
                    "Invalid parameter set for log level parameter. Please see above for issue on paramater name %s"
                    % self.name
                )
        if failures:
            failureString = ""
            for failure in failures:
                failureString += failure + "\n"
                failureString += failures[failure] + "\n"
            logger.error(
                "Detected %s errors setting up environment variable parameters:\n%s"
                % (len(failures), failureString)
            )
            print("Parameter errors detected for %s:\n%s" % (self.name, failureString))
        if failures:
            return False
        else:
            return True

    def validValidationListValues(self):
        for item in self.validationList:
            if not self.fitsTypeRequirement(item):
                return False
        return True

    def hasValidTypeRequirement(self):
        validTypesForList = [int, float, str]
        validTypesForSingle = validTypesForList.copy()
        validTypesForSingle.append(bool)
        if isinstance(self.typeRequirement, type):
            if self.typeRequirement in validTypesForSingle:
                return True
            else:
                return False
        else:
            if not isinstance(self.typeRequirement, collections.Iterable):
                return False
            else:
                for item in self.typeRequirement:
                    if not isinstance(item, type):
                        return False
                    if item not in validTypesForList:
                        return False
        return True

    def fitsTypeRequirement(self, value):
        if isinstance(self.typeRequirement, type):
            return isinstance(value, self.typeRequirement)
        else:
            return type(value) in self.typeRequirement

    def hasValidEnvironmentVariableName(self):
        envVar = self.name.upper()
        validCharacters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_"
        if envVar[0].isdigit():
            return False
        for character in envVar:
            if character not in validCharacters:
                return False
        return True

    def setEnvironmentValue(self):
        try:
            value = os.environ[self.environmentVariableName]
        except KeyError:
            if self.required:
                logger.error(
                    "Environment variable parameter %s was set as required, but no value was passed for it."
                    % self.name
                )
            value = self.value
            usingEnvironment = False
        else:
            usingEnvironment = True
        if value is None:
            return False
        if self.logLevel:
            value = value.upper()
        if self.typeRequirement == bool:
            value = self.setBooleanValue(value)
        else:
            value = self.setType(value)
        if self.usingValidationBounds:
            if value < self.lowerBound or value > self.upperBound:
                logger.warning(
                    "Got an out of bounds parameter set: %s set at %s. LowerBound: %s UpperBound: %s"
                    % (self.name, value, self.lowerBound, self.upperBound)
                )
        if self.usingValidationList:
            if value not in self.validationList:
                logger.warning(
                    "Got a parameter being set that is not on the validation list: %s set at %s. ValidationList: %s"
                    % (self.name, value, self.validationList)
                )
        if self.logLevel:
            value = self.setLogLevel(value)
        self.value = value
        return usingEnvironment

    def setLogLevel(self, value):
        valueTable = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
            "CRITICAL": logging.CRITICAL,
        }
        return valueTable[value]

    def setBooleanValue(self, value: str):
        if not value:
            return False
        elif self.value in ["FALSE", "false", "False", "0", 0]:
            return False
        else:
            return True

    def setType(self, value):
        if isinstance(self.typeRequirement, type):
            try:
                return self.typeRequirement(value)
            except Exception as err:
                logMessage = (
                    "Attempting to cast environment variable %s with value %s to type %s resulted in an exception as follows: \n%s"
                    % (self.name, value, self.typeRequirement, err)
                )
                logger.exception(logMessage)
                raise ArgumentTypeValidationError()
        if isinstance(self.typeRequirement, type):
            allowedTypes = [self.typeRequirement]
        else:
            allowedTypes = self.typeRequirement
        allowedTypes = [self.typeRequirement]
        for valueType in typeHeirarchy:
            if valueType in allowedTypes:
                try:
                    return valueType(value)
                except:  # Using a universal catch here, since I'm using this as a test and expect it to fail regularly
                    continue
        logger.error(
            "Unable to cast environment variable parameter %s to one of its required types: %s.  Env variable value: %s"
            % (self.name, self.typeRequirement, value)
        )
        raise ArgumentTypeValidationError(
            "Unable to cast environment variable parameter %s to one of its required types: %s.  Env variable value: %s"
            % (self.name, self.typeRequirement, value)
        )

    def formArgument(self):
        if not self.isArgument:
            return ""
        if self.value is None:
            return ""
        if self.typeRequirement == bool:
            return self.parseBooleanArg()
        if self.positionalArg:
            return str(self.value)
        if isinstance(self.value, str) and self.value.startswith("="):
            return "%s%s" % (self.flag, self.value)
        else:
            return "%s %s" % (self.flag, self.value)

    def parseBooleanArg(self):
        if not self.value:
            return ""
        elif self.value in ["FALSE", "false", "False", "0", 0]:
            return ""
        else:
            return self.flag

    def overview(self):
        returnDict = {
            "name": self.name,
            "type": type(self.value),
            "value": self.value,
            "flag": self.flag,
        }
        return str(returnDict)

    def __str__(self):
        return str(self.value)

    def __eq__(self, other):
        return self.value == other

    def __bool__(self):
        if isinstance(self.value, bool):
            return self.value
        else:
            return self.value is not None


class ParameterSideLoad(EnvVariable):

    def setEnvironmentValue(self):
        return False


class EnvParameters(object):
    """
    How to use:
    Initialize an empty parameter set with myInstance = EnvParameters()
    Add values that check environment variable using the following syntax:
    myInstance.addParameter(name, type, [optional values])
    add values in directly using the side load method when additional logic at checking time is required:
    myInstance.sideloadParameter(name, value, [optional values])
    """

    def __init__(self):
        self.parameters = {}
        self.variableNames = set()
        self.flags = set()

    def addParameter(
        self,
        name: str,
        typeRequirement: [type, list, tuple],
        default=None,
        flag: [int, str] = None,
        validationList: list = None,
        lowerBound: [int, float] = None,
        upperBound: [int, float] = None,
        expectedFile: bool = False,
        createdFile: bool = False,
        expectedDirectory: bool = False,
        createdDirectory: bool = False,
        logLevel: bool = False,
        required: bool = False,
        externalValidation: bool = False,
    ):
        parameter = EnvVariable(
            name,
            typeRequirement,
            default,
            flag,
            validationList,
            lowerBound,
            upperBound,
            expectedFile,
            createdFile,
            expectedDirectory,
            createdDirectory,
            logLevel,
            required,
            externalValidation,
        )
        if parameter.environmentVariableName not in self.variableNames:
            self.variableNames.add(parameter.environmentVariableName)
        else:
            logger.critical(
                "Environment variable name collision for %s"
                % parameter.environmentVariableName
            )
            raise ArgumentValueException(
                "Environment variable name collision for %s"
                % parameter.environmentVariableName
            )
        if parameter.isArgument:
            if parameter.flag not in self.flags:
                self.flags.add(parameter.flag)
            else:
                logger.critical(
                    "Environment variable argument flag collision for %s"
                    % parameter.flag
                )
                raise ArgumentValueException(
                    "Environment variable flag collision for %s" % parameter.flag
                )
        self.parameters[parameter.name] = parameter

    def sideLoadParameter(
        self,
        name: str,
        value,
        flag: [int, str] = None,
        expectedFile: bool = False,
        createdFile: bool = False,
        expectedDirectory: bool = False,
        createdDirectory: bool = False,
    ):
        parameter = ParameterSideLoad(
            name,
            type(value),
            default=value,
            validationList=[value],
            expectedFile=expectedFile,
            createdFile=createdFile,
            expectedDirectory=expectedDirectory,
            createdDirectory=createdDirectory,
        )
        if not parameter.environmentVariableName in self.variableNames:
            self.variableNames.add(parameter.environmentVariableName)
        else:
            logger.critical(
                "Environment variable name collision for %s on side load"
                % parameter.environmentVariableName
            )
            raise ArgumentValueException(
                "Environment variable name collision for %s"
                % parameter.environmentVariableName
            )
        if parameter.isArgument:
            if not parameter.flag in self.flags:
                self.flags.add(parameter.flag)
            else:
                logger.critical(
                    "Environment variable argument flag collision for %s"
                    % parameter.flag
                )
                raise ArgumentValueException(
                    "Environment variable flag collision for %s" % parameter.flag
                )
        self.parameters[parameter.name] = parameter

    def buildFlaggedArgumentString(self):
        flaggedArgs = []
        for key in self.parameters:
            parameter = self.parameters[key]
            if parameter.isArgument and not parameter.positionalArg:
                flaggedArgs.append(parameter.formArgument())
        return " ".join(flaggedArgs)

    def buildPositionalArgumentStrings(self):
        import operator

        prependArgs = []
        appendArgs = []
        for key in self.parameters:
            parameter = self.parameters[key]
            if parameter.isArgument and parameter.positionalArg:
                if parameter.flag >= 0:
                    prependArgs.append(parameter)
                else:
                    appendArgs.append(parameter)
        if not (prependArgs or appendArgs):
            return ("", "")
        if prependArgs:
            prependArgs.sort(key=operator.attrgetter("flag"))
            prependArgs = [arg.formArgument() for arg in prependArgs]
        if appendArgs:
            appendArgs.sort(key=operator.attrgetter("flag"))
            appendArgs = [arg.formArgument() for arg in appendArgs]
        prependArgString = " ".join(prependArgs)
        appendArgString = " ".join(appendArgs)
        return (prependArgString, appendArgString)

    def buildArgString(self):
        beginning, end = self.buildPositionalArgumentStrings()
        middle = self.buildFlaggedArgumentString()
        argList = [item for item in [beginning, middle, end] if item]
        return " ".join(argList)

    def checkCreatedFileStructures(self):
        for key in self.parameters:
            parameter = self.parameters[key]
            if parameter.createdDirectory:
                parameter.createDirectory()
        for key in self.parameters:
            parameter = self.parameters[key]
            if parameter.createdFile:
                parameter.validateCreatedFile()

    def __getattr__(self, item):
        if item in self.parameters:
            return self.parameters[item]
        else:
            if not isinstance(item, str):
                print(list(self.parameters.keys()))
                raise AttributeError(
                    "No parameter %s was found in the parameter set" % item
                )
            for key in self.parameters:
                if not isinstance(key, str):
                    continue
                keylower = key.lower()
                itemlower = item.lower()
                if keylower == itemlower:
                    return self.parameters[key]
        print(list(self.parameters.keys()))
        raise AttributeError("No parameter %s was found in the parameter set" % item)


class ArgumentTypeValidationError(Exception):
    pass


class ArgumentValueException(Exception):
    pass


class EnvironmentVariableParameterException(Exception):
    pass


def assertionFails(bool_: bool):
    try:
        assert bool_, "Critical assertion failed."
    except AssertionError:
        return True
    else:
        return False


if __name__ == "__main__":
    test = EnvParameters()
    test.addParameter("first", str, default="The first", validationList=["The first"])
    test.sideLoadParameter("sideload", "The side loaded one")
    test.addParameter(
        "second", str, "The second one", validationList=["The second one"]
    )
