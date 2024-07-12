import multiprocessing
import multiprocessing.pool
import logging

logger = logging.getLogger(__name__)


def calculateAvailableCores():
    import multiprocessing

    return max([multiprocessing.cpu_count() - 1, 1])


def calculateChunkSize(length, workers: int):
    return -1 * ((-1 * length) // workers)


class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess


class Deadpool(multiprocessing.pool.Pool):  # Deadpool has no class
    def __init__(self, *args, **kwargs):
        kwargs["context"] = NoDaemonContext()
        super(Deadpool, self).__init__(*args, **kwargs)


def parallelProcessRunner(
    processor,
    itemsToProcess,
    coreLimit: int = 0,
    filterFunction=False,
    totalSizeEstimate=None,
    coresPerProcess=1,
    nonDaemonic=False,
):
    logger.debug("Running import statements")
    import multiprocessing
    import inspect
    import collections

    logger.debug("Making assertions")
    assert callable(processor), "Processor must be a callable function/method"
    assert (
        len(inspect.signature(processor).parameters) == 1
    ), "Processor function must take one argument"
    assert isinstance(
        itemsToProcess, collections.Iterable
    ), "Items to process must be an iterable of some kind"
    assert coresPerProcess > 0, "Cores per process must be a positive integer"
    logger.debug("Calculating cores available")
    coreLimit = max([0, coreLimit])
    if not coreLimit:
        coreLimit = calculateAvailableCores()
    coreLimit = coreLimit // coresPerProcess
    logger.info("Using %s cores" % coreLimit)
    if nonDaemonic:
        logger.debug(
            "Setting pool using nonDaemonic processes. This can cause issues if not done with care."
        )
        workers = Deadpool(coreLimit)
    else:
        logger.debug("Setting up the process pool.")
        workers = multiprocessing.Pool(coreLimit)
    try:
        length = len(itemsToProcess)
        chunkSize = calculateChunkSize(length, coreLimit)
        logger.info(
            "Starting multiprocessing of %s objects in chunks of %s"
            % (length, chunkSize)
        )
        mapper = workers.map
    except TypeError:
        if not totalSizeEstimate:
            totalSizeEstimate = 50000
            logger.info(
                "Using default total size estimate of 50,000 because none was given."
            )
        chunkSize = calculateChunkSize(totalSizeEstimate, coreLimit)
        logger.info("Starting multiprocessing in chunks of %s" % chunkSize)
        mapper = workers.imap
    if not filterFunction:
        logger.debug("Returning results")
        return mapper(processor, itemsToProcess, chunkSize)
    else:
        logger.debug("Returning filtered results")
        results = mapper(processor, itemsToProcess, chunkSize)
        return [result for result in results if result]


if __name__ == "__main__":  # absolutely necessary for Windows machines
    import datetime

    class RandomDNASequenceMaker(object):

        def __init__(self, length):
            self.length = length
            self.bases = list("ATGC")

        def makeSequence(self, throwAway=0):
            import random

            sequence = []
            for i in range(self.length):
                sequence.append(random.choice(self.bases))
            return "".join(sequence)

    seqGen = RandomDNASequenceMaker(100)
    numberOfSequences = list(range(100000))
    start = datetime.datetime.now()
    sequenceCollector = []
    for i in numberOfSequences:
        sequenceCollector.append(seqGen.makeSequence(i))
    end = datetime.datetime.now()
    print("Serial process completed in %s" % (end - start))
    start = datetime.datetime.now()
    multiResult = parallelProcessRunner(seqGen.makeSequence, numberOfSequences)
    end = datetime.datetime.now()
    print("Parallel process completed in %s" % (end - start))
