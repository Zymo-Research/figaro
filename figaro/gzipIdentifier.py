import os
import gzip
import binascii


def isGzipped(path: str):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            "Unable to determine if file %s is gzipped because that file does not exist."
            % path
        )
    file = open(path, "rb")
    firstTwoBytes = file.read(2)
    file.close()
    if not binascii.hexlify(firstTwoBytes) == b"1f8b":
        return False
    try:
        file = gzip.open(path, "rb")
        tenBytes = file.read(10)
    except OSError:
        return False
    return True
