import sys

def MSG(msg_str):
    """
    Write a helpful progress message to stderr

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    sys.stderr.write("[gwas]: %s\n"%msg_str.strip())


def ERROR(msg_str):
    """
    Write an error message to stderr then quit
    with exit code 1

    Arguments
    ---------
    msg_str : str
       Message to write
    """
    MSG("ERROR: " + msg_str)
    sys.exit(1)