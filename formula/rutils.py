""" Utilities for working with R """
import shlex
from subprocess import Popen, PIPE


def r_cmd(cmd):
    """ Execute command `cmd` through R, return stdout, stderr
    """
    proc = Popen(shlex.split('R --slave'),
                 stdin = PIPE,
                 stdout = PIPE,
                 stderr = PIPE)
    return proc.communicate(cmd)


def checkfor_r():
    """ Return True if R can be called without error """
    try:
        out, err = r_cmd('')
    except OSError:
        return False
    return (out, err) == ('', '')


def r_check_library(libname):
    """ Return True if R library `libname` is installed
    """
    out, err = r_cmd("library('%s')" % libname)
    if len(out) != 0:
        return False
    if "Error" in err:
        return False
    return True
