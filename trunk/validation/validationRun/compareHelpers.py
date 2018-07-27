from subprocess import Popen
from shlex import split

def compareFiles(filea, fileb):
    child = Popen(split("diff " + filea + " " + fileb))
    child.communicate()[0]
    returnValue = child.returncode
    return returnValue