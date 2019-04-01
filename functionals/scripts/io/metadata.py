# External Modules
from typing     import Tuple as T
from os         import stat
from pwd        import getpwuid
from os.path    import getmtime
from datetime   import date

################################################################################
def metadata(stordir : str) -> T[str,date]:
    """
    Takes a path to a DFT calculation and extracts information from runtime.json
    Also keeps a record of where the information was taken from.
    """
    usr    = getpwuid(stat(stordir).st_uid).pw_name
    tstamp = date.fromtimestamp(getmtime(stordir))

    # WE DON'T SAVE ANYTHING TO IDENTIFY WHETHER JOB IS VIBJOB!!!
    return (usr,tstamp) # type: ignore

if __name__=='__main__':
    import sys
    print(metadata(sys.argv[1]))
