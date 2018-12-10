from typing     import List as L
from subprocess import getstatusoutput


################################################################################
def get_stray_gpaw(rootpath : str)-> L[str]:
    """
    Searches for GPAW calculation folders.
    """
    ######################
    # Initialize Variables
    #---------------------
    logfiles  = []
    badext     = ["sh", "traj", "json", "gpw", "py", "err", "xyz"]
    badextstr  = " ".join(['-not -name "*.%s"'%x for x in badext])

    ##############
    # Main program
    #-------------
    cmd = """for f in $(find {0}  -type f {1})
            do
                if $(head -2 $f | tail -1 | grep -q "___ ___ ___ _ _ _"); then
                    if $(grep -q "Free energy" $f); then
                        echo $f
                    fi
                fi
        done""".format(rootpath,badextstr)
    exit,gpawcheck = getstatusoutput(cmd)
    assert exit == 0,'Failure in get_gpaw_logfile bash execution: '+gpawcheck
    for logfile in gpawcheck.split('\n'):
        if logfile:
            dir = logfile[:logfile.rfind('/')]
            logfiles.append(logfile)

    return logfiles

if __name__=='__main__':
    import sys
    print(get_stray_gpaw(sys.argv[1]))
