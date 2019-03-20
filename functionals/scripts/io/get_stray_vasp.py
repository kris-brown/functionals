from typing import List as L
from subprocess import getstatusoutput

def get_stray_vasp(root:str)->L[str]:
    '''Looks for completed vasp jobs'''

    cmd = '''for f in $(find {}  -type d)
        do
            oz="$f/OSZICAR"
            out="$f/OUTCAR"
            inc="$f/INCAR"

            if [[ -e $out && $(grep "General timing" $out) ]]; then

                step="$(wc -l < $oz)"
                max="$(sed -n -e 's/NELM = //p' $inc)"

            if [[ $step -lt $max ]]; then

               echo $f

             fi
        fi
    done
    '''.format(root)

    exit,vaspcheck = getstatusoutput(cmd)
    assert exit == 0,'Failure in get_vasp_logfile bash execution: '+vaspcheck
    return [fi for fi in vaspcheck.split('\n') if fi]
