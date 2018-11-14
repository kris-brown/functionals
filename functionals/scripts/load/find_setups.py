from typing import List
import re
############################
def find_setups(log:str)->List[str]:
    '''Finds the checksum values for all setups used in a GPAW calculation'''

    s = r'((id: )(([a-z]|\d)+)(\s))'

    groups = re.findall(s,log)

    return [g[2] for g in groups]

if __name__=='__main__':
    import sys
    with open(sys.argv[1],'r') as f:
        print(find_setups(f.read()))
