from json import dumps


def parse_contribs_vasp(outcar: str) -> str:
    '''
    Reads outcar, gets a JSON structured (64,1,1) data
    '''
    start = outcar.find('BEEF xc energy contributions')
    lines = outcar[start:].split('\n')[1:65]
    contribs = [float(line.split(':')[1]) for line in lines]
    return dumps(contribs)
