from json import dumps

def parse_contribs_vasp(outcar : str)->str:
    '''
    Reads outcar, gets a JSON structured (64x5,1,1) data
    '''
    start    = outcar.find('BEEF xc energy contributions')
    lines    = outcar[start:].split('\n')[1:321]
    contribs = [float(line.split(':')[1]) for line in lines]
    # assuming that a1 isn't varied
    # assert contribs[:64] == contribs[64:128]
    return dumps(contribs[:64])
