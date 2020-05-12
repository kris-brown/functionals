from typing import Tuple as T


def parse_eigenval(fi: str) -> T[bool, float]:
    '''Parse Vasp EIGENVAL - check if all bands have integer occupation'''
    magmom = 0.
    lines = fi.split('\n')
    splitlines = list(filter(None, [x.split() for x in reversed(lines)]))
    if len(splitlines[0]) == 3:
        return True, 0.  # spin unpolarized!

    ints = True
    mags = [(float(l[-2]), float(l[-1])) for l in splitlines
            if len(l) == 5 and str.isdigit(l[0])]

    for up, down in mags:
        nonints = [abs(round(a)-a) > 0.01 for a in [up, down]]
        magmom += abs(up - down)
        if any(nonints):
            ints = False  # flip flag if at any point this is true

    return ints, magmom  # final band (reverse order)


if __name__ == '__main__':
    # DOESN'T GET THE SAME AS grepping MAG in OUTCAR for bulks
    for fi in ['atoms/pbe/Au']:  # 'atoms/pbe/W',
        with open('/Users/ksb/scp_tmp/vauto/%s/EIGENVAL' % fi, 'r') as f:
            print(parse_eigenval(f.read()))
