from typing import List as L, Tuple as T


def pop_surf() -> T[L[str], L[str], L[float], L[float]]:
    import os
    import re
    root = '/Users/ksb/scp_tmp/vauto/surf'
    mats, xcs, ots, hos = [], [], [], []

    for xc in os.listdir(root):
        for mat in os.listdir(os.path.join(root, xc)):

            def geteng(site: str) -> float:
                out = os.path.join(root, xc, mat, site, 'OUTCAR')
                with open(out, 'r') as f:
                    outcar = f.read()
                assert 'General timing' in outcar
                pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
                match = re.findall(pat, outcar)
                assert match
                return float(match[-1])
            try:
                oe, he = map(geteng, ['ontop', 'fcc'])
                mats.append(mat)
                xcs.append(xc)
                ots.append(oe)
                hos.append(he)
            except (AssertionError, FileNotFoundError):
                pass
    return mats, xcs, ots, hos


if __name__ == '__main__':
    for x in pop_surf():
        print(x)
