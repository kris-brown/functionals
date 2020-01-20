from typing import Tuple as T, Optional as O


def parse_csv(root: str, mat: str, volrat: float
              ) -> T[O[float], O[float], O[float], O[float]]:
    '''From the material name, assemble experimental data.'''
    from csv import reader
    import ase.units as units
    # if any([x in mat for x in ['Os', 'Re', 'Pb', 'Ir', 'La', "Pt", "Au"]]):
    #     print('Skipping expt data for ', mat)
    #     return None, None, None, None

    ce = bm = lat = mag = None

    with open(root % 'rungs_tran', 'r') as f:
        r = reader(f)
        for mat_, lat_, bm_, ce_ in r:
            if mat_ == mat:
                lat, bm, ce = map(float, [lat_, bm_, ce_])

    with open(root % 'errorestimate_lejaeghere', 'r') as f:
        r = reader(f)
        for mat_, ce_, bm_, _ in r:
            if mat_ == mat and ce is None:
                ce, bm = map(float, [ce_, bm_])
                ce *= units.kJ/units.mol

    with open(root % 'cohesive_guillermet', 'r') as f:
        r = reader(f)
        mRyd_to_eV = 0.013605691455111256
        for mat_, ce_ in r:
            if mat_ == mat and ce is None:
                ce = mRyd_to_eV * float(ce_)

    with open(root % 'sol58lp_54coh', 'r') as f:
        r = reader(f)
        for mat_, _, _, mag_, lp_, _, debye, ce_, bm_, corr, unit in r:
            if mat_ == mat:

                if mag_:
                    mag = float(mag_)

                if ce is None and ce_:
                    if unit == 'kcal/mol':
                        ce = float(ce_) * units.kcal/units.mol
                    elif unit == 'kJ/mol':
                        ce = float(ce_) * units.kJ/units.mol
                    else:
                        raise ValueError
                    if corr == 'False':
                        ce += (9./8.)*units.kB*float(debye)
                if bm is None and bm_:
                    bm = float(bm_)
                if lat is None and lp_:
                    lat = float(lp_)
    if False:  # add keld data
        with open(root % 'misc', 'r') as f:
            r = reader(f)
            for mat_, prop, val, unit, _, _, _ in r:
                if mat_ == mat:
                    if prop == 'lattice':
                        assert unit == 'A'
                        lat = float(val)
                    elif prop == 'cohesive':
                        assert unit == 'eV'
                        ce = float(val)
                    elif prop == 'magmom':
                        assert unit == 'bohr'
                        mag = float(val)

    return ce, bm, lat, mag
