from typing import Tuple as T, Optional as O
def parse_csv(root:str, mat:str, volrat: float
             ) -> T[O[float],O[float],O[float],O[float],O[float]]:
    '''From the material name, assemble experimental data.'''
    from csv import reader
    import ase.units as units# type: ignore
    ce = bm = lat = vol = mag = None

    with open(root%'rungs_tran','r') as f:
        r = reader(f)
        for mat_, lat_, bm_, ce_ in r:
            if mat_ == mat:
                lat, bm, ce = map(float, [lat_, bm_, ce_])
                vol = (lat * float(volrat))**3

    with open(root%'errorestimate_lejaeghere','r') as f:
        r = reader(f)
        kJ_mol_to_eV = 1.0364e-2
        for mat_, ce_, bm_, vol_ in r:
            if mat_ == mat and ce is None:
                ce, bm, vol = map(float, [ce_, bm_, vol_])
                ce = kJ_mol_to_eV * ce
                lat = vol**(1/3)/float(volrat)

    with open(root%'cohesive_guillermet','r') as f:
        r = reader(f)
        mRyd_to_eV = 0.013605691455111256
        for mat_, ce_ in r:
            if mat_ == mat and ce is None:
                ce = mRyd_to_eV * float(ce_)

    with open(root%'sol58lp_54coh','r') as f:
        r = reader(f)
        keys = ['cohesive energy','structure','debye temperature','lattice parameter','bulk modulus','magmom']
        for mat_, _, _, mag_, lp_, _, debye, ce_, bm_, corr, unit in r:
            if mat_ == mat:

                if mag_:
                    mag = float(mag_)

                if ce is None and ce_:
                    if unit == 'kcal/mol':
                        ce = float(ce_) * units.kcal/units.mol
                    elif unit == 'kJ/mol':
                        ce = float(ce_) * units.kcal/units.mol
                    else: raise ValueError
                    if corr == 'False':
                        ce += (9./8.)*units.kB*float(debye)
                if bm is None and bm_:
                    bm = float(bm_) 
                if lat is None and lp_:
                    lat = float(lp_)
                    vol = (lat * float(volrat))**3

    return ce, bm, lat, vol, mag
