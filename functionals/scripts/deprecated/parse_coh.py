from typing import Tuple as T, List as L, Optional as O

def parse_coh(pth: str, name: str) -> T[O[float], O[float], O[float], O[float]]:
    from csv import DictReader
    from ase.units import kB,kJ,kcal,mol # type: ignore
    assert 'hcp' not in name
    name = name.split('_')[0]
    keys = ['cohesive energy','structure','debye temperature','lattice parameter','bulk modulus','magmom']
    with open(pth,'r') as f:
        reader = DictReader(f)
        for row in reader:
            if row['solid_data'] == name:
                ce,struct,db_temp,l,bm,mag = map(row.get,keys)
                if ce:
                    if row['Ecoh unit']=='kcal/mol':
                        ce *= kcal/mol #junits.kcal/units.mol #unirs.kJ/units.mol
                    elif row['Ecoh unit']=='kJ/mol':
                        ce *= kJ/mol
                    else:
                        raise ValueError(row)

                    if row['Ecoh is phonon-corrected']=='False':
                        ce += (9./8.)*kB*float(db_temp) # type: ignore

                return ce, bm, l, mag # type: ignore
    raise ValueError(name+' not found!')


'''
    elif systype=='solid':
        struct = row[1]
        if struct=='hcp':
            return
        if struct in ['fcc', 'bcc']: atomspercell = 1.
        elif struct in ['diamond', 'hcp', 'rocksalt', 'zincblende', 'cesiumchloride']: atomspercell = 2.
        else:
            raise RuntimeError, 'unknown crystal structure ' + struct
        if struct in ['fcc', 'diamond', 'rocksalt', 'zincblende']:
            vol = float(row[4])**3 / 4.
        elif struct in ['bcc']:
            vol = float(row[4])**3 / 2.
        elif struct in ['cesiumchloride']:
            vol = float(row[4])**3
        solid = row[0]+'_'+struct
        no_coh = (solid in ['FeAl_cesiumchloride', 'BAs_zincblende', 'CoAl_cesiumchloride', 'NiAl_cesiumchloride'])
        if not no_coh:
           zpe_corrected = row[-2]
           e_coh = float(row[-4])
           unit = row[-1]
           if unit=='kcal/mol':
               e_coh *= kcpmol
           elif unit=='kJ/mol':
               e_coh *= kjpmol
           else:
               raise RuntimeError, 'unknown unit ' + unit

           if not zpe_corrected:
               debyeT = float(row[-5])
               e_coh += (9./8.)*units.kB*debyeT
           refdict[solid+':coh'] = [atomspercell, e_coh]

        if row[-3] is not None and row[-3].strip()!='':
            bulkmod = float(row[-3]) * units.GPa
            refdict[solid+':bulkmod'] = [bulkmod]
        refdict[solid+':latt'] = [vol]

'''
