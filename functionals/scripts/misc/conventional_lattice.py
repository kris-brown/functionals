def conventional_lattice(struct : str, vol_pa : float)->float:
    '''
    Struct = AFLOW prototype
    vol    = Volume of calculated structure at minimum energy
    n      = # of atoms

    Returns the lattice constant for a conventional unit cell

    '''
    nick_dict = {'AB_1_a_b_225'       : ('rocksalt',        'cubic',8),
                 'AB_1_a_c_216'       : ('zincblende',      'cubic',8),
                 'AB_1_a_b_221'       : ('cesium chloride', 'cubic',2),
                 'A_1_a_225'          : ('fcc',             'cubic',4),
                 'A_1_a_229'          : ('bcc',             'cubic',2),
                 'A_2_c_194'          : ('hcp',             'hexagonal',2),
                 'A_2_a_227'          : ('diamond',         'cubic',8),
                 'AB3_1_a_d_221'      : ('anti-ReO3',       'cubic',0),
                 'A2B3_8_ad_e_206'    : ('antibixbyite',    'cubic',0),
                 'AB3C_cP5_221_a_c_b' : ('perovskite',      'cubic',0)
                } # rutile?

    corrdict = {'cubic':1,'hexagonal':1/1.633}
    nick,kind,n_unit = nick_dict[struct] # error if we have a new structure to consider
    assert n_unit, "Need to figure out how many "

    vol_ = float(vol_pa * n_unit) # get a cell vol for conventional # of atoms
    correction = corrdict[kind]

    return (vol_ * correction)**(1/3)
