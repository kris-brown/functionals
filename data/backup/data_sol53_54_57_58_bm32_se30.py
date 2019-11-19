from typing import Optional as O
from ase.atoms import Atoms, string2symbols   # type: ignore
from ase.units import kB, kJ, kcal, mol   # type: ignore
from ase.io import write    # type: ignore
import numpy as np   # type: ignore

# Sol57LC
sol57_lp    = ['Li_bcc','Na_bcc','K_bcc','Rb_bcc','Ca_fcc','Sr_fcc','Ba_bcc',
               'V_bcc','Nb_bcc','Ta_bcc','Mo_bcc','W_bcc','Fe_bcc',
               'Rh_fcc','Ir_fcc','Ni_fcc','Pd_fcc','Pt_fcc','Cu_fcc',
               'Ag_fcc','Au_fcc','Al_fcc',
               'C_diamond','Si_diamond','Ge_diamond','Sn_diamond',
               'LiH_b1','LiF_b1','LiCl_b1','NaF_b1','NaCl_b1','MgO_b1','MgS_b1',
               'CaO_b1','TiC_b1','TiN_b1','ZrC_b1','ZrN_b1','VC_b1','VN_b1',
               'NbC_b1','NbN_b1','FeAl_b2','CoAl_b2','NiAl_b2','BN_b3',
               'BP_b3','BAs_b3','AlN_b3','AlP_b3','AlAs_b3','GaN_b3',
               'GaP_b3','GaAs_b3','InP_b3','InAs_b3','SiC_b3']
assert len(sol57_lp) == 57
# Sol58LC = Sol57LC + Pb_fcc
sol58_lp    = ['Li_bcc','Na_bcc','K_bcc','Rb_bcc','Ca_fcc','Sr_fcc','Ba_bcc',
               'V_bcc','Nb_bcc','Ta_bcc','Mo_bcc','W_bcc','Fe_bcc',
               'Rh_fcc','Ir_fcc','Ni_fcc','Pd_fcc','Pt_fcc','Cu_fcc',
               'Ag_fcc','Au_fcc','Al_fcc','Pb_fcc',
               'C_diamond','Si_diamond','Ge_diamond','Sn_diamond',
               'LiH_b1','LiF_b1','LiCl_b1','NaF_b1','NaCl_b1','MgO_b1','MgS_b1',
               'CaO_b1','TiC_b1','TiN_b1','ZrC_b1','ZrN_b1','VC_b1','VN_b1',
               'NbC_b1','NbN_b1','FeAl_b2','CoAl_b2','NiAl_b2','BN_b3',
               'BP_b3','BAs_b3','AlN_b3','AlP_b3','AlAs_b3','GaN_b3',
               'GaP_b3','GaAs_b3','InP_b3','InAs_b3','SiC_b3']
assert len(sol58_lp) == 58
# Sol53Ec
sol53_coh    = ['Li_bcc','Na_bcc','K_bcc','Rb_bcc','Ca_fcc','Sr_fcc','Ba_bcc',
                'V_bcc','Nb_bcc','Ta_bcc','Mo_bcc','W_bcc','Fe_bcc',
                'Rh_fcc','Ir_fcc','Ni_fcc','Pd_fcc','Pt_fcc','Cu_fcc',
                'Ag_fcc','Au_fcc','Al_fcc',
                'C_diamond','Si_diamond','Ge_diamond','Sn_diamond',
                'LiH_b1','LiF_b1','LiCl_b1','NaF_b1','NaCl_b1','MgO_b1','MgS_b1',
                'CaO_b1','TiC_b1','TiN_b1','ZrC_b1','ZrN_b1','VC_b1','VN_b1',
                'NbC_b1','NbN_b1','BN_b3','BP_b3','AlN_b3','AlP_b3','AlAs_b3',
                'GaN_b3','GaP_b3','GaAs_b3','InP_b3','InAs_b3','SiC_b3']
assert len(sol53_coh) == 53
# Sol54Ec
sol54_coh    = ['Li_bcc','Na_bcc','K_bcc','Rb_bcc','Ca_fcc','Sr_fcc','Ba_bcc',
                'V_bcc','Nb_bcc','Ta_bcc','Mo_bcc','W_bcc','Fe_bcc',
                'Rh_fcc','Ir_fcc','Ni_fcc','Pd_fcc','Pt_fcc','Cu_fcc',
                'Ag_fcc','Au_fcc','Al_fcc','Pb_fcc',
                'C_diamond','Si_diamond','Ge_diamond','Sn_diamond',
                'LiH_b1','LiF_b1','LiCl_b1','NaF_b1','NaCl_b1','MgO_b1','MgS_b1',
                'CaO_b1','TiC_b1','TiN_b1','ZrC_b1','ZrN_b1','VC_b1','VN_b1',
                'NbC_b1','NbN_b1','BN_b3','BP_b3','AlN_b3','AlP_b3','AlAs_b3',
                'GaN_b3','GaP_b3','GaAs_b3','InP_b3','InAs_b3','SiC_b3']
assert len(sol54_coh) == 54
# BM32
bm32        = ['Li_bcc','Na_bcc','Ca_fcc','Sr_fcc','Ba_bcc',
               'Rh_fcc','Pd_fcc','Cu_fcc','Ag_fcc','Al_fcc',
               'C_diamond','Si_diamond','Ge_diamond','Sn_diamond',
               'LiH_b1','LiF_b1','LiCl_b1','NaF_b1','NaCl_b1','MgO_b1',
               'BN_b3','BP_b3','BAs_b3','AlN_b3','AlP_b3','AlAs_b3',
               'GaN_b3','GaP_b3','GaAs_b3','InP_b3','InAs_b3','SiC_b3']
assert len(bm32) == 32
# hcp
hcp_solids_10 = ['Cd_hcp', 'Co_hcp', 'Os_hcp', 'Ru_hcp', 'Zn_hcp',
                 'Zr_hcp', 'Sc_hcp', 'Be_hcp', 'Mg_hcp', 'Ti_hcp']
assert len(hcp_solids_10) == 10
# SE30
se21 = ['Li110','Na110','K110','Rb110','Ba110','Ca111','Sr111',
        'Nb110','Ta110','Mo110','W110','Fe110','Al111','Ni111',
        'Cu111','Rh111','Pd111','Ag111','Ir111','Pt111','Au111']
se9  = ['Mg0001','Zn0001','Cd0001','Sc0001','Ti0001',
        'Co0001','Zr0001','Ru0001','Os0001']
se30 = se21+se9
assert len(se30) == 30
#
#
#
#==========================================================
# SURFACE ENERGIES
#==========================================================
# see Vitos, Ruban, Skriver, Kollar, Surf Sci 411 (1998) 186-202
# experimental surface energies in J/m**2
#
#
#
surf_data = {
# (110) facets of bcc metals
'Li110': {
    'solid': 'Li_bcc',
    'facet': '110',
    'E_surf': [0.522, 0.525]},
'Na110': {
    'solid': 'Na_bcc',
    'facet': '110',
    'E_surf': [0.261, 0.260]},
'K110': {
    'solid': 'K_bcc',
    'facet': '110',
    'E_surf': [0.145, 0.130]},
'Rb110': {
    'solid': 'Rb_bcc',
    'facet': '110',
    'E_surf': [0.117, 0.110]},
'Cs110': {
    'solid': 'Cs_bcc',
    'facet': '110',
    'E_surf': [0.095, 0.095]},
'Ba110': {
    'solid': 'Ba_bcc',
    'facet': '110',
    'E_surf': [0.380, 0.370]},
'V110': {
    'solid': 'V_bcc',
    'facet': '110',
    'E_surf': [2.622, 2.550]},
'Fe110': {
    'solid': 'Fe_bcc',
    'facet': '110',
    'E_surf': [2.417, 2.475]},
'Nb110': {
    'solid': 'Nb_bcc',
    'facet': '110',
    'E_surf': [2.655, 2.700]},
'Ta110': {
    'solid': 'Ta_bcc',
    'facet': '110',
    'E_surf': [2.902, 3.150]},
'Mo110': {
    'solid': 'Mo_bcc',
    'facet': '110',
    'E_surf': [2.907, 3.000]},
'Ta110': {
    'solid': 'Ta_bcc',
    'facet': '110',
    'E_surf': [2.902, 3.150]},
'W110': {
    'solid': 'W_bcc',
    'facet': '110',
    'E_surf': [3.265, 3.675]},
# (0001) facets of hcp metals
'Mg0001': {
    'solid': 'Mg_hcp',
    'facet': '0001',
    'E_surf': [0.785, 0.760]},
'Zn0001': {
    'solid': 'Zn_hcp',
    'facet': '0001',
    'E_surf': [0.993, 0.990]},
'Cd0001': {
    'solid': 'Cd_hcp',
    'facet': '0001',
    'E_surf': [0.762, 0.740]},
'Sc0001': {
    'solid': 'Sc_hcp',
    'facet': '0001',
    'E_surf': [1.275]},
'Ti0001': {
    'solid': 'Ti_hcp',
    'facet': '0001',
    'E_surf': [1.989, 2.100]},
'Co0001': {
    'solid': 'Co_hcp',
    'facet': '0001',
    'E_surf': [2.522, 2.550]},
'Zr0001': {
    'solid': 'Zr_hcp',
    'facet': '0001',
    'E_surf': [1.909, 2.000]},
'Ru0001': {
    'solid': 'Ru_hcp',
    'facet': '0001',
    'E_surf': [3.043, 3.050]},
'La0001': {
    'solid': 'La_hcp',
    'facet': '0001',
    'E_surf': [1.020]},
'Os0001': {
    'solid': 'Os_hcp',
    'facet': '0001',
    'E_surf': [3.439, 3.450]},
# (111) facets of fcc metals
'Ca111': {
    'solid': 'Ca_fcc',
    'facet': '111',
    'E_surf': [0.502, 0.490]},
'Sr111': {
    'solid': 'Sr_fcc',
    'facet': '111',
    'E_surf': [0.419, 0.410]},
'Al111': {
    'solid': 'Al_fcc',
    'facet': '111',
    'E_surf': [1.143, 1.160]},
'Pb111': {
    'solid': 'Pb_fcc',
    'facet': '111',
    'E_surf': [0.593, 0.600]},
'Ni111': {
    'solid': 'Ni_fcc',
    'facet': '111',
    'E_surf': [2.380, 2.450]},
'Cu111': {
    'solid': 'Cu_fcc',
    'facet': '111',
    'E_surf': [1.790, 1.825]},
'Rh111': {
    'solid': 'Rh_fcc',
    'facet': '111',
    'E_surf': [2.659, 2.700]},
'Pb111': {
    'solid': 'Pb_fcc',
    'facet': '111',
    'E_surf': [0.593, 0.600]},
'Pd111': {
    'solid': 'Pd_fcc',
    'facet': '111',
    'E_surf': [2.003, 2.050]},
'Ag111': {
    'solid': 'Ag_fcc',
    'facet': '111',
    'E_surf': [1.246, 1.250]},
'Ir111': {
    'solid': 'Ir_fcc',
    'facet': '111',
    'E_surf': [3.048, 3.000]},
'Pt111': {
    'solid': 'Pt_fcc',
    'facet': '111',
    'E_surf': [2.489, 2.475]},
'Au111': {
    'solid': 'Au_fcc',
    'facet': '111',
    'E_surf': [1.506, 1.500]},
} # type: dict

#Surface energy functions
def in_surf_data(name:str)->None:
    if name not in surf_data:
        raise KeyError('System %s not in database.' % name)

def get_esurf_bulk_data(name:str)->tuple:
    in_surf_data(name)
    solid = surf_data[name]['solid']
    lp = get_solid_lattice_parameter(solid)
    facet = surf_data[name]['facet']
    struct = get_solid_crystal_structure(solid)
    symb = get_solid_symbols(solid)
    mm = get_solid_magmom(solid)
    if facet is '111':
        assert struct is 'fcc'
    elif facet is '110':
        assert struct is 'bcc'
    elif facet is '0001':
        assert struct is 'hcp'
    else:
        raise ValueError()
    return solid, symb, facet, struct, lp, mm

def get_surface_area(facet:str, lp:float)->float:
    if facet is '111':
        A = lp**2. * np.sqrt(3.)/4.
    elif facet is '110':
        A = lp**2. / np.sqrt(2.)
    elif facet is '0001':
        A = lp**2. * np.sqrt(3.)/2.
    return A

def convert_to_Jm2(A:float, e:float)->float:
    A *= 1.e-20 # m**2
    J = kJ/1000.
    e /= (J*A)
    return e

def get_exp_surface_energy(name:str)->tuple:
    in_surf_data(name)
    e = surf_data[name]['E_surf']
    return e, np.mean(e)

def getarg(sys:str, systems:str)->str:
    k = np.where(systems==sys)[0]
    if len(k) < 1:
        raise KeyError('%s not in keylist' % sys)
    elif len(k) > 1:
        raise KeyError('multiple %s in keylist' % sys)
    assert len(k) == 1
    return k[0]

def get_dft_latpar(solid:str, xc:str)->float:
    in_bulk_data(solid)
    from data_handler_045 import Data_Handler   # type: ignore
    data = Data_Handler()
    if solid in sol58_lp:
        db = 'solids_1_lp'
    elif solid in hcp_solids_10:
        db = 'hcp10_lp'
    else:
        raise KeyError('solid %s not found for xc=%s' % (solid,xc))
    a, b, c, d = data.get_accurate_ref_data(db, xc)
    arg = getarg(solid, b)
    return c[arg]
#
#
#
#===============================================================
# Solids_1 lattice constants, cohesive energies, and bulk moduli
#===============================================================
#
#
#
# References:
# 0) CRC, 64th edition, 1983-1984
# 1) Hao, Fang, Sun, Csonka, Philipsen, Perdew, PRB 85, 014111 (2012)
# 2) Schimka, Harl, Kresse, JCP 134, 024116 (2011)
# 3) Kittel, Introduction to Solid State Physics, 8th edition (2005)
# 4) CODATA, http://www.codata.org/resources/databases/key1.html (01/29/2013)
# 5a) Haglund: 3d metals, PRB (1991)
# 5b) Haglund: 4d metals, PRB (1992)
# 5c) Haglund: 5d metals, PRB (1993)
# 6) Sun, Marsman, Csonka, Ruzsinszky, Hao, Kim, Kresse, Perdew, PRB 84, 035117 (2011)
#
# In general: zero-point phonon corrected lattice constants,
#             cohesive energies, and bulk moduli.
#
bulk_data = {
# Cubic pure solids
'Li_bcc': {
    'name': 'Li_bcc',
    'lattice parameter': 3.451, #ref 1
    'pbe_lp': 3.429,            #ref 1 (dft-code=BAND)
    'cohesive energy': 37.70,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 344.,  #ref 3
    'bulk modulus': 13.9,       #ref 2 (phonon-corrected)
    'symbols': 'Li',
    'latex name': "Li",
    'structure':'bcc',
    'magmom': None},
'Na_bcc': {
    'name': 'Na_bcc',
    'lattice parameter': 4.207, #ref 1
    'pbe_lp': 4.197,            #ref 1 (dft-code=BAND)
    'cohesive energy': 25.76,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 158.,  #ref 3
    'bulk modulus': 7.7,        #ref 2 (phonon-corrected)
    'symbols': 'Na',
    'latex name': "Na",
    'structure':'bcc',
    'magmom': None},
'K_bcc': {
    'name': 'K_bcc',
    'lattice parameter': 5.211, #ref 1
    'pbe_lp': 5.281,            #ref 1 (dft-code=BAND)
    'cohesive energy': 21.48,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 91.,   #ref 3
    'bulk modulus': None,
    'symbols': 'K',
    'latex name': "K",
    'structure':'bcc',
    'magmom': None},
'Rb_bcc': {
    'name': 'Rb_bcc',
    'lattice parameter': 5.580, #ref 1
    'pbe_lp': 5.665,            #ref 1 (dft-code=BAND)
    'cohesive energy': 19.64,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 56.,   #ref 3
    'bulk modulus': None,
    'symbols': 'Rb',
    'latex name': "Rb",
    'structure':'bcc',
    'magmom': None},
'Ca_fcc': {
    'name': 'Ca_fcc',
    'lattice parameter': 5.555, #ref 1
    'pbe_lp': 5.521,            #ref 1 (dft-code=BAND)
    'cohesive energy': 42.39,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 230.,  #ref 3
    'bulk modulus': 18.7,       #ref 6 (phonon-corrected)
    'symbols': 'Ca',
    'latex name': "Ca",
    'structure':'fcc',
    'magmom': None},
'Sr_fcc': {
    'name': 'Sr_fcc',
    'lattice parameter': 6.042, #ref 1
    'pbe_lp': 6.013,            #ref 1 (dft-code=BAND)
    'cohesive energy': 39.3,    #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 147.,  #ref 3
    'bulk modulus': 12.5,       #ref 6 (phonon-corrected)
    'symbols': 'Sr',
    'latex name': "Sr",
    'structure':'fcc',
    'magmom': None},
'Ba_bcc': {
    'name': 'Ba_bcc',
    'lattice parameter': 5.004, #ref 1
    'pbe_lp': 5.022,            #ref 1 (dft-code=BAND)
    'cohesive energy': 43.2,    #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 110.,  #ref 3
    'bulk modulus': 9.4,        #ref 6 (phonon-corrected)
    'symbols': 'Ba',
    'latex name': "Ba",
    'structure':'bcc',
    'magmom': None},
'V_bcc': {
    'name': 'V_bcc',
    'lattice parameter': 3.024, #ref 1
    'pbe_lp': 2.997,            #ref 1 (dft-code=BAND)
    'cohesive energy': 122.12,  #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 380.,  #ref 3
    'bulk modulus': None,
    'symbols': 'V',
    'latex name': "V",
    'structure':'bcc',
    'magmom': None},
'Nb_bcc': {
    'name': 'Nb_bcc',
    'lattice parameter': 3.293, #ref 1
    'pbe_lp': 3.310,            #ref 1 (dft-code=BAND)
    'cohesive energy': 172.758, #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 275.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Nb',
    'latex name': "Nb",
    'structure':'bcc',
    'magmom': None},
'Ta_bcc': {
    'name': 'Ta_bcc',
    'lattice parameter': 3.299, #ref 1
    'pbe_lp': 3.347,            #ref 1 (dft-code=BAND)
    'cohesive energy': 186.765, #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 240.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Ta',
    'latex name': "Ta",
    'structure':'bcc',
    'magmom': None},
'Mo_bcc': {
    'name': 'Mo_bcc',
    'lattice parameter': 3.141, #ref 1
    'pbe_lp': 3.161,            #ref 1 (dft-code=BAND)
    'cohesive energy': 156.92,  #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 450.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Mo',
    'latex name': "Mo",
    'structure':'bcc',
    'magmom': None},
'W_bcc': {
    'name': 'W_bcc',
    'lattice parameter': 3.161, #ref 1
    'pbe_lp': 3.170,            #ref 1 (dft-code=BAND)
    'cohesive energy': 202.70,  #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 400.,  #ref 3
    'bulk modulus': None,
    'symbols': 'W',
    'latex name': "W",
    'structure':'bcc',
    'magmom': None},
'Fe_bcc': {
    'name': 'Fe_bcc',
    'lattice parameter': 2.855, #ref 1
    'pbe_lp': 2.834,            #ref 1 (dft-code=BAND)
    'cohesive energy': 98.94,   #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 470.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Fe',
    'latex name': "Fe",
    'structure':'bcc',
    'magmom': 2.22},
'Rh_fcc': {
    'name': 'Rh_fcc',
    'lattice parameter': 3.793, #ref 1
    'pbe_lp': 3.827,            #ref 1 (dft-code=BAND)
    'cohesive energy': 132.79,   #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 480.,  #ref 3
    'bulk modulus': 272.1,      #ref 2 (phonon-corrected)
    'symbols': 'Rh',
    'latex name': "Rh",
    'structure':'fcc',
    'magmom': None},
'Ir_fcc': {
    'name': 'Ir_fcc',
    'lattice parameter': 3.832, #ref 1
    'pbe_lp': 3.872,            #ref 1 (dft-code=BAND)
    'cohesive energy': 158.78,  #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 420.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Ir',
    'latex name': "Ir",
    'structure':'fcc',
    'magmom': None},
'Ni_fcc': {
    'name': 'Ni_fcc',
    'lattice parameter': 3.509, #ref 1
    'pbe_lp': 3.517,            #ref 1 (dft-code=BAND)
    'cohesive energy': 102.213, #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 450.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Ni',
    'latex name': "Ni",
    'structure':'fcc',
    'magmom': 0.64},
'Pd_fcc': {
    'name': 'Pd_fcc',
    'lattice parameter': 3.876, #ref 1
    'pbe_lp': 3.932,            #ref 1 (dft-code=BAND)
    'cohesive energy': 90.20,   #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 274.,  #ref 3
    'bulk modulus': 198.1,      #ref 2 (phonon-corrected)
    'symbols': 'Pd',
    'latex name': "Pd",
    'structure':'fcc',
    'magmom': None},
'Pt_fcc': {
    'name': 'Pt_fcc',
    'lattice parameter': 3.913, #ref 1
    'pbe_lp': 3.971,            #ref 1 (dft-code=BAND)
    'cohesive energy': 134.9,   #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 240.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Pt',
    'latex name': "Pt",
    'structure':'fcc',
    'magmom': None},
'Cu_fcc': {
    'name': 'Cu_fcc',
    'lattice parameter': 3.595, #ref 1
    'pbe_lp': 3.630,            #ref 1 (dft-code=BAND)
    'cohesive energy': 80.36,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 343.,  #ref 3
    'bulk modulus': 145.0,      #ref 2 (phonon-corrected)
    'symbols': 'Cu',
    'latex name': "Cu",
    'structure':'fcc',
    'magmom': None},
'Ag_fcc': {
    'name': 'Ag_fcc',
    'lattice parameter': 4.063, #ref 1
    'pbe_lp': 4.150,            #ref 1 (dft-code=BAND)
    'cohesive energy': 68.0,    #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 225.,  #ref 3
    'bulk modulus': 110.8,      #ref 2 (phonon-corrected)
    'symbols': 'Ag',
    'latex name': "Ag",
    'structure':'fcc',
    'magmom': None},
'Au_fcc': {
    'name': 'Au_fcc',
    'lattice parameter': 4.061, #ref 1
    'pbe_lp': 4.147,            #ref 1 (dft-code=BAND)
    'cohesive energy': 87.46,   #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 165.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Au',
    'latex name': "Au",
    'structure':'fcc',
    'magmom': None},
'Al_fcc': {
    'name': 'Al_fcc',
    'lattice parameter': 4.019, #ref 1
    'pbe_lp': 4.037,            #ref 1 (dft-code=BAND)
    'cohesive energy': 78.3,    #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 428.,  #ref 3
    'bulk modulus': 82.0,       #ref 2 (phonon-corrected)
    'symbols': 'Al',
    'latex name': "Al",
    'structure':'fcc',
    'magmom': None},
'Pb_fcc': {
    'name': 'Pb_fcc',
    'lattice parameter': 4.912, #ref 1
    'pbe_lp': 5.040,            #ref 1 (dft-code=BAND)
    'cohesive energy': 46.81,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 105.,  #ref 3
    'bulk modulus': None,
    'symbols': 'Pb',
    'latex name': "Pb",
    'structure':'fcc',
    'magmom': None},
'C_diamond': {
    'name': 'C_diamond',
    'lattice parameter': 3.555, #ref 1
    'pbe_lp': 3.571,            #ref 1 (dft-code=BAND)
    'cohesive energy': 170.0,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 2230., #ref 3
    'bulk modulus': 454.7,      #ref 2 (phonon-corrected)
    'symbols': 'C',
    'latex name': "C",
    'structure':'diamond',
    'magmom': None},
'Si_diamond': {
    'name': 'Si_diamond',
    'lattice parameter': 5.422, #ref 1
    'pbe_lp': 5.468,            #ref 1 (dft-code=BAND)
    'cohesive energy': 106.52,  #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 645.,  #ref 3
    'bulk modulus': 100.8,      #ref 2 (phonon-corrected)
    'symbols': 'Si',
    'latex name': "Si",
    'structure':'diamond',
    'magmom': None},
'Ge_diamond': {
    'name': 'Ge_diamond',
    'lattice parameter': 5.644, #ref 1
    'pbe_lp': 5.764,            #ref 1 (dft-code=BAND)
    'cohesive energy': 88.25,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 374.,  #ref 3
    'bulk modulus': 77.3,       #ref 2 (phonon-corrected)
    'symbols': 'Ge',
    'latex name': "Ge",
    'structure':'diamond',
    'magmom': None},
'Sn_diamond': {
    'name': 'Sn_diamond',
    'lattice parameter': 6.476, #ref 1
    'pbe_lp': 6.659,            #ref 1 (dft-code=BAND)
    'cohesive energy': 72.01,   #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 200.,  #ref 3
    'bulk modulus': 55.5,       #ref 2 (phonon-corrected)
    'symbols': 'Sn',
    'latex name': "Sn",
    'structure':'diamond',
    'magmom': None},
# Cubic compounds
'LiH_b1': {
    'name': 'LiH_b1',
    'lattice parameter': 3.979, #ref 2
    'pbe_lp': 4.006,            #ref 2 (dft-code=VASP)
    'cohesive energy': 240.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 40.1,       #ref 2 (phonon-corrected)
    'symbols': 'LiH',
    'latex name': "LiH",
    'structure':'rocksalt',
    'magmom': None},
'LiF_b1': {
    'name': 'LiF_b1',
    'lattice parameter': 3.974, #ref 1
    'pbe_lp': 4.064,            #ref 1 (dft-code=BAND)
    'cohesive energy': 430.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 76.3,       #ref 2 (phonon-corrected)
    'symbols': 'LiF',
    'latex name': "LiF",
    'structure':'rocksalt',
    'magmom': None},
'LiCl_b1': {
    'name': 'LiCl_b1',
    'lattice parameter': 5.072, #ref 1
    'pbe_lp': 5.147,            #ref 1 (dft-code=BAND)
    'cohesive energy': 346.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 38.7,       #ref 2 (phonon-corrected)
    'symbols': 'LiCl',
    'latex name': "LiCl",
    'structure':'rocksalt',
    'magmom': None},
'NaF_b1': {
    'name': 'NaF_b1',
    'lattice parameter': 4.570, #ref 1
    'pbe_lp': 4.700,            #ref 1 (dft-code=BAND)
    'cohesive energy': 383.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 53.1,       #ref 2 (phonon-corrected)
    'symbols': 'NaF',
    'latex name': "NaF",
    'structure':'rocksalt',
    'magmom': None},
'NaCl_b1': {
    'name': 'NaCl_b1',
    'lattice parameter': 5.565, #ref 1
    'pbe_lp': 5.695,            #ref 1 (dft-code=BAND)
    'cohesive energy': 322.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 27.6,       #ref 2 (phonon-corrected)
    'symbols': 'NaCl',
    'latex name': "NaCl",
    'structure':'rocksalt',
    'magmom': None},
'MgO_b1': {
    'name': 'MgO_b1',
    'lattice parameter': 4.188, #ref 1
    'pbe_lp': 4.255,            #ref 1 (dft-code=BAND)
    'cohesive energy': 502.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 169.8,       #ref 2 (phonon-corrected)
    'symbols': 'MgO',
    'latex name': "MgO",
    'structure':'rocksalt',
    'magmom': None},
'MgS_b1': {
    'name': 'MgS_b1',
    'lattice parameter': 5.188, #ref 1
    'pbe_lp': 5.228,            #ref 1 (dft-code=BAND)
    'cohesive energy': 91.52,   #ref 0
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 650.,  #ref 1
    'bulk modulus': None,
    'symbols': 'MgS',
    'latex name': "MgS",
    'structure':'rocksalt',
    'magmom': None},
'CaO_b1': {
    'name': 'CaO_b1',
    'lattice parameter': 4.781, #ref 1
    'pbe_lp': 4.832,            #ref 1 (dft-code=BAND)
    'cohesive energy': 126.16,  #ref 4
    'Ecoh unit': 'kcal/mol',
    'Ecoh is phonon-corrected': False,
    'debye temperature': 648.,  #ref 1
    'bulk modulus': None,
    'symbols': 'CaO',
    'latex name': "CaO",
    'structure':'rocksalt',
    'magmom': None},
'TiC_b1': {
    'name': 'TiC_b1',
    'lattice parameter': 4.318, #ref 1
    'pbe_lp': 4.332,            #ref 1 (dft-code=BAND)
    'cohesive energy': 690.5,   #ref 5a
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'TiC',
    'latex name': "TiC",
    'structure':'rocksalt',
    'magmom': None},
'TiN_b1': {
    'name': 'TiN_b1',
    'lattice parameter': 4.226, #ref 1
    'pbe_lp': 4.247,            #ref 1 (dft-code=BAND)
    'cohesive energy': 645.9,   #ref 5a
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'TiN',
    'latex name': "TiN",
    'structure':'rocksalt',
    'magmom': None},
'ZrC_b1': {
    'name': 'ZrC_b1',
    'lattice parameter': 4.687, #ref 1
    'pbe_lp': 4.708,            #ref 1 (dft-code=BAND)
    'cohesive energy': 765.3,   #ref 5b
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'ZrC',
    'latex name': "ZrC",
    'structure':'rocksalt',
    'magmom': None},
'ZrN_b1': {
    'name': 'ZrN_b1',
    'lattice parameter': 4.576, #ref 1
    'pbe_lp': 4.594,            #ref 1 (dft-code=BAND)
    'cohesive energy': 726.0,   #ref 5b
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'ZrN',
    'latex name': "ZrN",
    'structure':'rocksalt',
    'magmom': None},
'VC_b1': {
    'name': 'VC_b1',
    'lattice parameter': 4.148, #ref 1
    'pbe_lp': 4.154,            #ref 1 (dft-code=BAND)
    'cohesive energy': 669.5,   #ref 5a
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'VC',
    'latex name': "VC",
    'structure':'rocksalt',
    'magmom': None},
'VN_b1': {
    'name': 'VN_b1',
    'lattice parameter': 4.130, #ref 1
    'pbe_lp': 4.116,            #ref 1 (dft-code=BAND)
    'cohesive energy': 602.6,   #ref 5a
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'VN',
    'latex name': "VN",
    'structure':'rocksalt',
    'magmom': None},
'NbC_b1': {
    'name': 'NbC_b1',
    'lattice parameter': 4.461, #ref 1
    'pbe_lp': 4.484,            #ref 1 (dft-code=BAND)
    'cohesive energy': 796.8,   #ref 5b
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'NbC',
    'latex name': "NbC",
    'structure':'rocksalt',
    'magmom': None},
'NbN_b1': {
    'name': 'NbN_b1',
    'lattice parameter': 4.383, #ref 1
    'pbe_lp': 4.422,            #ref 1 (dft-code=BAND)
    'cohesive energy': 723.3,   #ref 5a
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'NbN',
    'latex name': "NbN",
    'structure':'rocksalt',
    'magmom': None},
'FeAl_b2': {
    'name': 'FeAl_b2',
    'lattice parameter': 2.881, #ref 1
    'pbe_lp': 2.868,            #ref 1 (dft-code=BAND)
    'cohesive energy': None,
    'bulk modulus': None,
    'symbols': 'FeAl',
    'latex name': "FeAl",
    'structure':'cesiumchloride',
    'magmom': 0.35},
'CoAl_b2': {
    'name': 'CoAl_b2',
    'lattice parameter': 2.854, #ref 1
    'pbe_lp': 2.851,            #ref 1 (dft-code=BAND)
    'cohesive energy': None,
    'bulk modulus': None,
    'symbols': 'CoAl',
    'latex name': "CoAl",
    'structure':'cesiumchloride',
    'magmom': None},
'NiAl_b2': {
    'name': 'NiAl_b2',
    'lattice parameter': 2.881, #ref 1
    'pbe_lp': 2.892,            #ref 1 (dft-code=BAND)
    'cohesive energy': None,
    'bulk modulus': None,
    'symbols': 'NiAl',
    'latex name': "NiAl",
    'structure':'cesiumchloride',
    'magmom': None},
'BN_b3': {
    'name': 'BN_b3',
    'lattice parameter': 3.594, #ref 1
    'pbe_lp': 3.624,            #ref 1 (dft-code=BAND)
    'cohesive energy': 652.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 410.2,      #ref 2 (phonon-corrected)
    'symbols': 'BN',
    'latex name': "BN",
    'structure':'zincblende',
    'magmom': None},
'BP_b3': {
    'name': 'BP_b3',
    'lattice parameter': 4.527, #ref 1
    'pbe_lp': 4.548,            #ref 1 (dft-code=BAND)
    'cohesive energy': 496.,    #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 168.0,      #ref 2 (phonon-corrected)
    'symbols': 'BP',
    'latex name': "BP",
    'structure':'zincblende',
    'magmom': None},
'BAs_b3': {
    'name': 'BAs_b3',
    'lattice parameter': 4.764, #ref 1
    'pbe_lp': 4.809,            #ref 1 (dft-code=BAND)
    'cohesive energy': None,
    'bulk modulus': 151.1,      #ref 2 (phonon-corrected)
    'symbols': 'BAs',
    'latex name': "BAs",
    'structure':'zincblende',
    'magmom': None},
'AlN_b3': {
    'name': 'AlN_b3',
    'lattice parameter': 4.368, #ref 2
    'pbe_lp': 4.402,            #ref 2 (dft-code=VASP)
    'cohesive energy': 564.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 206.0,      #ref 2 (phonon-corrected)
    'symbols': 'AlN',
    'latex name': "AlN",
    'structure':'zincblende',
    'magmom': None},
'AlP_b3': {
    'name': 'AlP_b3',
    'lattice parameter': 5.450, #ref 1
    'pbe_lp': 5.504,            #ref 1 (dft-code=BAND)
    'cohesive energy': 417.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 87.4,       #ref 2 (phonon-corrected)
    'symbols': 'AlP',
    'latex name': "AlP",
    'structure':'zincblende',
    'magmom': None},
'AlAs_b3': {
    'name': 'AlAs_b3',
    'lattice parameter': 5.649, #ref 1
    'pbe_lp': 5.728,            #ref 1 (dft-code=BAND)
    'cohesive energy': 369.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 75.0,       #ref 2 (phonon-corrected)
    'symbols': 'AlAs',
    'latex name': "AlAs",
    'structure':'zincblende',
    'magmom': None},
'GaN_b3': {
    'name': 'GaN_b3',
    'lattice parameter': 4.523, #ref 1
    'pbe_lp': 4.549,            #ref 1 (dft-code=BAND)
    'cohesive energy': 439.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 213.7,       #ref 2 (phonon-corrected)
    'symbols': 'GaN',
    'latex name': "GaN",
    'structure':'zincblende',
    'magmom': None},
'GaP_b3': {
    'name': 'GaP_b3',
    'lattice parameter': 5.441, #ref 1
    'pbe_lp': 5.506,            #ref 1 (dft-code=BAND)
    'cohesive energy': 348.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 89.6,       #ref 2 (phonon-corrected)
    'symbols': 'GaP',
    'latex name': "GaP",
    'structure':'zincblende',
    'magmom': None},
'GaAs_b3': {
    'name': 'GaAs_b3',
    'lattice parameter': 5.641, #ref 1
    'pbe_lp': 5.751,            #ref 1 (dft-code=BAND)
    'cohesive energy': 322.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 76.7,       #ref 2 (phonon-corrected)
    'symbols': 'GaAs',
    'latex name': "GaAs",
    'structure':'zincblende',
    'magmom': None},
'InP_b3': {
    'name': 'InP_b3',
    'lattice parameter': 5.858, #ref 1
    'pbe_lp': 5.963,            #ref 1 (dft-code=BAND)
    'cohesive energy': 335.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 72.0,       #ref 2 (phonon-corrected)
    'symbols': 'InP',
    'latex name': "InP",
    'structure':'zincblende',
    'magmom': None},
'InAs_b3': {
    'name': 'InAs_b3',
    'lattice parameter': 6.048, #ref 1
    'pbe_lp': 6.188,            #ref 1 (dft-code=BAND)
    'cohesive energy': 297.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 58.6,       #ref 2 (phonon-corrected)
    'symbols': 'InAs',
    'latex name': "InAs",
    'structure':'zincblende',
    'magmom': None},
'SiC_b3': {
    'name': 'SiC_b3',
    'lattice parameter': 4.348, #ref 1
    'pbe_lp': 4.378,            #ref 1 (dft-code=BAND)
    'cohesive energy': 625.0,   #ref 2
    'Ecoh unit': 'kJ/mol',
    'Ecoh is phonon-corrected': True,
    'debye temperature': None,
    'bulk modulus': 229.1,      #ref 2 (phonon-corrected)
    'symbols': 'SiC',
    'latex name': "SiC",
    'structure':'zincblende',
    'magmom': None},
# hcp pure solids
'Cd_hcp': {
    'name': "Cd_hcp",
    'lattice parameter': 2.979,
    'lattice parameter 2': 5.620,
    'cohesive energy': 1.16,
    'bulk modulus': 46.7, # GPa
    'symbols': 'Cd',
    'latex name': "Cd",
    'structure':'hcp',
    'magmom': None},
'Co_hcp': {
    'name': "Co_hcp",
    'lattice parameter': 2.507,
    'lattice parameter 2': 4.069,
    'cohesive energy': 4.39,
    'bulk modulus': 191.4, # GPa
    'symbols': 'Co',
    'latex name': "Co",
    'structure':'hcp',
    'magmom': 1.72},
'Os_hcp': {
    'name': "Os_hcp",
    'lattice parameter': 2.734,
    'lattice parameter 2': 4.392,
    'cohesive energy': 8.17,
    'bulk modulus': 418., # GPa
    'symbols': 'Os',
    'latex name': "Os",
    'structure':'hcp',
    'magmom': None},
'Ru_hcp': {
    'name': "Ru_hcp",
    'lattice parameter': 2.706,
    'lattice parameter 2': 4.282,
    'cohesive energy': 6.74,
    'bulk modulus': 320.8, # GPa
    'symbols': 'Ru',
    'latex name': "Ru",
    'structure':'hcp',
    'magmom': None},
'Zn_hcp': {
    'name': "Zn_hcp",
    'lattice parameter': 2.665,
    'lattice parameter 2': 4.947,
    'cohesive energy': 1.35,
    'bulk modulus': 59.8, # GPa
    'symbols': 'Zn',
    'latex name': "Zn",
    'structure':'hcp',
    'magmom': None},
'Ti_hcp': {
    'name': "Ti_hcp",
    'lattice parameter': 2.951,
    'lattice parameter 2': 4.684,
    'cohesive energy': 4.85,
    'bulk modulus': 105.1, # GPa
    'symbols': 'Ti',
    'latex name': "Ti",
    'structure':'hcp',
    'magmom': None},
'Zr_hcp': {
    'name': "Zr_hcp",
    'lattice parameter': 3.232,
    'lattice parameter 2': 5.148,
    'cohesive energy': 6.25,
    'bulk modulus': 83.3, # GPa
    'symbols': 'Zr',
    'latex name': "Zr",
    'structure':'hcp',
    'magmom': None},
'Sc_hcp': {
    'name': "Sc_hcp",
    'lattice parameter': 3.309,
    'lattice parameter 2': 5.268,
    'cohesive energy': 3.90,
    'bulk modulus': 43.5, # GPa
    'symbols': 'Sc',
    'latex name': "Sc",
    'structure':'hcp',
    'magmom': None},
'Be_hcp': {
    'name': "Be_hcp",
    'lattice parameter': 2.286,
    'lattice parameter 2': 3.585,
    'cohesive energy': 3.32,
    'bulk modulus': 100.3, # GPa
    'symbols': 'Be',
    'latex name': "Be",
    'structure':'hcp',
    'magmom': None},
'Mg_hcp': {
    'name': "Mg_hcp",
    'lattice parameter': 3.209,
    'lattice parameter 2': 5.211,
    'cohesive energy': 1.51,
    'bulk modulus': 35.4, # GPa
    'symbols': 'Mg',
    'latex name': "Mg",
    'structure':'hcp',
    'magmom': None},
} # type: dict
# bulk solid functions
def in_bulk_data(name:str)->None:
    if name not in bulk_data:
        raise KeyError('System %s not in database.' % name)

def get_solid_lattice_parameter(name:str)->float:
    in_bulk_data(name)
    return bulk_data[name]['lattice parameter']

def get_hcp_covera(name:str)->float:
    in_bulk_data(name)
    if name in hcp_solids_10:
        a = bulk_data[name]['lattice parameter']
        c = bulk_data[name]['lattice parameter 2']
        return c/a
    else:
        return 0.

def get_solid_pbe_lp(name:str)->O[float]:
    in_bulk_data(name)
    if name in sol58_lp:
        return bulk_data[name]['pbe_lp']
    else:
        return None

def get_solid_cohesive_energy(name:str)->float:
    #in_bulk_data(name)
    #assert name in sol58_lp
    d = bulk_data[name]
    e = d['cohesive energy']
    if e == None:
        return e
    else:
        # units
        unit = d['Ecoh unit']
        if unit == 'kJ/mol':
            e *= kJ/mol
        elif unit == 'kcal/mol':
            e *= kcal/mol
        else:
            raise ValueError()
        corrected = d['Ecoh is phonon-corrected']
        if corrected == False:
            theta = d['debye temperature']
            assert theta != None
            e += (9./8.)*kB*theta
        return e

def get_solid_bulk_modulus(name:str)->O[float]:
    in_bulk_data(name)
    if name in bm32:
        return bulk_data[name]['bulk modulus']
    else:
        return None

def get_solid_magmom(name:str)->float:
    in_bulk_data(name)
    return bulk_data[name]['magmom']

def get_solid_crystal_structure(name:str)->str:
    in_bulk_data(name)
    struct = bulk_data[name]['structure']
    assert struct in ['fcc', 'bcc', 'diamond', 'hcp',
                      'rocksalt', 'cesiumchloride', 'zincblende']
    return struct

def get_solid_symbols(name:str)->str:
    in_bulk_data(name)
    return bulk_data[name]['symbols']

def get_solids_latex_name(name:str)->str:
    in_bulk_data(name)
    d = bulk_data[name]
    struct = d['structure']
    if struct == 'diamond':
        struct = 'dia'
    return d['latex name'] + ' ('+struct+')'

def get_solids_common_name(name:str)->str:
    return get_solids_latex_name(name)

def setup_bulk(solid:str, offset:float=0., set_lp:bool=False)->Atoms:
    from ase.build import bulk   # type: ignore
    in_bulk_data(solid)
    symbol = get_solid_symbols(solid)
    s = get_solid_crystal_structure(solid)
    if set_lp:
        lp = offset
    else:
        lp = get_solid_lattice_parameter(solid) + offset
    m = get_solid_magmom(solid)
    if s == 'hcp':
        cov = get_hcp_covera(solid)
        atoms = bulk(symbol, s, a=lp, covera=cov)
    else:
        atoms = bulk(symbol, s, a=lp,orthorhombic = True )
    atoms.set_pbc(True)
    if m != None:
        mm = np.zeros(len(atoms))
        mm[:] = m
        atoms.set_initial_magnetic_moments(mm)
    return atoms

def write_bulks()->None:
    for k,v in bulk_data.items():
        a = setup_bulk(k)
        write('structures/'+k+'.traj',a)

def setup_atom(symbol:str, vac:float=3.0, t:float=0.0, odd:float=0.0)->Atoms:
    from ase import Atoms   # type: ignore
    from ase.data import chemical_symbols, atomic_numbers, ground_state_magnetic_moments   # type: ignore
    assert symbol in chemical_symbols
    Z = atomic_numbers[symbol]
    m = ground_state_magnetic_moments[Z]
    atoms = Atoms(symbol)
    atoms.center(vacuum=vac)
    atoms.translate([t,t,t])
    atoms.cell[1][1] += odd
    atoms.cell[2][2] += 2*odd
    atoms.set_initial_magnetic_moments([m])
    return atoms

######################
# atomic constituents
######################
atoms_ = np.array([])
for aa in sol53_coh:
    symbs = get_solid_symbols(aa)
    atoms_ = np.append(atoms_, string2symbols(symbs))
sol53_coh_atoms = np.unique(atoms_)
atoms_ = np.array([])
for aa in sol54_coh:
    symbs = get_solid_symbols(aa)
    atoms_ = np.append(atoms_, string2symbols(symbs))
sol54_coh_atoms = np.unique(atoms_)
atoms_ = np.array([])
for aa in sol58_lp:
    symbs = get_solid_symbols(aa)
    atoms_ = np.append(atoms_, string2symbols(symbs))
sol58_atoms = np.unique(atoms_)
del atoms_, aa, symbs

######################
# complete set of solids
######################
solids_tot = np.append(sol58_lp, sol54_coh)
solids_tot = np.append(solids_tot, bm32)
solids_tot = np.unique(solids_tot)
np.unique(sol58_lp)
assert len(solids_tot) == len(np.unique(sol58_lp)) == 58
solids_tot_atoms = np.unique(sol54_coh_atoms)
assert len(solids_tot_atoms) == len(sol54_coh_atoms)
