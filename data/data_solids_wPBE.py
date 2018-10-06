from ase.atoms import Atoms, string2symbols   # type: ignore
from ase.units import kB, kJ   # type: ignore
import numpy as np   # type: ignore

qubic_solids_26 = ['Li_bcc','Na_bcc','K_bcc','Rb_bcc',
                   'Ca_fcc','Sr_fcc','Ba_bcc',
                   'Nb_bcc','Ta_bcc',
                   'Mo_bcc','W_bcc','Fe_bcc',
                   'Rh_fcc','Ir_fcc',
                   'Ni_fcc','Pd_fcc','Pt_fcc',
                   'Cu_fcc','Ag_fcc','Au_fcc',
                   'Al_fcc','Pb_fcc',              #KRIS: WHY WERE Pb_fcc and Au_fcc commented out???
                   'C_diamond','Si_diamond','Ge_diamond','Sn_diamond']
qubic_solids_27 = np.append(qubic_solids_26, 'V_bcc')

hcp_solids_10 = ['Cd_hcp', 'Co_hcp', 'Os_hcp', 'Ru_hcp', 'Zn_hcp',
                 'Zr_hcp', 'Sc_hcp', 'Be_hcp', 'Mg_hcp', 'Ti_hcp']


qubic_compounds_28 = ['LiH_b1', 'LiF_b1', 'LiCl_b1', 'NaF_b1', 'NaCl_b1',
                      'MgO_b1', 'MgS_b1', 'CaO_b1', 'TiC_b1', 'TiN_b1',
                      'ZrC_b1', 'ZrN_b1', 'NbC_b1', 'NbN_b1',
                      'FeAl_b2', 'CoAl_b2', 'NiAl_b2', 'BN_b3', 'BP_b3',
                      'AlN_b3', 'AlP_b3', 'AlAs_b3', 'GaN_b3',
                      'GaP_b3', 'GaAs_b3', 'InP_b3', 'InAs_b3', 'SiC_b3']
qubic_compounds_30 = np.append(qubic_compounds_28, ['VC_b1', 'VN_b1'])

qubic_compounds_31 = np.append(qubic_compounds_30, ['BAs_b3'])

pearson_compounds_15 = ['KBr_b1', 'CaSe_b1', 'SeAs_b1', 'LiI_b1',
                        'CsF_b1', 'CsI_b2', 'AgF_b1', 'AgCl_b1', 'AgBr_b1',
                        'CaS_b1', 'BaO_b1', 'BaSe_b1', 'CdO_b1', 'MnO_b1',
                        'MnS_b1'] # 'RbI_b1',

haglund_compounds_31 = ['ScC_b1', 'MnC_b1', 'FeC_b1', 'CoC_b1', 'NiC_b1',
                        'ScN_b1', 'CrN_b1', 'MnN_b1', 'CoN_b1', 'NiN_b1',
                        'MoC_b1', 'CrC_b1', 'RuC_b1', 'RhC_b1', 'PdC_b1',
                        'MoN_b1', 'RuN_b1', 'RhN_b1', 'PdN_b1',
                        'LaC_b1', 'TaC_b1', 'WC_b1', 'OsC_b1', 'IrC_b1',  'PtC_b1',
                        'LaN_b1', 'TaN_b1', 'WN_b1', 'OsN_b1', 'IrN_b1',  'PtN_b1'] #'FeN_b1',
qubic_compounds_46 = np.append(pearson_compounds_15, haglund_compounds_31)

solids_1 = np.append(qubic_solids_26, qubic_compounds_28)
solids_2 = np.append(hcp_solids_10, qubic_compounds_46)
solids_0 = np.append(solids_1, solids_2)

sl20 = ['Li_bcc', 'Na_bcc', 'Ca_fcc', 'Sr_fcc', 'Ba_bcc',
        'Al_fcc', 'Cu_fcc', 'Rh_fcc', 'Pd_fcc', 'Ag_fcc',
        'C_diamond', 'Si_diamond', 'Ge_diamond',
        'SiC_b3', 'GaAs_b3',
        'LiF_b1', 'LiCl_b1', 'NaF_b1', 'NaCl_b1', 'MgO_b1']
assert len(sl20) == 20

# ZPAE (semi-empirical) corrected lattice parameters from Haas, Tran, Blaha, PRB 79, 085104 (2009)
# and ZPAE (phonon) corrected lattice parameters + Debye temps + bulk moduli (GPa) from Hao et al., PRB 85, 014111 (2012)
# cohesive energies (at 0 K) and Debye temperatures from Kittel and Hao et al.
data = {
# Cubic pure solids
'Li_bcc': {
    'name': 'Li_bcc',
    'lattice parameter': 3.451,
    'pbe_lp': 3.429,
    'cohesive energy': 1.63,
    'debye temperature': 344.,
    'bulk modulus': 13.,
    'symbols': 'Li',
    'latex name': "Li",
    'structure':'bcc',
    'magmom': None},
'Na_bcc': {
    'name': 'Na_bcc',
    'lattice parameter': 4.207,
    'pbe_lp': 4.197,
    'cohesive energy': 1.113,
    'bulk modulus': 7.5,
    'debye temperature': 158.,
    'symbols': 'Na',
    'latex name': "Na",
    'structure':'bcc',
    'magmom': None},
'K_bcc': {
    'name': 'K_bcc',
    'lattice parameter': 5.211,
    'pbe_lp': 5.281,
    'cohesive energy': 0.934,
    'bulk modulus': 3.7,
    'debye temperature': 91.,
    'symbols': 'K',
    'latex name': "K",
    'structure':'bcc',
    'magmom': None},
'Rb_bcc': {
    'name': 'Rb_bcc',
    'lattice parameter': 5.580,
    'pbe_lp': 5.665,
    'cohesive energy': 0.852,
    'bulk modulus': 3.1,
    'debye temperature': 56.,
    'symbols': 'Rb',
    'latex name': "Rb",
    'structure':'bcc',
    'magmom': None},
'Ca_fcc': {
    'name': 'Ca_fcc',
    'lattice parameter': 5.555,
    'pbe_lp': 5.521,
    'cohesive energy': 1.84,
    'bulk modulus': 18.36,
    'debye temperature': 230.,
    'symbols': 'Ca',
    'latex name': "Ca",
    'structure':'fcc',
    'magmom': None},
'Sr_fcc': {
    'name': 'Sr_fcc',
    'lattice parameter': 6.042,
    'pbe_lp': 6.013,
    'cohesive energy': 1.72,
    'bulk modulus': 12.36,
    'debye temperature': 147.,
    'symbols': 'Sr',
    'latex name': "Sr",
    'structure':'fcc',
    'magmom': None},
'Ba_bcc': {
    'name': 'Ba_bcc',
    'lattice parameter': 5.004,
    'pbe_lp': 5.022,
    'cohesive energy': 1.90,
    'bulk modulus': 9.3,
    'debye temperature': 110.,
    'symbols': 'Ba',
    'latex name': "Ba",
    'structure':'bcc',
    'magmom': None},
'V_bcc': {
    'name': 'V_bcc',
    'lattice parameter': 3.024,
    'pbe_lp': 2.997,
    'cohesive energy': 5.31,
    'bulk modulus': 161.9,
    'debye temperature': 380.,
    'symbols': 'V',
    'latex name': "V",
    'structure':'bcc',
    'magmom': None},
'Nb_bcc': {
    'name': 'Nb_bcc',
    'lattice parameter': 3.293,
    'pbe_lp': 3.310,
    'cohesive energy': 7.57,
    'bulk modulus': 170.2,
    'debye temperature': 275.,
    'symbols': 'Nb',
    'latex name': "Nb",
    'structure':'bcc',
    'magmom': None},
'Ta_bcc': {
    'name': 'Ta_bcc',
    'lattice parameter': 3.299,
    'pbe_lp': 3.347,
    'cohesive energy': 8.10,
    'bulk modulus': 194.,
    'debye temperature': 240.,
    'symbols': 'Ta',
    'latex name': "Ta",
    'structure':'bcc',
    'magmom': None},
'Mo_bcc': {
    'name': 'Mo_bcc',
    'lattice parameter': 3.141,
    'pbe_lp': 3.161,
    'cohesive energy': 6.82,
    'bulk modulus': 272.5,
    'debye temperature': 450.,
    'symbols': 'Mo',
    'latex name': "Mo",
    'structure':'bcc',
    'magmom': None},
'W_bcc': {
    'name': 'W_bcc',
    'lattice parameter': 3.161,
    'pbe_lp': 3.170,
    'cohesive energy': 8.90,
    'bulk modulus': 314.,
    'debye temperature': 400.,
    'symbols': 'W',
    'latex name': "W",
    'structure':'bcc',
    'magmom': None},
'Fe_bcc': {
    'name': 'Fe_bcc',
    'lattice parameter': 2.855,
    'pbe_lp': 2.834,
    'cohesive energy': 4.28,
    'bulk modulus': 168.3,
    'debye temperature': 470.,
    'symbols': 'Fe',
    'latex name': "Fe",
    'structure':'bcc',
    'magmom': 2.22},
'Rh_fcc': {
    'name': 'Rh_fcc',
    'lattice parameter': 3.793,
    'pbe_lp': 3.827,
    'cohesive energy': 5.75,
    'bulk modulus': 269.,
    'debye temperature': 480.,
    'symbols': 'Rh',
    'latex name': "Rh",
    'structure':'fcc',
    'magmom': None},
'Ir_fcc': {
    'name': 'Ir_fcc',
    'lattice parameter': 3.832,
    'pbe_lp': 3.872,
    'cohesive energy': 6.94,
    'bulk modulus': 355.,
    'debye temperature': 420.,
    'symbols': 'Ir',
    'latex name': "Ir",
    'structure':'fcc',
    'magmom': None},
'Ni_fcc': {
    'name': 'Ni_fcc',
    'lattice parameter': 3.509,
    'pbe_lp': 3.517,
    'cohesive energy': 4.44,
    'bulk modulus': 188.,
    'debye temperature': 450.,
    'symbols': 'Ni',
    'latex name': "Ni",
    'structure':'fcc',
    'magmom': 0.64},
'Pd_fcc': {
    'name': 'Pd_fcc',
    'lattice parameter': 3.876,
    'pbe_lp': 3.932,
    'cohesive energy': 3.89,
    'bulk modulus': 195.,
    'debye temperature': 274.,
    'symbols': 'Pd',
    'latex name': "Pd",
    'structure':'fcc',
    'magmom': None},
'Pt_fcc': {
    'name': 'Pt_fcc',
    'lattice parameter': 3.913,
    'pbe_lp': 3.971,
    'cohesive energy': 5.84,
    'bulk modulus': 278.3,
    'debye temperature': 240.,
    'symbols': 'Pt',
    'latex name': "Pt",
    'structure':'fcc',
    'magmom': None},
'Cu_fcc': {
    'name': 'Cu_fcc',
    'lattice parameter': 3.595,
    'pbe_lp': 3.630,
    'cohesive energy': 3.49,
    'bulk modulus': 142.,
    'debye temperature': 343.,
    'symbols': 'Cu',
    'latex name': "Cu",
    'structure':'fcc',
    'magmom': None},
'Ag_fcc': {
    'name': 'Ag_fcc',
    'lattice parameter': 4.063,
    'pbe_lp': 4.150,
    'cohesive energy': 2.95,
    'bulk modulus': 109.,
    'debye temperature': 225.,
    'symbols': 'Ag',
    'latex name': "Ag",
    'structure':'fcc',
    'magmom': None},
'Au_fcc': {
    'name': 'Au_fcc',
    'lattice parameter': 4.061,
    'pbe_lp': 4.147,
    'cohesive energy': 3.81,
    'bulk modulus': 180.,
    'debye temperature': 165.,
    'symbols': 'Au',
    'latex name': "Au",
    'structure':'fcc',
    'magmom': None},
'Al_fcc': {
    'name': 'Al_fcc',
    'lattice parameter': 4.019,
    'pbe_lp': 4.037,
    'cohesive energy': 3.39,
    'bulk modulus': 79.4,
    'debye temperature': 428.,
    'symbols': 'Al',
    'latex name': "Al",
    'structure':'fcc',
    'magmom': None},
'C_diamond': {
    'name': 'C_diamond',
    'lattice parameter': 3.555,
    'pbe_lp': 3.571,
    'cohesive energy': 7.37,
    'bulk modulus': 443.,
    'debye temperature': 2230.,
    'symbols': 'C',
    'latex name': "C (diamond)",
    'structure':'diamond',
    'magmom': None},
'Si_diamond': {
    'name': 'Si_diamond',
    'lattice parameter': 5.422,
    'pbe_lp': 5.468,
    'cohesive energy': 4.63,
    'bulk modulus': 99.2,
    'debye temperature': 645.,
    'symbols': 'Si',
    'latex name': "Si (diamond)",
    'structure':'diamond',
    'magmom': None},
'Ge_diamond': {
    'name': 'Ge_diamond',
    'lattice parameter': 5.644,
    'pbe_lp': 5.764,
    'cohesive energy': 3.85,
    'bulk modulus': 75.8,
    'debye temperature': 374.,
    'latex name': "Ge",
    'symbols': 'Ge',
    'structure':'diamond',
    'magmom': None},
'Sn_diamond': {
    'name': 'Sn_fcc',
    'lattice parameter': 6.476,
    'pbe_lp': 6.659,
    'cohesive energy': 3.14,
    'bulk modulus': 53.1,
    'debye temperature': 200.,
    'symbols': 'Sn',
    'latex name': "Sn",
    'structure':'diamond',
    'magmom': None},
'Pb_fcc': {
    'name': 'Pb_fcc',
    'lattice parameter': 4.912,
    'pbe_lp': 5.040,
    'cohesive energy': 2.03,
    'bulk modulus': 48.8,
    'debye temperature': 105.,
    'symbols': 'Pb',
    'latex name': "Pb",
    'structure':'fcc',
    'magmom': None},
# hcp pure solids
'Cd_hcp': {
    'name': "Cd_hcp",
    'lattice parameter': 2.979,
    'lattice parameter 2': 5.620,
    'revtpss_lp': 2.95057912237,
    'cohesive energy': 1.16,
    'bulk modulus kittel': 46.7, # GPa
    'symbols': 'Cd',
    'latex name': "Cd",
    'structure':'hcp',
    'magmom': None},
'Co_hcp': {
    'name': "Co_hcp",
    'lattice parameter': 2.507,
    'lattice parameter 2': 4.069,
    'revtpss_lp':  2.45918897125,
    'cohesive energy': 4.39,
    'bulk modulus kittel': 191.4, # GPa
    'symbols': 'Co',
    'latex name': "Co",
    'structure':'hcp',
    'magmom': 1.72},
    #'gpaw_conv': [1.e-5, 1.e-6, 1.e-4]},
'Os_hcp': {
    'name': "Os_hcp",
    'lattice parameter': 2.734,
    'lattice parameter 2': 4.392,
    'revtpss_lp': 2.73044397786,
    'cohesive energy': 8.17,
    'bulk modulus kittel': 418., # GPa
    'symbols': 'Os',
    'latex name': "Os",
    'structure':'hcp',
    'magmom': None},
'Ru_hcp': {
    'name': "Ru_hcp",
    'lattice parameter': 2.706,
    'lattice parameter 2': 4.282,
    'revtpss_lp':  2.7191530774,
    'cohesive energy': 6.74,
    'bulk modulus kittel': 320.8, # GPa
    'symbols': 'Ru',
    'latex name': "Ru",
    'structure':'hcp',
    'magmom': None},
'Zn_hcp': {
    'name': "Zn_hcp",
    'lattice parameter': 2.665,
    'lattice parameter 2': 4.947,
    'revtpss_lp': 2.59872461224,
    'cohesive energy': 1.35,
    'bulk modulus kittel': 59.8, # GPa
    'symbols': 'Zn',
    'latex name': "Zn",
    'structure':'hcp',
    'magmom': None},
'Ti_hcp': {
    'name': "Ti_hcp",
    'lattice parameter': 2.951,
    'lattice parameter 2': 4.684,
    'revtpss_lp': 2.92208023815,
    'cohesive energy': 4.85,
    'bulk modulus kittel': 105.1, # GPa
    'symbols': 'Ti',
    'latex name': "Ti",
    'structure':'hcp',
    'magmom': None,
    'gpaw_conv': [1.e-6, 1.e-4, 1.e-6]},
'Zr_hcp': {
    'name': "Zr_hcp",
    'lattice parameter': 3.232,
    'lattice parameter 2': 5.148,
    'revtpss_lp': 3.23504460776,
    'cohesive energy': 6.25,
    'bulk modulus kittel': 83.3, # GPa
    'symbols': 'Zr',
    'latex name': "Zr",
    'structure':'hcp',
    'magmom': None},
'Sc_hcp': {
    'name': "Sc_hcp",
    'lattice parameter': 3.309,
    'lattice parameter 2': 5.268,
    'revtpss_lp': 3.30710158266,
    'cohesive energy': 3.90,
    'bulk modulus kittel': 43.5, # GPa
    'symbols': 'Sc',
    'latex name': "Sc",
    'structure':'hcp',
    'magmom': None},
'Be_hcp': {
    'name': "Be_hcp",
    'lattice parameter': 2.286,
    'lattice parameter 2': 3.585,
    'revtpss_lp': 2.27878331698,
    'cohesive energy': 3.32,
    'bulk modulus kittel': 100.3, # GPa
    'symbols': 'Be',
    'latex name': "Be",
    'structure':'hcp',
    'magmom': None},
'Mg_hcp': {
    'name': "Mg_hcp",
    'lattice parameter': 3.209,
    'lattice parameter 2': 5.211,
    'revtpss_lp': 3.17456515719,
    'cohesive energy': 1.51,
    'bulk modulus kittel': 35.4, # GPa
    'symbols': 'Mg',
    'latex name': "Mg",
    'structure':'hcp',
    'magmom': None},
# compounds
'LiH_b1': {
    'name': 'LiH_b1',
    'lattice parameter': 3.979,
    'pbe_lp': 4.006,
    'cohesive energy': 2.49,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'LiH',
    'latex name': "LiH",
    'structure':'rocksalt',
    'magmom': None},
'LiF_b1': {
    'name': 'LiF_b1',
    'lattice parameter': 3.974,
    'pbe_lp': 4.064,
    'cohesive energy': 4.46,
    'debye temperature': 732.,
    'bulk modulus': 69.8,
    'symbols': 'LiF',
    'latex name': "LiF",
    'structure':'rocksalt',
    'magmom': None},
'LiCl_b1': {
    'name': 'LiCl_b1',
    'lattice parameter': 5.072,
    'pbe_lp': 5.147,
    'cohesive energy': 3.59,
    'debye temperature': 422.,
    'bulk modulus': 35.4,
    'symbols': 'LiCl',
    'latex name': "LiCl",
    'structure':'rocksalt',
    'magmom': None},
'NaF_b1': {
    'name': 'NaF_b1',
    'lattice parameter': 4.570,
    'pbe_lp': 4.700,
    'cohesive energy': 3.97,
    'debye temperature': 492.,
    'bulk modulus': 51.4,
    'symbols': 'NaF',
    'latex name': "NaF",
    'structure':'rocksalt',
    'magmom': None},
'NaCl_b1': {
    'name': 'NaCl_b1',
    'lattice parameter': 5.565,
    'pbe_lp': 5.695,
    'cohesive energy': 3.34,
    'debye temperature': 321.,
    'bulk modulus': 26.6,
    'symbols': 'NaCl',
    'latex name': "NaCl",
    'structure':'rocksalt',
    'magmom': None},
'MgO_b1': {
    'name': 'MgO_b1',
    'lattice parameter': 4.188,
    'pbe_lp': 4.255,
    'cohesive energy': 5.20,
    'debye temperature': 946.,
    'bulk modulus': 165.,
    'symbols': 'MgO',
    'latex name': "MgO",
    'structure':'rocksalt',
    'magmom': None},
'MgS_b1': {
    'name': 'MgS_b1',
    'lattice parameter': 5.188,
    'pbe_lp': 5.228,
    'cohesive energy': 4.07,
    'debye temperature': 650.,
    'bulk modulus': 80.,
    'symbols': 'MgS',
    'latex name': "MgS",
    'structure':'rocksalt',
    'magmom': None},
'CaO_b1': {
    'name': 'CaO_b1',
    'lattice parameter': 4.781,
    'pbe_lp': 4.832,
    'cohesive energy': 5.57,
    'debye temperature': 648.,
    'bulk modulus': 115.,
    'symbols': 'CaO',
    'latex name': "CaO",
    'structure':'rocksalt',
    'magmom': None},
'TiC_b1': {
    'name': 'TiC_b1',
    'lattice parameter': 4.318,
    'pbe_lp': 4.332,
    'cohesive energy': 7.16,
    'debye temperature': 940.,
    'bulk modulus': 233.,
    'symbols': 'TiC',
    'latex name': "TiC",
    'structure':'rocksalt',
    'magmom': None},
'TiN_b1': {
    'name': 'TiN_b1',
    'lattice parameter': 4.226,
    'pbe_lp': 4.247,
    'cohesive energy': 6.69,
    'debye temperature': 757.,
    'bulk modulus': 277.,
    'symbols': 'TiN',
    'latex name': "TiN",
    'structure':'rocksalt',
    'magmom': None},
'ZrC_b1': {
    'name': 'ZrC_b1',
    'lattice parameter': 4.687,
    'pbe_lp': 4.708,
    'cohesive energy': 7.93,
    'debye temperature': 700.,
    'bulk modulus': 207.,
    'symbols': 'ZrC',
    'latex name': "ZrC",
    'structure':'rocksalt',
    'magmom': None},
'ZrN_b1': {
    'name': 'ZrN_b1',
    'lattice parameter': 4.576,
    'pbe_lp': 4.594,
    'cohesive energy': 7.52,
    'debye temperature': 700.,
    'bulk modulus': 240.,
    'symbols': 'ZrN',
    'latex name': "ZrN",
    'structure':'rocksalt',
    'magmom': None},
'VC_b1': {
    'name': 'VC_b1',
    'lattice parameter': 4.148,
    'pbe_lp': 4.154,
    'cohesive energy': 6.94,
    'debye temperature': 971.,
    'bulk modulus': 303.,
    'symbols': 'VC',
    'latex name': "VC",
    'structure':'rocksalt',
    'magmom': None,
    'gpaw_conv': [1.e-6, 1.e-3, 1.e-6]},
'VN_b1': {
    'name': 'VN_b1',
    'lattice parameter': 4.130,
    'pbe_lp': 4.116,
    'cohesive energy': 6.25,
    'debye temperature': 755.,
    'bulk modulus': 268.,
    'symbols': 'VN',
    'latex name': "VN",
    'structure':'rocksalt',
    'magmom': None},
'NbC_b1': {
    'name': 'NbC_b1',
    'lattice parameter': 4.461,
    'pbe_lp': 4.484,
    'cohesive energy': 8.26,
    'debye temperature': 761.,
    'bulk modulus': 300.,
    'symbols': 'NbC',
    'latex name': "NbC",
    'structure':'rocksalt',
    'magmom': None},
'NbN_b1': {
    'name': 'NbN_b1',
    'lattice parameter': 4.383,
    'pbe_lp': 4.422,
    'cohesive energy': 7.50,
    'debye temperature': 730.,
    'bulk modulus': 287.,
    'symbols': 'NbN',
    'latex name': "NbN",
    'structure':'rocksalt',
    'magmom': None},
'FeAl_b2': {
    'name': 'FeAl_b2',
    'lattice parameter': 2.881,
    'pbe_lp': 2.868,
    'cohesive energy': 4.23,
    'debye temperature': 514.,
    'bulk modulus': 136.1,
    'symbols': 'FeAl',
    'latex name': "FeAl",
    'structure':'cesiumchloride',
    'magmom': 0.35},
'CoAl_b2': {
    'name': 'CoAl_b2',
    'lattice parameter': 2.854,
    'pbe_lp': 2.851,
    'cohesive energy': 4.35,
    'debye temperature': 500.,
    'bulk modulus': 162.,
    'symbols': 'CoAl',
    'latex name': "CoAl",
    'structure':'cesiumchloride',
    'magmom': None},
'NiAl_b2': {
    'name': 'NiAl_b2',
    'lattice parameter': 2.881,
    'pbe_lp': 2.892,
    'cohesive energy': 4.46,
    'debye temperature': 402.,
    'bulk modulus': 166.,
    'symbols': 'NiAl',
    'latex name': "NiAl",
    'structure':'cesiumchloride',
    'magmom': None},
'BN_b3': {
    'name': 'BN_b3',
    'lattice parameter': 3.594,
    'pbe_lp': 3.624,
    'cohesive energy': 6.76,
    'debye temperature': 1700.,
    'bulk modulus': 372.3,
    'symbols': 'BN',
    'latex name': "BN",
    'structure':'zincblende',
    'magmom': None},
'BP_b3': {
    'name': 'BP_b3',
    'lattice parameter': 4.527,
    'pbe_lp': 4.548,
    'cohesive energy': 5.14,
    'debye temperature': 985.,
    'bulk modulus': 152.,
    'symbols': 'BP',
    'latex name': "BP",
    'structure':'zincblende',
    'magmom': None},
'BAs_b3': {
    'name': 'BAs_b3',
    'lattice parameter': 4.764,
    'pbe_lp': 4.809,
    'cohesive energy': None,
    'debye temperature': 800.,
    'bulk modulus': 138.,
    'symbols': 'BAs',
    'latex name': "BAs",
    'structure':'zincblende',
    'magmom': None},
'AlN_b3': {
    'name': 'AlN_b3',
    'lattice parameter': 4.368,
    'pbe_lp': 4.402,
    'cohesive energy': 5.85,
    'debye temperature': None,
    'bulk modulus': None,
    'symbols': 'AlN',
    'latex name': "AlN",
    'structure':'zincblende',
    'magmom': None},
'AlP_b3': {
    'name': 'AlP_b3',
    'lattice parameter': 5.450,
    'pbe_lp': 5.504,
    'cohesive energy': 4.32,
    'debye temperature': 588.,
    'bulk modulus': 74.4,
    'symbols': 'AlP',
    'latex name': "AlP",
    'structure':'zincblende',
    'magmom': None},
'AlAs_b3': {
    'name': 'AlAs_b3',
    'lattice parameter': 5.649,
    'pbe_lp': 5.728,
    'cohesive energy': 3.82,
    'debye temperature': 292.,
    'bulk modulus': 74.,
    'symbols': 'AlAs',
    'latex name': "AlAs",
    'structure':'zincblende',
    'magmom': None},
'GaN_b3': {
    'name': 'GaN_b3',
    'lattice parameter': 4.523,
    'pbe_lp': 4.549,
    'cohesive energy': 4.55,
    'debye temperature': 600.,
    'bulk modulus': 195.,
    'symbols': 'GaN',
    'latex name': "GaN",
    'structure':'zincblende',
    'magmom': None},
'GaP_b3': {
    'name': 'GaP_b3',
    'lattice parameter': 5.441,
    'pbe_lp': 5.506,
    'cohesive energy': 3.61,
    'debye temperature': 445.,
    'bulk modulus': 87.4,
    'symbols': 'GaP',
    'latex name': "GaP",
    'structure':'zincblende',
    'magmom': None},
'GaAs_b3': {
    'name': 'GaAs_b3',
    'lattice parameter': 5.641,
    'pbe_lp': 5.751,
    'cohesive energy': 3.34,
    'debye temperature': 344.,
    'bulk modulus': 75.6,
    'symbols': 'GaAs',
    'latex name': "GaAs",
    'structure':'zincblende',
    'magmom': None},
'InP_b3': {
    'name': 'InP_b3',
    'lattice parameter': 5.858,
    'pbe_lp': 5.963,
    'cohesive energy': 3.47,
    'debye temperature': 321.,
    'bulk modulus': 71.,
    'symbols': 'InP',
    'latex name': "InP",
    'structure':'zincblende',
    'magmom': None},
'InAs_b3': {
    'name': 'InAs_b3',
    'lattice parameter': 6.048,
    'pbe_lp': 6.188,
    'cohesive energy': 3.08,
    'debye temperature': 247.,
    'bulk modulus': 58.,
    'symbols': 'InAs',
    'latex name': "InAs",
    'structure':'zincblende',
    'magmom': None},
'SiC_b3': {
    'name': 'SiC_b3',
    'lattice parameter': 4.348,
    'pbe_lp': 4.378,
    'cohesive energy': 6.48,
    'debye temperature': 1232.,
    'bulk modulus': 225.,
    'symbols': 'SiC',
    'latex name': "SiC",
    'structure':'zincblende',
    'magmom': None},
# Pearson
'KBr_b1': {
    'name': 'KBr_b1',
    'lattice parameter': 6.60,
    'revtpss_lp': 6.75894864163,
    'cohesive energy': 3.08,
    'symbols': 'KBr',
    'latex name': "KBr",
    'structure':'rocksalt',
    'magmom': None},
'CaSe_b1': {
    'name': 'CaSe_b1',
    'lattice parameter': 5.80,
    'revtpss_lp': 5.98316194871,
    'cohesive energy': 4.01,
    'symbols': 'CaSe',
    'latex name': "CaSe",
    'structure':'rocksalt',
    'magmom': None},
'SeAs_b1': {
    'name': 'SeAs_b1',
    'lattice parameter': 5.48,
    'revtpss_lp': 5.55199610297,
    'cohesive energy': 2.46, #4.92,
    'symbols': 'SeAs',
    'latex name': "SeAs",
    'structure':'rocksalt',
    'magmom': None},
'RbI_b1': {
    'name': 'RbI_b1',
    'lattice parameter': 7.00,
    'revtpss_lp': None,
    'cohesive energy': 2.71,
    'symbols': 'RbI',
    'latex name': "RbI",
    'structure':'rocksalt',
    'magmom': None},
'LiI_b1': {
    'name': 'LiI_b1',
    'lattice parameter': 6.00,
    'revtpss_lp': 5.98205708599,
    'cohesive energy': 2.78,
    'symbols': 'LiI',
    'latex name': "LiI",
    'structure':'rocksalt',
    'magmom': None},
'CsF_b1': {
    'name': 'CsF_b1',
    'lattice parameter': 6.02,
    'revtpss_lp': 6.21937874326,
    'cohesive energy': 3.66,
    'symbols': 'CsF',
    'latex name': "CsF",
    'structure':'rocksalt',
    'magmom': None},
'CsI_b2': {
    'name': 'CsI_b2',
    'lattice parameter': 4.42,
    'revtpss_lp': 4.6376073499,
    'cohesive energy': 2.75,
    'symbols': 'CsI',
    'latex name': "CsI",
    'structure':'cesiumchloride',
    'magmom': None},
'AgF_b1': {
    'name': 'AgF_b1',
    'lattice parameter': 4.92,
    'revtpss_lp': 4.95954689488,
    'cohesive energy': 2.95,
    'symbols': 'AgF',
    'latex name': "AgF",
    'structure':'rocksalt',
    'magmom': None},
'AgCl_b1': {
    'name': 'AgCl_b1',
    'lattice parameter': 5.56,
    'revtpss_lp': 5.55845581698,
    'cohesive energy': 2.75,
    'symbols': 'AgCl',
    'latex name': "AgCl",
    'structure':'rocksalt',
    'magmom': None},
'AgBr_b1': {
    'name': 'AgBr_b1',
    'lattice parameter': 5.78,
    'revtpss_lp': 5.7885713515,
    'cohesive energy': 2.58,
    'symbols': 'AgBr',
    'latex name': "AgBr",
    'structure':'rocksalt',
    'magmom': None},
'CaS_b1': {
    'name': 'CaS_b1',
    'lattice parameter': 5.70,
    'revtpss_lp': 5.7588414509,
    'cohesive energy': 4.81,
    'symbols': 'CaS',
    'latex name': "CaS",
    'structure':'rocksalt',
    'magmom': None},
'BaO_b1': {
    'name': 'BaO_b1',
    'lattice parameter': 5.52,
    'revtpss_lp': 5.60522427797,
    'cohesive energy': 5.10,
    'symbols': 'BaO',
    'latex name': "BaO",
    'structure':'rocksalt',
    'magmom': None},
'BaSe_b1': {
    'name': 'BaSe_b1',
    'lattice parameter': 6.60,
    'revtpss_lp': 6.67681841117,
    'cohesive energy': 4.03,
    'symbols': 'BaSe',
    'latex name': "BaSe",
    'structure':'rocksalt',
    'magmom': None},
'CdO_b1': {
    'name': 'CdO_b1',
    'lattice parameter': 4.70,
    'revtpss_lp': 4.75039242719,
    'cohesive energy': 3.21,
    'symbols': 'CdO',
    'latex name': "CdO",
    'structure':'rocksalt',
    'magmom': None},
'MnO_b1': {
    'name': 'MnO_b1',
    'lattice parameter': 4.44,
    'revtpss_lp': 4.4520326675,
    'cohesive energy': 4.75,
    'symbols': 'MnO',
    'latex name': "MnO",
    'structure':'rocksalt',
    'magmom': 2.4},
'MnS_b1': {
    'name': 'MnS_b1',
    'lattice parameter': 5.22,
    'revtpss_lp': 5.13039706051,
    'cohesive energy': 4.01,
    'symbols': 'MnS',
    'latex name': "MnS",
    'structure':'rocksalt',
    'magmom': 2.4},
# Haglund
'ScC_b1': {
    'name': 'ScC_b1',
    'lattice parameter': 4.72,
    'revtpss_lp': 4.69537210811,
    'cohesive energy': 6.37,
    'symbols': 'ScC',
    'latex name': "ScC",
    'structure':'rocksalt',
    'magmom': None},
'CrC_b1': { # slow gpaw convergence
    'name': 'CrC_b1',
    'lattice parameter': 4.12,
    'revtpss_lp': 4.12680477339,
    'cohesive energy': 5.80,
    'symbols': 'CrC',
    'latex name': "CrC",
    'structure':'rocksalt',
    'magmom': 0.6},
    #'gpaw_conv': [1.e-5, 1.e-7, 1.e-4]},
'MnC_b1': {
    'name': 'MnC_b1',
    'lattice parameter': 4.12,
    'revtpss_lp': 4.11651335966,
    'cohesive energy': 5.14,
    'symbols': 'MnC',
    'latex name': "MnC",
    'structure':'rocksalt',
    'magmom': 1.2},
    #'gpaw_conv': [1.e-5, 1.e-6, 1.e-4]},
'FeC_b1': {
    'name': 'FeC_b1',
    'lattice parameter': 4.09,
    'revtpss_lp': 3.98592429783,
    'cohesive energy': 5.67,
    'symbols': 'FeC',
    'latex name': "FeC",
    'structure':'rocksalt',
    'magmom': None},
'CoC_b1': {
    'name': 'CoC_b1',
    'lattice parameter': 4.05,
    'revtpss_lp': 3.98558621985,
    'cohesive energy': 5.69,
    'symbols': 'CoC',
    'latex name': "CoC",
    'structure':'rocksalt',
    'magmom': None},
'NiC_b1': {
    'name': 'NiC_b1',
    'lattice parameter': 3.99,
    'revtpss_lp': 4.03667215219,
    'cohesive energy': 5.65,
    'symbols': 'NiC',
    'latex name': "NiC",
    'structure':'rocksalt',
    'magmom': None},
'ScN_b1': {
    'name': 'ScN_b1',
    'lattice parameter': 4.51,
    'revtpss_lp': 4.52753210905,
    'cohesive energy': 6.72,
    'symbols': 'ScN',
    'latex name': "ScN",
    'structure':'rocksalt',
    'magmom': None},
'CrN_b1': {
    'name': 'CrN_b1',
    'lattice parameter': 4.15,
    'revtpss_lp': 4.18609030986,
    'cohesive energy': 5.14,
    'symbols': 'CrN',
    'latex name': "CrN",
    'structure':'rocksalt',
    'magmom': 1.3},
'MnN_b1': {
    'name': 'MnN_b1',
    'lattice parameter': 4.20,
    'revtpss_lp': 4.15891205836,
    'cohesive energy': 4.08,
    'symbols': 'MnN',
    'latex name': "MnN",
    'structure':'rocksalt',
    'magmom': 1.6},
    #'gpaw_conv': [1.e-5, 1.e-7, 1.e-4]},
'FeN_b1': {
    'name': 'FeN_b1',
    'lattice parameter': 4.13,
    'revtpss_lp': 4.12,
    'cohesive energy': 4.59,
    'symbols': 'FeN',
    'latex name': "FeN",
    'structure':'rocksalt',
    'magmom': 1.3},
    #'gpaw_conv': [1.e-5, 1.e-4, 1.e-7]},
'CoN_b1': {
    'name': 'CoN_b1',
    'lattice parameter': 4.10,
    'revtpss_lp': 3.99392851055,
    'cohesive energy': 4.53,
    'symbols': 'CoN',
    'latex name': "CoN",
    'structure':'rocksalt',
    'magmom': None},
'NiN_b1': {
    'name': 'NiN_b1',
    'lattice parameter': 4.10,
    'revtpss_lp': 4.03914234139,
    'cohesive energy': 4.48,
    'symbols': 'NiN',
    'latex name': "NiN",
    'structure':'rocksalt',
    'magmom': None},
'MoC_b1': {
    'name': 'MoC_b1',
    'lattice parameter': 4.278,
    'revtpss_lp': 4.39395801553,
    'cohesive energy': 7.22,
    'symbols': 'MoC',
    'latex name': "MoC",
    'structure':'rocksalt',
    'magmom': None},
'RuC_b1': {
    'name': 'RuC_b1',
    'lattice parameter': 4.129,
    'revtpss_lp': 4.33092551761,
    'cohesive energy': 6.73,
    'symbols': 'RuC',
    'latex name': "RuC",
    'structure':'rocksalt',
    'magmom': None},
'RhC_b1': {
    'name': 'RhC_b1',
    'lattice parameter': 4.145,
    'revtpss_lp': 4.35180356284,
    'cohesive energy': 6.23,
    'symbols': 'RhC',
    'latex name': "RhC",
    'structure':'rocksalt',
    'magmom': None},
'PdC_b1': {
    'name': 'PdC_b1',
    'lattice parameter': 4.221,
    'revtpss_lp': 4.4101528149,
    'cohesive energy': 5.36,
    'symbols': 'PdC',
    'latex name': "PdC",
    'structure':'rocksalt',
    'magmom': None},
'MoN_b1': {
    'name': 'MoN_b1',
    'lattice parameter': 4.214,
    'revtpss_lp': 4.37920656361,
    'cohesive energy': 6.20,
    'symbols': 'MoN',
    'latex name': "MoN",
    'structure':'rocksalt',
    'magmom': None},
'RuN_b1': {
    'name': 'RuN_b1',
    'lattice parameter': 4.058,
    'revtpss_lp': 4.34894810722,
    'cohesive energy': 5.02,
    'symbols': 'RuN',
    'latex name': "RuN",
    'structure':'rocksalt',
    'magmom': None},
'RhN_b1': {
    'name': 'RhN_b1',
    'lattice parameter': 4.082,
    'revtpss_lp': 4.35633367368,
    'cohesive energy': 4.78,
    'symbols': 'RhN',
    'latex name': "RhN",
    'structure':'rocksalt',
    'magmom': None},
'PdN_b1': {
    'name': 'PdN_b1',
    'lattice parameter': 4.145,
    'revtpss_lp': 4.41700841872,
    'cohesive energy': 4.03,
    'symbols': 'PdN',
    'latex name': "PdN",
    'structure':'rocksalt',
    'magmom': None},
'LaC_b1': {
    'name': 'LaC_b1',
    'lattice parameter': 5.429,
    'revtpss_lp': 5.50510598304,
    'cohesive energy': 5.74,
    'symbols': 'LaC',
    'latex name': "LaC",
    'structure':'rocksalt',
    'magmom': None},
'TaC_b1': {
    'name': 'TaC_b1',
    'lattice parameter': 4.457,
    'revtpss_lp': 4.48669679096,
    'cohesive energy': 8.56,
    'symbols': 'TaC',
    'latex name': "TaC",
    'structure':'rocksalt',
    'magmom': None},
'WC_b1': {
    'name': 'WC_b1',
    'lattice parameter': 4.266,
    'revtpss_lp': 4.40260619124,
    'cohesive energy': 8.25,
    'symbols': 'WC',
    'latex name': "WC",
    'structure':'rocksalt',
    'magmom': None},
'OsC_b1': {
    'name': 'OsC_b1',
    'lattice parameter': 4.176,
    'revtpss_lp': 4.35787411212,
    'cohesive energy': 7.36,
    'symbols': 'OsC',
    'latex name': "OsC",
    'structure':'rocksalt',
    'magmom': None},
'IrC_b1': {
    'name': 'IrC_b1',
    'lattice parameter': 4.129,
    'revtpss_lp': 4.39419683105,
    'cohesive energy': 6.84,
    'symbols': 'IrC',
    'latex name': "IrC",
    'structure':'rocksalt',
    'magmom': None},
'PtC_b1': {
    'name': 'PtC_b1',
    'lattice parameter': 4.206,
    'revtpss_lp': 4.45741920763,
    'cohesive energy': 6.34,
    'symbols': 'PtC',
    'latex name': "PtC",
    'structure':'rocksalt',
    'magmom': None},
'LaN_b1': {
    'name': 'LaN_b1',
    'lattice parameter': 5.305,
    'revtpss_lp': 5.35277831785,
    'cohesive energy': 6.27,
    'symbols': 'LaN',
    'latex name': "LaN",
    'structure':'rocksalt',
    'magmom': None},
'TaN_b1': {
    'name': 'TaN_b1',
    'lattice parameter': 4.340,
    'revtpss_lp': 4.44437817766,
    'cohesive energy': 7.63,
    'symbols': 'TaN',
    'latex name': "TaN",
    'structure':'rocksalt',
    'magmom': None},
'WN_b1': {
    'name': 'WN_b1',
    'lattice parameter': 4.202,
    'revtpss_lp': 4.39478109563,
    'cohesive energy': 7.01,
    'symbols': 'WN',
    'latex name': "WN",
    'structure':'rocksalt',
    'magmom': None},
'OsN_b1': {
    'name': 'OsN_b1',
    'lattice parameter': 4.058,
    'revtpss_lp': 4.37261212422,
    'cohesive energy': 5.61,
    'symbols': 'OsN',
    'latex name': "OsN",
    'structure':'rocksalt',
    'magmom': None},
'IrN_b1': {
    'name': 'IrN_b1',
    'lattice parameter': 4.074,
    'revtpss_lp': 4.40354481322,
    'cohesive energy': 5.13,
    'symbols': 'IrN',
    'latex name': "IrN",
    'structure':'rocksalt',
    'magmom': None},
'PtN_b1': {
    'name': 'PtN_b1',
    'lattice parameter': 4.137,
    'revtpss_lp': 4.48127067537,
    'cohesive energy': 4.63,
    'symbols': 'PtN',
    'latex name': "PtN",
    'structure':'rocksalt',
    'magmom': None},
}

def in_data(name):
    if name not in data:
        raise KeyError('System %s not in database.' % name)

def get_solid_lattice_parameter(name):
    in_data(name)
    return data[name]['lattice parameter']

def get_hcp_covera(name):
    in_data(name)
    if name in hcp_solids_10:
        a = data[name]['lattice parameter']
        c = data[name]['lattice parameter 2']
        return c/a
    else:
        return 0.

def get_solid_pbe_lp(name):
    in_data(name)
    if name in solids_1:
        return data[name]['pbe_lp']
    else:
        return None

def get_solid_revtpss_lp(name):
    in_data(name)
    if name in solids_2:
        return data[name]['revtpss_lp']
    else:
        return None

def lattice_parameter(name):
    in_data(name)
    assert name in solids_0, name
    if name in solids_1:
        return get_solid_lattice_parameter(name)
    elif name in solids_2:
        return get_solid_revtpss_lp(name) - 0.043 #KRIS: WHY????
    else:
        kdkdkd

def get_solid_cohesive_energy(name, ZPVE_corr=True):
    in_data(name)
    e = data[name]['cohesive energy']
    if ZPVE_corr == True and name in qubic_solids_27:
        e += (9./8.)*kB*data[name]['debye temperature']
    return e

def get_solid_bulk_modulus(name):
    in_data(name)
    if name in solids_1:
        return data[name]['bulk modulus']
    else:
        return None

def get_solid_magmom(name):
    in_data(name)
    return data[name]['magmom']

def get_solid_crystal_structure(name):
    in_data(name)
    struct = data[name]['structure']
    assert struct in ['fcc', 'bcc', 'diamond', 'hcp',
                      'rocksalt', 'cesiumchloride', 'zincblende']
    return struct

def get_solid_symbols(name):
    in_data(name)
    return data[name]['symbols']

def get_solids_latex_name(name):
    in_data(name)
    d = data[name]
    struct = d['structure']
    if struct == 'diamond':
        struct = 'dia'
    return d['latex name'] + ' ('+struct+')'

def get_solids_common_name(name):
    in_data(name)
    return get_solids_latex_name(name)

def get_solid_bulk_gpaw_convergence(name, xc=None):
    in_data(name)
    d = data[name]
    std = [1.e-6, 1.e-7, 1.e-6] # energy, eigenstates, density
    if name == 'LiH_b1':
        if xc in ['TPSS', 'revTPSS']:
            return [1.e-5, 1.e-4, 1.e-6]
        else:
            return std
    elif 'gpaw_conv' not in d:
        return std
    else:
        return d['gpaw_conv']

# atomic constituents
atoms_ = np.array([])
for aaa in solids_1:
    symbs = get_solid_symbols(aaa)
    atoms_ = np.append(atoms_, string2symbols(symbs))
solids_1_unique_atoms = np.unique(atoms_)
del atoms_, aaa, symbs

atoms_ = np.array([])
for aaa in solids_2:
    symbs = get_solid_symbols(aaa)
    atoms_ = np.append(atoms_, string2symbols(symbs))
solids_2_unique_atoms = np.unique(atoms_)
del atoms_, aaa, symbs

solids_0_unique_atoms = np.append(solids_1_unique_atoms, solids_2_unique_atoms)
solids_0_unique_atoms = np.unique(solids_0_unique_atoms)

getKey={'Li-bcc':'Li_bcc'
        ,'IrC-rocksalt': 'IrC_b1'
        , 'NaF-rocksalt': 'NaF_b1'
        , 'NiC-rocksalt': 'NiC_b1'
        , 'AlN-zincblende': 'AlN_b3'
        , 'Sn-diamond': 'Sn_diamond'
        , 'GaN-zincblende': 'GaN_b3'
        , 'TaC-rocksalt': 'TaC_b1'
        , 'FeAl-cscl': 'FeAl_b2'
        , 'CoN-rocksalt': 'CoN_b1'
        , 'FeN-rocksalt': 'FeN_b1'
        , 'CrC-rocksalt': 'CrC_b1'
        , 'Fe-bcc': 'Fe_bcc'
        , 'VN-rocksalt': 'VN_b1'
        , 'BaSe-rocksalt': 'BaSe_b1'
        , 'AlAs-zincblende': 'AlAs_b3'
        , 'Os-hcp': 'Os_hcp', 'VC-rocksalt': 'VC_b1'
        , 'GaP-zincblende': 'GaP_b3'
        , 'RbI-rocksalt': 'RbI_b1'
        , 'C-diamond': 'C_diamond'
        , 'Sc-hcp': 'Sc_hcp'
        , 'Ru-hcp': 'Ru_hcp'
        , 'TiC-rocksalt': 'TiC_b1'
        , 'Co-hcp': 'Co_hcp'
        , 'MgO-rocksalt': 'MgO_b1'
        , 'Rh-fcc': 'Rh_fcc'
        , 'CsI-cscl': 'CsI_b2'
        , 'Zr-hcp': 'Zr_hcp'
        , 'Cd-hcp': 'Cd_hcp'
        , 'Ni-fcc': 'Ni_fcc'
        , 'Nb-bcc': 'Nb_bcc'
        , 'WN-rocksalt': 'WN_b1'
        , 'NaCl-rocksalt': 'NaCl_b1'
        , 'LiH-rocksalt': 'LiH_b1'
        , 'LaN-rocksalt': 'LaN_b1'
        , 'MgS-rocksalt': 'MgS_b1'
        , 'CdO-rocksalt': 'CdO_b1'
        , 'ZrC-rocksalt': 'ZrC_b1'
        , 'NbN-rocksalt': 'NbN_b1'
        , 'MnC-rocksalt': 'MnC_b1'
        , 'NbC-rocksalt': 'NbC_b1'
        , 'InAs-zincblende': 'InAs_b3'
        , 'K-bcc': 'K_bcc'
        , 'PtC-rocksalt': 'PtC_b1'
        , 'PdN-rocksalt': 'PdN_b1'
        , 'LiF-rocksalt': 'LiF_b1'
        , 'Si-diamond': 'Si_diamond'
        , 'BN-zincblende': 'BN_b3'
        , 'BaO-rocksalt': 'BaO_b1'
        , 'PtN-rocksalt': 'PtN_b1'
        , 'Ba-bcc': 'Ba_bcc'
        , 'Ti-hcp': 'Ti_hcp'
        , 'V-bcc': 'V_bcc'
        , 'Au-fcc': 'Au_fcc'
        , 'Pt-fcc': 'Pt_fcc'
        , 'Zn-hcp': 'Zn_hcp'
        , 'Rb-bcc': 'Rb_bcc'
        , 'CaS-rocksalt': 'CaS_b1'
        , 'CoC-rocksalt': 'CoC_b1'
        , 'Sr-fcc': 'Sr_fcc'
        , 'RuC-rocksalt': 'RuC_b1'
        , 'GaAs-zincblende': 'GaAs_b3'
        , 'IrN-rocksalt': 'IrN_b1'
        , 'MoC-rocksalt': 'MoC_b1'
        , 'OsC-rocksalt': 'OsC_b1'
        , 'InP-zincblende': 'InP_b3'
        , 'BAs-rocksalt': 'BAs_b3'
        , 'MnN-rocksalt': 'MnN_b1'
        , 'Ca-fcc': 'Ca_fcc'
        , 'NiAl-cscl': 'NiAl_b2'
        , 'TiN-rocksalt': 'TiN_b1'
        , 'Pd-fcc': 'Pd_fcc'
        , 'RuN-rocksalt': 'RuN_b1'
        , 'Na-bcc': 'Na_bcc'
        , 'CrN-rocksalt': 'CrN_b1'
        , 'Al-fcc': 'Al_fcc'
        , 'RhC-rocksalt': 'RhC_b1'
        , 'Cu-fcc': 'Cu_fcc'
        , 'AgCl-rocksalt': 'AgCl_b1'
        , 'ZrN-rocksalt': 'ZrN_b1'
        , 'LaC-rocksalt': 'LaC_b1'
        , 'Mo-bcc': 'Mo_bcc'
        , 'Ag-fcc': 'Ag_fcc'
        , 'Mg-hcp': 'Mg_hcp'
        , 'AgBr-rocksalt': 'AgBr_b1'
        , 'ScC-rocksalt': 'ScC_b1'
        , 'Ge-diamond': 'Ge_diamond'
        , 'AlP-zincblende': 'AlP_b3'
        , 'MnS-rocksalt': 'MnS_b1'
        , 'W-bcc': 'W_bcc'
        , 'WC-rocksalt': 'WC_b1'
        , 'MnO-rocksalt': 'MnO_b1'
        , 'TaN-rocksalt': 'TaN_b1'
        , 'ScN-rocksalt': 'ScN_b1'
        , 'FeC-rocksalt': 'FeC_b1'
        , 'OsN-rocksalt': 'OsN_b1'
        , 'Pb-fcc': 'Pb_fcc'
        , 'Be-hcp': 'Be_hcp'
        , 'Ta-bcc': 'Ta_bcc'
        , 'CoAl-cscl': 'CoAl_b2'
        , 'Ir-fcc': 'Ir_fcc'
        , 'CaO-rocksalt': 'CaO_b1'
        , 'PdC-rocksalt': 'PdC_b1'
        , 'RhN-rocksalt': 'RhN_b1'
        , 'MoN-rocksalt': 'MoN_b1'
        , 'NiN-rocksalt': 'NiN_b1'
        , 'AgF-rocksalt': 'AgF_b1'
        , 'LiCl-rocksalt': 'LiCl_b1'
        , 'CsF-rocksalt': 'CsF_b1'
        , 'BP-zincblende': 'BP_b3'
        , 'SiC-zincblende': 'SiC_b3'}
