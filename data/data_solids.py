
from ase.atoms import Atoms, string2symbols   # type: ignore
from ase.units import kB, kJ, kcal, mol   # type: ignore
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
#SMS 03/2016
#Sol27LC
sol27_lp = sol57_lp[0:26]
sol27_lp.append('Pb_fcc')
assert len(sol27_lp) == 27
data_info = {
            'sol27_lp':sol27_lp,
            'slab':[],
            'molecule':[],
            'ads':[],
            'dbh24':[],
            'bulk':[]
            }
#SMS 03/2016

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
