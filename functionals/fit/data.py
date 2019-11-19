# External
from typing import (Any,
                    Set      as S,
                    Dict     as D,
                    List     as L,
                    Tuple    as T,
                    Optional as O,
                    Callable as C)

from numpy        import array,empty,vstack # type: ignore
from random       import shuffle

'''Data Preprocessing for nonlinear fitting.'''

################################################################################
bmdata = ['Ag_fcc', 'AlAs_b3', 'AlP_b3', 'Al_fcc', 'Au_fcc', 'BAs_b3', 'BN_b3', 'BP_b3', 'Ba_bcc', 'Be_hcp', 'C_diamond', 'CaO_b1', 'Ca_fcc', 'Cd_hcp', 'CoAl_b2', 'Co_hcp', 'Cu_fcc', 'FeAl_b2', 'Fe_bcc', 'GaAs_b3', 'GaN_b3', 'GaP_b3', 'Ge_diamond', 'InAs_b3', 'InP_b3', 'Ir_fcc', 'K_bcc', 'LiCl_b1', 'LiF_b1', 'Li_bcc', 'MgO_b1', 'MgS_b1', 'Mg_hcp', 'Mo_bcc', 'NaCl_b1', 'NaF_b1', 'Na_bcc', 'NbC_b1', 'NbN_b1', 'Nb_bcc', 'NiAl_b2', 'Ni_fcc', 'Os_hcp', 'Pb_fcc', 'Pd_fcc', 'Pt_fcc', 'Rb_bcc', 'Rh_fcc', 'Ru_hcp', 'Sc_hcp', 'SiC_b3', 'Si_diamond', 'Sn_diamond', 'Sr_fcc', 'Ta_bcc', 'TiC_b1', 'TiN_b1', 'Ti_hcp', 'VC_b1', 'VN_b1', 'V_bcc', 'W_bcc', 'Zn_hcp', 'ZrC_b1', 'ZrN_b1', 'Zr_hcp']

# Note that RbI and FeAl have lattice constant data but no BM data, and we use
lcdata = ['AgBr_b1', 'AgCl_b1', 'AgF_b1', 'Ag_fcc', 'AlAs_b3', 'AlN_b3', 'AlP_b3', 'Al_fcc', 'Au_fcc', 'BAs_b3', 'BN_b3', 'BP_b3', 'BaO_b1', 'BaSe_b1', 'Ba_bcc', 'Be_hcp', 'C_diamond', 'CaO_b1', 'CaS_b1', 'CaSe_b1', 'Ca_fcc', 'CdO_b1', 'Cd_hcp', 'CoAl_b2', 'CoC_b1', 'CoN_b1', 'Co_hcp', 'CrC_b1', 'CrN_b1', 'CsF_b1', 'CsI_b2', 'Cu_fcc', 'FeAl_b2', 'FeC_b1', 'FeN_b1', 'Fe_bcc', 'GaAs_b3', 'GaN_b3', 'GaP_b3', 'Ge_diamond', 'InAs_b3', 'InP_b3', 'IrC_b1', 'IrN_b1', 'Ir_fcc', 'KBr_b1', 'K_bcc', 'LaC_b1', 'LaN_b1', 'LiCl_b1', 'LiF_b1', 'LiH_b1', 'LiI_b1', 'Li_bcc', 'MgO_b1', 'MgS_b1', 'Mg_hcp', 'MnC_b1', 'MnN_b1', 'MnO_b1', 'MnS_b1', 'MoC_b1', 'MoN_b1', 'Mo_bcc', 'NaCl_b1', 'NaF_b1', 'Na_bcc', 'NbC_b1', 'NbN_b1', 'Nb_bcc', 'NiAl_b2', 'NiC_b1', 'NiN_b1', 'Ni_fcc', 'OsC_b1', 'OsN_b1', 'Os_hcp', 'Pb_fcc', 'PdC_b1', 'PdN_b1', 'Pd_fcc', 'PtC_b1', 'PtN_b1', 'Pt_fcc', 'RbI_b1', 'Rb_bcc', 'RhC_b1', 'RhN_b1', 'Rh_fcc', 'RuC_b1', 'RuN_b1', 'Ru_hcp', 'ScC_b1', 'ScN_b1', 'Sc_hcp', 'SeAs_b1', 'SiC_b3', 'Si_diamond', 'Sn_diamond', 'Sr_fcc', 'TaC_b1', 'TaN_b1', 'Ta_bcc', 'TiC_b1', 'TiN_b1', 'Ti_hcp', 'VC_b1', 'VN_b1', 'V_bcc', 'WC_b1', 'WN_b1', 'W_bcc', 'Zn_hcp', 'ZrC_b1', 'ZrN_b1', 'Zr_hcp']

cedata = ['AgBr_b1', 'AgCl_b1', 'AgF_b1', 'Ag_fcc', 'AlAs_b3', 'AlN_b3', 'AlP_b3', 'Al_fcc', 'Au_fcc', 'BN_b3', 'BP_b3', 'BaO_b1', 'BaSe_b1', 'Ba_bcc', 'Be_hcp', 'C_diamond', 'CaO_b1', 'CaS_b1', 'CaSe_b1', 'Ca_fcc', 'CdO_b1', 'Cd_hcp', 'CoAl_b2', 'CoC_b1', 'CoN_b1', 'Co_hcp', 'CrC_b1', 'CrN_b1', 'CsF_b1', 'CsI_b2', 'Cu_fcc', 'FeAl_b2', 'FeC_b1', 'FeN_b1', 'Fe_bcc', 'GaAs_b3', 'GaN_b3', 'GaP_b3', 'Ge_diamond', 'InAs_b3', 'InP_b3', 'IrC_b1', 'IrN_b1', 'Ir_fcc', 'KBr_b1', 'K_bcc', 'LaC_b1', 'LaN_b1', 'LiCl_b1', 'LiF_b1', 'LiH_b1', 'LiI_b1', 'Li_bcc', 'MgO_b1', 'MgS_b1', 'Mg_hcp', 'MnC_b1', 'MnN_b1', 'MnO_b1', 'MnS_b1', 'MoC_b1', 'MoN_b1', 'Mo_bcc', 'NaCl_b1', 'NaF_b1', 'Na_bcc', 'NbC_b1', 'NbN_b1', 'Nb_bcc', 'NiAl_b2', 'NiC_b1', 'NiN_b1', 'Ni_fcc', 'OsC_b1', 'OsN_b1', 'Os_hcp', 'Pb_fcc', 'PdC_b1', 'PdN_b1', 'Pd_fcc', 'PtC_b1', 'PtN_b1', 'Pt_fcc', 'RbI_b1', 'Rb_bcc', 'RhC_b1', 'RhN_b1', 'Rh_fcc', 'RuC_b1', 'RuN_b1', 'Ru_hcp', 'ScC_b1', 'ScN_b1', 'Sc_hcp', 'SeAs_b1', 'SiC_b3', 'Si_diamond', 'Sn_diamond', 'Sr_fcc', 'TaC_b1', 'TaN_b1', 'Ta_bcc', 'TiC_b1', 'TiN_b1', 'Ti_hcp', 'VC_b1', 'VN_b1', 'V_bcc', 'WC_b1', 'WN_b1', 'W_bcc', 'Zn_hcp', 'ZrC_b1', 'ZrN_b1', 'Zr_hcp']

mags = ['MnC_b1','MnO_b1','FeN_b1','Fe_bcc','CrC_b1','CrN_b1','Ni_fcc','FeAl_b2','Co_hcp','MnS_b1','MnN_b1']

bmsplit = [['Ag_fcc', 'Rh_fcc', 'Ta_bcc', 'Ti_hcp',     'VC_b1',  'BAs_b3',  'CaO_b1'],  # 0
           ['Al_fcc', 'Sr_fcc', 'V_bcc',  'Zn_hcp',     'ZrC_b1', 'GaAs_b3', 'MgO_b1'],  # 1
           ['Au_fcc', 'Ba_bcc', 'W_bcc',  'Zr_hcp',     'N_b1',   'InAs_b3', 'LiF_b1'],  # 2
           ['Ca_fcc', 'Fe_bcc', 'Be_hcp', 'C_diamond',  'BN_b3',  'AlP_b3',  'NaF_b1'],  # 3
           ['Cu_fcc', 'K_bcc',  'Cd_hcp', 'Ge_diamond', 'GaN_b3', 'BP_b3',   'LiCl_b1'], # 4
           ['Ir_fcc', 'Li_bcc', 'Co_hcp', 'Si_diamond', 'NbN_b1', 'GaP_b3',  'NaCl_b1'], # 5
           ['Ni_fcc', 'Mo_bcc', 'Mg_hcp', 'Sn_diamond', 'TiN_b1', 'InP_b3'],             # 6
           ['Pb_fcc', 'Na_bcc', 'Os_hcp', 'NbC_b1',     'VN_b1',  'CoAl_b2'],            # 7
           ['Pd_fcc', 'Nb_bcc', 'Ru_hcp', 'SiC_b3',     'ZrN_b1', 'FeAl_b2'],            # 8
           ['Pt_fcc', 'Rb_bcc', 'Sc_hcp', 'TiC_b1',     'AlAs_b3','NiAl_b2'],            # 9
           ]

new = [x for x in lcdata if x not in bmdata]
# Things not in bmdata
newsplit = [ ['CoC_b1', 'PtC_b1', 'IrN_b1', 'ScN_b1', 'LiH_b1'],  # 0
             ['CrC_b1', 'RhC_b1', 'LaN_b1', 'TaN_b1', 'CsI_b2'], # 1
             ['FeC_b1', 'RuC_b1', 'MnN_b1', 'WN_b1',  'LiI_b1'],  # 2
             ['IrC_b1', 'ScC_b1', 'MoN_b1', 'SeAs_b1','RbI_b1'], # 3
             ['LaC_b1', 'TaC_b1', 'NiN_b1', 'BaO_b1', 'CaS_b1'], # 4
             ['MnC_b1', 'WC_b1',  'OsN_b1', 'CdO_b1', 'MnS_b1'], # 5
             ['MoC_b1', 'AlN_b3', 'PdN_b1', 'MnO_b1', 'AgBr_b1'], # 6
             ['NiC_b1', 'CoN_b1', 'PtN_b1', 'AgF_b1', 'BaSe_b1'], # 7
             ['OsC_b1', 'CrN_b1', 'RhN_b1', 'CsF_b1', 'CaSe_b1'],# 8
             ['PdC_b1', 'FeN_b1', 'RuN_b1', 'AgCl_b1', 'KBr_b1'],# 9
             ]
cesplit = [n+bmsplit[(i+3)%10] for i,n in enumerate(newsplit)]
lcsplit = [n+bmsplit[(i+6)%10] for i,n in enumerate(reversed(newsplit))]

class Datum(object):
    '''Something to be fit: a CE/BM/LC'''
    def __init__(self, mat : str, kind : str, vec : L[float], offset : float, target : float)->None:
        assert kind in ['ce', 'bm', 'lc']
        self.mat    = mat
        self.kind   = kind
        self.vec    = vec
        self.offset = offset
        self.target = target

    def __eq__(self, other:object) -> bool:
        return False if not isinstance(other,Datum) else \
             vars(self) == vars(other)

    def __hash__(self) -> int:
        return hash((self.mat,self.kind))

    def err(self, x : array, vol : bool = False) -> float:
        '''
        Compute error for this data point:
         - for lc, compute either vol or lattice constant
         '''
        y = array(self.vec) @ x + self.offset
        if vol or self.kind != 'lc': return y - self.target
        else:
            hcp     = 'hcp' in self.mat
            factor  = 1./1.4142 if hcp else 1
            new_y   = (max(y,0)*factor)**(1/3)
            new_tar = (self.target*factor)**(1/3)
            return new_y - new_tar

class Data(object):
    def __init__(self, data : S[Datum],full : bool = True) -> None:
        self.ce = sorted([d for d in data if d.kind == 'ce'], key = lambda x: x.mat)
        self.bm = sorted([d for d in data if d.kind == 'bm'], key = lambda x: x.mat)
        self.lc = sorted([d for d in data if d.kind == 'lc'], key = lambda x: x.mat)

        assert self.ce and self.bm and self.lc

        if full: # for now, full dataset is fixed:
            missing_ce = set(cedata) - set([x.mat for x in self.ce])
            missing_bm = set(bmdata) - set([x.mat for x in self.bm])
            missing_lc = set(lcdata) - set([x.mat for x in self.lc])
            print('SKIPPING CHECK')
            if False:
                assert not missing_ce, missing_ce
                assert not missing_lc, missing_lc
                assert not missing_bm, missing_bm

    def __eq__(self, other : object) -> bool:
        return False if not isinstance(other,Data) else \
             all([getattr(self,x)==getattr(other,x) for x in ['ce','bm','lc']])

    def __len__(self)->int: return len(self.ce) + len(self.bm) + len(self.lc)

    @property
    def data(self)->S[Datum]:
        return set(self.ce) | set(self.bm) | set(self.lc)

    def mse(self,x:array,key:str)->float:
        '''Mean squared error'''
        assert key in ['ce','bm','lc','vol']
        if key in ['ce','bm','lc']:
            return sum(d.err(x)**2 for d in getattr(self,key))/len(getattr(self,key))
        elif key == 'vol':
            return sum(d.err(x,vol=True)**2 for d in self.lc)/len(self.lc)
        else: raise ValueError()


    def xy(self, ce : float, bm : float, lc : float, mag:float)->T[array,array]:
        '''Matrices for fast computation of cost function, weighted by scale.'''
        X,Y = empty((0,64)),[]
        for d in self.data:
            scale = dict(ce=ce,bm=bm,lc=lc*(mag if d.mat in mags else 1))
            s = scale[d.kind]
            assert s >= 0
            if s > 0:
                X = vstack((X,array([d.vec])/s))
                Y.append((d.target - d.offset)/s)
        return X,array(Y)

    @classmethod
    def from_list(cls,xs:list,full:bool=False) -> 'Data':
        return cls({Datum(**x) for x in xs},full=full)

    def to_list(self)->L[dict]: return [vars(d) for d in self.data]

    def split(self, i : int) -> T['Data','Data']:
        '''Return a 90/10 split of the data for cross validation.'''
        assert i in range(10)
        test =  {x for x in self.ce if x.mat in cesplit[i]} \
              | {x for x in self.bm if x.mat in bmsplit[i]} \
              | {x for x in self.lc if x.mat in lcsplit[i]}
        testd  = Data(test,full=False)
        traind = Data(self.data - test,full=False)
        return traind,testd
