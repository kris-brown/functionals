# External
from typing import Tuple as T
from json   import loads,dumps
from ast    import literal_eval
import numpy as np # type: ignore

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, AND, Env, defaultEnv, Import,
                     Literal as Lit,  EQ, JPath, Constraint, LEFT, One, SUM,
                     LIKE, CONCAT, LT, ABS, GROUP_CONCAT, NOT, NULL, Const)
################################################################################
################################################################################
################################################################################
# Constants
allspecies = \
    ['Li_bcc','Na_bcc','K_bcc','Rb_bcc',
                       'Ca_fcc','Sr_fcc','Ba_bcc',
                       'Nb_bcc','Ta_bcc',
                       'Mo_bcc','W_bcc','Fe_bcc',
                       'Rh_fcc','Ir_fcc',
                       'Ni_fcc','Pd_fcc','Pt_fcc',
                       'Cu_fcc','Ag_fcc','Au_fcc',
                       'Al_fcc','Pb_fcc',
                       'C_diamond','Si_diamond','Ge_diamond','Sn_diamond'] +\
    ['Cd_hcp', 'Co_hcp', 'Os_hcp', 'Ru_hcp', 'Zn_hcp',
                 'Zr_hcp', 'Sc_hcp', 'Be_hcp', 'Mg_hcp', 'Ti_hcp'] +\
    ['LiH_b1', 'LiF_b1', 'LiCl_b1', 'NaF_b1', 'NaCl_b1',
                          'MgO_b1', 'MgS_b1', 'CaO_b1', 'TiC_b1', 'TiN_b1',
                          'ZrC_b1', 'ZrN_b1', 'NbC_b1', 'NbN_b1',
                          'FeAl_b2', 'CoAl_b2', 'NiAl_b2', 'BN_b3', 'BP_b3',
                          'AlN_b3', 'AlP_b3', 'AlAs_b3', 'GaN_b3',
                          'GaP_b3', 'GaAs_b3', 'InP_b3', 'InAs_b3', 'SiC_b3'] +\
    ['KBr_b1', 'CaSe_b1', 'SeAs_b1', 'LiI_b1',
                            'CsF_b1', 'CsI_b2', 'AgF_b1', 'AgCl_b1', 'AgBr_b1',
                            'CaS_b1', 'BaO_b1', 'BaSe_b1', 'CdO_b1', 'MnO_b1',
                            'MnS_b1', 'RbI_b1',] + \
   ['ScC_b1', 'MnC_b1', 'FeC_b1', 'CoC_b1', 'NiC_b1',
                        'ScN_b1', 'CrN_b1', 'MnN_b1', 'CoN_b1', 'NiN_b1',
                        'MoC_b1', 'CrC_b1', 'RuC_b1', 'RhC_b1', 'PdC_b1',
                        'MoN_b1', 'RuN_b1', 'RhN_b1', 'PdN_b1',
                        'LaC_b1', 'TaC_b1', 'WC_b1', 'OsC_b1', 'IrC_b1',  'PtC_b1',
                        'LaN_b1', 'TaN_b1', 'WN_b1', 'OsN_b1', 'IrN_b1',  'PtN_b1', 'FeN_b1']
def analysis(mod : Model) -> None:
    # Extract tables
    tabs = ['Functional','Beef','Bulks','Atoms','Refs','Job','Calc']
    Fx,Beef,Bulks,Atoms,Refs,Job,Calc = map(mod.get,tabs)

    refs__atoms, refs__bulks,job__calc,bulks__job,atoms__job = \
        map(mod.get_rel,[Refs.r('atoms'),Refs.r('bulks'),Job.r('calc'),
                         Bulks.r('job'),Atoms.r('job')])
    ########################################################################
    ########################################################################
    ########################################################################
    ibq = Query(dict(f = Fx.id(),
                     b = LEFT(Fx['data'](), One) |EQ| Lit('[')))
    isbeef =                                                                    \
        Gen(name    = 'isbeef',
            desc    = 'Determines whether fx was BEEF style',
            query   = ibq,
            actions = [Fx(functional=ibq['f'],beef=ibq['b'])])
    ########################################################################
    bcols = ['data']+ ['a1'+x for x in '12345']+ ['msb']

    bq    = Query(dict(f = Fx.id(), b  =Fx['data']()),
                  constr = Fx['beef']())

    def extract_beef(j:str)->T[str,float,float,float,float,float,float]:
        try:
            a,b,c,d,e,f,g = loads(j)
            return dumps(a),float(b),float(c),float(d),float(e),float(f),float(g)
        except Exception as e:
            print(e); import pdb;pdb.set_trace(); assert False
    bpb = PyBlock(extract_beef, args = [bq['b']], outnames = bcols)

    beef =                                                                      \
        Gen(name    = 'pop_beef',
            desc    = 'Populates BEEF table',
            query   = bq,
            funcs   = [bpb],
            actions = [Beef(insert=True,
                            functional= bq['f'],
                            **{x:bpb[x] for x in bcols})])
    ############################################################################
    name  = CONCAT(Lit('%'),Atoms['name'](),Lit('%'))
    p1,p2 = [JPath('calc',[job__calc,x]) for x in [bulks__job,atoms__job]]

    erq = Query(exprs = dict(b = Bulks.id(), a = Atoms.id(),
                             c = Bulks['composition'](),
                             n = Atoms['num']()),
                basis  = [Bulks,Atoms],
                constr = EQ(Calc.id(p1),Calc.id(p2))
                            |AND| LIKE(Bulks['name'](),name)
                            |AND| Atoms['int_occupation']()
                            |AND| LT(ABS(Atoms['true_mag']()-Atoms['mag']()),Lit(0.05)))

    erpb = PyBlock(lambda c,n: literal_eval(c).get(n,[]),
                   env  = Env(Import('ast','literal_eval')),
                   args = [erq['c'],erq['n']])

    exptref =                                                                   \
        Gen(name    = 'refs',
            desc    = 'Populates expt_refs mapping table',
            query   = erq,
            funcs   = [erpb],
            actions = [Refs(insert = True,
                            bulks  = erq['b'],
                            atoms  = erq['a'],
                            num    = erpb['out'])])

    ########################################################################
    # Paths
    refpth =  JPath("refs", [refs__bulks])
    apth   =  JPath('atoms',[refs__atoms, refs__bulks])

    # Attributes
    Num    = Refs['num'](refpth)
    Eng    = Atoms['energy'](apth)
    Norm   = SUM(Eng * Num)

    eq = Query(exprs    = dict(b = Bulks.id(),
                               e = Bulks['energies'](),
                               r = Norm),
                aggcols = [Bulks.id()],
                basis   = [Bulks])

    epb = PyBlock(lambda engs,refeng: loads(engs)[0] - float(refeng),
                  args  = [eq[x] for x in 'er'])

    eform =                                                                     \
        Gen(name    = 'eform',
            desc    = 'Difference between a relaxed energy the sum of reference energies ',
            query   = eq,
            funcs   = [epb],
            actions = [Bulks(bulks  = eq['b'],
                             eform  = epb['out'])])

    ########################################################################

    ceq = Query(exprs    = dict(b = Bulks.id(),
                                ce = Bulks['eform']()/Bulks['n_atoms']()))
    ce =                                                                     \
        Gen(name    = 'ce',
            desc    = 'Difference between a relaxed energy (per atom) and '     \
                      'the sum of reference energies (divided by # of atoms)',
            query   = ceq,
            actions = [Bulks(bulks  = ceq['b'],
                             ce     = ceq['ce'])])

    ########################################################################
    bmq = Query(exprs = dict(b = Bulks.id(),
                             v = Bulks['volumes'](),
                             e = Bulks['energies']()))

    doc         = '\nConvert eV/A³ to GPa or GJ/m³'
    bmpb1,bmpb2 = [PyBlock(lambda x: loads(x),args = [bmq[x]]) for x in 've']
    bmpb3 = PyBlock(lambda xs,ys: 2 * xs[0] * np.polyfit(xs,ys,2)[0] * (10**-9) * (1.602 * 10**-19) * (10**10)**3,
                   args = [x['out'] for x in [bmpb1,bmpb2]])
    bm =                                                                        \
        Gen(name    = 'bulkmod',
            desc    = 'Calculates bulk modulus given the 5 minimum energy jobs'+doc,
            funcs   = [bmpb1,bmpb2,bmpb3],
            query   = bmq,
            actions = [Bulks(bulks=bmq['b'],bulkmod=bmpb3['out'])])
    ############################################################################

    all  = ','.join(sorted(set(allspecies)))
    bpth = JPath('bulks',[bulks__job,job__calc])

    dq = Query(exprs   = dict(c = Calc.id(),
                              t = GROUP_CONCAT(Bulks['name'](bpth))),
               basis   = [Calc],
               constr  = NOT(NULL(Bulks['ce'](bpth))),
               aggcols = [Calc.id()])

    dqpb = PyBlock(lambda s,all: ','.join(set(all.split(','))-set(s.split(','))),
                   args = [dq['t'],Const(all)])
    dqpb2 = PyBlock(lambda s: s=='',args=[dqpb['out']])
    done =                                                                      \
        Gen(name    = 'done',
            desc    = 'Determines if calculator is ready for fitting (has all data)',
            query   = dq,
            funcs   = [dqpb,dqpb2],
            actions = [Calc(calc=dq['c'],missing=dqpb['out'],done=dqpb2['out'])])

    ############################################################################


    dq2 = Query(exprs   = dict(c = Calc.id(),
                              t = GROUP_CONCAT(Bulks['name'](bpth))),
               basis   = [Calc],
               constr  = Bulks['success'](bpth),
               aggcols = [Calc.id()])

    dqpb3 = PyBlock(lambda s,all: ','.join(set(all.split(','))-set(s.split(','))),
                   args = [dq['t'],Const(all)])
    done2 =                                                                      \
        Gen(name    = 'done2',
            desc    = 'Determines if calculator is ready for fitting (has all data)',
            query   = dq2,
            funcs   = [dqpb3],
            actions = [Calc(calc=dq2['c'],missing_bulk=dqpb3['out'])])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [isbeef,beef,exptref,eform,ce,bm,done,done2]

    mod.add(gens)
