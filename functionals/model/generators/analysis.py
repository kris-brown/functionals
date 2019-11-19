# External
from typing import Tuple as T
from json   import loads,dumps
from ast    import literal_eval
import numpy as np # type: ignore
from scipy.stats import linregress         # type: ignore

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, AND, Env, defaultEnv, Import,
                  Literal as Lit,  EQ, JPath, Constraint, LEFT, One, SUM, GT,
                  LIKE, CONCAT, LT, ABS, GROUP_CONCAT, NOT, NULL, Const, MAX,
                  RIGHT, AVG, R2, COUNT, CASE, LIKE)

from functionals.scripts.load.true_mag import true_mag

################################################################################
################################################################################
################################################################################
# Constants
allspecies = \
    ['Li_bcc','Na_bcc','K_bcc','Rb_bcc', 'Ca_fcc','Sr_fcc','Ba_bcc',
       'Nb_bcc','Ta_bcc', 'Mo_bcc','W_bcc','Fe_bcc',
       'Rh_fcc','Ir_fcc','Ni_fcc','Pd_fcc','Pt_fcc',
       'Cu_fcc','Ag_fcc','Au_fcc',  'Al_fcc','Pb_fcc',
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
    tabs = ['Functional','Bulks','Atoms','Refs','Calc']
    Fx,Bulks,Atoms,Refs,Calc = map(mod.get,tabs)

    refs__atoms, refs__bulks,bulks__calc,atoms__calc,\
    calc__functional = \
        map(mod.get_rel,[Refs.r('atoms'),Refs.r('bulks'),Bulks.r('calc'),
                          Atoms.r('calc'), Calc.r('functional')])
    ########################################################################
    ########################################################################
    ########################################################################
    tmagq = Query(dict(i=Atoms.id(),n=Atoms['name']()))
    tmagpb = PyBlock(true_mag,args=[tmagq['n']])
    tmag = Gen(name='truemag',
               desc = 'populates atoms.true_mag',
               query=tmagq, funcs= [tmagpb], actions=[Atoms(atoms=tmagq['i'],true_mag=tmagpb['out'])])
    ########################################################################
    # xvq = Query(dict(b = Bulks.id(),
    #                  n = RIGHT(Bulks['name'](),Lit(3)),
    #                  l = Bulks['expt_l']()))
    #
    # xvpb=PyBlock(lambda n,l: float(l)**3 * (1.4142 if n=='hcp' else 1.),
    #              args = [xvq[x] for x in 'nl'])
    #
    # xvol = Gen(name='xvol',
    #            desc = 'Computes the experimental volume, given a lattice '
    #                   'constant + crystal structure',
    #            query = xvq,
    #            funcs = [xvpb],
    #            actions = [Bulks(bulks=xvq['b'],expt_vol=xvpb['out'])])
    ########################################################################
    ibq = Query(exprs = dict(f = Fx.id(),
                             b = EQ(LEFT(Fx['data'](), One), Lit('['))))
    isbeef =                                                                    \
        Gen(name    = 'isbeef',
            desc    = 'Determines whether fx was BEEF style by checking if functional is a word or an array',
            query   = ibq,
            actions = [Fx(functional=ibq['f'],beef=ibq['b'])])

    ############################################################################
    name  = CONCAT([Lit('%,'),Atoms['num'](),Lit(',%')])

    bcalc,acalc = [JPath('calc',[x]) for x in [bulks__calc,atoms__calc]]

    erq = Query(exprs = dict(b = Bulks.id(), a = Atoms.id(),
                             c = Bulks['composition'](),
                             n = Atoms['num']()),
                basis  = [Bulks,Atoms],
                constr = AND([EQ(Calc.id(bcalc),Calc.id(acalc)),
                              LIKE(Bulks['elems'](),name),
                              Atoms['int_occupation'](),
                              LT(ABS(Atoms['true_mag']()-Atoms['mag']()),Lit(0.05))]))

    erpb = PyBlock(lambda c,n: literal_eval(c).get(n,[]),
                   env  = Env([Import('ast',['literal_eval'])]),
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

    eq = Query(exprs    = dict(b = Bulks.id(),
                               e = Bulks['eng']() - SUM(Refs['energy'](refpth))),
                aggcols = [Bulks.id()],
                basis   = [Bulks],
                aconstr = EQ(SUM(Refs['num'](refpth)), Bulks['n_atoms']()))


    eform =                                                                     \
        Gen(name    = 'eform',
            desc    = 'Difference between a relaxed energy the sum of reference energies ',
            query   = eq, actions = [Bulks(bulks  = eq['b'], eform  = eq['e'])])

    ########################################################################

    ceq = Query(exprs    = dict(b = Bulks.id(),
                                ce = Lit(-1)*Bulks['eform']()/Bulks['n_atoms']()))
    ce =                                                                     \
        Gen(name    = 'ce',
            desc    = 'Difference between a relaxed energy (per atom) and '     \
                      'the sum of reference energies (divided by # of atoms)',
            query   = ceq,
            actions = [Bulks(bulks  = ceq['b'],
                             ce     = ceq['ce'])])

    ########################################################################

    all  = ','.join(sorted(set(allspecies)))
    bpth = JPath('bulks',[bulks__calc])

    dq = Query(exprs   = dict(c = Calc.id(),
                              t = GROUP_CONCAT(Bulks['name'](bpth))),
               basis   = [Calc],
               constr  = AND([NOT(NULL(Bulks['ce'](bpth))), NOT(Bulks['irregular'](bpth))]),
               aggcols = [Calc.id()])

    dqpb = PyBlock(lambda s,all: ','.join(set(all.split(','))-set(s.split(','))),
                   args = [dq['t'],Const(all)])

    dqpb2 = PyBlock(lambda s: s=='',args=[dqpb['out']])

    done =                                                                      \
        Gen(name    = 'done',
            desc    = 'Determines if calculator is ready for fitting (has all data)',
            query   = dq,
            funcs   = [dqpb,dqpb2],
            actions = [Calc(calc=dq['c'], missing = dqpb['out'],done=dqpb2['out'])])

    ############################################################################
    def missbulk(s:str,all:str)->T[int,str]:
        miss = set(all.split(','))-set(s.split(','))
        return len(miss),','.join(miss)

    dq2 = Query(exprs   = dict(c = Calc.id(),
                              t = GROUP_CONCAT(Bulks['name'](bpth))),
               basis   = [Calc],
               constr  = Bulks['success'](bpth),
               aggcols = [Calc.id()])

    dqpb3 = PyBlock(missbulk,
                   args = [dq2['t'],Const(all)],
                   outnames = ['n','s'])
    done2 =                                                                      \
        Gen(name    = 'done2',
            desc    = 'Determines if calculator is ready for fitting (has all data)',
            query   = dq2,
            funcs   = [dqpb3],
            actions = [Calc(calc         = dq2['c'],
                            n_missing    = dqpb3['n'],
                            missing_bulk = dqpb3['s'])])

    ########################################################################
    refatom = JPath('atoms', [refs__atoms])

    rcq = Query(exprs = dict(r = Refs.id(),
                             c = Atoms['contribs'](refatom),
                             n = Refs['num'](),
                             e = Refs['num']()*Atoms['energy'](refatom)),
                basis  = [Refs],
                opt_attr = [Atoms['contribs'](refatom)])

    rcpb = PyBlock(lambda c,n: dumps((n*np.array(loads(c))).tolist()) if c else None,
                   args = [rcq[x] for x in 'cn'])

    refcontribs =                                                               \
        Gen(name    = 'refcontribs',
            desc    = 'Takes atomic contribs and multiplies by stoichiometry',
            query   = rcq,
            funcs   = [rcpb],
            actions = [Refs(refs     = rcq['r'],
                            energy   = rcq['e'],
                            contribs = rcpb['out'])])


    ########################################################################
    cols = dict(ce=('ce','expt_ce'),bm=('bulkmod','expt_bm'),lat=('lattice','expt_l'))
    cbp  = JPath('bulks',[bulks__calc])
    msegens = []
    for k,(x_,y_) in cols.items():
        x = Bulks.get(x_)(cbp); y = Bulks.get(y_)(cbp)
        mseq = Query(exprs   = dict(c = Calc.id(), m = AVG((x-y)**Lit(2))),
                     basis   = [Calc],
                     #constr  = Calc['done'](),
                     aggcols = [Calc.id()])
        msegens.append(Gen(name='mse_'+k,
                           query=mseq,
                           actions=[Calc(calc=mseq['c'],
                                          **{'mse_'+k:mseq['m']})]))
    msec,mseb,msel = msegens

    ############################################################################
    aqcols = ['pw','fx']
    acq = Query(dict(a=Atoms.id(), **{x:Atoms[y]() for x,y in zip('pf',aqcols)}))
    icalc = Calc(insert=True, pw=acq['p'],
                 functional=Fx(insert=True, data=acq['f']))
    atom_calc =                                                                 \
        Gen(name = 'atom_calc',
            desc = 'populate atoms.calc',
            query = acq,
            actions = [Atoms(atoms=acq['a'],calc=icalc)]
            )
    ############################################################################

    annq = Query(dict(a=Atoms.id(), s = Atoms['stordir']()))

    def annfunc(pth:str) -> T[str, int]:
        from ase.data import chemical_symbols # type: ignore
        symb = pth.split('/')[-1]
        return symb, chemical_symbols.index(symb)

    annpb = PyBlock(annfunc, args=[annq['s']], outnames = ['n','m'])
    atom_namenum =                                                              \
        Gen(name = 'atom_namenum',
            desc = 'populate atoms.name and atoms.num',
            query = annq,
            funcs = [annpb],
            actions = [Atoms(atoms=annq['a'], name=annpb['n'], num=annpb['m'])])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [isbeef,exptref,eform,ce,done,done2,refcontribs,
            msec,mseb,msel, atom_namenum, atom_calc, tmag]

    mod.add(gens)
