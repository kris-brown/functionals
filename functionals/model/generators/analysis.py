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
                  RIGHT, AVG, R2, COUNT, CASE, OR, LIKE)

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



from functionals.scripts.fit.a_ce import a_ce
from functionals.scripts.fit.a_bm import a_bm

def analysis(mod : Model) -> None:
    # Extract tables
    tabs = ['Functional','Beef','Bulks','Atoms','Refs','Job','Calc']
    Fx,Beef,Bulks,Atoms,Refs,Job,Calc = map(mod.get,tabs)

    refs__atoms, refs__bulks,job__calc,bulks__job,atoms__job, beef__functional,\
    calc__functional = \
        map(mod.get_rel,[Refs.r('atoms'),Refs.r('bulks'),Job.r('calc'),
                         Bulks.r('job'),Atoms.r('job'),Beef.r('functional'),
                         Calc.r('functional')])
    ########################################################################
    ########################################################################
    ########################################################################
    xvq = Query(dict(b = Bulks.id(),
                     n = RIGHT(Bulks['name'](),Lit(3)),
                     l = Bulks['expt_l']()))

    xvpb=PyBlock(lambda n,l: float(l)**2 * (1.633 if n=='hcp' else float(l)),
                 args = [xvq[x] for x in 'nl'])

    xvol = Gen(name='xvol',
               desc = 'Computes the experimental volume, given a lattice '
                      'constant + crystal structure',
               query = xvq,
               funcs = [xvpb],
               actions = [Bulks(bulks=xvq['b'],expt_vol=xvpb['out'])])
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
    bcalc,acalc = [JPath('calc',[job__calc,x]) for x in [bulks__job,atoms__job]]

    erq = Query(exprs = dict(b = Bulks.id(), a = Atoms.id(),
                             c = Bulks['composition'](),
                             n = Atoms['num']()),
                basis  = [Bulks,Atoms],
                constr = EQ(Calc.id(bcalc),Calc.id(acalc))
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

    eq = Query(exprs    = dict(b = Bulks.id(),
                               e = Bulks['energies'](),
                               r = SUM(Refs['energy'](refpth))),
                aggcols = [Bulks.id()],
                basis   = [Bulks],
                aconstr = SUM(Refs['num'](refpth)) |EQ| Bulks['n_atoms']())

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
                                ce = Lit(-1)*Bulks['eform']()/Bulks['n_atoms']()))
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
    uconv       = (10**-9) * (1.602 * 10**-19) * (10**10)**3
    doc         = '\nConvert eV/A³ to GPa or GJ/m³'
    bmpb1,bmpb2 = [PyBlock(lambda x: loads(x),args = [bmq[x]]) for x in 've']
    bmpb3 = PyBlock(lambda xs,ys,conv: 2 * xs[0] * np.polyfit(xs,ys,2)[0] * conv,
                   args = [bmpb1['out'],bmpb2['out'],Const(uconv)])
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
               #constr  = NOT(NULL(Bulks['ce'](bpth))),
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
                   args = [dq['t'],Const(all)],
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
    beefpth = JPath('beef',[beef__functional,calc__functional,job__calc,bulks__job])
    aceq = Query(exprs=dict(b = Bulks.id(),
                            r = GROUP_CONCAT(Refs['energy'](refpth)),
                            c = GROUP_CONCAT(Refs['contribs'](refpth),delim='$'),
                            x = Bulks['contribs'](),
                            e = Bulks['energies'](),
                            n = Bulks['n_atoms'](),
                            f = MAX(Beef['data'](beefpth))),
                 aggcols = [Bulks.id()],
                 basis   = [Bulks])

    acepb = PyBlock(a_ce,
                    env  = defaultEnv + Env(Import('functools','reduce')),
                    args = [aceq[x] for x in 'rcexnf'],
                    outnames = ['a','b'])

    ace =                                                                       \
        Gen(name    = 'ace',
            desc    = 'Computes A matrix vectors to yield cohesive energy',
            query   =  aceq,
            funcs   = [acepb],
            actions = [Bulks(bulks=aceq['b'],a_ce=acepb['a'],b_ce=acepb['b'])])

    ########################################################################

    abmq = Query(exprs=dict(b = Bulks.id(),
                            x = Bulks['expt_bm'](),
                            c = Bulks['contribs'](),
                            e = Bulks['energies'](),
                            v = Bulks['volumes'](),
                            f = Beef['data'](beefpth)),
                 basis   = [Bulks])

    bmout = ['a_bm','b_bm','a_l','b_l']
    abmpb = PyBlock(a_bm,
                    env      = defaultEnv + Env(Import('numpy.linalg','inv')),
                    args     = [abmq[x] for x in 'evcfx'],
                    outnames = bmout)

    abm =                                                                       \
        Gen(name    = 'abm',
            desc    = 'Computes A matrix vectors to yield cohesive energy',
            query   =  abmq,
            funcs   = [abmpb],
            actions = [Bulks(bulks=aceq['b'],**{x:abmpb[x] for x in bmout})])

    ########################################################################
    diff = Bulks['eosbm']() - Bulks['bulkmod']()

    iq = Query(dict(b=Bulks.id(),i = GT(ABS(diff)/Bulks['expt_bm'](),
                                          Lit(0.1))))
    irreg =                                                                     \
        Gen(name = 'irreg',
            desc = 'Sets Bulks.irregular to TRUE if ASE eos bulkmod differs from discrete by too much',
            query= iq,
            actions=[Bulks(bulks=iq['b'],irregular=iq['i'])])


    ########################################################################
    cols = dict(ce=('ce','expt_ce'),bm=('bulkmod','expt_bm'),lat=('lattice','expt_l'))
    cbp  = JPath('bulks',[bulks__job,job__calc])
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
    ########################################################################
    r2gens = []

    def fx(x:str,y:str)->float:
        x_,y_ = [loads('[%s]'%z) for z in [x,y]]
        #print(sorted(list(zip(y_,x_))))
        out = linregress(x_,y_)[2]**2
        return out

    for k,(x_,y_) in cols.items():
        x = Bulks.get(x_)(cbp); y = Bulks.get(y_)(cbp)
        r2q = Query(exprs   = dict(c = Calc.id(),
                                   x = GROUP_CONCAT(x),
                                   y = GROUP_CONCAT(y)),
                     basis   = [Calc],
                     #constr  = Fx['beef'](JPath('functional',[calc__functional])), #delete this
                     aggcols = [Calc.id()])
        r2 = PyBlock(fx,
            env = Env(Import('scipy.stats','linregress'))+defaultEnv,
            args = [r2q['x'],r2q['y']])

        r2gens.append(Gen(name='r2_'+k,
                           query=r2q, funcs=[r2],
                           actions=[Calc(calc=mseq['c'],
                                          **{'r2_'+k:r2['out']})]))
    r2c,r2b,r2l = r2gens
    ########################################################################
    cd   = dict(Hydride=['H'],III=['Al'],IV=['C','Se','As'],V=['N','V','P'],
                VI=['S','O'],VII=['I','F','Br'])

    case = CASE({OR(*[LIKE(Bulks['name'](),
                           Lit('_%{}%'.format(v))) for v in vs]
                   ) : Lit(k)
                    for k,vs in cd.items()},
                Lit('Metal'))

    aq = Query(dict(b=Bulks.id(),k=case))

    alloy = Gen(
        name    = 'alloy',
        desc    = 'Populates bulks.alloy',
        query   = aq,
        actions = [Bulks(bulks=aq['b'],alloy=aq['k'])]
    )



    ########################################################################
    ########################################################################
    ########################################################################

    gens = [alloy,xvol,isbeef,beef,exptref,eform,ce,bm,done,done2,refcontribs,ace,abm,
            irreg,msec,mseb,msel,r2c,r2b,r2l]

    mod.add(gens)
