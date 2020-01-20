# External
from typing import Tuple as T
from json import loads, dumps
from ast import literal_eval
import numpy as np

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, AND, Env, Import,
                   Literal as Lit,  EQ, JPath, LEFT, One, SUM,
                   CONCAT, AVG, LIKE, GROUP_CONCAT, ABS)

from functionals.scripts.load.true_mag import true_mag

#########################################################################
#########################################################################
#########################################################################


def analysis(mod: Model) -> None:
    # Extract tables
    tabs = ['Functional', 'Bulks', 'Atoms', 'Refs', 'Calc', 'Surf']
    Fx, Bulks, Atoms, Refs, Calc, Surf = map(mod.get, tabs)

    refs__atoms, refs__bulks, bulks__calc, atoms__calc, calc__functional = \
        map(mod.get_rel, [Refs.r('atoms'), Refs.r('bulks'), Bulks.r('calc'),
                          Atoms.r('calc'), Calc.r('functional')])
    ########################################################################
    ########################################################################

    def aefun(intocc: bool, true_mag: float, mag: float, fx: str) -> str:
        if fx is None:
            return 'unconverged'
        elif not intocc:
            return 'non-integer magmom'
        elif (true_mag - mag) > 0.05:
            return 'wrong magmom'
        else:
            return ''
    int_occupation, tm, mag, fx = [Atoms[x]() for x in [
        'int_occupation', 'true_mag', 'mag', 'fx'
    ]]
    aeq = Query(dict(a=Atoms.id(), i=int_occupation,
                     t=tm, m=mag, f=fx),
                opt_attr=[int_occupation, tm, mag, fx])
    aepb = PyBlock(aefun, args=[aeq[x] for x in 'itmf'])
    aerror = Gen(
        name='atomerror',
        desc='detects anything wrong with atomic job',
        query=aeq, funcs=[aepb],
        actions=[Atoms(atoms=aeq['a'], error=aepb['out'])]
    )
    ########################################################################
    tmagq = Query(dict(i=Atoms.id(), n=Atoms['name']()))
    tmagpb = PyBlock(true_mag, args=[tmagq['n']])
    tmag = Gen(name='truemag',
               desc='populates atoms.true_mag',
               query=tmagq, funcs=[tmagpb],
               actions=[Atoms(atoms=tmagq['i'], true_mag=tmagpb['out'])])

    ########################################################################
    xvq = Query(dict(b=Bulks.id(),
                     v=(Bulks['expt_lat']()/Bulks['volrat']())**Lit(3)))
    xvol = Gen('xvol', query=xvq, actions=[Bulks(bulks=xvq['b'],
                                                 expt_vol=xvq['v'])])
    ########################################################################
    cnq = Query(dict(b=Bulks.id(), n=Bulks['stordir']()))
    cnpb = PyBlock(lambda x: x.split('/')[-2], args=[cnq['n']])
    calcname = Gen('calcname', query=cnq, funcs=[cnpb],
                   actions=[Bulks(bulks=cnq['b'], calcname=cnpb['out'])])
    ########################################################################
    ibq = Query(exprs=dict(f=Fx.id(),
                           b=EQ(LEFT(Fx['data'](), One), Lit('['))))
    isbeef =                                                                  \
        Gen(name='isbeef',
            desc='Check whether fx was BEEF - if fx is a word or array',
            query=ibq,
            actions=[Fx(functional=ibq['f'], beef=ibq['b'])])

    #####################################################################
    name = CONCAT([Lit('%,'), Atoms['num'](), Lit(',%')])

    bcalc, acalc = [JPath('calc', [x]) for x in [bulks__calc, atoms__calc]]

    erq = Query(exprs=dict(b=Bulks.id(), a=Atoms.id(),
                           c=Bulks['composition'](),
                           n=Atoms['num']()),
                basis=[Bulks, Atoms],
                constr=AND([
                    EQ(Calc.id(bcalc), Calc.id(acalc)),
                    LIKE(Bulks['elems'](), name),
                    EQ(Atoms['error'](), Lit(''))]))

    erpb = PyBlock(lambda c, n: literal_eval(c).get(n, []),
                   env=Env([Import('ast', ['literal_eval'])]),
                   args=[erq['c'], erq['n']])

    exptref =                                                                \
        Gen(name='refs',
            desc='Populates expt_refs mapping table',
            query=erq, funcs=[erpb],
            actions=[Refs(insert=True,
                          bulks=erq['b'],
                          atoms=erq['a'],
                          num=erpb['out'])])

    ########################################################################
    # Paths
    refpth = JPath("refs", [refs__bulks])

    eq = Query(exprs=dict(b=Bulks.id(),
                          e=Bulks['energy']() - SUM(Refs['energy'](refpth))),
               aggcols=[Bulks.id()],
               basis=[Bulks],
               aconstr=EQ(SUM(Refs['num'](refpth)), Bulks['n_atoms']()))

    eform =                                                                  \
        Gen(name='eform',
            desc='Diff between a relaxed energy the sum of reference engs ',
            query=eq, actions=[Bulks(bulks=eq['b'], eform=eq['e'])])

    ########################################################################

    ceq = Query(exprs=dict(b=Bulks.id(),
                           ce=Lit(-1)*Bulks['eform']()/Bulks['n_atoms']()))
    ce =                                                                     \
        Gen(name='ce',
            desc='Difference between a relaxed energy (per atom) and '
                 'the sum of reference energies (divided by # of atoms)',
            query=ceq,
            actions=[Bulks(bulks=ceq['b'], ce=ceq['ce'])])

    ########################################################################

    bpth = JPath('bulks', [bulks__calc])
    com = Lit(',')
    bname, bce, bbm = [Bulks[x](bpth) for x in ['name', 'ce', 'bm']]
    catargs = [bname, com, bce, com, bbm]
    dq = Query(exprs=dict(c=Calc.id(), w=Calc['allmat'](),
                          h=GROUP_CONCAT(CONCAT(catargs), delim=' ')),
               basis=[Calc], aggcols=[Calc.id()], opt_attr=[bce, bbm])

    def dfun(have_: str, want_: str) -> T[str, int, bool]:
        wce, wbml = [set(x.split()) for x in want_.split('|')]

        for name, ce, bm in [x.split(',') for x in have_.split()]:
            if ce and name in wce:
                wce.remove(name)
            if bm and name in wbml:
                wbml.remove(name)
        missing = [x+'_ce' for x in wce] + [x+'_bml' for x in wbml]
        n = len(missing)
        return ' '.join(sorted(missing)), n, n == 0

    dcols = ['missing', 'n_missing', 'done']
    dqpb = PyBlock(dfun, args=[dq['h'], dq['w']], outnames=dcols)

    done =                                                                  \
        Gen(name='done',
            desc='Determines if calculator is ready for fitting (has data)',
            query=dq, funcs=[dqpb],
            actions=[Calc(calc=dq['c'], **{x: dqpb[x] for x in dcols})])

    ########################################################################
    refatom = JPath('atoms', [refs__atoms])

    rcq = Query(exprs=dict(r=Refs.id(),
                           c=Atoms['contribs'](refatom),
                           n=Refs['num'](),
                           e=Refs['num']()*Atoms['energy'](refatom)),
                basis=[Refs],
                opt_attr=[Atoms['contribs'](refatom)])

    rcpb = PyBlock(
        lambda c, n: dumps((n*np.array(loads(c))).tolist()) if c else None,
        args=[rcq[x] for x in 'cn'])

    refcontribs =                                                   \
        Gen(name='refcontribs',
            desc='Takes atomic contribs and multiplies by stoichiometry',
            query=rcq,
            funcs=[rcpb],
            actions=[Refs(refs=rcq['r'],
                          energy=rcq['e'],
                          contribs=rcpb['out'])])

    ########################################################################
    cols = dict(ce=('ce', 'expt_ce'), bm=('bm', 'expt_bm'),
                lat=('lat', 'expt_lat'),  mag=('mag', 'expt_mag'))
    cbp = JPath('bulks', [bulks__calc])
    msegens = []
    for k, (x_, y_) in cols.items():
        x = Bulks.get(x_)(cbp)
        y = Bulks.get(y_)(cbp)
        mseq = Query(exprs=dict(c=Calc.id(), m=AVG(ABS(x-y))),
                     basis=[Calc],
                     # constr=Calc['done'](),
                     aggcols=[Calc.id()])
        msegens.append(Gen(name='mae_'+k,
                           query=mseq,
                           actions=[Calc(calc=mseq['c'],
                                         **{'mae_'+k: mseq['m']})]))
    msec, mseb, msel, msem = msegens

    #####################################################################
    aqcols = ['pw', 'fx']
    z = zip('pf', aqcols)
    acq = Query(dict(a=Atoms.id(), **{x: Atoms[y]() for x, y in z}))
    icalc = Calc(insert=True, pw=acq['p'],
                 functional=Fx(insert=True, data=acq['f']))
    atom_calc =                                                          \
        Gen(name='atom_calc',
            desc='populate atoms.calc',
            query=acq,
            actions=[Atoms(atoms=acq['a'], calc=icalc)])
    #####################################################################

    annq = Query(dict(a=Atoms.id(), s=Atoms['stordir']()))

    def annfunc(pth: str) -> T[str, int]:
        from ase.data import chemical_symbols
        symb = pth.split('/')[-1]
        return symb, chemical_symbols.index(symb)

    annpb = PyBlock(annfunc, args=[annq['s']], outnames=['n', 'm'])
    atom_namenum =                                                       \
        Gen(name='atom_namenum',
            desc='populate atoms.name and atoms.num',
            query=annq,
            funcs=[annpb],
            actions=[Atoms(atoms=annq['a'], name=annpb['n'], num=annpb['m'])])

    ########################################################################
    cnp = JPath('functional', [calc__functional])

    cnq = Query(dict(c=Calc.id(), d=Fx['data'](cnp)), basis=[Calc])

    def cnfunc(data: str) -> str:
        import os
        import json
        root = '/Users/ksb/functionals/data/beefs'
        if data[0] == '[':
            data_ = json.loads(data)
            for beef in os.listdir(root):
                with open(os.path.join(root, beef), 'r') as f:
                    if json.load(f) == data_:
                        return beef
            import pdb
            pdb.set_trace()
            assert False
            return ''
        else:
            return data

    cnpb = PyBlock(cnfunc, args=[cnq['d']])
    calcname2 =                                                           \
        Gen(name='calcname2',
            desc='Figures out the name of a BEEF calculator from its coefs',
            query=cnq, funcs=[cnpb],
            actions=[Calc(calc=cnq['c'], name=cnpb['out'])])
    ########################################################################
    seq = Query(dict(s=Surf.id(), m=Surf['mat'](
    ), o=Surf['ontop'](), h=Surf['hollow']()))

    def sef(m: str, o: float, h: float) -> float:
        if m == 'Pd':
            return h-o
        else:
            return o-h
    sepb = PyBlock(sef, args=[seq[x] for x in 'moh'])
    surfan = \
        Gen(name='surferr',
            query=seq, funcs=[sepb], tags=['surf'],
            actions=[Surf(surf=seq['s'], err=sepb['out'])])
    ########################################################################
    ########################################################################
    ########################################################################

    gens = [isbeef, exptref, eform, ce, refcontribs, xvol, msem,
            calcname, msec, mseb, msel, atom_namenum, atom_calc,
            tmag, calcname2, done, aerror, surfan]

    mod.add(gens)
