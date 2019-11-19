# External Modules
from typing import Any, Type

# Internal Modules
from dbgen import (Model, Obj, Attr, Rel, Int, Varchar, Text, Decimal, Boolean,
                   Date, PathEQ, Path, RawView)

#################
# INPUT SUMMARY #
#################

################################################################################
mets = ['ce','bm','lat']
funmetrics = ['mse_%s'%(x) for x in mets]

calc = Obj(
    name  = 'calc',
    desc  = 'Calculator details',
    attrs = [Attr('pw',                identifying= True, desc = 'Planewave cutoff, eV'),
             Attr('done',   Boolean(),            desc = 'Whether or not all calculations have been finished'),
             Attr('missing',Text(),               desc = 'Materials not yet finished (bulk or atom related)'),
             Attr('n_missing',Int(),              desc = 'Materials not yet finished (bulk or atom related)'),
             Attr('missing_bulk',Text(),          desc = 'Materials for which not enough bulk calculations have finished'),
             Attr('fitdata',Text(),               desc = 'For a BEEF calc, the input data to fitting')]
           +[Attr(m, Decimal(), desc = 'Only fills out if results in for all materials')
                for m in funmetrics],
    fks  = [Rel('functional',identifying=True)])

################################################################################

functional = Obj(
    name  = 'functional',
    desc  = 'DFT functional',
    attrs = [Attr('data',Text(),identifying=True), Attr('beef',Boolean())])

################################################################################

#########
# JOBS  #
#########

################################################################################

atoms = Obj(
    name  = 'atoms',
    desc  = 'DFT calculation',
    attrs = [Attr('stordir',       Varchar(), identifying=True),
             Attr('runtime',       Int(),      desc='Appx runtime, hours'),
             Attr('pw',                        desc = 'Planewave cutoff, eV'),
             Attr('econv',   Decimal(),            desc = 'EDIFF in the INCAR'),
             Attr('name',          Varchar()),
             Attr('fx',            Text()),
             Attr('num',           Int(),      desc = 'Atomic number'),
             Attr('true_mag',      Int(),      desc = 'Magnetic moment'),
             Attr('mag',           Decimal(),  desc = 'Magnetic moment'),
             Attr('int_occupation',Boolean(),  desc = 'Whether or not all orbitals have integer electron occupation'),
             Attr('energy',        Decimal(),  desc = 'eV'),
             Attr('contribs',      Text(),     desc = 'list of five 64-element exchange contribs, if beef calculation')],
    fks  = [Rel('calc')]
    )

################################################################################
bulkjob = Obj(
    name ='bulkjob',
    desc = 'Individual singlepoints for a bulk.',
    attrs = [Attr('stordir',        Varchar(), identifying=True),
             Attr('parent',         Varchar()),
             Attr('fx',             Text()),
             Attr('contribs',       Text()),
             Attr('int_occupation', Boolean(),  desc = 'Whether or not all orbitals have integer electron occupation'),
             Attr('energy',         Decimal(),  desc = 'Energy'),
             Attr('mag',            Decimal(),  desc = 'magnetic moment'),
             Attr('runtime',        Int(),      desc='Appx runtime, hours'),
             Attr('pw',                         desc = 'Planewave cutoff, eV'),
             Attr('econv',   Decimal(),            desc = 'EDIFF in the INCAR'),
             Attr('strain',         Int(),      desc='Relative strain, percent.'),
             Attr('volume',         Decimal(),  desc = 'A^3'),
             Attr('lattice',        Decimal(),  desc = 'A'),
             ])

bulks = Obj(
   name  = 'bulks',
   desc  = 'Set of singlepoints on a particular material under strain',
   attrs = [Attr('stordir',     Varchar(), identifying=True),
            Attr('name',       Varchar(), desc = 'Species nickname'),
            Attr('composition',Varchar(), desc = 'Species composition'),
            Attr('elems',Varchar(), desc = 'common separated list of atomic nums'),
            Attr('n_atoms',               desc = 'Number of atoms in unit cell'),
            Attr('n_elems',               desc = 'Number of distinct chemical species'),
            Attr('alloy',      Varchar(), desc = 'What type of alloy it is, if any'),
            # Figuring out if job is complete
            Attr('allvol',      Text('long'),   desc = 'volumes  of all jobs'),
            Attr('alleng',      Text('long'),   desc = 'energies of all jobs'),
            Attr('strains',    Text(),    desc = 'List of strains that have completed'),
            Attr('morejobs',   Int(),     desc = '-1 = needs more jobs at lower strain, +1 higher strain'),
            Attr('success',    Boolean(), desc = 'We have enough calculations to do an EOS'),

            # Analysis Results from minimum 5 jobs
            Attr('volumes',     Text('long'),   desc = 'volumes of the 5 most optimal jobs'),
            Attr('energies',    Text('long'),   desc = 'energies of the 5 most optimal jobs'),
            Attr('contribs',    Text('long'),   desc = 'xc contribs of 5 most optimal jobs'),
            Attr('eng',         Decimal(),      desc = 'Minimum energy according to quadratic fit'),
            Attr('mag',         Decimal(),      desc = 'Magmom of the minimum energy job'),
            Attr('bulkmod',     Decimal(),      desc = 'Bulk modulus, GPa, from 5 optimal jobs + stencil'),
            Attr('bulkmod_lstsq', Decimal(),    desc = 'Bulk modulus, GPa, from quadratic fit of 5 optimal jobs'),
            Attr('lstsq_abc',   Text(),         desc = 'JSON of quadratic fit of 5 optimal jobs'),
            Attr('lattice',     Decimal(),      desc = 'Conventional unit cell lattice from quad fit'),
            Attr('volrat',      Decimal(),      desc = 'lattice/vol^1/3'),
            Attr('volume',      Decimal(),      desc = 'Volume from quad fit'),
            Attr('eform',       Decimal(),      desc = 'Formation energy, eV, given ENG'),
            Attr('ce',          Decimal(),      desc = 'Cohesive  energy, eV'),
            # Analysis from ase.eos
            Attr('eosbm',       Decimal(),      desc='Bulk modulus, GPa, from ASE EOS'),
            Attr('eos_diff',   Decimal(),      desc='Ratio of discrete BM to EOS BM'),
            Attr('irregular',   Boolean(),      desc='Whether discrete BM differs significantly from EOS BM'),

            # Experimental results
            Attr('expt_ce',    Decimal(15,3),   desc = 'Experimental cohesive energy, eV'),
            Attr('expt_bm',    Decimal(15,3),   desc = 'Experimental bulk modulus, GPa'),
            Attr('expt_l',     Decimal(15,6),   desc = 'Experimental lattice parameter of reference stoichiometry, A'),
            Attr('expt_vol',   Decimal(),       desc = 'Experimental volume of reference stoichiometry, A^3'),
            Attr('expt_mag',   Decimal(),       desc = 'Experimental magnetic moment, bohr'),

            # Fitting inputs
            Attr('ab_ce', Text(), desc = 'a vectors, when dotted w/ BEEF coef + offset gives formation energy in eV'),
            Attr('ab_bm', Text(), desc = 'a vectors, when dotted w/ BEEF coef + offset gives bulk modulus in GPa'),
            Attr('ab_vol',  Text(), desc = 'a vectors, when dotted w/ BEEF coef + offset gives conventional lattice constant in A'),
            ],
    fks = [Rel('calc')])


################################################################################

expt_refs = Obj(
    name  = 'refs',
    desc  = 'Mapping table between bulk and relevant atomic calculations',
    attrs = [Attr('num',     Int(),    desc='Stoichiometric coefficient'),
             Attr('energy',  Decimal(),desc='Energy of atom, weighted by stoichiometry'),
             Attr('contribs',Text(),   desc='Exchange contributions of atom, weighted by stoichiometry')],
    fks   = [Rel('bulks', identifying=True),
             Rel('atoms', identifying=True)])

################################################################################

###########
# FITTING #
###########

################################################################################
fitparams = Obj(
    name = 'fitparams',
    desc = 'Input parameters to a fit',
    attrs = [Attr('ce_scale', Decimal(),  identifying= True, desc = 'One cost unit'),
             Attr('bm_scale', Decimal(),  identifying= True, desc = 'One cost unit'),
             Attr('lc_scale', Decimal(),  identifying= True, desc = 'One cost unit'),
             Attr('mag_scale', Decimal(),  identifying= True, desc = 'Extra weight on magnetic material LC'),
             Attr('consts',   Text(),     identifying= True, desc = 'Sorted string list of linear constraints to be included'),
             Attr('reg',      Decimal(),  identifying= True, desc = 'Regularization penalty') ])

fit = Obj(
    name  = 'fit',
    desc  = 'A constrained fit to some subset of cohesive/BM/lattice data',
    attrs = [Attr('name',Varchar(),identifying=True),
             Attr('pth',Varchar(),   desc='Location of fit job'),
             Attr('done',            desc='0 = not done at all, 1 = main fit done, 2 = cv done'),
             Attr('cvioldict',Text(),desc='Dict of constraints and their violation at "optimum"'),
             # Analysis of results
             Attr('plt', Boolean(),desc='Whether or not plots were made'),]
            +[Attr(m, Decimal(), desc = 'Only fills out if results in for all materials')
                 for m in funmetrics]
            +[Attr('err_'+x,Decimal()) for x in mets],
            #+[Attr('r2_'+x,Decimal()) for x in mets],
    fks = [Rel('fitparams', identifying= True), Rel('calc')])

################################################################################

# Views #

err = RawView('errs','''
SELECT CASE WHEN beef THEN 'beef' ELSE data END,
       mse_ce,mse_bm,mse_lat
    FROM calc JOIN functional ON functional=functional_id
UNION
SELECT name,mse_ce,mse_bm,mse_lat from fit
''')

################################################################################
################################################################################
################################################################################

objs = [calc, functional, bulkjob, atoms, bulks, fitparams, fit, expt_refs, ]


def new_model() -> Model:
    m = Model('functionals',objlist=objs, viewlist=[err])
    return m
