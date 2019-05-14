# External Modules
from typing import Any, Type

# Internal Modules
from dbgen import Model, Obj, Attr, Rel, Int, Varchar, Text, Decimal, Boolean, Date, PathEQ, Path



#################
# INPUT SUMMARY #
#################

################################################################################

funmetrics = ['%s_%s'%(x,y) for x in ['mse','r2'] for y in ['ce','bm','lat']]

calc = Obj(
    name  = 'calc',
    desc  = 'Calculator details',
    attrs = [Attr('pw',                id = True, desc = 'Planewave cutoff, eV'),
             Attr('econv',  Decimal(), id = True, desc = 'Energy convergance criterion'),
             Attr('done',   Boolean(),            desc = 'Whether or not all calculations have been finished'),
             Attr('missing',Text(),               desc = 'Materials not yet finished (bulk or atom related)'),
             Attr('n_missing',Int(),              desc = 'Materials not yet finished (bulk or atom related)'),
             Attr('missing_bulk',Text(),          desc = 'Materials for which not enough bulk calculations have finished')]
           +[Attr(m, Decimal(), desc = 'Only fills out if results in for all materials')
                for m in funmetrics])

calc_rels = [Rel('functional','calc',id=True)]

################################################################################

functional = Obj(
    name  = 'functional',
    desc  = 'DFT functional',
    attrs = [Attr('data',Text(),id=True), Attr('beef',Boolean())])

################################################################################
ams = ['a1'+x for x in '12345']+ ['msb']

beeffunc = Obj(
    name  = 'beef',
    desc  = 'BEEF functional',
    attrs = [Attr('data',Text(),id=True)]
           +[Attr(x, Decimal(), id=True) for x in ams])

bf_rels = [Rel('functional','beef',id=True)]

################################################################################

#########
# JOBS  #
#########

################################################################################

job = Obj(
    name = 'job',
    desc = 'either an atom singlepoint or a series of strained bulks',
    attrs = [Attr('stordir',     Varchar(), id=True),
             Attr('runtime',     Int(),desc='Appx runtime, hours')])


job_rels = [Rel('calc',  'job'),]

################################################################################

atoms = Obj(
    name = 'atoms',
    desc  = 'DFT calculation',
    attrs = [Attr('name',          Varchar()),
             Attr('num',           Int(),      desc = 'Atomic number'),
             Attr('true_mag',      Int(),      desc = 'Magnetic moment'),
             Attr('mag',           Decimal(),  desc = 'Magnetic moment'),
             Attr('int_occupation',Boolean(),  desc = 'Whether or not all orbitals have integer electron occupation'),
             Attr('energy',        Decimal(),  desc = 'eV'),
             Attr('contribs',      Text(),     desc = 'list of five 64-element exchange contribs, if beef calculation')])

atoms_rels = [Rel('job','atoms',id=True)] # Rel('potcar','atoms'),

################################################################################

bulks = Obj(
   name  = 'bulks',
   desc  = 'Set of singlepoints on a particular material under strain',
   attrs = [Attr('name',       Varchar(), desc = 'Species nickname'),
            Attr('n_atoms',               desc = 'Number of atoms in unit cell'),
            Attr('n_elems',               desc = 'Number of distinct chemical species'),
            Attr('alloy',      Varchar(), desc = 'What type of alloy it is, if any'),
            # Figuring out if job is complete
            Attr('strain_low', Int(),     desc = 'Lowest strain with a directory for calculation'),
            Attr('strain_hi',  Int(),     desc = 'Highest strain with a directory for calculation'),
            Attr('incomplete', Text(),    desc = 'List of strains that have completed'),
            Attr('morejobs',   Int(),     desc = '-1 = needs more jobs at lower strain, +1 more at higher strain'),
            Attr('success',    Boolean(), desc = 'We have enough calculations to do an EOS'),

            # Properties when enough jobs have completed
            Attr('volumes',     Text('long'),   desc = 'volumes of the 5 most optimal jobs'),
            Attr('energies',    Text('long'),   desc = 'energies of the 5 most optimal jobs'),
            Attr('contribs',    Text('long'),   desc = 'xc contribs of 5 most optimal jobs'),
            Attr('img',         Text('long'),   desc = 'base64 encoded image of EOS fit'),

            Attr('composition',     Varchar(),   desc='Species composition'),

            # Analysis Results from minimum 5 jobs
            Attr('eform',       Decimal(),      desc = 'Formation energy, eV'),
            Attr('ce',          Decimal(),      desc = 'Cohesive  energy, eV'),
            Attr('bulkmod',     Decimal(),      desc = 'Bulk modulus, GPa'),
            Attr('lattice',     Decimal(),      desc = 'Conventional unit cell lattice (optimized)'),
            # Analysis from ase.eos
            Attr('eosbm',       Decimal(),      desc='Bulk modulus, GPa'),
            Attr('irregular',   Boolean(),      desc='Whether discrete BM differs significantly from EOS BM'),

            # Experimental results
            Attr('expt_ce',    Decimal(15,3),   desc = 'Experimental cohesive energy, eV'),
            Attr('expt_bm',    Decimal(15,3),   desc = 'Experimental bulk modulus, GPa'),
            Attr('expt_l',     Decimal(15,6),   desc = 'Experimental lattice parameter of reference stoichiometry, A'),
            Attr('expt_vol',     Decimal(),   desc = 'Experimental volume of reference stoichiometry, A^3'),
            Attr('expt_mag',   Decimal(),   desc = 'Experimental magnetic moment, bohr'),

            # Fitting inputs
            Attr('a_ce', Text(), desc = '5 vectors (for each a1 value), when dotted w/ BEEF coef + offset gives formation energy in eV'),
            Attr('a_bm', Text(), desc = '5 vectors (for each a1 value), when dotted w/ BEEF coef + offset gives bulk modulus in GPa'),
            Attr('a_l',  Text(), desc = '5 vectors (for each a1 value), when dotted w/ BEEF coef + offset gives conventional lattice constant in A'),

            Attr('b_ce', Text(), desc = '5 offsets (for each a1 value), eV'),
            Attr('b_bm', Text(), desc = '5 offsets (for each a1 value), GPa'),
            Attr('b_l',  Text(), desc = '5 offsets (for each a1 value), A'),

            ])

bulks_rels = [Rel('job','bulks',id=True)]

################################################################################

################################################################################

expt_refs = Obj(
    name = 'refs',
    desc = 'Mapping table between bulk and relevant atomic calculations',
    attrs = [Attr('num',     Int(), desc='Stoichiometric coefficient'),
             Attr('energy',     Decimal(),desc='Energy of atom, weighted by stoichiometry'),
             Attr('contribs',Text(),desc='Exchange contributions of atom, weighted by stoichiometry')])

ea_rels = [Rel('bulks', 'refs', id=True),
           Rel('atoms', 'refs', id=True)]

################################################################################

###########
# FITTING #
###########

################################################################################
fitparams = Obj(
    name = 'fitparams',
    desc = 'Input parameters to a fit',
    attrs = [Attr('bm_weight',    Decimal(),  id = True, desc = 'Relative weight of bulk-modulus data to cohesive energy'),
             Attr('lat_weight',   Decimal(),  id = True, desc = 'Relative weight of lattice data to cohesive energy'),
             Attr('consts',       Text(),     id = True, desc = 'Ordered (by importance) list linear constraints to be included'),
             Attr('reg',          Decimal(),  id = True, desc = 'Regularization penalty') ])

fit = Obj(
    name  = 'fit',
    desc  = 'A constrained fit to some subset of cohesive/BM/lattice data',
    attrs = [Attr('name',Varchar(),id=True),
             Attr('pth',Varchar(),desc='Location of fit job'),
             # Result params
             Attr('timestamp',  Date(),       desc = 'Timestamp of fitting'),

             # Analysis of results
             Attr('opt',        Varchar(),   desc = '5 indices referring to optimal r2 c_viol tradeoff'),
             Attr('score',      Decimal(),   desc = "Arbitrary combination of R2's and c_viol"),
             Attr('result',     Text(),      desc = '5 Flattened NxN fitted coefficients'),
             Attr('decaycosts', Text(),      desc = 'JSONd list of 5 avg R2 values'),
             Attr('c_viol',     Text(),      desc = 'Constraint violations'),
             ]+[Attr(m,Decimal(20,6)) for m in funmetrics])

fit_rels = [Rel('fitparams','fit', id = True), Rel('calc','fit')]

################################################################################

# const = Obj(
#     name  = 'const',
#     desc  ='Constrain Fx(s,a)',
#     attrs = [Attr('const_name',  Varchar(),  id  = True),
#              Attr('description', Text(),     desc='Description of constraint'),
#              Attr('val',         Decimal(),  desc='Value of Fx(s,alpha)'),
#              Attr('kind',        Varchar(),  desc='GT/LT/EQ'),
#              Attr('s',           Decimal(),  desc='If valid only for particular s, else None'),
#              Attr('alpha',       Decimal(),  desc='If valid only for particular alpha, else None')])
#

################################################################################
################################################################################
################################################################################

objs = [calc,functional,beeffunc,job,atoms,bulks,expt_refs,
        fitparams,fit]

rels =  calc_rels + bf_rels + job_rels + atoms_rels + bulks_rels + ea_rels \
        + fit_rels


def new_model() -> Model:
    m = Model('functionals')
    m.add(objs); m.add(rels);
    return m
