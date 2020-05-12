
# Internal Modules
from dbgen import (Model, Obj, Attr, Rel, Int, Varchar, Text, Decimal, Boolean,
                   RawView)

#################
# INPUT SUMMARY #
#################

###############################################################################
mets = ['ce', 'bm', 'lat', 'vol', 'mag']
funmetrics = [y+'mae_%s' % (x) for x in mets for y in ['', 'rel']]
errmetrics = [y+'err_%s' % (x) for x in mets for y in ['', 'rel']]

ms = [Attr(m, Decimal(), desc='Only True if results in for all mats')
      for m in funmetrics]
calc = Obj(
    name='calc',
    desc='Calculator details',
    attrs=[Attr('allmat', Text(), desc='List of all (bulk) mats with data'),
           Attr('done', Boolean(),
                desc='Whether or not all calculations have been finished'),
           Attr('missing', Text(),
                desc='Materials not yet finished (bulk or atom related)'),
           Attr('n_missing', Int(),
                desc='Materials not yet finished (bulk or atom related)'),
           Attr('fitdata', Text(),
                desc='For a BEEF calc, the input data to fitting'),
           Attr('name', Text(), identifying=True, desc='Functional nickname'),
           Attr('data', Text(), desc='BEEF coefs, if any'),
           Attr('beef', Boolean(), desc='True if BEEF style')] + ms)
# fks=[Rel('functional', identifying=True)])

##############################################################################

# functional = Obj(
#     name='functional',
#     desc='DFT functional',
#     attrs=[Attr('data', Text(), identifying=True), Attr('beef', Boolean())])

##############################################################################

#########
# JOBS  #
#########

##############################################################################

atoms = Obj(
    name='atoms',
    desc='DFT calculation',
    attrs=[
        Attr('stordir', Varchar(), identifying=True),
        # Attr('runtime', Int(), desc='Appx runtime, hours'),
        # Attr('pw', desc='Planewave cutoff, eV'),
        # Attr('econv', Decimal(), desc='EDIFF in the INCAR'),
        Attr('name', Varchar()),
        Attr('calcname', Varchar()),
        Attr('err', Text(), desc='No problems detected if empty string'
             ' - Else: econv, mag diff, intocc, maxstep, unconverged'),
        Attr('num', Int(), desc='Atomic number'),
        Attr('true_mag', Int(), desc='Magnetic moment'),
        Attr('energy', Decimal(), desc='eV'),
        Attr('contribs', Text(), desc='JSONd list xc contribs')],
    fks=[Rel('calc')]
)

bulks = Obj(
    name='bulks',
    desc='Set of singlepoints on a particular material under strain',
    attrs=[
        Attr('stordir', Varchar(), identifying=True),
        Attr('name', Varchar(), desc='Species nickname'),
        Attr('err', Varchar(), desc='Error in latopt'),
        Attr('calcname', Varchar(), desc='Calc nickname'),
        Attr('composition', Varchar(), desc='Species composition'),
        Attr('elems', Varchar(), desc='common separated list of atomic nums'),
        Attr('struct', Varchar(), desc='fcc/bcc/rocksalt/etc.'),
        Attr('n_atoms', desc='Number of atoms in unit cell'),
        Attr('n_elems', desc='Number of distinct chemical species'),
        Attr('alloy', Varchar(), desc='What type of alloy it is, if any'),

        # Analysis Results from minimum 5 jobs
        Attr('success', Boolean(), desc='DFT data for each experiment data'),
        # Attr('optstrains', Text(), desc='List of completed opt strains'),
        # Attr('energy', Decimal(),desc='Minimum energy according to quadfit'),
        Attr('mag', Decimal(), desc='Magmom of the minimum energy job'),
        Attr('contrib', Text('long'), desc='xc contribs of 5 optimal jobs'),
        Attr('eng', Decimal(), desc='unit cell lattice from quadfit'),
        Attr('lat', Decimal(), desc='unit cell lattice from quadfit'),
        Attr('vol', Decimal(), desc='Volume from quad fit'),
        Attr('volrat', Decimal(), desc='lattice/vol^1/3'),

        Attr('volumes', Text('long'), desc='volumes of the 5 optimal jobs'),
        Attr('energies', Text('long'), desc='energies of the 5 optimal jobs'),
        Attr('contribs', Text('long'), desc='xc contribs of 5 optimal jobs'),
        Attr('bm', Decimal(), desc='Bulk modulus, GPa, from 5 jobs + stencil'),
        Attr('bulkmod_lstsq', Decimal(), desc='Bulkmod, GPa, from quadfit'),
        Attr('lstsq_abc', Text(), desc='JSON of quad fit of 5 optimal jobs'),
        Attr('eform', Decimal(), desc='Formation energy, eV, given ENG'),
        Attr('ce', Decimal(), desc='Cohesive  energy, eV'),
        # Analysis from ase.eos
        Attr('eosbm', Decimal(), desc='Bulk modulus, GPa, from ASE EOS'),
        Attr('eos_diff', Decimal(), desc='Ratio of discrete BM to EOS BM'),
        Attr('irregular', Boolean(), desc='discrete BM differs from EOS'),

        # Experimental results
        Attr('expt_ce', Decimal(15, 3), desc='Expt cohesive energy, eV'),
        Attr('expt_bm', Decimal(15, 3), desc='Expt bulk modulus, GPa'),
        Attr('expt_lat', Decimal(15, 6),
             desc='Experimental lattice parameter, A'),
        Attr('expt_vol', Decimal(), desc='Experimental volume, A^3'),
        Attr('expt_mag', Decimal(), desc='Experimental magnetic moment, bohr'),
        Attr('hasdata', Boolean(), desc='We have any data for this'),

        # Relative error

        # Fitting inputs
        Attr('ab_ce', Text(), desc='a vector+offset to be dotted w/ '\
             'BEEF coef + offset gives formation energy in eV'),
        Attr('ab_bm', Text(), desc='a vector+offset to be dotted w/ '\
             'BEEF coef + offset gives bulk modulus in GPa'),
        Attr('ab_vol', Text(), desc='a vector+offset to be dotted w/ '\
             'BEEF coef + offset gives conventional lattice in A')] +
    [Attr(x, Decimal()) for x in errmetrics],
    fks=[Rel('calc')])


##############################################################################

expt_refs = Obj(
    name='refs',
    desc='Mapping table between bulk and relevant atomic calculations',
    attrs=[Attr('num', Int(), desc='Stoichiometric coefficient'),
           Attr('energy', Decimal(), desc='Atom energy, weighted by stoich'),
           Attr('contribs', Text(),
                desc='Exchange contributions of atom, weighted by stoich')],
    fks=[Rel('bulks', identifying=True),
         Rel('atoms', identifying=True)])

##############################################################################

###########
# FITTING #
###########

fitparams = Obj(
    name='fitparams',
    desc='Input parameters to a fit',
    attrs=[Attr('ce_scale', Decimal(), identifying=True, desc='One cost unit'),
           Attr('bm_scale', Decimal(), identifying=True, desc='One cost unit'),
           Attr('lc_scale', Decimal(), identifying=True, desc='One cost unit'),
           Attr('consts', Text(), identifying=True,
                desc='Special weights or points to use'),
           Attr('abdata', Text(), desc='computed from consts')])

# cons = Obj(
#     name='const', desc='CONSTRAINT',
#     attrs=[Attr('name', Varchar(), identifying=True),
#            Attr('abdata', Text()), Attr('kind', Varchar()),
#            Attr('points', Text()),
#            Attr('val', Decimal()), Attr('func', Varchar())]
# )

fit = Obj(
    name='fit',
    desc='A constrained fit to some subset of cohesive/BM/lattice data',
    attrs=[Attr('x', Text(), desc='Coefficients of optimal point in the traj'),
           Attr('cv', Text(), desc='Summary of cross validation data'),
           Attr('err', Text(), desc='Scipy error message'),
           Attr("step", Int(), identifying=True,
                desc='Optimization iteration step')]
    + [Attr(m, Decimal(), desc='') for m in funmetrics],
    fks=[Rel('fitparams', identifying=True), Rel('calc', identifying=True)])

surf = Obj(
    name='surf',
    desc='2x2x4 surfaces with CO in either on top or hollow position',
    attrs=[Attr('mat', Varchar(), identifying=True, desc='Which material'),
           Attr('xc', Varchar(), identifying=True, desc='Which fx'),
           Attr('ontop', Decimal(), desc='On top energy'),
           Attr('hollow', Decimal(), desc='Hollow energy'),
           Attr('err', Decimal(), desc='Agreement with expt relative error')]
)
##############################################################################
# Views #

err = RawView('errs', '''
SELECT name,mae_ce,mae_bm,mae_lat,mae_mag,
       relmae_ce,relmae_bm,relmae_lat,relmae_mag,
       true as isfx
FROM calc
UNION
SELECT fit_id::text,mae_ce,mae_bm,mae_lat,NULL,
       relmae_ce,relmae_bm,relmae_lat, NULL,
       false as isfx
FROM fit
''')

##############################################################################
##############################################################################
##############################################################################

objs = [calc, atoms, bulks, fitparams, fit, expt_refs,
        surf]


def new_model() -> Model:
    m = Model('functionals', objlist=objs, viewlist=[err])
    return m
