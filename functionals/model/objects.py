# External Modules
from typing import Any, Type

# Internal Modules
from dbgen import Model, Obj, Attr, Rel, Int, Varchar, Text, Decimal, Boolean, Date, PathEQ, Path

################################################################################

elem = Obj(
    name  = 'element',
    desc  = 'chemical element',
    attrs = [Attr('atomic_number',                      desc = '# of protons', id = True),
             Attr('symbol',                  Varchar(), desc = 'E.g. He, K, Li'),
             Attr('atomic_weight',           Decimal(), desc = 'atomic_units'),
             Attr('name',                    Varchar()),
             Attr('atomic_radius',                       desc = 'Angstrom'),
             Attr('phase',                   Varchar(),  desc = 'Phase of matter'),
             Attr('group_id',                            desc = 'column in periodic table'),
             Attr('period',                              desc = 'row in periodic table'),
             Attr('evaporation_heat',        Decimal(),  desc = 'kJ/mol'),
             Attr('fusion_heat',             Decimal(),  desc = 'kJ/mol'),
             Attr('melting_point',           Decimal(),  desc = 'K'),
             Attr('is_radioactive',          Boolean()),
             Attr('lattice_struct',          Varchar(),  desc = 'e.g. HEX, SC'),
             Attr('econf',                   Varchar(),  desc = 'electron configuration'),
             Attr('heat_of_formation',       Decimal(),  desc = 'kJ/mol'),
             Attr('electron_affinity',       Decimal(),  desc = 'eV'),
             Attr('boiling_point',           Decimal(),  desc = 'K'),
             Attr('proton_affinity',         Decimal(),  desc = 'kJ/mol'),
             Attr('en_pauling',              Decimal(),  desc = 'electronegativity'),
             Attr('pointgroup',              Varchar()),
             Attr('spacegroup'),
             Attr('metallic_radius',),
             Attr('vdw_radius',              Decimal()),
             Attr('density',                 Decimal()),
             Attr('en_allen',                Decimal()),
             Attr('en_ghosh',                Decimal()),
             Attr('covalent_radius_bragg',   Decimal()),
             Attr('covalent_radius_slater',  Decimal()),
             Attr('geochemical_class',       Varchar()),
             Attr('abundance_crust',         Decimal()),
             Attr('abundance_sea',           Decimal()),
             Attr('atomic_volume',           Decimal()),
             Attr('lattice_constant',        Decimal()),
             Attr('dipole_polarizability',   Decimal()),
             Attr('thermal_conductivity',    Decimal()),
             Attr('gas_basicity',            Decimal()),
            ])

################################################################################

potcar = Obj(
    name  = 'potcar',
    desc  = 'POTCAR pseudopotential file',
    attrs = [Attr('titel', Varchar(),  desc='title? [sic]', id = True),
              Attr('lultra',Boolean(),desc='use ultrasoft pp'),
              Attr('iunscr',Int(),      desc='unscreen: 0-lin, 1-nonlin, 2-no'),
              Attr('rpacor',Decimal(),  desc='partial core radius'),
              Attr('pomass',Decimal(),  desc='mass'),
              Attr('zval',  Decimal(),  desc='valence'),
              Attr('rcore', Decimal(),  desc='outmost cutoff radius'),
              Attr('rwigs', Decimal(),  desc='wigner-seitz radius (au A)'),
              Attr('enmax', Decimal()),
              Attr('enmin', Decimal()),
              Attr('lcor',  Boolean(),desc='correct aug charges'),
              Attr('lpaw',  Boolean(),desc='paw pp'),
              Attr('eaug',  Decimal()),
              Attr('rmax',  Decimal(),  desc='core radius for proj-oper'),
              Attr('raug',  Decimal(),  desc='factor for augmentation sphere'),
              Attr('rdep',  Decimal(),  desc='radius for radial grids'),
              Attr('rdept', Decimal(),  desc='core radius for aug-charge')])

################################################################################
funmetrics = ['%s_%s'%(x,y) for x in ['mse','r2'] for y in ['ce','bm','lat']]

calc = Obj(
    name  = 'calc',
    desc  = 'Calculator details',
    attrs = [Attr('pw',                id = True, desc = 'Planewave cutoff, eV'),
             Attr('econv',  Decimal(),            desc = 'Energy convergance criterion'),
             Attr('a11',    Decimal(), id = True, desc=''),
             Attr('a12',    Decimal(), id = True, desc=''),
             Attr('a13',    Decimal(), id = True, desc=''),
             Attr('a14',    Decimal(), id = True, desc=''),
             Attr('a15',    Decimal(), id = True, desc=''),
             Attr('msb',    Decimal(), id = True, desc='')]
           +[Attr(m,Decimal()) for m in funmetrics])

calc_rels = [Rel('functional','calc',id=True)]

functional = Obj('functional',desc='DFT functional',
                 attrs = [Attr('data',Text(),id=True),Attr('beef',Boolean())
                          ])

beeffunc = Obj('beef',desc='BEEF functional')
bf_rels = [Rel('functional','beef',id=True)]

################################################################################

cellinit = [Attr(x,Decimal(),id=True) for x in [a+b for a in 'abc' for b in '123']]
noninit  = ['a','b','c','surface_area','volume']

cell = Obj(
    name  = 'cell',
    desc  = 'Periodic cells defined by three vectors (all units Angstrom)',
    attrs = cellinit + [Attr(x,Decimal()) for x in noninit])

################################################################################


species = Obj(
    name  = 'species',
    desc  = """
        Abstraction of a struct which throws away position info.

        Like Pure_struct, but contains stoich information and special considerations
        for molecules which are defined by pointgroup rather than spacegroup

        Composition is a python dictionary mapping atomic number to NORMALIZED stoich
        """,
    attrs = [Attr('composition', Varchar(),id=True,  desc='Stringified python dict (ordered, normalized)'),
             Attr('symmetry',    Varchar(),id=True,  desc='''For molecules, pointgroup;
                                                            for bulks, Prototype name;
                                                            for surfaces, underlying prototype+facet'''),
             Attr('nickname',    Varchar(),          desc='Human-readable name, ought (but need not) be unique'),
             Attr('n_elems',                         desc='# of distinct chemical elements in species'),
             Attr('n_atoms',                         desc='Total # of atoms in *normalized* stoich')])

################################################################################

s_d = Obj(
    name  = 'species_dataset',
    desc  = 'Datasets that contain information about chemical species',
    attrs = [Attr('dataset_name', Varchar(), id = True)])

################################################################################


s_d_e = Obj(
    name  = 'species_dataset_element',
    desc  = 'Datum about a particular species',
    attrs = [Attr('property',Varchar(),id=True),
             Attr('value',Text()),
             Attr('datatype',Varchar())])

sde_rels = [Rel('species',        'species_dataset_element', id = True),
            Rel('species_dataset','species_dataset_element', id = True)]

################################################################################


species_comp = Obj(
    name  = 'species_comp',
    desc  = 'Mapping table between species and element to show composition',
    attrs = [Attr('num',desc='Number of this element in lowest integer terms')])

sc_rels = [Rel('species','species_comp',id=True),
           Rel('element','species_comp',id=True)]

################################################################################


pure = Obj(
    name  = 'pure_struct',
    desc  = 'Structure abstraction based on AFLOW prototypes refined by Ankit Jain',
    attrs = [Attr('prototype',Varchar(),id=True),
             Attr('nickname',Varchar())])

################################################################################


struct = Obj(
    name  = 'struct',
    desc  = 'Chemical structure defined in periodic cell',
    attrs = [Attr('raw',        Text(),id=True,desc='JSON encoding of ASE atoms object'),
             Attr('system_type',Varchar(),desc='One of: bulk, molecule, surface'),
             Attr('composition_norm',Text()),
             Attr('n_atoms'),
             Attr('n_elems'),
             Attr('composition',    Text()),
             Attr('metal_comp',     Text()),
             Attr('str_symbols',    Text()),
             Attr('str_constraints',Text())])

struct_rels = [Rel('cell',       'struct'),
               Rel('species',    'struct'),
               Rel('pure_struct','struct')]

################################################################################
incar = Obj(
    name = 'incar',
    desc  = 'VASP input file',
    attrs = [Attr('encut',  Decimal(), id =True, desc='PW cutoff, eV'),
             Attr('sigma',  Decimal(), id = True,desc='Fermi smearing, eV'),
             Attr('metagga',Varchar(), id = True, desc=''),
             Attr('gga',    Varchar(), id = True, desc=''),
             Attr('prec',   Varchar(), id = True, desc=''),
             Attr('ediff',  Decimal(), id = True, desc=''),
             Attr('algo',   Varchar(), id = True, desc=''),
             Attr('ismear', Int(),     id = True, desc=''),
             Attr('npar',   Int(),     id = True, desc=''),
             Attr('nelm',   Int(),     id = True, desc=''),
             Attr('ispin',  Int(),     id = True, desc=''),
             Attr('ibrion', Int(),     id = True, desc=''),
             Attr('lcharg', Boolean(),  id = True, desc=''),
             Attr('lbeefens',Boolean(), id = True, desc=''),
             Attr('addgrid', Boolean(), id = True, desc=''),
             Attr('lasph',  Boolean(),  id = True, desc=''),
             Attr('lwave',  Boolean(), default=False, id = True, desc=''),
             Attr('a11',    Decimal(), id = True, desc=''),
             Attr('a12',    Decimal(), id = True, desc=''),
             Attr('a13',    Decimal(), id = True, desc=''),
             Attr('a14',    Decimal(), id = True, desc=''),
             Attr('a15',    Decimal(), id = True, desc=''),
             Attr('msb',    Decimal(), id = True, desc=''),
             Attr('magmom', Int(),     id = True, desc=''),
             ])

job = Obj(
    name = 'job',
    desc  = 'DFT calculation',
    attrs = [#Attr('logfile',      Varchar(),id=True,desc='Path to primary log file'),
           Attr('stordir',      Varchar(), id=True,     desc='Contains logfile'),
           Attr('log',          Text('long'),   desc='Content of primary log file'),
           Attr('timestamp',    Date(),         desc='Unix timestamp on log file'),
           Attr('user',         Varchar(),      desc='Who owns the directory'),
           #Attr('distribs',     Text('long'),   desc='3 Serialized arrays: s, alpha, and density'),
           Attr('int_occupation',Boolean(),   desc='Whether or not all orbitals have integer electron occupation '
                                                    '(something that ought be true for single atom calcs)'),
           Attr('spinpol',      Boolean(),    desc='Whether calculation was spin polarized'),
           Attr('energy',       Decimal(),      desc='eV')]+[
           Attr('contribs',     Text(),   desc='list of five 64-element exchange contribs, if beef calculation')] + [
           Attr('kptden_'+x,     Decimal(),     desc='K point density')
             for x in 'xyz'] + [
           Attr('k'+x,           Decimal(),     desc='# of K points')
             for x in 'xyz'] )

job_rels = [Rel('calc',  'job'),
            Rel('struct','job'),
            Rel('incar', 'job')]

atom = Obj(
    name  = 'atom',
    desc  = 'An atom, considered within a specific chemical structure',
    attrs = [Attr('ind',id=True,             desc = 'ASE atom index'),
            Attr('number',                  desc = 'Atomic number'),
            Attr('x',                       desc = 'position'),
            Attr('y',                       desc = 'position'),
            Attr('z',                       desc = 'position'),
            Attr('constrained',Boolean(), desc = 'Whether or not there was a FixAtoms constraint'),
            Attr('magmom',                  desc = 'Units: Bohr'),
            Attr('tag',                     desc = 'ASE atom tag')])

atom_rels = [Rel('struct', 'atom',id=True),
             Rel('element','atom'),
             Rel('potcar','atom')]

################################################################################

expt = Obj(
    name = 'expt',
   desc  = 'Set of single points on a particular material (with a particular calc)',
   attrs = [Attr('n_atoms',     id=True,        desc='Value for all jobs in this experiment'),
            Attr('energy_pa',   Decimal(),      desc='Per atom, eV'),
            Attr('bulkmod',     Decimal(),      desc='Bulk modulus, GPa'),
            Attr('volume_pa',   Decimal(),      desc='Per atom, A^3'),
            Attr('all_vols',    Text('long'),   desc='every volume of the related bulk jobs'),
            Attr('volumes',     Text('long'),   desc='volumes of the 5 most optimal jobs'),
            Attr('energies',    Text('long'),   desc='energies of the 5 most optimal jobs'),
            Attr('img',         Text('long'),   desc='base64 encoded image of EOS fit'),
            Attr('lattice',     Decimal(),      desc='Conventional unit cell lattice (optimized)'),
            Attr('n_data',                      desc='Number of aggregated data points'),
            Attr('min_gap',     Decimal(),      desc='Absolute difference between best singlepoint volume_pa and the fitted optimum'),
            Attr('complete',    Boolean(),    desc='Whether enough jobs are done to do fitting (at least 5 w/ clear minimum)'),

            Attr('name',            Varchar(),   desc='Species nickname'),
                Attr('eform',       Decimal(),      desc='Per atom, eV'),
                Attr('coefs',           Text('long'),desc='Calc Coefs'),
                Attr('composition',     Varchar(),   desc='Species composition'),
                #Attr('x_bulk',       Text('long'),  desc='Serialized matrix with 5 rows, each corresponding to the exchange energy contribs of the best 5 bulk jobs'),
                #Attr('dx_bulk_atom',    Text('long'),desc='Serialized vector of the difference in exchange energy contributions between (optimum) bulk and (stoichiometrically-weighted) atoms'),
                Attr('atomic_contribs', Text('long'),desc='Serialized dict of relevant reference info'),
                Attr('atomic_energies', Text(),      desc='Serialized dict of relevant reference info'),
                Attr('bulk_contribs',   Text('long'),desc='Best job xc contribs'),
                Attr('bulk_energy',     Decimal(),   desc='Best job energy'),
                Attr('expt_cohesive_energy',Decimal(),desc='Experimental cohesive energy'),
                Attr('expt_bm',         Decimal(),   desc='Experimental bulk modulus'),
                Attr('expt_volume',     Decimal(),   desc='Experimental volume of reference stoichiometry'),
                # Attr('energy_vector',   Text('long'),desc="JSON'd vector of 5 energies"),
                # Attr('volume_vector',   Text('long'),desc="JSON'd vector of 5 volumes"),
                # Attr('contrib_vector',  Text('long'),desc="JSON'd 5x64 matrix with exchange contributions"),
                Attr('bulk_ratio',                   desc='Ratio of bulk system to size of normalized species'),
                Attr('contribs',                 Text(),         desc='exchange contribs of the 5 most optimal jobs')])


expt_rels = [Rel('species','expt',id = True),
             Rel('calc',   'expt',id = True),
             Rel('best_job','expt',   'bulk_job')]

################################################################################


bj = Obj(
    name = 'bulk_job',
    desc  = 'A subset of jobs which have a many-one relationship linking jobs to an experiment',
    attrs =[Attr('dv',      Decimal(),  desc = "Difference in volume from 'minimum' of the expt it belongs to"),
             Attr('gap',     Decimal(),  desc = 'Abs(dv)'),
             Attr('near_min',Boolean(),desc = 'Whether this job is in the bottom 5 data points')])

bj_rels = [Rel('job',     'bulk_job',id=True),
           Rel('expt',    'bulk_job')]

################################################################################

ref = Obj(
    name = 'reference',desc='A single calculation that gives the energy of an isolated atom (w/ a calc)',
    attrs = [Attr('energy', Decimal(), desc='Energy copied over from job table'),
             Attr('contribs', Text(),  desc = 'BEEF Contributions, copied over')])

ref_rels = [Rel('job',    'reference', id = True),
            Rel('calc',   'reference'),
            Rel('element','reference')]


################################################################################

expt_refs = Obj(
    name = 'expt_refs',
    desc = 'Mapping table between bulk expt and relevant atomic calculations')

ea_rels = [Rel('expt',      'expt_refs', id=True),
           Rel('reference', 'expt_refs', id=True)]

################################################################################
fitparams = Obj(
    name = 'fitparams',
    desc = 'Input parameters to a fit',
    attrs = [Attr('bm_weight',    Decimal(),  id = True, desc = 'Relative weight of bulk-modulus data to cohesive energy'),
             Attr('lat_weight',   Decimal(),  id = True, desc = 'Relative weight of lattice data to cohesive energy'),
             Attr('constden',                 id = True, desc = 'Number of s or alpha points linearly constrained between logscale(-2,2)'),
             Attr('consts',       Text(),     id = True, desc = 'Ordered (by importance) list linear constraints to be included'),
             Attr('dataconst',    Text(),     id = True, desc = 'Regex on what data should be included'),
             Attr('reg',        Decimal(),    id=True,   desc = 'Regularization penalty')
             ])

fit = Obj(
    name  = 'fit',
    desc  = 'A constrained fit to some subset of cohesive/BM/lattice data',
    attrs = [Attr('name',Varchar(),id=True), Attr('pth',Varchar(),desc='Location of fit job'),
             Attr('a1', Decimal(), id=True,desc='Which value of a1 contribs was used'),
             Attr('msb',Decimal(),desc='Which value of msb was used'),
             Attr('decay',                  desc = '1-5'),
             # Result params
             Attr('nsteps',                  desc = 'Number of steps'),
             Attr('steps',      Text('long'), desc = 'JSON of trajectory data'),
             Attr('timestamp',  Date(),      desc = 'Timestamp of fitting'),
             Attr('runtime',    Decimal(),   desc = 'Duration of fitting, s'),
             Attr('score',      Decimal(),   desc = "Arbitrary combination of R2's and c_viol"),
             Attr('result',     Text(),      desc = 'Flattened NxN fitted coefficients'),
             Attr('c_viol',     Decimal(),   desc = 'Constraint violation'),
             Attr('lda_viol',   Decimal(),   desc = 'Degree to which LDA Limit was violated'),
             Attr('h_viol',     Decimal(),   desc = 'Degree to which hydrogen energy was violated'),
             ]+[Attr(m,Decimal()) for m in funmetrics])

fit_rels = [Rel('fitparams','fit', id = True),Rel('calc','fit')]



const = Obj(
    name = 'const',
    desc='Constrain Fx(s,a)',
    attrs=[Attr('const_name',  Varchar(),id=True),
            Attr('description', Text(),     desc='Description of constraint'),
            Attr('val',         Decimal(),  desc='Value of Fx(s,alpha)'),
            Attr('kind',        Varchar(),  desc='GT/LT/EQ'),
            Attr('s',           Decimal(),  desc='If valid only for particular s, else None'),
            Attr('alpha',       Decimal(),  desc='If valid only for particular alpha, else None')])

################################################################################
# fit_const = Obj('fit_const',desc='Mapping table to denote which constraints were used in a given fit',
#                 attrs=[Attr('const_weight',Decimal())])
#
# fc_rels = [Rel('fit','fit_const',id=True),Rel('const','fit_const',id=True)]
# nl_const = Obj('nonlin_const',desc='Nonlinear constraints',
#                attrs = [Attr('nlconst_name', Varchar(),id=True),
#                         Attr('description',  Text(),    desc='Description of nonlinear constraint'),
#                         Attr('f',            Text(),    desc='Source code for nonlin function to be minimized'),
#                         Attr('df',           Text(),    desc='Source code for derivative of f'),
#                         Attr('hess',         Text(),    desc='Source code for hessian of  f'),
#                         Attr('lb',           Decimal(), desc='Lower bound'),
#                         Attr('ub',           Decimal(), desc='Upper bound')])
# fnlc = Obj('fit_nonlin_const',desc='Mapping table to denote which constraints were used in a given fit',
#            attrs = [Attr('nl_const_weight',Decimal())])
#
# fnlc_rels = [Rel('fit',         'fit_nonlin_const',id=True),
#              Rel('nonlin_const','fit_nonlin_const',id=True)]
################################################################################
# fit_step = Obj('fit_step',desc='A single iteration in a fitting process',
#                 attrs = [Attr('niter',id=True,   desc='Iteration number'),
#                          Attr('cost',  Decimal(),desc='Objective function cost'),
#                          Attr('c_viol',Decimal(),desc='Constraint cost')])
#
# fs_rels = [Rel('fit','fit_step',id=True)]
# fit_data = Obj('fit_expt',
#                desc = 'Mapping table specifying which  was used in a given fit',)
#
# fd_rels = [Rel('fitparams', 'fit_expt',id=True),
#            Rel('expt',      'fit_expt',id=True)]
################################################################################

################################################################################
################################################################################
################################################################################

objs = [elem, potcar, calc, cell, species, s_d, s_d_e, incar,fitparams,
        species_comp, pure, struct, job,  atom, expt, bj, ref,
        fit,  const,  expt_refs,
        functional, beeffunc ] # nl_const,fit_const, fnlc, fit_data,fit_step,

rels =  sde_rels+ sc_rels + struct_rels + job_rels + atom_rels + fit_rels\
        + expt_rels + bj_rels + ref_rels  +   ea_rels \
         + calc_rels + bf_rels # + fc_rels + fnlc_rels+ fd_rels +fs_rels +

peqs = [PathEQ(Path(attr=ref['energy']),
               Path(rels=[ref.r('job')],attr=job['energy'])),
        PathEQ(Path(attr=expt.id),
               Path(rels=[expt.r('best_job'),bj.r('expt')],attr=expt.id))]

def new_model() -> Model:
    m = Model('functionals')
    m.add(objs); m.add(rels); m.add(peqs)
    return m
