# External Modules
from typing import Any, Type

# Internal Modules
from dbgen2 import Model, Obj, Attr, Rel, Int, Varchar, Text, Decimal, Date

################################################################################
# Add objects and relations
#-------------------------
globals = Obj('globals',
              desc  = 'Properties of entire database. There can be max one row',
              attrs = [Attr('all_data',Text('long'),
                            desc = 'All cohesive data'),
                       Attr('all_constraints',Text('long'),
                            desc = 'All linear constraints'),
                        Attr('all_nlconstraints',Text('long'),
                            desc = 'All nonlinear constraints')])
elem = Obj('element',desc = 'chemical element',
           attrs = [Attr('atomic_number', desc = '# of protons', id = True),
                    Attr('symbol', Varchar(), desc = 'E.g. He, K, Li'),
                    Attr('atomic_weight', Decimal(), desc = 'atomic_units'),
                    Attr('name', Varchar()),
                    Attr('atomic_radius',  desc = 'Angstrom'),
                    Attr('phase',  Varchar(),desc = 'Phase of matter'),
                    Attr('group_id', desc = 'column in periodic table'),
                    Attr('period', desc = 'row in periodic table'),
                    Attr('evaporation_heat', Decimal(), desc = 'kJ/mol'),
                    Attr('fusion_heat', Decimal(), desc = 'kJ/mol'),
                    Attr('melting_point', Decimal(), desc = 'K'),
                    Attr('is_radioactive', Int('tiny')),
                    Attr('lattice_struct', Varchar(), desc = 'e.g. HEX, SC'),
                    Attr('econf', Varchar(),desc = 'electron configuration'),
                    Attr('heat_of_formation', Decimal(), desc = 'kJ/mol'),
                    Attr('electron_affinity', Decimal(), desc = 'eV'),
                    Attr('boiling_point', Decimal(), desc = 'K'),
                    Attr('proton_affinity', Decimal(), desc = 'kJ/mol'),
                    Attr('en_pauling', Decimal(), desc = 'electronegativity'),
                    Attr('pointgroup',Varchar()),
                    Attr('spacegroup'),
                    Attr('metallic_radius',),
                    Attr('vdw_radius',Decimal()),
                    Attr('density',Decimal()),
                    Attr('en_allen',Decimal()),
                    Attr('en_ghosh',Decimal()),
                    Attr('covalent_radius_bragg',Decimal()),
                    Attr('covalent_radius_slater',Decimal()),
                    Attr('geochemical_class',Varchar()),
                    Attr('abundance_crust',Decimal()),
                    Attr('abundance_sea',Decimal()),
                    Attr('atomic_volume',Decimal()),
                    Attr('lattice_constant',Decimal()),
                    Attr('dipole_polarizability',Decimal()),
                    Attr('thermal_conductivity',Decimal()),
                    Attr('gas_basicity',Decimal()),
                    ])

set_fam = Obj('setup_family',desc = "Class of Setups that are the same 'type' ",
              attrs = [Attr('name',Varchar(),id=True),
                      Attr('kind',Varchar(),id=True),
                      Attr('xc',Varchar(),id=True)])

setup = Obj('setup',desc='Pseudopotential',
            attrs=[Attr('checksum',Varchar(),desc='MD5 hash of file',id=True),
                   Attr('val',desc='Number of valence electrons')])

set_rels = [Rel('element','setup'), Rel('setup_family','setup')]
calc = Obj('calc',desc='Calculator details',
           attrs = [Attr('beef',Int('tiny'),id=True,desc='Whether or not functional is BEEF-style, (i.e. whether or not the coefs attr is meaningful)'),
                    Attr('coefs',Text(),id=True,desc='JSON dumped 8x8 matrix of exchange coefficients which defines a functional as a sum of bases'),
                    Attr('xc',Varchar(),id=True,desc='Name of functional'),
                    Attr('pw',id=True,desc='Planewave cutoff, eV'),
                    Attr('econv',Decimal(),desc='Energy convergance criterion'),])

cellinit = [Attr(x,Decimal(),id=True) for x in [a+b for a in 'abc' for b in '123']]
noninit  = ['a','b','c','surface_area','volume']
cell = Obj('cell', desc = 'Periodic cells defined by three vectors (all units Angstrom)',
           attrs = cellinit +
                   [Attr(x,Decimal()) for x in noninit])

species = Obj('species',desc = """
Abstraction of a struct which throws away position info.

Like Pure_struct, but contains stoich information and special considerations
for molecules which are defined by pointgroup rather than spacegroup

Composition is a python dictionary mapping atomic number to NORMALIZED stoich
""", attrs = [
Attr('composition',Varchar(),id=True,desc='Stringified python dict (ordered, normalized)'),
Attr('symmetry',Varchar(),id=True,desc='''For molecules, pointgroup;
for bulks, Prototype name;
for surfaces, underlying prototype+facet'''),
Attr('nickname',Varchar(),desc='Human-readable name, ought (but need not) be unique'),
Attr('n_elems',desc='# of distinct chemical elements in species'),
Attr('n_atoms',desc='Total # of atoms in *normalized* stoich')])

s_d = Obj('species_dataset',desc='Datasets that contain information about chemical species',
          attrs = [Attr('dataset_name', Varchar(), id = True)])

s_d_e = Obj('species_dataset_element',desc='Datum about a particular species',
            attrs = [Attr('property',Varchar(),id=True),
                     Attr('value',Text()),
                     Attr('datatype',Varchar())])

sde_rels = [Rel('species','species_dataset_element', id = True),
            Rel('species_dataset','species_dataset_element', id = True)]

species_comp = Obj('species_comp',desc= 'Mapping table between species and element to show composition',
                    attrs = [Attr('num',desc='Number of this element in lowest integer terms')])
sc_rels = [Rel('species','species_comp',id=True),
            Rel('element','species_comp',id=True)]

pure = Obj('pure_struct',desc='Structure abstraction based on AFLOW prototypes refined by Ankit Jain',
           attrs =[Attr('prototype',Varchar(),id=True),
                   Attr('spacegroup'),
                   Attr('free',desc='# of free params'),
                   Attr('nickname',Varchar())])

struct = Obj('struct',desc='Chemical structure defined in periodic cell',
             attrs =[Attr('raw',Text(),id=True,desc='JSON encoding of ASE atoms object'),
                     Attr('system_type',Varchar(),desc='One of: bulk, molecule, surface'),
                     Attr('composition_norm',Text()),
                     Attr('n_atoms'),
                     Attr('n_elems'),
                     Attr('composition',Text()),
                     Attr('metal_comp',Text()),
                     Attr('str_symbols',Text()),
                     Attr('str_constraints',Text())])

struct_rels = [Rel('cell','struct'),Rel('species','struct'),Rel('pure_struct','struct')]

job = Obj('job',desc='DFT calculation',
          attrs = [Attr('logfile',Varchar(),id=True,desc='Path to primary log file'),
                   Attr('stordir',Varchar(),desc='Contains logfile'),
                   Attr('log',Text('long'),desc='Content of primary log file'),
                   Attr('timestamp',Date(),desc='Unix timestamp on log file'),
                   Attr('user',Varchar(),desc='Who owns the directory'),
                   Attr('kx',desc='K points'),
                   Attr('ky',desc='K points'),
                   Attr('kz',desc='K points'),
                   Attr('kptden_x',Decimal(),desc='K point density'),
                   Attr('kptden_y',Decimal(),desc='K point density'),
                   Attr('kptden_z',Decimal(),desc='K point density'),
                   Attr('spinpol',Int('tiny'),desc='Whether calculation was spin polarized'),
                   Attr('has_contribs',Int('tiny'),desc='Whether xc_contribs.txt exists'),
                   Attr('contribs',Text(),desc='output exchange contribs, if beef calculation'),
                   Attr('energy',Decimal(),desc='eV')] )

job_rels = [Rel('calc','job'),Rel('struct','job')]

jobset = Obj('job_setup',desc='mapping table between Setup and Job')
jobsetrel = [Rel('job','job_setup',id=True),Rel('setup','job_setup',id=True)]

atom = Obj('atom',desc='An atom, considered within a specific chemical structure',
           attrs = [Attr('ind',id=True,desc='ASE atom index'),
                    Attr('number',desc='Atomic number'),
                    Attr('x',desc='position'),
                    Attr('y',desc='position'),
                    Attr('z',desc='position'),
                    Attr('constrained',Int('tiny'),desc='Whether or not there was a FixAtoms constraint'),
                    Attr('magmom',desc='Units: Bohr'),
                    Attr('tag',desc='ASE atom tag')])
atom_rels = [Rel('struct','atom',id=True),Rel('element','atom')]

expt = Obj('expt',desc='Set of single points on a particular material (with a particular calc)',
           attrs = [Attr('n_atoms',id=True,desc='Value for all jobs in this experiment'),
                    Attr('energy_pa',Decimal(),desc='Per atom, eV'),
                    Attr('bulkmod',Decimal(),desc='Bulk modulus, GPa'),
                    Attr('volume_pa',Decimal(),desc='Per atom, A^3'),
                    Attr('all_vols',Text('long'),desc='every volume of the related bulk jobs'),
                    Attr('volumes',Text('long'),desc='volumes of the 5 most optimal jobs'),
                    Attr('energies',Text('long'),desc='energies of the 5 most optimal jobs'),
                    Attr('contribs',Text('long'),desc='exchange contribs of the 5 most optimal jobs'),
                    Attr('img',Text('long'),desc='base64 encoded image of EOS fit'),
                    Attr('eform',Decimal(),desc='Per atom, eV'),
                    Attr('lattice',Decimal(),desc='Conventional unit cell lattice (optimized)'),
                    Attr('n_data',desc='Number of aggregated data points'),
                    Attr('min_gap',Decimal(),desc='Absolute difference between best singlepoint volume_pa and the fitted optimum')])

expt_rels = [Rel('species','expt',id=True),Rel('calc','expt',id=True)]

bj = Obj('bulk_job',desc='A subset of jobs which have a many-one relationship linking jobs to an experiment',
         attrs =[Attr('dv',Decimal(),desc="Difference in volume from 'minimum' of the expt it belongs to"),
                 Attr('gap',Decimal(),desc='Abs(dv)'),
                 Attr('near_min',Int('tiny'),desc='Whether this job is in the bottom 5 data points')])
bj_rels = [Rel('job','bulk_job',id=True),Rel('expt','bulk_job')]


ref = Obj('reference',desc='A single calculation that gives the energy of an isolated atom (w/ a calc)',
          attrs=[Attr('energy',Decimal())])
ref_rels = [Rel('job','reference',id=True),
            Rel('calc','reference'),Rel('element','reference')]

dft_data = Obj('dft_data',desc='',
               attrs = [Attr('name',Varchar(),desc='Species nickname'),
                        Attr('coefs',Text('long'),desc='Calc Coefs'),
                        Attr('composition',Varchar(),desc='Species composition'),
                        Attr('atomic_contribs',Text('long'),desc='Serialized dict of relevant reference info'),
                        Attr('atomic_energies',Text(),desc='Serialized dict of relevant reference info'),
                        Attr('bulk_contribs',Text('long'),desc='Best job xc contribs'),
                        Attr('bulk_energy',Decimal(),desc='Best job energy'),
                        Attr('expt_cohesive_energy',Decimal(),desc='Experimental cohesive energy'),
                        Attr('expt_bm',Decimal(),desc='Experimental bulk modulus'),
                        Attr('expt_volume',Decimal(),desc='Experimental volume of reference stoichiometry'),
                        Attr('energy_vector',Text('long'),desc="JSON'd vector of 5 energies"),
                        Attr('volume_vector',Text('long'),desc="JSON'd vector of 5 volumes"),
                        Attr('contrib_vector',Text('long'),desc="JSON'd 5x64 matrix with exchange contributions"),
                        Attr('bulk_ratio',desc='Ratio of bulk system to size of normalized species'),

                        ])
dft_rels = [Rel('expt','dft_data',id=True),Rel('best_job','dft_data','job')]

fit = Obj('fit',desc = 'A constrained fit to some subset of cohesive/BM/lattice data',
          attrs = [Attr('name',Varchar(),id=True),
                   # input params
                   Attr('bm_weight',Decimal(),desc='Relative weight of bulk-modulus data to cohesive energy'),
                   Attr('lat_weight',Decimal(),desc='Relative weight of lattice data to cohesive energy'),
                   Attr('consts',Text(),desc='SQL const on what linear constraints should be included'),
                   Attr('nlconsts',Text(),desc='SQL const on what linear constraints should be included'),
                   Attr('dataconst',Text(),desc='SQL const on what data should be included'),
                   Attr('basis',desc='Size of fitted functional'),
                   Attr('initfit',Int('tiny'),desc='If true: initialize with lstsq fit w/o constraints (else with 0)'),
                   Attr('bound',Decimal(),desc='Range over which to search for coefficients'),
                   Attr('maxiter',desc='Stop nonlinear fitting after this step'),
                   Attr('constden',desc='Number of s or alpha points linearly constrained between logscale(-2,2)'),
                   # Intermediate computation
                   Attr('raw_data',Text('long'),desc='Data to be fed to fitting script'),
                   Attr('raw_const',Text('long'),desc='Data to be fed to fitting script'),
                   Attr('raw_nlconst',Text('long'),desc='Data to be fed to fitting script'),
                   Attr('n_const',default=0,desc='number of constraints'),
                   Attr('n_data',default=0,desc='Number of data points'),
                   Attr('n_nlconst',default=0,desc='Number of nonlinear constraints'),
                   # Result params
                   Attr('timestamp',Date(),desc='Timestamp of fitting'),
                   Attr('runtime',Decimal(),desc='Duration of fitting, s'),
                   Attr('r2_ce',Decimal(),desc='R2 fit of cohesive energies'),
                   Attr('r2_bm',Decimal(),desc='R2 fit of bulk moduli'),
                   Attr('r2_lat',Decimal(),desc='R2 fit of lattice constants'),
                   Attr('c_viol',Decimal(),desc='Constraint violation'),
                   Attr('score',Decimal(),desc="Arbitrary combination of R2's and c_viol"),
                   Attr('result',Text(),desc='Flattened NxN fitted coefficients'),
                   Attr('fullresult',Text(),desc='Flattened 8x8 fitted coefficients'),
                   Attr('log',Text('long'),desc='Output from fitting'),
                   Attr('err',Text(),desc='Error during scipy.minimize()'),
                   Attr('lda_viol',Decimal(),desc='Degree to which LDA Limit was violated'),
                   Attr('beefdist',Decimal(),desc='Score for cartesian distance btw output and BEEF')])

fit_step = Obj('fit_step',desc='A single iteration in a fitting process',
                attrs = [Attr('niter',id=True,desc='Iteration number'),
                         Attr('cost',Decimal(),desc='Objective function cost'),
                        Attr('c_viol',Decimal(),desc='Constraint cost')])

fs_rels = [Rel('fit','fit_step',id=True)]

fit_data = Obj('fit_data',desc='Mapping table specifying which Cohesive_data was used in a given fit')
fd_rels = [Rel('fit','fit_data',id=True),Rel('dft_data','fit_data',id=True)]

const = Obj('const',desc='Constrain Fx(s,a)',
             attrs=[Attr('const_name',Varchar(),id=True),
                    Attr('description',Text(),desc='Description of constraint'),
                    Attr('val',Decimal(),desc='Value of Fx(s,alpha)'),
                    Attr('kind',Varchar(),desc='GT/LT/EQ'),
                    Attr('s',Decimal(),desc='If valid only for particular s, else None'),
                    Attr('alpha',Decimal(),desc='If valid only for particular alpha, else None'),
                    Attr('vec',Text(),desc='Raw 64-float JSON, ignore s/a if defined')])

fit_const = Obj('fit_const',desc='Mapping table to denote which constraints were used in a given fit',
                attrs=[Attr('const_weight',Decimal())])
fc_rels = [Rel('fit','fit_const',id=True),Rel('const','fit_const',id=True)]

nl_const = Obj('nonlin_const',desc='Nonlinear constraints',
               attrs = [Attr('nlconst_name',Varchar(),id=True),
                        Attr('description',Text(),desc='Description of nonlinear constraint'),
                        Attr('f',Text(),desc='Source code for nonlin function to be minimized'),
                        Attr('df',Text(),desc='Source code for derivative of f'),
                        Attr('hess',Text(),desc='Source code for hessian of  f'),
                        Attr('lb',Decimal(),desc='Lower bound'),
                        Attr('ub',Decimal(),desc='Upper bound')])
fnlc = Obj('fit_nonlin_const',desc='Mapping table to denote which constraints were used in a given fit',
           attrs = [Attr('nl_const_weight',Decimal())])
fnlc_rels = [Rel('fit','fit_nonlin_const',id=True),
             Rel('nonlin_const','fit_nonlin_const',id=True)]

objs = [globals, elem, set_fam, setup, calc, cell, species, s_d, s_d_e,
        species_comp, pure, struct, job, jobset, atom, expt, bj, ref, dft_data,
        fit, fit_step,fit_data, const, fit_const, nl_const, fnlc ]
rels = set_rels + sde_rels + sc_rels + struct_rels + job_rels+ jobsetrel + atom_rels \
        + expt_rels + bj_rels + ref_rels + dft_rels + fs_rels + fd_rels \
        + fc_rels + fnlc_rels

def new_model() -> Model:
    m = Model('functionals')
    m.add(objs); m.add(rels)
    return m
