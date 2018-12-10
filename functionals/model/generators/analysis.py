# External Modules
from json   import loads

# Internal Modules
from dbgen import (Model, Generator as Gen, Expr, AS, AND, GROUP_CONCAT, SUM, COUNT,
                   Literal, NOT, NULL, CONCAT, ABS, MIN,
                   SimpleFunc, Unpack, defaultEnv, Env,Import)

from functionals.scripts.load.get_kptden             import get_kptden
from functionals.scripts.load.get_kpts_gpaw          import get_kpts_gpaw
from functionals.scripts.misc.eos_func               import eos_func
from functionals.scripts.misc.conventional_lattice   import conventional_lattice

################################################################################
################################################################################
################################################################################

def make_volume(val:str,sym:str)->float:
    sg = int(sym.split('_')[-1])
    if sg < 168: raise ValueError
    elif sg < 195:
        return float(val)**2 * 1.633
    else:
        return float(val)**3

def get_near_min(dv:float,all_vols:str)->int:
    '''
    Given a deviation from optimum volume (and a list of all volumes), return
    whether or not the given deviation is in the following subset:
        [..., opt - 2, opt - 1, opt, opt + 1, opt + 2, ...]
    '''

    vols  = sorted(loads(all_vols))
    n     = len(vols)
    m     = vols.index(float(dv))
    if n % 2 == 1:
        min_n = (n-1)//2
    else:
        mid   = n//2
        lo,hi = mid - 1, mid + 1
        min_n = lo if abs(vols[lo]) < abs(vols[hi]) else hi

    return int(m in range(min_n - 2, min_n + 3))

################################################################################
################################################################################
################################################################################

def analysis(mod : Model) -> None:
    # Extract tables
    #---------------
    tabs = ['job','atom','element','struct','calc','cell',
            'species','expt','bulk_job','reference','species_dataset',
            'species_comp','species_dataset_element','dft_data']

    Job, Atom, Element, Struct, Calc, Cell, Species, Expt, Bulk_job,\
    Reference, Species_dataset,Species_comp,Species_dataset_element,\
    Dft_data = map(mod.get, tabs) # type: ignore

    # Abbreviations
    #--------------
    D   = Dft_data
    SDE = Species_dataset_element
    BJ  = Bulk_job.Job

    ########################################################################

    bulkexpt =                                                                  \
        Gen(name = 'bulkexpt',
            desc = 'All pairs of (bulk) species + calc., if mBEEF calc'\
                    'Links the bulk_job table to Expt table, too',
            targets = [Bulk_job.Expt, Expt],
            basis   = 'Bulk_job',
            input   = [Struct.n_atoms,
                 *[e(Job.Calc)   for e in Calc.select(Calc._inits)],
                 *[e(Job.Struct) for e in Species.select(Species._inits)]],
            constr = Calc.xc == 'mBEEF')
    ########################################################################
    aggs = ['contribs','energies','volumes']

    agg_dict = dict(contribs = Job.contribs(Bulk_job.Expt),
                    energies = Job.energy(Bulk_job.Expt),
                    volumes  = Cell.volume(Bulk_job.Expt))

    pop_aggs =                                                                  \
        Gen(name    = 'pop_aggs',
            desc    = 'Store information about closest 5 singlepoints near optimum',
            targets = Expt.select(aggs) ,
            input   = [CONCAT('[',GROUP_CONCAT(y),']') |AS| x
                            for x,y in agg_dict.items()],
            links   = Bulk_job.Expt,
            constr  = Bulk_job.near_min,
            agg     = 'expt')

    ########################################################################

    all_vols =                                                                  \
        Gen(name    = 'all_vols',
            desc    = 'Store all volumes',
            targets = Expt.all_vols,
            input   = CONCAT('[',GROUP_CONCAT(Bulk_job.dv),']')
                            |AS| 'all_volumes',
            links   = Bulk_job.Expt,
            agg     = 'Expt' )
    ########################################################################

    near_min =                                                                  \
        Gen(name = 'near_min',
            desc = '?',
            targets = Bulk_job.near_min,
            input   = [Bulk_job.dv,Expt.all_vols],
            func    = SimpleFunc(get_near_min,['dv','all_vols'],['near_min']))

    ########################################################################
    Num = Species_comp.num
    Eng  = Reference.energy(multi=[[Reference.Calc], [Reference.Element]])
    eform =                                                                 \
        Gen(name    = 'eform',
            desc    = 'Difference between a relaxed energy (per atom) and '\
                      'the sum of reference energies (divided by # of atoms)',
            targets = Expt.eform,
            input   = (Expt.energy_pa - SUM(Eng* Num)/ SUM(Num)) |AS| 'eform',
            basis   =  'Expt',
            constr  =  NOT(NULL(Job.contribs(Reference.Job))), #Edge case - what if job has an energy but x_contribs failed?
            agg     = 'expt',
            aggconst= COUNT(Literal(1)) == Species.n_elems,
            links   = [Species_comp.Species,
                     Reference.Calc,
                     Reference.Element])

    ########################################################################

    ks    = ['kx','ky','kz','kptden_x','kptden_y','kptden_z']

    k_env = defaultEnv + Env(ase=['Atoms'])
    p_env = defaultEnv + Env([Import('dbgen.core.parsing',objs=['parse_line'])])

    kpts =                                                                  \
        Gen(name    = 'kpts',
            desc    = 'Kpoints and kpoint density for each job',
            targets = Job.select(ks),
            input   = [Job.log]+Cell.select(Cell._inits),
            func    = SimpleFunc((get_kptden, k_env),
                                   inputs = ks[:3]+Cell._inits,
                                   outputs= ks[-3:])
                    + SimpleFunc((get_kpts_gpaw,p_env),['log'],['kpts'])
                    + Unpack('kpts',['kx','ky','kz']))

    ########################################################################

    eos_cols = ['energy_pa', 'bulkmod', 'volume_pa', 'img', 'n_data']
    eos_env  = Env([Import('string',objs=['ascii_lowercase']),
                    Import('random',objs=['choices']),
                    Import('os',objs=['environ','remove']),
                    Import('os.path',objs=['join']),
                    Import('base64',objs=['b64encode']),
                    Import('ase.eos',objs=['EquationOfState']),
                    Import('numpy',objs=['polyfit','poly1d','mean','std','array'])
                    ]) + defaultEnv
    eos =                                                                   \
        Gen(name    = 'eos',
            desc    = 'Uses ASE equation of state over aggregated volumes '
                      'and energies for an experiment',
            targets = Expt.select(eos_cols),
            input   = [GROUP_CONCAT(Cell.volume) |AS| 'v',
                       GROUP_CONCAT(Job.energy)  |AS| 'e',
                       Expt.n_atoms,
                       COUNT(Literal(1)) |AS| 'n_data'],
            func    = SimpleFunc((eos_func, eos_env),
                                   inputs  = ['v','e','n_atoms'],
                                   outputs = ['volume_pa','energy_pa','bulkmod','img']),
            basis   = 'Bulk_job',
            agg     = 'Expt')

    ########################################################################

    lattice =                                                                   \
        Gen(name    = 'lattice',
            desc    = 'Computes lattice constant of conventional unit cell',
            targets =  Expt.lattice,
            input   = [Species.symmetry,Expt.volume_pa],
            func    = conventional_lattice)

    ########################################################################

    dcols = ['coefs','name','composition',
             'energy_vector','volume_vector','contrib_vector']
    cohesive =                                                              \
        Gen(name    = 'cohesive',
            desc    = '?',
            targets = [D,*D.select(dcols)],
            input   = [Calc.coefs,
                         Species.nickname |AS| 'name',
                         Species.composition,
                         Expt.energies |AS| 'energy_vector',
                         Expt.volumes  |AS| 'volume_vector',
                         Expt.contribs |AS| 'contrib_vector'],
            basis   = 'expt',
            constr  = (Calc.xc == 'mBEEF')
                        |AND| NOT(NULL(Expt.eform))
                        |AND| (Expt.n_data > 5))

    ########################################################################

    sde_dict = {'cohesive energy' : D.expt_cohesive_energy,
                'bulk modulus'    : D.expt_bm}

    cohesive_target,bm_target = [
        Gen(name    = x,
            desc    = 'Copying %s info from Keld data into dataset'%x,
            targets = y,
            input   = SDE.value,
            links   = SDE.Species,
            constr  = SDE.property == x)
        for x,y in sde_dict.items()]

    ########################################################################

    vol_target =                                                                \
        Gen(name    = 'vol_target',
            desc    = '?',
            targets = D.expt_volume ,
            input   = [SDE.value, Species.symmetry],
            links   = SDE.Species,
            func    = SimpleFunc(make_volume,
                                   ['value','symmetry'],
                                   ['expt_volume']),
            constr  = (SDE.property == 'lattice parameter'))


    ########################################################################

    def make_dict(attr : str, name : str) -> Expr:
        '''
        Makes a serialized python dictionary mapping element_ids to Job
        properties of their reference calculation
        '''
        elemattr = Element.atomic_number(Species_comp.Element)
        jobattr  = Job._get(attr)(multi=[[Reference.Element],[Reference.Calc]])
        group    = GROUP_CONCAT(CONCAT(elemattr,':',jobattr))
        return CONCAT('{',group,'}') |AS| name

    cohesive_elems =                                                        \
        Gen(name    = 'cohesive_elems',
            desc    = '?',
            targets = [D.atomic_contribs,D.atomic_energies],
            input   = [make_dict('contribs','atomic_contribs'),
                         make_dict('energy','atomic_energies')],
            links   = [Species_comp.Species,
                         Reference.Element,
                         Reference.Calc],
            agg     = D._name)
    ########################################################################

    vol_pa = Cell.volume / Struct.n_atoms

    gap =                                                                       \
        Gen(name    = 'gap',
            desc    = 'Computes deviation from optimum for each Bulk_job '
                      'singlepoint',
            targets = [Bulk_job.dv,Bulk_job.gap],
            input   = [vol_pa - Expt.volume_pa      |AS| 'dv',
                       ABS(vol_pa - Expt.volume_pa) |AS| 'gap'])

    ########################################################################

    mingap =                                                                    \
        Gen(name    = 'mingap',
            desc    = 'Stores a value under "Expt" unique to the best corresponding Bulk_job' ,
            targets = Expt.min_gap ,
            input   = MIN(Bulk_job.gap) |AS| 'min_gap',
            links   = Bulk_job.Expt,
            agg     = 'expt')


    ########################################################################

    ratio = Struct.n_atoms(BJ) / Species.n_atoms

    cohesive_optjob =                                                       \
        Gen(name    = 'cohesive_optjob',
            desc    = 'Populates FK using the data populated by generators '
                       '"gap" and "mingap"',
            targets = [D.Job, D.bulk_energy, D.bulk_contribs, D.bulk_ratio],
            input   = [Job.energy(BJ)   |AS| 'bulk_energy',
                         Job.contribs(BJ) |AS| 'bulk_contribs',
                         ratio            |AS| 'bulk_ratio'],
            links   = Bulk_job.Expt,
            constr  = Bulk_job.gap == Expt.min_gap)


    gens = [bulkexpt, pop_aggs, all_vols, near_min, eform, kpts, eos, lattice,
            cohesive, cohesive_target, bm_target, vol_target, cohesive_elems,
            gap, mingap, cohesive_optjob,]

    mod.add(gens)
