# External Modules
from typing import Type
from json   import loads

# Internal Modules
from dbgen import (Model,CONST, DESC, INPUT, FUNC, GET,
                   BASIS, LINKS, AGG, AGGCONST,
                   Expr, AS, AND, GROUP_CONCAT, SUM, COUNT,
                   Literal, NOT, NULL, CONCAT, ABS, MIN,
                   SimpleFunc, Unpack)

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

def analysis(mod:Type[Model]) -> None:
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

    with mod:
        ########################################################################

        bulkexpt =                                                              \
            ((Bulk_job.Expt, Expt, ) ==
                GET /BASIS/ Bulk_job
                    /INPUT/ [Struct.n_atoms,
                             *[e(Job.Calc)   for e in Calc.select(Calc._inits)],
                             *[e(Job.Struct) for e in Species.select(Species._inits)]]
                    /CONST/ (Calc.xc == 'mBEEF')
                    /DESC/ 'All pairs of (bulk) species + calc., if mBEEF calc'\
                           'Links the bulk_job table to Expt table, too')
        ########################################################################
        aggs = ['contribs','energies','volumes']

        agg_dict = dict(contribs = Job.contribs(Bulk_job.Expt),
                        energies = Job.energy(Bulk_job.Expt),
                        volumes  = Cell.volume(Bulk_job.Expt))

        pop_aggs =                                                              \
            (Expt.select(aggs) ==
                GET /INPUT/ [CONCAT('[',GROUP_CONCAT(y),']') |AS| x
                                for x,y in agg_dict.items()]
                    /LINKS/ Bulk_job.Expt
                    /DESC/ 'Store information about closest 5 singlepoints near optimum'
                    /CONST/ Bulk_job.near_min
                    /AGG/   Expt )
        ########################################################################

        all_vols =                                                              \
            (Expt.all_vols ==
                GET /INPUT/ (CONCAT('[',GROUP_CONCAT(Bulk_job.dv),']')
                                |AS| 'all_volumes')
                    /LINKS/ Bulk_job.Expt
                    /DESC/ 'Store all volumes'
                    /AGG/   Expt )
        ########################################################################

        near_min =                                                              \
            (Bulk_job.near_min ==
                GET /INPUT/ [Bulk_job.dv,Expt.all_vols]
                    /FUNC/ SimpleFunc(get_near_min,['dv','all_vols'],['near_min'])
                    )

        ########################################################################

        eform =                                                                 \
            (Expt.eform ==
                GET /INPUT/ ((Expt.energy_pa -
                                SUM(Reference.energy(multi=[[Reference.Calc],
                                                            [Reference.Element]])
                                    * Species_comp.num)
                                    / SUM(Species_comp.num))
                              |AS| 'eform')
                    /BASIS/ Expt
                    /CONST/ NOT(NULL(Job.contribs(Reference.Job))) #Edge case - what if job has an energy but x_contribs failed?
                    /AGG/   Expt
                    /AGGCONST/ (COUNT(Literal(1))==Species.n_elems)
                    /LINKS/ [Species_comp.Species,
                             Reference.Calc,
                             Reference.Element]
                    /DESC/  'Difference between a relaxed energy (per atom) and '\
                            'the sum of reference energies (divided by # of atoms)')

        ########################################################################

        ks   = ['kx','ky','kz','kptden_x','kptden_y','kptden_z']

        kpts =                                                                  \
            (Job.select(ks) ==
                GET /INPUT/ ([Job.log]+Cell.select(Cell._inits))
                    /FUNC/  (SimpleFunc(get_kptden,
                                       inputs = ks[:3]+Cell._inits,
                                       outputs= ks[-3:])
                           + SimpleFunc(get_kpts_gpaw,['log'],['kpts'])
                           + Unpack('kpts',['kx','ky','kz']))
                    /DESC/ 'Kpoints and kpoint density for each job')

        ########################################################################

        eos_cols = ['energy_pa', 'bulkmod', 'volume_pa', 'img', 'n_data']

        eos =                                                                   \
            (Expt.select(eos_cols) ==
                GET /INPUT/ [GROUP_CONCAT(Cell.volume) |AS| 'v',
                             GROUP_CONCAT(Job.energy)  |AS| 'e',
                             Expt.n_atoms,
                             COUNT(Literal(1)) |AS| 'n_data']
                    /FUNC/  SimpleFunc(eos_func,
                                       inputs  = ['v','e','n_atoms'],
                                       outputs = ['volume_pa','energy_pa','bulkmod','img'])
                    /BASIS/ Bulk_job
                    /AGG/   Expt
                    /DESC/ 'Uses ASE equation of state over aggregated volumes '
                           'and energies for an experiment')

        ########################################################################

        lattice =                                                               \
            (Expt.lattice ==
                GET /INPUT/ [Species.symmetry,Expt.volume_pa]
                    /FUNC/  conventional_lattice
                    /DESC/ 'Computes lattice constant of conventional unit cell')

        ########################################################################

        dcols = ['coefs','name','composition',
                 'energy_vector','volume_vector','contrib_vector']
        cohesive =                                                              \
            ((D,D.select(dcols)) ==
                GET /INPUT/ [Calc.coefs,
                             Species.nickname |AS| 'name',
                             Species.composition,
                             Expt.energies |AS| 'energy_vector',
                             Expt.volumes  |AS| 'volume_vector',
                             Expt.contribs |AS| 'contrib_vector']
                    /BASIS/ Expt
                    /CONST/ (      (Calc.xc == 'mBEEF')
                             |AND| NOT(NULL(Expt.eform))
                             |AND| (Expt.n_data > 5)))

        ########################################################################

        sde_dict = {'cohesive energy' : D.expt_cohesive_energy,
                    'bulk modulus'    : D.expt_bm}

        cohesive_target,bm_target = [
            y == GET /INPUT/ SDE.value
                     /LINKS/ SDE.Species
                     /DESC/ ('Copying %s info from Keld data into dataset'%x)
                     /CONST/ ((SDE.property == x))
            for x,y in sde_dict.items()]

        ########################################################################

        vol_target =                                                            \
            (D.expt_volume ==
                GET /INPUT/ [SDE.value, Species.symmetry]
                    /LINKS/ SDE.Species
                    /FUNC/  SimpleFunc(make_volume,
                                       ['value','symmetry'],
                                       ['expt_volume'])
                    /CONST/ ((SDE.property == 'lattice parameter')))


        ########################################################################

        def make_dict(attr : str, name : str) -> Expr:
            '''
            Makes a serialized python dictionary mapping element_ids to Job
            properties of their reference calculation
            '''
            elemattr= Element.atomic_number(Species_comp.Element)
            jobattr = Job._get(attr)(multi=[[Reference.Element],[Reference.Calc]])
            group   = GROUP_CONCAT(CONCAT(elemattr,':',jobattr))
            return CONCAT('{',group,'}') |AS| name

        cohesive_elems =                                                        \
            ((D.atomic_contribs,D.atomic_energies) ==
                GET /INPUT/ [make_dict('contribs','atomic_contribs'),
                             make_dict('energy','atomic_energies')]
                    /LINKS/ [Species_comp.Species,
                             Reference.Element,
                             Reference.Calc]
                    /AGG/ D)
        ########################################################################

        vol_pa = Cell.volume / Struct.n_atoms

        gap =                                                                   \
            ((Bulk_job.dv,Bulk_job.gap) ==
                GET /INPUT/ [vol_pa - Expt.volume_pa      |AS| 'dv',
                             ABS(vol_pa - Expt.volume_pa) |AS| 'gap']
                    /DESC/ 'Computes deviation from optimum for each Bulk_job '
                            'singlepoint')

        ########################################################################

        mingap =                                                                \
            (Expt.min_gap == GET /INPUT/ (MIN(Bulk_job.gap) |AS| 'min_gap')
                                 /LINKS/ Bulk_job.Expt
                                 /AGG/   Expt
                                 /DESC/ 'Stores a value under "Expt" unique to '
                                         'the best corresponding Bulk_job ')

        ########################################################################

        ratio = Struct.n_atoms(BJ) / Species.n_atoms

        cohesive_optjob =                                                       \
            ((D.Job, D.bulk_energy, D.bulk_contribs, D.bulk_ratio) ==
                GET /INPUT/ [Job.energy(BJ)   |AS| 'bulk_energy',
                             Job.contribs(BJ) |AS| 'bulk_contribs',
                             ratio            |AS| 'bulk_ratio']
                    /LINKS/ Bulk_job.Expt
                    /CONST/ (Bulk_job.gap == Expt.min_gap)
                    /DESC/  'Populates FK using the data populated by generators '
                            '"gap" and "mingap"')
