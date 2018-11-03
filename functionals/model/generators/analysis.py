# External Modules
from typing import Any, Type, TYPE_CHECKING

# Internal Modules
if TYPE_CHECKING:
    from dbgen.support.model     import Model

from dbgen.support.get      import (CONST, DESC, INPUT, FUNC, GET, TAGS,
                                     BASIS, LINKS, AGG, AGGCONST)
from dbgen.support.expr     import (Expr, IN, AS, AND, GROUP_CONCAT, SUM, COUNT,
                                     Literal, NOT, NULL, CONCAT, ABS, MIN)
from dbgen.support.funclike import SimpleFunc, Unpack, SimplePipe

from functionals.scripts.get_kptden             import get_kptden
from functionals.scripts.get_kpts_gpaw          import  get_kpts_gpaw
from functionals.scripts.eos_func               import eos_func
from functionals.scripts.conventional_lattice   import conventional_lattice

################################################################################
def analysis(mod:Type['Model']) -> None:
    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell',
            'species','expt','bulk_job','reference','species_dataset',
            'species_comp','species_dataset_element','cohesive_data']

    Job, Atom, Element, Struct, Calc, Cell, Species, Expt, Bulk_job,\
    Reference, Species_dataset,Species_comp,Species_dataset_element,\
    Cohesive_data = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################
        bulkexpt =                                                              \
            ((Bulk_job.Expt, Expt) ==
                GET /BASIS/ Bulk_job
                    /INPUT/ [Struct.n_atoms,
                             *[e(Job.Calc)   for e in Calc.select(Calc._inits)],
                             *[e(Job.Struct) for e in Species.select(Species._inits)]]
                    /DESC/ 'All pairs of (bulk) species + calc. '\
                           'Links the bulk_job table to Expt table, too')
        ########################################################################
        eform =                                                                 \
            (Expt.eform ==
                GET /INPUT/ ((Expt.energy_pa -
                                SUM(Reference.energy(multi=[[Reference.Calc],[Reference.Element]])
                                    * Species_comp.num)
                                    / SUM(Species_comp.num))
                              |AS| 'eform')
                    /BASIS/ Expt
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
        eos =                                                                   \
            (Expt.select('energy_pa', 'bulkmod', 'volume_pa', 'img', 'n_data') ==
                GET /INPUT/ [GROUP_CONCAT(Cell.volume) |AS| 'v',
                             GROUP_CONCAT(Job.energy)  |AS| 'e',
                             Expt.n_atoms,
                             COUNT(Literal(1)) |AS| 'n_data']
                    /FUNC/  SimpleFunc(eos_func,
                                       inputs  = ['v','e','n_atoms'],
                                       outputs = ['volume_pa','energy_pa','bulkmod','img'])
                    /BASIS/ Bulk_job
                    /AGG/   Expt
                    /DESC/ 'Uses ASE equation of state over aggregated volumes '\
                           'and energies for an experiment')

        ########################################################################
        lattice =                                                               \
            (Expt.lattice == GET /INPUT/ [Species.symmetry,Expt.volume_pa]
                                 /FUNC/ conventional_lattice
                                 /DESC/ 'Computes lattice constant of conventional unit cell')
        ########################################################################
        CD  = Cohesive_data
        SDE = Species_dataset_element
        cohesive =                                                              \
            ((CD,CD.coefs,CD.name,CD.composition) ==
                GET /INPUT/ [Calc.coefs,
                             Species.nickname |AS| 'name',
                             Species.composition]
                    /BASIS/ Expt
                    /CONST/ (      (Calc.xc == 'mBEEF')
                             |AND| NOT(NULL(Expt.eform))
                             |AND| (Expt.n_data > 5)))

        ########################################################################

        cohesive_target =                                                       \
            (CD.target ==
                GET /INPUT/ SDE.value
                    /LINKS/ SDE.Species
                    /CONST/ ((SDE.property == 'cohesive energy')))
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
            ((CD.atomic_contribs,CD.atomic_energies) ==
                GET /INPUT/ [make_dict('contribs','atomic_contribs'),
                             make_dict('energy','atomic_energies')]
                    /LINKS/ [Species_comp.Species,
                             Reference.Element,
                             Reference.Calc]
                    /AGG/ CD)
        ########################################################################
        ########################################################################
        vol_pa = Cell.volume / Struct.n_atoms
        gap =                                                                   \
            (Bulk_job.gap ==
                GET /INPUT/ (ABS(vol_pa - Expt.volume_pa) |AS| 'gap')
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
        BJ = Bulk_job.Job
        ratio = Struct.n_atoms(BJ) / Species.n_atoms

        cohesive_optjob =                                                       \
            ((CD.Job, CD.bulk_energy, CD.bulk_contribs, CD.bulk_ratio) ==
                GET /INPUT/ [Job.energy(BJ)   |AS| 'bulk_energy',
                             Job.contribs(BJ) |AS| 'bulk_contribs',
                             ratio            |AS| 'bulk_ratio']
                    /LINKS/ Bulk_job.Expt
                    /CONST/ (Bulk_job.gap == Expt.min_gap)
                    /DESC/  'Populates FK using the data populated by generators '
                            '"gap" and "mingap"')

        ########################################################################
        # This takes time and really isn't used for anything? nice to have later when not nuking so frequently
        # sg =                                                                    \
        #     (Struct.sg ==
        #         GET /INPUT/ Struct.raw
        #             /FUNC/  SimplePipe([json_to_traj, get_spacegroup],['raw'],['sg'])
        #             /CONST/ (Struct.system_type |IN| ["bulk","surface"])
        #             /DESC/  'Computes spacegroup of a bulk')
