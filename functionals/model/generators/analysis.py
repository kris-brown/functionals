# External Modules
from typing import Any, Type, TYPE_CHECKING

# Internal Modules
if TYPE_CHECKING:
    from dbgen.support.model     import Model

from dbgen.support.get      import (CONST, DESC, INPUT, FUNC, GET, TAGS,
                                     BASIS, LINKS, AGG, AGGCONST)
from dbgen.support.expr     import (IN, AS, AND, GROUP_CONCAT, SUM, COUNT,
                                     Literal)
from dbgen.support.funclike import SimpleFunc, Unpack, SimplePipe

from functionals.scripts.get_kptden        import get_kptden
from functionals.scripts.get_kpts_gpaw     import  get_kpts_gpaw
from functionals.scripts.eos_func          import eos_func
################################################################################
def analysis(mod:Type['Model'])->None:
    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell',
            'species','expt','bulk_job','reference',
            'species_comp','species_dataset',
            'species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Species, Expt, Bulk_job,\
    Reference, Species_comp, Species_dataset,         \
    Species_dataset_element = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################
        eform =                                                                 \
            (Expt.eform ==
                GET /INPUT/ ((Expt.energy -
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
            ((Expt.energy, Expt.bulkmod, Expt.volume, Expt.img) ==
                GET /INPUT/ [GROUP_CONCAT(Cell.volume) |AS| 'v',
                             GROUP_CONCAT(Job.energy)  |AS| 'e',
                             Struct.n_atoms]
                    /FUNC/  SimpleFunc(eos_func,['v','e','n_atoms'],
                                       ['volume','energy','bulkmod','img'])
                    /BASIS/ Bulk_job
                    /AGG/   Expt
                    /DESC/ 'Uses ASE equation of state over aggregated volumes '\
                           'and energies for an experiment')
        ########################################################################
        # This takes time and really isn't used for anything? nice to have later when not nuking so frequently
        # sg =                                                                    \
        #     (Struct.sg ==
        #         GET /INPUT/ Struct.raw
        #             /FUNC/  SimplePipe([json_to_traj, get_spacegroup],['raw'],['sg'])
        #             /CONST/ (Struct.system_type |IN| ["bulk","surface"])
        #             /DESC/  'Computes spacegroup of a bulk')
