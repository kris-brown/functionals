# External Modules
from typing     import Any, Type, Tuple, List, Callable as C
from re         import search
from ase.data   import chemical_symbols # type: ignore
from ast        import literal_eval

# Internal Modules
from dbgen import (Model, CONST, DESC, INPUT, FUNC, CONSTS,
                    GET, TAGS, BASIS, LINKS, IO, OPTION, AGG,
                    AGGCONST, SimpleFunc, PyBlock, noIndex,
                    Unpack, SimplePipe, AS,AND)

from functionals.scripts.find_setups       import find_setups
from functionals.scripts.anytraj           import  anytraj
from functionals.scripts.metadata          import  metadata
from functionals.scripts.get_atoms         import  get_atoms
from functionals.scripts.get_cell          import  get_cell
from functionals.scripts.get_system_type   import get_system_type
from functionals.scripts.get_pointgroup    import get_pointgroup
from functionals.scripts.get_spacegroup    import get_spacegroup
from functionals.scripts.json_to_traj      import json_to_traj
from functionals.scripts.get_bulk          import get_bulk
from functionals.scripts.get_pure_struct   import get_pure_struct
from functionals.scripts.cell_info         import cell_info
from functionals.scripts.countatm          import countatm

##############################################################################
nick_dict = {'AB_1_a_b_225'       : 'rocksalt',
             'AB_1_a_c_216'       : 'zincblende',
             'AB_1_a_b_221'       : 'cesium chloride',
             'A_1_a_225'          : 'fcc',
             'A_1_a_229'          : 'bcc',
             'A_2_c_194'          : 'hcp',
             'A_2_a_227'          : 'diamond',
             'AB3_1_a_d_221'      : 'anti-ReO3',
             'A2B3_8_ad_e_206'    : 'antibixbyite',
             'AB3C_cP5_221_a_c_b' : 'perovskite'
            } # rutile?

def species_nick(comp:str, sym:str)->str:
    '''Create a nickname for a species'''
    elems = ''.join([chemical_symbols[int(e)]+(str(num) if num>1 else '')
                        for e, num in literal_eval(comp).items()])
    return elems+'_'+nick_dict.get(sym,'')

def eng(s:str)->float:
    '''Parse the result free energy from a GPAW logfile'''
    pat = r'Free energy:\s+([-+]?\d+\.\d+)'
    match = search(pat, s); assert match
    return float(match.groups()[0])
##############################################################################

def load(mod:Type['Model'])->None:
    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell','pure_struct',
            'species','bulk_job','reference','setup','setup_family',
            'job_setup','species_comp','species_dataset',
            'species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Pure_struct, Species, Bulk_job,\
    Reference, Setup, Setup_family, Job_setup, Species_comp, Species_dataset,         \
    Species_dataset_element = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################

        jobsetup =                                                              \
            ((Job_setup, Setup) ==
                GET /INPUT/ Job.log
                    /BASIS/ Job
                    /FUNC/ SimpleFunc(find_setups,
                                    ['log'],['checksum'])
                    /DESC/ 'Populate mapping table between Job and Setup')
        ########################################################################

        psnick =                                                                \
            (Pure_struct.nickname ==
                GET /INPUT/ Pure_struct.prototype
                    /FUNC/  (lambda x: nick_dict.get(x))
                    /TAGS/  ['species']
                    /DESC/  'Use a dict to map nicknames to some Pure structs')

        ########################################################################

        struct =                                                               \
            ((Job.Struct, Struct) ==
                GET /INPUT/ Job.stordir
                    /FUNC/  anytraj
                    /DESC/  '''Assumes there's only one traj in the directory''')

        ########################################################################

        jobeng = Job.energy == GET /INPUT/ Job.log /FUNC/ eng

        jobmetadata =                                                           \
            ((Job.user, Job.timestamp) ==
                GET /INPUT/ Job.stordir
                    /FUNC/ metadata
                    /DESC/ 'Scrape metadata about job')

        ########################################################################

        atomdetails = ['number','x','y','z','constrained','magmom','tag']
        atom =                                                                  \
            ((Atom, Atom.Element, Atom.select(atomdetails)) ==
                GET /INPUT/ Struct.raw
                    /FUNC/  SimpleFunc(get_atoms,
                                        ['raw'],
                                        atomdetails + ['ind','atomic_number'])
                    /DESC/  'Initializes atoms from knowing n_atoms')

        ########################################################################

        elems =                                                                 \
            (Element ==
                GET /FUNC/ (lambda:list(range(1,100)))
                    /DESC/ 'Seeds the Element table by providing atomic numbers')


        ########################################################################

        cell =                                                                  \
            ((Cell, Struct.Cell)
                == GET /INPUT/ Struct.raw
                       /FUNC/  get_cell
                       /DESC/  'Populate cells from Structs')

        ########################################################################

        celinfo =                                                               \
            (Cell.select('surface_area','volume','a','b','c') ==
                GET /INPUT/ Cell.select(Cell._inits)
                    /FUNC/  cell_info
                    /DESC/  'Basic geometric formulas applied to 3x3 cell '
                            'representation')

        ########################################################################
        systype =                                                               \
            (Struct.system_type ==
                GET /INPUT/ Struct.raw
                    /FUNC/ SimplePipe([json_to_traj, get_system_type],
                                      inputs  = ['raw'],
                                      outputs = ['system_type'])

                    /DESC/ 'Apply function to identify whether system is '
                             'bulk/mol/surf')

        ########################################################################
        ps   =                                                                  \
            ((Pure_struct, Struct.Pure_struct) ==
                GET /INPUT/ Struct.raw
                    /FUNC/ SimplePipe([json_to_traj, get_bulk, get_pure_struct],
                                            inputs  = ['raw'],
                                            outputs = ['prototype'])
                    /TAGS/  ['long','parallel']
                    /CONST/ (Struct.system_type == 'bulk')
                    /DESC/ 'Determine the "pure structure" code for each Bulk')
        ########################################################################

        spnick =                                                                \
            (Species.nickname ==
                GET /INPUT/ [Species.composition, Species.symmetry]
                    /FUNC/ SimpleFunc(species_nick,
                                      ['composition','symmetry'],
                                      ['nickname'])
                    /DESC/ "Determine species' nicknames by analyzing their "\
                           "composition and symmetry")
        ########################################################################
        def getNum(d:str)->Tuple[List[int], List[int]]:
            elems, counts = zip(*literal_eval(d).items())
            return list(elems), list(counts)

        sp_comp =                                                               \
            (Species_comp == GET /INPUT/ Species.composition
                                 /FUNC/ SimpleFunc(getNum,
                                                   ['composition'],
                                                   ['atomic_number','num'])
                                  /DESC/ 'Populates species composition table')
        ########################################################################
        pop_species =                                                           \
            ((Species, Struct.Species, Species.n_elems) ==
                GET /INPUT/ [Struct.composition_norm |AS| 'composition',
                             Pure_struct.prototype   |AS| "symmetry",
                             Struct.n_elems]
                    /DESC/ 'Populate species from Struct and Purestruct tables')
        ########################################################################
        def sum_values(x:str)->int:
            return sum(dict(literal_eval(x)).values())

        species_natoms =                                                        \
            (Species.n_atoms ==
                GET /INPUT/ Species.composition
                    /FUNC/ sum_values
                    /DESC/ 'Total number of atoms in the most reduced stoichiometric ratios')
        ########################################################################
        refr =                                                                   \
            ((Reference, Reference.Calc, Reference.Element, Reference.energy) ==
                GET /INPUT/ [Calc.xc,
                             Element.atomic_number(Atom.Element),
                             Job.energy]
                    /BASIS/ Job
                    /LINKS/ Atom.Struct
                    /CONST/ ((Struct.system_type == 'molecule')
                             |AND| (Struct.n_atoms == 1))
                    /DESC/ 'Collects calculations of isolated atoms')

        ########################################################################
        bulkjob =                                                               \
            ((Bulk_job) ==
                GET /BASIS/ Job
                    /CONST/ (Struct.system_type == 'bulk')
                    /DESC/  "Subset of jobs that are bulk calculations")
        ########################################################################
        spin     = lambda x: int('Spin-polarized calculation' in x)

        spinpol =                                                               \
            (Job.spinpol ==
                GET /INPUT/ Job.log
                    /FUNC/ spin
                    /DESC/ 'Populates Job.spinpol')
        ########################################################################
        compcols = Struct.select('n_atoms','n_elems','composition',
                         'composition_norm','metal_comp','str_symbols',
                         'str_constraints')
        countatoms =                                                            \
            (compcols ==
                GET /INPUT/ Struct.raw
                    /FUNC/ countatm
                    /DESC/ 'Get stoichiometric properties of Structure from raw')
