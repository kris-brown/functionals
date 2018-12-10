# External Modules
from typing     import Tuple as T, List as L
from re         import search
from ase.data   import chemical_symbols # type: ignore
from ast        import literal_eval

# Internal Modules
from dbgen import (Model, Generator as Gen, SimpleFunc,SimplePipe, AS,AND,
                   Env, defaultEnv, Import)

from functionals.scripts.io.anytraj             import anytraj
from functionals.scripts.load.find_setups       import find_setups
from functionals.scripts.io.metadata            import metadata
from functionals.scripts.atoms.get_atoms        import get_atoms
from functionals.scripts.atoms.get_cell         import get_cell
from functionals.scripts.atoms.get_system_type  import get_system_type
from functionals.scripts.atoms.json_to_traj     import json_to_traj
from functionals.scripts.atoms.get_bulk         import get_bulk
from functionals.scripts.atoms.get_pure_struct  import get_pure_struct
from functionals.scripts.atoms.cell_info        import cell_info
from functionals.scripts.atoms.countatm         import countatm

from catalysis_model.scripts.Pure.Atoms.traj_to_json import traj_to_json

##############################################################################
##############################################################################
##############################################################################

def species_nick(comp:str, sym:str)->str:
    '''Create a nickname for a species'''

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

    elems = ''.join([chemical_symbols[int(e)]+(str(num) if num>1 else '')
                        for e, num in literal_eval(comp).items()])
    return elems+'_'+nick_dict.get(sym,'')

eng_env = Env(re=['search'])
def eng(s:str)->float:
    '''Parse the result free energy from a GPAW logfile'''
    pat = r'Free energy:\s+([-+]?\d+\.\d+)'
    match = search(pat, s); assert match
    return float(match.groups()[0])

##############################################################################
##############################################################################
##############################################################################

def load(mod:Model)->None:
    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell','pure_struct',
            'species','bulk_job','reference','setup','setup_family',
            'job_setup','species_comp','species_dataset',
            'species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Pure_struct, Species, Bulk_job,\
    Reference, Setup, Setup_family, Job_setup, Species_comp, Species_dataset,         \
    Species_dataset_element = map(mod.get, tabs) # type: ignore

    ########################################################################
    fs_env = defaultEnv + Env(re=['findall'])
    jobsetup =                                                                  \
        Gen(name    = 'jobsetup',
            desc    = 'Populate mapping table between Job and Setup',
            targets = [Job_setup, Setup],
            input   =  Job.log,
            basis   = 'job',
            func    =  SimpleFunc((find_setups,fs_env),
                                ['log'],['checksum']))
    ########################################################################

    psnick =                                                                    \
        Gen(name    = 'psnick',
            desc    = 'Use a dict to map nicknames to some Pure structs',
            targets = Pure_struct.nickname,
            input   = Pure_struct.prototype,
            func    = (lambda x: {'AB_1_a_b_225'       : 'rocksalt',
                                    'AB_1_a_c_216'       : 'zincblende',
                                    'AB_1_a_b_221'       : 'cesium chloride',
                                    'A_1_a_225'          : 'fcc',
                                    'A_1_a_229'          : 'bcc',
                                    'A_2_c_194'          : 'hcp',
                                    'A_2_a_227'          : 'diamond',
                                    'AB3_1_a_d_221'      : 'anti-ReO3',
                                    'A2B3_8_ad_e_206'    : 'antibixbyite',
                                    'AB3C_cP5_221_a_c_b' : 'perovskite'
                            }.get(x)),
            tags  =  'species')

    ########################################################################
    anytraj_env = Env([Import('ase.io',objs=['read']),
                       Import('glob',objs=['glob']),
                       Import('os.path',objs=['getsize'])])
    struct =                                                                    \
        Gen(name    = 'struct',
            desc    = "Assumes there's only one traj in the directory",
            targets = [Job.Struct, Struct],
            input   = Job.stordir,
            func    = SimplePipe([(anytraj, anytraj_env),traj_to_json],
                                    ['stordir'],['raw']))

    ########################################################################

    jobeng = Gen(name    = 'jobeng',
                 desc    = '?',
                 targets = Job.energy,
                 input   = Job.log,
                 func    = (eng, eng_env))

    ########################################################################

    md_env = Env([Import('os',     objs=['stat']),
                  Import('pwd',    objs=['getpwuid']),
                  Import('os.path',objs=['getmtime'])]) + defaultEnv
    jobmetadata =                                                           \
        Gen(name    = 'jobmetadata',
            desc    = 'Scrape metadata about job',
            targets = [Job.user, Job.timestamp],
            input   = Job.stordir,
            func    = (metadata, md_env))

    ########################################################################

    atomdetails = ['number','x','y','z','constrained','magmom','tag']
    atom =                                                                  \
        Gen(name    = 'atom',
            desc    = 'Initializes atoms from knowing n_atoms',
            targets = [Atom, Atom.Element, *Atom.select(atomdetails)],
            input   = Struct.raw,
            func    = SimpleFunc(get_atoms,
                                    ['raw'],
                                    atomdetails + ['ind','atomic_number']))

    ########################################################################

    elems =                                                                 \
        Gen(name    = 'elems',
            desc    = 'Seeds the Element table by providing atomic numbers',
            targets = Element,
            func    = (lambda:list(range(1,100))))


    ########################################################################

    cell =                                                                  \
        Gen(name    = 'cell',
            desc    = 'Populate cells from Structs',
            targets = [Cell, Struct.Cell],
            input   = Struct.raw,
            func    = get_cell)

    ########################################################################

    celinfo =                                                               \
        Gen(name    = 'celinfo',
            desc    = 'Basic geometric formulas applied to 3x3 cell '
                       'representation',
            targets = Cell.select('surface_area','volume','a','b','c'),
            input   = Cell.select(Cell._inits),
            func    = cell_info)

    ########################################################################
    jtt_env = Env([Import('ase.constraints',objs=['FixAtoms']),
                   Import('ase',objs=['Atoms'])]) + defaultEnv

    jtt = (json_to_traj,jtt_env)
    systype =                                                               \
        Gen(name = 'systype',
            desc = 'Apply function to identify whether system is '
                     'bulk/mol/surf',
            targets = Struct.system_type,
            input   =  Struct.raw,
            func    =  SimplePipe([jtt, get_system_type],
                                  inputs  = ['raw'],
                                  outputs = ['system_type']))

    ########################################################################
    gb_env = Env([Import('bulk_enumerator.bulk', objs = ['BULK']),
                  Import('ase.io', objs=['write']),
                  Import('os',     objs=['remove','environ']),
                  Import('os.path',objs=['join']),
                  Import('random', objs=['choices']),
                  Import('string', objs=['ascii_lowercase'])])
    ps   =                                                                  \
        Gen(name    = 'ps',
            desc    = 'Determine the "pure structure" code for each Bulk',
            targets = [Pure_struct, Struct.Pure_struct],
            input   = Struct.raw,
            func    = SimplePipe([jtt, (get_bulk, gb_env), get_pure_struct],
                                        inputs  = ['raw'],
                                        outputs = ['prototype']),
            tags    =  ['long','parallel'],
            constr  = Struct.system_type == 'bulk')
    ########################################################################
    le_env = defaultEnv + Env(ast=['literal_eval'])
    cs_env = Env([Import('ase.data',objs=['chemical_symbols'])])

    spnick =                                                                \
        Gen(name    = 'spnick',
            desc    = "Determine species' nicknames by analyzing their "\
                       "composition and symmetry" ,
            targets = Species.nickname,
            input   = [Species.composition, Species.symmetry],
            func    = SimpleFunc((species_nick,le_env + cs_env),
                                  ['composition','symmetry'],
                                  ['nickname']))
    ########################################################################
    def getNum(d:str) -> T[L[int], L[int]]:
        elems, counts = zip(*literal_eval(d).items())
        return list(elems), list(counts)


    sp_comp =                                                                  \
        Gen(name    = 'sp_comp',
            desc    = 'Populates species composition table',
            targets = Species_comp ,
            input   = Species.composition,
            func    = SimpleFunc((getNum,le_env),
                                   ['composition'],
                                   ['atomic_number','num']))
    ########################################################################
    pop_species =                                                               \
        Gen(name    = 'pop_species',
            desc    = 'Populate species from Struct and Purestruct tables',
            targets = [Species, Struct.Species, Species.n_elems],
            input   = [Struct.composition_norm |AS| 'composition',
                         Pure_struct.prototype   |AS| "symmetry",
                         Struct.n_elems])
    ########################################################################
    def sum_values(x:str)->int:
        return sum(dict(literal_eval(x)).values())

    species_natoms =                                                            \
        Gen(name    = 'species_natoms',
            desc    = 'Total number of atoms in the most reduced stoichiometric ratios',
            targets = Species.n_atoms,
            input   = Species.composition,
            func    = (sum_values, le_env))
    ########################################################################
    refr =                                                                    \
        Gen(name    = 'refr',
            desc    = 'Collects calculations of isolated atoms',
            targets = [Reference, Reference.Calc, Reference.Element, Reference.energy],
            input   = [Calc.xc,
                         Element.atomic_number(Atom.Element),
                         Job.energy],
            basis   = 'job',
            links   =  Atom.Struct,
            constr   =  ((Struct.system_type == 'molecule')
                         |AND| (Struct.n_atoms == 1)) )

    ########################################################################
    bulkjob =                                                               \
        Gen(name    = 'bulkjob',
            desc    = "Subset of jobs that are bulk calculations",
            targets = Bulk_job,
            basis   = 'job',
            constr  = (Struct.system_type == 'bulk'))
    ########################################################################
    spin     = lambda x: int('Spin-polarized calculation' in x)

    spinpol =                                                               \
        Gen(name    = 'spinpol',
            desc    = 'Populates Job.spinpol',
            targets = Job.spinpol,
            input   = Job.log,
            func    = spin)
    ########################################################################
    compcols = Struct.select('n_atoms','n_elems','composition',
                     'composition_norm','metal_comp','str_symbols',
                     'str_constraints')
    ca_env = Env([Import('dbgen.core.lists',objs=['normalize_list','nub'])])+defaultEnv
    countatoms =                                                            \
        Gen(name    = 'countatoms',
            desc    = 'Get stoichiometric properties of Structure from raw',
            targets = compcols,
            input   =  Struct.raw,
            func    = (countatm, ca_env))


    gens = [jobsetup,psnick,struct,jobeng,jobmetadata,atom,elems,cell,celinfo,
            systype,ps,spnick,sp_comp,pop_species,species_natoms,refr,bulkjob,
            spinpol, countatoms]

    mod.add(gens)
