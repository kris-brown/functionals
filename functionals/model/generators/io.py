# External Modules
from typing     import Callable as C
from os         import environ
from os.path    import exists, join

# Internal Modules
from dbgen import (Model, Generator as Gen, SimpleFunc, Env, Import, defaultEnv)

from functionals.scripts.io.parse_setup       import parse_setup
from functionals.scripts.io.get_stray_gpaw    import get_stray_gpaw
from functionals.scripts.io.parse_mendeleev   import parse_mendeleev
from functionals.scripts.io.parse_keld        import parse_keld
from functionals.scripts.load.parse_pw_gpaw   import parse_pw_gpaw
from functionals.scripts.load.parse_xc_gpaw   import parse_xc_gpaw
from functionals.scripts.load.get_econv_gpaw  import get_econv_gpaw

##############################################################################
elempath = environ['ELEMPATH']
keldpath = environ['KELDPATH']
logpth   = '/Users/ksb/scp_tmp/auto' #
psppth   = '/Users/ksb/scp_tmp/norm_conserving_setups'

pthEnv  = Env([Import('os.path',objs=['join','exists'])])

def readfile(pth:str)->str:
    with open(pth,'r') as f: return f.read()

def getcoef(stor:str, bf:str)->str:
    def coefpath(x:str)->str:
        return join(x,'BEEFoftheDay.txt')

    def readfile(pth:str)->str:
        with open(pth,'r') as f: return f.read()

    if not bf:
        return ''
    else:
        return (readfile(coefpath(stor))
                    if exists(coefpath(stor))
                        else readfile('/Users/ksb/functionals/data/beef.json'))

parse_env = defaultEnv + Env([Import('dbgen.core.parsing',objs=['parse_line'])])

calc_plan =  SimpleFunc((parse_pw_gpaw,parse_env), ['log'],['pw'])              \
           + SimpleFunc((parse_xc_gpaw,parse_env), ['log'],['xc'])              \
           + SimpleFunc((get_econv_gpaw,parse_env),['log'],['econv'])           \
           + SimpleFunc((getcoef,pthEnv),['stordir','beef'],['coefs'])          \
           + SimpleFunc(lambda x: int('beef' in x.lower()),['xc'],['beef'])

elemcols = ['symbol', 'name', 'atomic_weight','atomic_radius'
          , 'phase','evaporation_heat', 'pointgroup','spacegroup'
          , 'melting_point', 'metallic_radius', 'vdw_radius'
          , 'density', 'en_allen' , 'is_radioactive'
          , 'lattice_struct' , 'fusion_heat'
          , 'econf', 'period', 'covalent_radius_bragg'
          , 'geochemical_class', 'abundance_crust', 'heat_of_formation'
          , 'electron_affinity', 'atomic_volume',  'boiling_point'
          , 'proton_affinity', 'covalent_radius_slater'
          , 'lattice_constant', 'dipole_polarizability'
          , 'en_ghosh', 'thermal_conductivity', 'group_id', 'en_pauling'
          , 'gas_basicity'
          ,'abundance_sea']
#################################################################################
def io(mod:Model) -> None:

    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell','pure_struct',
            'species','expt','bulk_job','reference','setup','setup_family',
            'job_setup','species_comp','species_dataset',
            'species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Pure_struct, Species, Expt, Bulk_job,\
    Reference, Setup, Setup_family, Job_setup, Species_comp, Species_dataset,    \
    Species_dataset_element = map(mod.get, tabs) # type: ignore

    ########################################################################
    gl_env = defaultEnv + Env(subprocess = ['getstatusoutput'])

    getlogs  =                                                              \
        Gen(name    = 'getlogs',
            desc    = 'Scrape folder recursively for any GPAW jobs',
            targets = Job,
            func    = SimpleFunc((get_stray_gpaw, gl_env),
                                      ['root'],['logfile']),
            consts  =  dict(root=logpth))
    ########################################################################
    getsd    = Gen(name    = 'getsd',
                   desc    = 'Get folder from logfile path',
                   targets = Job.stordir,
                   input   = Job.logfile,
                   func    = (lambda x: x[:x.rfind('/')]))
    ########################################################################
    datagpaw = Gen(name    = 'datagpaw',
                   desc    =  'Read log file from disk',
                   targets = Job.log,
                   input   = Job.logfile,
                   func    =  readfile)
    ########################################################################
    su_env = defaultEnv + Env([Import('gzip',objs=[('open','gopen')]),
                               Import('hashlib',objs=['md5']),
                               Import('os',objs=['listdir']),
                               Import('os.path',objs=['isdir','join']),
                               Import('xml.etree.ElementTree',objs=['fromstring']),
                               Import('ase.data',objs=['chemical_symbols'])
                               ])
    setups =                                                                    \
        Gen(name    = 'setups',
            desc    = 'Populate Setup table from a path containing setups',
            targets = [Setup, Setup_family, Setup.Setup_family, Setup.Element,
                        Setup.val],
            func    =  SimpleFunc((parse_setup, su_env),
                                  ['path'],
                                  ['checksum','xc','kind','name',
                                   'atomic_number','val']),
            consts  = {'path' : psppth})
    ########################################################################

    elemzinfo =                                                                 \
        Gen(name    = 'elemzinfo',
            desc    = '''Read a JSON file containing element reference data
                        Requires providing path to this file as a constant''',
            targets = Element.select(elemcols) ,
            input   = Element.atomic_number,
            func    =  SimpleFunc(parse_mendeleev,
                                      inputs  = ['atomic_number','pth'],
                                      outputs = elemcols),
            consts  = {'pth':elempath})


    ########################################################################
    # things that get populated by kelddata
    keldpop  = [Species_dataset, Species_dataset_element,
                Species_dataset_element.value,
                Species_dataset_element.datatype,
                Species, Species.n_elems]
    keldenv  = Env([Import('ase.data',objs=['chemical_symbols']),
                    Import('re',objs=['findall'])]) + defaultEnv
    kelddata =                                                              \
        Gen(name    = 'kelddata',
            desc    = 'Parses a .py Keld data file',
            targets = keldpop,
            func    = SimpleFunc((parse_keld,keldenv),
                                  inputs =['keld_path'],
                                  outputs=['dataset_name','composition',
                                           'symmetry','atomic_number','num','n_elems',
                                           'property','value','datatype']),
            consts  = {'keld_path':keldpath})
    ########################################################################
    calc =                                                                      \
        Gen(name    = 'calc',
            desc    = 'Populate calc table + F.K. from Relax_job',
            targets = [Job.Calc, Calc],
            input   =  Job.select('stordir','log'),
            basis   = 'job',
            func    = calc_plan)
    ########################################################################
    def xx(x:str)->bool:
        return exists(join(x,'xccontribs.txt'))
    has_contribs =                                                              \
        Gen(name    = 'has_contribs',
            desc    = '?',
            targets = Job.has_contribs,
            input   = Job.stordir,
            func    = (xx,pthEnv),
            constr  = (Calc.xc == 'mBEEF'))

    def get_c(x:str)->str:
        with open(join(x,'xccontribs.txt'),'r') as f:
            return f.read()
    contribs =                                                              \
        Gen(name    = 'contribs',
            desc    = '?',
            targets = Job.contribs ,
            input   = Job.stordir,
            func    =  (get_c,pthEnv),
            constr  = Job.has_contribs)

    gens = [getlogs,getsd,setups,datagpaw,kelddata,elemzinfo,calc,has_contribs,contribs]
    mod.add(gens)
