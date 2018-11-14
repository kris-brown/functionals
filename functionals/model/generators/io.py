# External Modules
from typing     import Any, Type, Tuple, List, TYPE_CHECKING, Callable as C
from os         import environ
from os.path    import exists, join

# Internal Modules
from dbgen import (Model, CONST, DESC, INPUT, FUNC, CONSTS,
                    GET, TAGS, BASIS, LINKS, IO, OPTION, AGG,
                    AGGCONST, AND, SimpleFunc, PyBlock, noIndex,
                    Unpack, SimplePipe)

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
logpth   = '/Users/ksb/scp_tmp/auto'
psppth   = '/Users/ksb/scp_tmp/norm_conserving_setups'

def readfile(pth:str)->str:
    with open(pth,'r') as f: return f.read()


beefpath = '/Users/ksb/functionals/data/beef.json'
coefpath = lambda x: join(x,'BEEFoftheDay.txt') # type: C[[str], str]
xcpath   = lambda x: join(x,'xccontribs.txt')  # type: C[[str], str]
getcoef  = lambda stor, bf: (readfile(coefpath(stor))
                            if exists(coefpath(stor))
                            else readfile(beefpath)) \
                             if bf else ''

calc_plan =  SimpleFunc(parse_pw_gpaw, ['log'],['pw'])                          \
           + SimpleFunc(parse_xc_gpaw, ['log'],['xc'])                          \
           + SimpleFunc(get_econv_gpaw,['log'],['econv'])                       \
           + SimpleFunc(getcoef,['stordir','beef'],['coefs'])                   \
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
def io(mod:Type['Model'])->None:

    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell','pure_struct',
            'species','expt','bulk_job','reference','setup','setup_family',
            'job_setup','species_comp','species_dataset',
            'species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Pure_struct, Species, Expt, Bulk_job,\
    Reference, Setup, Setup_family, Job_setup, Species_comp, Species_dataset,         \
    Species_dataset_element = map(mod.get, tabs) # type: ignore

    with mod:
        ########################################################################
        getlogs  =                                                              \
            (Job == GET /FUNC/ (lambda : get_stray_gpaw(logpth))
                        /DESC/ 'Scrape folder recursively for any GPAW jobs')
        ########################################################################
        getsd    = (Job.stordir == GET /INPUT/ Job.logfile
                                       /FUNC/ (lambda x: x[:x.rfind('/')])
                                       /DESC/ 'Get folder from logfile path')
        ########################################################################
        datagpaw = (Job.log == GET /INPUT/ Job.logfile
                                   /FUNC/ readfile
                                   /DESC/ 'Read log file from disk')
        ########################################################################
        setups =                                                                \
            ((Setup, Setup_family, Setup.Setup_family, Setup.Element,
              Setup.val) ==
                GET /FUNC/ SimpleFunc(parse_setup,
                                      ['path'],
                                      ['checksum','xc','kind','name',
                                       'atomic_number','val'])
                    /CONSTS/ {'path' : psppth}
                    /DESC/   'Populate Setup table from a path containing setups')
        ########################################################################

        elemzinfo =                                                             \
            (Element.select(elemcols) ==
                 GET /INPUT/ Element.atomic_number
                     /DESC/   '''Read a JSON file containing element reference data
                                 Requires providing path to this file as a constant'''
                     /FUNC/   SimpleFunc(parse_mendeleev,
                                          inputs  = ['atomic_number','pth'],
                                          outputs = elemcols)
                     /CONSTS/ {'pth':elempath})


        ########################################################################
        # things that get populated by kelddata
        keldpop  = (Species_dataset, Species_dataset_element,
                    Species_dataset_element.value,
                    Species_dataset_element.datatype,
                    Species, Species.n_elems)

        kelddata =                                                              \
            (keldpop ==
                GET /FUNC/ SimpleFunc(parse_keld, inputs=['keld_path'],
                                      outputs=['dataset_name','composition',
                                               'symmetry','atomic_number','num','n_elems',
                                               'property','value','datatype'])
                    /CONSTS/ {'keld_path':keldpath}
                    /DESC/  'Parses a .py Keld data file')
        ########################################################################
        calc =                                                                  \
            ((Job.Calc, Calc) ==
                GET /INPUT/  Job.select('stordir','log')
                    /BASIS/  Job
                    /FUNC/   calc_plan
                    /DESC/   'Populate calc table + F.K. from Relax_job')
        ########################################################################
        has_contribs =                                                          \
            (Job.has_contribs ==
                GET /INPUT/ Job.stordir
                    /FUNC/ (lambda x: exists(xcpath(x)))
                    /CONST/ (Calc.xc == 'mBEEF'))

        contribs =                                                              \
            (Job.contribs == GET /INPUT/ Job.stordir
                             /FUNC/ (lambda x: readfile(xcpath(x)))
                             /CONST/ Job.has_contribs)
