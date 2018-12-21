# External Modules
from typing     import Callable as C, Union as U, List as L
from os         import environ
from os.path    import exists, join

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Env, Query, Const, Import, defaultEnv, EQ, Arg)

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
logpth   = '/Users/ksb/scp_tmp/auto/auto_small' #
psppth   = '/Users/ksb/scp_tmp/norm_conserving_setups'

pthEnv  = Env(Import('os.path','join','exists'))

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

parse_env = defaultEnv + Env(Import('dbgen.utils.parsing','parse_line'))

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
    Reference, Setup, Setup_family, Job_setup, Species_comp, SD, SDE \
        = map(mod.get, tabs)

    ########################################################################
    gl_env = defaultEnv + Env(Import('subprocess','getstatusoutput'))
    glpb = PyBlock(get_stray_gpaw,
                   env  = gl_env,
                   args = [Const(logpth)])
    getlogs  =                                                              \
        Gen(name    = 'getlogs',
            desc    = 'Scrape folder recursively for any GPAW jobs',
            actions = [Job(insert  = True,
                           logfile = glpb['out'])],
            funcs    = [glpb])
    ########################################################################
    sdq  = Query(exprs={'j'   : Job.id,
                        'log' : Job['logfile']})

    sdpb = PyBlock(lambda x: x[:x.rfind('/')],
                   args = [sdq['log']])

    getsd    = Gen(name    = 'getsd',
                   desc    = 'Get folder from logfile path',
                   actions = [Job(job     = sdq['j'],
                                  stordir = sdpb['out'])],
                   query   = sdq,
                   funcs   = [sdpb])
    ########################################################################
    dgq  = Query(exprs={'j' : Job.id, 'log' : Job['logfile']})
    dgpb = PyBlock(readfile, args = [dgq['log']])
    datagpaw = Gen(name    = 'datagpaw',
                   desc    = 'Read log file from disk',
                   query   = dgq,
                   actions = [Job(job = dgq['j'], log = dgpb['out'])],
                   funcs   = [dgpb])
    ########################################################################
    su_env = defaultEnv + Env(Import('gzip',open='gopen'),
                              Import('hashlib','md5'),
                              Import('os','listdir'),
                              Import('os.path','isdir','join'),
                              Import('xml.etree.ElementTree','fromstring'),
                              Import('ase.data','chemical_symbols'))

    spb = PyBlock(parse_setup,
                  env      = su_env,
                  args     = [Const(psppth)],
                  outnames = ['check','xc','kind','name','z','val'])

    isetup = Setup_family(insert = True,
                          name   = spb['name'],
                          kind   = spb['kind'],
                          xc     = spb['xc'])

    ielem  = Element(insert        = True,
                     atomic_number = spb['z'])
    setups =                                                                    \
        Gen(name    = 'setups',
            desc    = 'Populate Setup table from a path containing setups',
            actions = [Setup(insert       = True,
                             setup_family = isetup ,
                             element      = ielem,
                             checksum     = spb['check'],
                             val          = spb['val'])],
            funcs    = [spb])
    ########################################################################
    eiq       = Query(exprs={'e':Element.id,
                             'z':Element['atomic_number']})

    eipb      = PyBlock(parse_mendeleev,
                        args     = [eiq['z'], Const(elempath)],
                        outnames = elemcols)
    elemzinfo =                                                                 \
        Gen(name    = 'elemzinfo',
            desc    = '''Read a JSON file containing element reference data
                        Requires providing path to this file as a constant''',
            actions = [Element(element = eiq['e'],
                               **{x:eipb[x] for x in elemcols})],
            query   = eiq,
            funcs   = [eipb])


    ########################################################################
    # things that get populated by kelddata
    keldenv  = Env(Import('ase.data','chemical_symbols'),
                   Import('re','findall')) + defaultEnv

    kpb = PyBlock(parse_keld,
                  env      = keldenv,
                  args     = [Const(keldpath)],
                  outnames = ['ds_name','comp',
                              'sym','UNUSEDVAL','UNUSEDVAL2','n_elems',
                              'prop','val','dt'])

    ispecies = Species(insert      = True,
                       composition = kpb['comp'],
                       symmetry    = kpb['sym'])

    isd  = SD(insert       = True,
              dataset_name = kpb['ds_name'])

    isde = SDE(insert          = True ,
              property        = kpb['prop'],
              value           = kpb['val'],
              datatype        = kpb['dt'],
              species_dataset = isd,
              species         = ispecies)

    kelddata =                                                              \
        Gen(name    = 'kelddata',
            desc    = 'Parses a .py Keld data file',
            actions = [isde],
            funcs   = [kpb])
    ########################################################################

    cq    = Query(exprs={'j':Job.id,'log':Job['log'],'sd':Job['stordir']})

    pw    = PyBlock(parse_pw_gpaw, env = parse_env,
                    args = [cq['log']])
    xc    = PyBlock(parse_xc_gpaw, env = parse_env,
                    args = [cq['log']])
    econv = PyBlock(get_econv_gpaw,env=parse_env,
                    args=[cq['log']])
    beef  = PyBlock(lambda x: int('beef' in x.lower()),
                    args=[xc['out']])
    coef  = PyBlock(getcoef,env=pthEnv,
                    args = [cq['sd'],beef['out']])

    calc =                                                                      \
        Gen(name    = 'calc',
            desc    = 'Populate calc table + F.K. from Relax_job',
            query   = cq,
            actions = [Job(job  = cq['j'],
                           calc = Calc(insert = True,
                                     pw       = pw['out'],
                                     xc       = xc['out'],
                                     coefs    = coef['out'],
                                     econv    = econv['out'],
                                     beef     = beef['out']))],
            funcs   = [pw,xc,econv,coef,beef])
    ########################################################################
    def xx(x:str)->bool:
        return exists(join(x,'xccontribs.txt'))

    hcq  = Query(exprs  = {'j'  : Job.id,
                           'sd' : Job['stordir']},
                 basis  = ['job'],
                 constr = EQ(Calc['xc'],'mBEEF'))

    hcpb = PyBlock(xx,env=pthEnv,args=[hcq['sd']])
    has_contribs =                                                              \
        Gen(name    = 'has_contribs',
            desc    = '?',
            actions = [Job(job          = hcq['j'],
                           has_contribs = hcpb['out'])],
            query   = hcq,
            funcs   = [hcpb])
    ########################################################################
    def get_c(x:str)->str:
        with open(join(x,'xccontribs.txt'),'r') as f:
            return f.read()

    ctrq  = Query(exprs  = {'j':Job.id,
                            'sd':Job['stordir']},
                  constr = Job['has_contribs'])

    ctrpb = PyBlock(get_c,
                    env      = pthEnv,
                    args     = [ctrq['sd']])

    contribs =                                                              \
        Gen(name    = 'contribs',
            desc    = '?',
            actions = [Job(job      = ctrq['j'],
                           contribs = ctrpb['out'])],
            query   = ctrq,
            funcs   = [ctrpb])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [getlogs,getsd,setups,datagpaw,kelddata,elemzinfo,calc,has_contribs,contribs]

    mod.add(gens)
