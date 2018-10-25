# External Modules
from typing import Any,Type,Tuple,TYPE_CHECKING,Callable as C
from os import environ
from os.path import exists,join
from json import loads

if TYPE_CHECKING:
    from dbgen.support.model     import Model

from dbgen.support.get      import CONST,DESC,INPUT,FUNC,CONSTS,GET,TAGS,BASIS,LINKS,IO
from dbgen.support.expr     import IN
from dbgen.support.funclike import SimpleFunc,PyBlock,noIndex,Unpack,SimplePipe
from dbgen.core.lists       import normalize_list,nub

from functionals.scripts.get_stray_gpaw  import  get_stray_gpaw
from functionals.scripts.anytraj_local   import  anytraj_local
from functionals.scripts.metadata        import  metadata
from functionals.scripts.get_atoms       import  get_atoms
from functionals.scripts.parse_mendeleev import  parse_mendeleev
from functionals.scripts.parse_pw_gpaw   import  parse_pw_gpaw
from functionals.scripts.parse_psp_gpaw  import  parse_psp_gpaw
from functionals.scripts.parse_xc_gpaw   import  parse_xc_gpaw
from functionals.scripts.get_kpts_gpaw   import  get_kpts_gpaw
from functionals.scripts.get_econv_gpaw  import  get_econv_gpaw
from functionals.scripts.get_cell        import  get_cell

from functionals.scripts.get_system_type   import get_system_type
from functionals.scripts.get_pointgroup    import get_pointgroup
from functionals.scripts.get_spacegroup    import get_spacegroup
from functionals.scripts.json_to_traj      import json_to_traj


###########################################################
elempath = environ['ELEMPATH']
logpth = '/Users/ksb/scp_tmp/auto'
def readfile(pth:str)->str:
    with open(pth,'r') as f:
        return f.read()

spin     = lambda x: int('Spin-polarized calculation' in x)
beefpath = '/Users/ksb/functionals/data/beef.json'
coefpath = lambda x: join(x,'coefs.json') # type: C[[str],str]
getcoef  = lambda stor,bf: (readfile(coefpath(stor))
                           if exists(coefpath(stor))
                           else readfile(beefpath)) \
                            if bf else ''

calc_plan =  SimpleFunc(parse_pw_gpaw,['log'],['pw'])           \
           + SimpleFunc(parse_psp_gpaw,['log'],['psp'])         \
           + SimpleFunc(parse_xc_gpaw,['log'],['xc'])           \
           + SimpleFunc(get_econv_gpaw,['log'],['econv'])       \
           + SimpleFunc(get_kpts_gpaw,['log'],['kpts'])         \
           + SimpleFunc(spin,['log'],['spinpol'])               \
           + SimpleFunc(getcoef,['stordir','beef'],['coefs'])   \
           + Unpack('kpts',['kx','ky','kz'])                    \
           + SimpleFunc(lambda x: int('beef' in x.lower()),['xc'],['beef'])

def cell_info(a0:float,a1:float,a2:float,
              b0:float,b1:float,b2:float,
              c0:float,c1:float,c2:float
              ) -> Tuple[float,float,float,float,float]:
    """ Basic geometry (inputs are 3 <x,y,z> vectors) """
    surface_area = float((a0*b2-a2*b1) * (a1*b2 - a2*b1)
                     +(a0*b2-a2*b0) * (a0*b2 - a2*b0)
                     +(a0*b1-a1*b0) * (a0*b1 - a1*b0))**0.5
    volume = a0*(b1*c2-b2*c1)\
            -a1*(b0*c2-b2*c0)\
            +a2*(b0*c1-b1*c0)
    a = float(a0*a0+a1*a1+a2*a2)**0.5
    b = float(b0*b0+b1*b1+b2*b2)**0.5
    c = float(c0*c0+c1*c1+c2*c2)**0.5
    return surface_area,volume,a,b,c

def countatm(raw:str)->Tuple[int,int,str,str,str,str,str]:
    """Stoichiometry analysis"""
    nonmets  = [1,2,6,7,8,9,10,16,17,18]
    nums     = sorted([a['number'] for a in loads(raw)['atomdata']])
    nnorm    = normalize_list(nums)
    metnums  = [n for n in nums if n not in nonmets]
    uniqnums = nub(nums)
    comp     = {n:nums.count(n)  for n in uniqnums}
    compnorm = {n:nnorm.count(n) for n in uniqnums}
    metcomp  = {n:metnums.count(n) for n in uniqnums if n not in nonmets}
    consts   = str([a['constrained'] for a in loads(raw)['atomdata']])
    return len(nums),len(uniqnums),str(comp),str(compnorm),str(metcomp),\
            str(nums),consts

############################################################################
############################################################################
############################################################################
############################################################################
def add_generators(mod:Type['Model'])->None:
    tabs = ['job','atom','element','struct','calc','cell','pure_struct']
    Job,Atom,Element,Struct,Calc,Cell,Pure_struct = map(mod.get,tabs) # type: ignore
    with mod:
        ########################################################################
        getlogs  = Job == GET /FUNC/ (lambda : get_stray_gpaw(logpth))
        getsd    = Job.stordir == GET /INPUT/ Job.logfile /FUNC/ (lambda x: x[:x.rfind('/')])
        ########################################################################
        datagpaw = Job.log == GET /INPUT/ Job.logfile /FUNC/ readfile

        ########################################################################
        anytraj =                                                               \
            ((Job.Struct,Struct) ==
                GET /INPUT/ Job.stordir /FUNC/ anytraj_local /DESC/  '''Get any trajectory''')
        ########################################################################

        jobmetadata =                                                           \
            ((Job.user,Job.timestamp) ==
                GET /INPUT/ Job.stordir /FUNC/ metadata /DESC/ 'Scrape metadata about job')

        ########################################################################
        atomdetails = ['number','x','y','z','constrained','magmom','tag']
        atom =                                                                  \
            ((Atom,Atom.Element,Atom.select(atomdetails)) ==
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
        elemcols = ['symbol','atomic_weight','lattice_struct','econf']

        elemzinfo =                                                             \
            (Element.select(elemcols) ==
                 GET /INPUT/ Element.atomic_number
                     /DESC/   '''Read a JSON file containing element reference data
                                 Requires providing path to this file as a constant'''
                     /FUNC/   PyBlock(parse_mendeleev,
                                      args = noIndex(['atomic_number','pth']))
                     /CONSTS/ {'pth':elempath})

        ########################################################################
        calc =                                                                  \
            ((Job.Calc,Calc) ==
                GET /INPUT/  Job.select('stordir','log')
                    /BASIS/  Job
                    /FUNC/   calc_plan
                    /DESC/   'Populate calc table + F.K. from Relax_job')

        ########################################################################
        cell =                                                                  \
            ((Cell,Struct.Cell)
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
                    /FUNC/ SimplePipe([json_to_traj,get_system_type],
                                      inputs  = ['raw'],
                                      outputs = ['system_type'])

                    /DESC/ 'Apply function to identify whether system is '
                             'bulk/mol/surf')


        ########################################################################
        ps   =                                                                  \
            ((Pure_struct,Struct.Pure_struct) ==
                GET /INPUT/ Struct.raw
                    /FUNC/ SimplePipe([json_to_traj,get_bulk,get_pure_struct],
                                            inputs  = ['raw'],
                                            outputs = ['prototype'])
                    /TAGS/ ['long','parallel']
                    /DESC/ 'Determine the "pure structure" code for each Bulk')


        ########################################################################
        sg =                                                                    \
            (Struct.sg ==
                GET /INPUT/ Struct.raw
                    /FUNC/  SimplePipe([json_to_traj,get_spacegroup],
                                        ['raw'],['sg'])
                    /CONST/ (Struct.system_type |IN| ["bulk","surface"])
                    /DESC/  'Computes spacegroup of a bulk')


        ########################################################################

        ########################################################################

        ########################################################################

        ########################################################################

        ########################################################################
