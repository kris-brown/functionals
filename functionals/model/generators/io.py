# External Modules
from typing     import Callable as C, Union as U, List as L, Tuple as T, Optional as O, Dict as D
from os         import environ
from os.path    import exists, join
from csv        import DictReader, reader
from json       import load, dumps
# Internal Modules
from dbgen import (Model, Gen, PyBlock, Env, Query, Const, Import, defaultEnv,
                    EQ, Literal as Lit, Arg, JPath, CONCAT, LIKE)

from functionals.scripts.io.parse_setup             import parse_setup
from functionals.scripts.io.parse_mendeleev         import parse_mendeleev
from functionals.scripts.io.parse_keld              import parse_keld
from functionals.scripts.io.parse_fit               import parse_fit
from functionals.scripts.io.get_stray_vasp          import get_stray_vasp
from functionals.scripts.load.parse_incar           import parse_incar
from functionals.scripts.load.parse_potcar          import parse_potcar
from functionals.scripts.load.parse_pw_vasp         import parse_pw_vasp
from functionals.scripts.load.parse_xc_vasp         import parse_xc_vasp
from functionals.scripts.load.parse_contribs_vasp   import parse_contribs_vasp

##############################################################################
# CONSTANTS #
 ##########
elempath = environ['ELEMPATH']
keldpath = environ['KELDPATH']
logpth   = '/Users/ksb/scp_tmp/vauto' #
psppth   = '/Users/ksb/vossj/vasp/potpaw_PBE'
fitpth   = '/Users/ksb/scp_tmp/fitresult'
pthEnv  = Env(Import('os.path','join','exists'))

a1msb = ['a11','a12','a13','a14','a15','msb']

def const() -> T[L[O[str]],L[O[str]],L[O[float]],L[O[float]],L[O[float]],L[O[str]],L[O[str]]]:

    with open('/Users/ksb/functionals/data/constraints.csv','r') as f:
        r = reader(f)
        next(r)
        name,desc,s,a,val,typ,vec = tuple(map(list,zip(*r)))
    s_,a_,val_  = [[float(str(x)) if x else None for x in xs] for xs in [s,a,val]]
    n_,d_,t_,v_ = [[str(x) if x else None for x in xs] for xs in [name,desc,typ,vec]]
    return n_,d_,s_,a_,val_,t_,v_

def readfile(pth:str)->str:
    with open(pth,'r') as f: return f.read()

parse_env = defaultEnv + Env(Import('dbgen.utils.parsing','parse_line'))

elemcols = ['symbol', 'name', 'atomic_weight','atomic_radius',
            'phase','evaporation_heat', 'pointgroup','spacegroup',
            'melting_point', 'metallic_radius', 'vdw_radius',
            'density', 'en_allen' , 'is_radioactive',
            'lattice_struct' , 'fusion_heat',
            'econf', 'period', 'covalent_radius_bragg',
            'geochemical_class', 'abundance_crust', 'heat_of_formation',
            'electron_affinity', 'atomic_volume',  'boiling_point',
            'proton_affinity', 'covalent_radius_slater',
            'lattice_constant', 'dipole_polarizability',
            'en_ghosh', 'thermal_conductivity', 'group_id', 'en_pauling',
            'gas_basicity', 'abundance_sea']
#################################################################################
def io(mod : Model) -> None:

    # Extract tables
    tabs = ['job', 'atom', 'element', 'struct', 'calc', 'cell', 'pure_struct',
            'species', 'expt', 'bulk_job', 'reference',
            'species_comp', 'species_dataset',
            'species_dataset_element', 'incar', 'potcar', 'functional', 'beef',
            'const','fitparams','fit']

    Job, Atom, Element, Struct, Calc, Cell, Pure_struct, Species, Expt, Bulk_job,\
    Reference, Species_comp, SD, SDE,Incar, Potcar,Functional,Beef, Cnst, \
    Fitparams, Fit \
        = map(mod.get, tabs)

    job__calc, job__incar, job__struct \
        = map(mod.get_rel,[Job.r('calc'),Job.r('incar'),Job.r('struct')])
    ########################################################################
    gl_env = defaultEnv + Env(Import('subprocess','getstatusoutput'))

    glpb = PyBlock(get_stray_vasp,
                   env  = gl_env,
                   args = [Const(logpth)])
    getlogs  =                                                                  \
        Gen(name    = 'getlogs',
            desc    = 'Scrape folder recursively for any VASP jobs',
            actions = [Job(insert  = True,
                           stordir = glpb['out'])],
            funcs    = [glpb])
    ########################################################################
    dvq  = Query(exprs={'j' : Job.id(), 'log' : CONCAT(Job['stordir'](),Lit('/OUTCAR'))})
    dvpb = PyBlock(readfile, args = [dvq['log']])
    datavasp = Gen(name    = 'datavasp',
                   desc    = 'Read log file from disk',
                   query   = dvq,
                   actions = [Job(job = dvq['j'], log = dvpb['out'])],
                   funcs   = [dvpb])
    ########################################################################
    su_env = defaultEnv + Env(Import('gzip',open='gopen'),
                              Import('hashlib','md5'),
                              Import('os','listdir'),
                              Import('os.path','isdir','join'),
                              Import('xml.etree.ElementTree','fromstring'),
                              Import('ase.data','chemical_symbols'))


    ########################################################################
    eiq       = Query(exprs={'e':Element.id(),
                             'z':Element['atomic_number']()})

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

    kelddata =                                                                  \
        Gen(name    = 'kelddata',
            desc    = 'Parses a .py Keld data file',
            actions = [isde],
            funcs   = [kpb])
    ########################################################################


    ipth = JPath(Incar, [job__incar])


    cq    = Query(exprs={'j'    : Job.id(),
                         'sd'   : Job['stordir'](),
                         'p'    : Incar['encut'](ipth),
                         'e'    : Incar['ediff'](ipth),
                         **{x:Incar[x](ipth) for x in a1msb}},
                  basis = [Job])

    xc    = PyBlock(parse_xc_vasp,
                    env  = parse_env + Env(Import('os.path','exists')),
                    args = [cq['sd']])

    def floatfun(*xs:float)->tuple:
        return tuple(map(float,xs))

    ffloat = PyBlock(floatfun,
                     args     = [cq[z] for z in ['p','e']+a1msb],
                     outnames = ['p','e'] + a1msb)

    calc =                                                                      \
        Gen(name    = 'calc',
            desc    = 'Populate calc table + F.K. from Relax_job',
            query   = cq,
            actions = [Job(job  = cq['j'],
                           calc = Calc(insert     = True,
                                       pw         = ffloat['p'],
                                       econv      = ffloat['e'],
                                       **{x:ffloat[x] for x in a1msb},
                                       functional = Functional(insert = True,
                                                               data   = xc['out'])))],
            funcs   = [xc,ffloat])

    bq = Query(exprs  = {'f':Functional.id()},
               constr = Functional['data']() |LIKE| Lit('[%'))
    fxbeef = Gen(name = 'beef',
                 desc = 'Populate BEEF table from Functionals table',
                 query = bq ,
                 actions = [Beef(insert=True, functional = bq['f'])])

    ########################################################################

    icols = ['encut','sigma','metagga','prec','ediff','algo','ismear','npar',
            'nelm','ispin','ibrion','lcharg','lbeefens','addgrid','lasph','lwave',
            'a11','a12','a13','a14','a15','msb','magmom']

    iq    = Query(exprs={'j'   : Job.id(),
                         'log' : CONCAT(Job['stordir'](), Lit('/INCAR'))})
    ipb   = PyBlock(parse_incar, args = [iq['log']], outnames = icols)
    incar = Gen(name    = 'incar',
                desc    = 'Read incar file from disk',
                query   = iq,
                actions = [Job(job   = iq['j'],
                               incar = Incar(insert=True,
                                             **{k:ipb[k] for k in icols}))],
                funcs   = [ipb])

    ########################################################################
    ctrq  = Query(exprs  = {'j':Job.id(),'l':Job['log']()})

    ctrpb = PyBlock(parse_contribs_vasp,args = [ctrq['l']])

    contribs =                                                                  \
        Gen(name    = 'contribs',
            desc    = '?',
            actions = [Job(job      = ctrq['j'],
                           contribs = ctrpb['out'])],
            query   = ctrq,
            funcs   = [ctrpb])

    ########################################################################
    js  = JPath(Struct,[job__struct])
    pcols = ['ind','titel','lultra','iunscr','rpacor','pomass','zval','rcore',
             'rwigs','enmax','enmin','lcor','lpaw','eaug','rmax','raug','rdep',
             'rdept']
    ppq = Query(exprs = dict(sd = CONCAT(Job['stordir'](),Lit('/POTCAR')),
                             s  = Struct.id(js)),
                basis = [Job])
    pppb = PyBlock(parse_potcar,
                   env  = defaultEnv + Env(Import('re','search')),
                   args = [ppq['sd']],
                   outnames = pcols)
    potcar =                                                                    \
        Gen(name    = 'potcar',
            desc    ='Finds potcar of every atom',
            query   = ppq,
            funcs   = [pppb],
            actions = [Atom(insert = True,
                            struct = ppq['s'],
                            ind    = pppb['ind'],
                            potcar = Potcar(insert = True,
                                            **{x:pppb[x] for x in pcols[1:]}))])

    ########################################################################

    def int_occupied(pth:str)->bool:
        '''Parse Vasp EIGENVAL - check if all bands have integer occupation'''
        with open(pth,'r') as f:
            lines = reversed(f.readlines())
        for l in lines:
            n,_,_,a,b = map(float,l.split())
            if int(a)!=a or int(b)!= b: return False
            if n == 1.0: return True
        raise ValueError('Should not reach this part of code...')

    iq    = Query(exprs={'j'   : Job.id(),
                         'log' : CONCAT(Job['stordir'](), Lit('/EIGENVAL'))})
    iopb  = PyBlock(int_occupied,args=[iq['log']])

    iogen =                                                                     \
        Gen(name    = 'int_occupied',
            desc    = 'populates',
            query   = iq,
            funcs   = [iopb],
            actions = [Job(job=iq['j'],int_occupation=iopb['out'])])


    ########################################################################

    constcols = ['const_name','description','s','alpha','val','kind','vec']

    csv_env   = defaultEnv + Env(Import('csv','reader'))
    pcpb      = PyBlock(const, env = csv_env, outnames = constcols)

    pop_constraint =                                                            \
        Gen(name    = 'pop_constraint',
            desc    = 'populate constraints',
            funcs   = [pcpb],
            tags    = ['fit'],
            actions = [Cnst(insert = True, **{x:pcpb[x] for x in constcols})])

################################################################################
    fitcols = ['constden','consts','reg','dataconst', 'bm_weight','lat_weight']

    pf_env = Env(Import('os','environ'), Import('csv','DictReader')) + defaultEnv

    def parse_fitparams() -> T[L[str],L[str],L[str],L[str],L[str],L[str]]:
        fitcols = ['constden','consts','reg','dataconst','bm_weight','lat_weight']
        vals = {c:[] for c in fitcols} # type: D[str,L[str]]
        with open(environ['FITPATH'],'r') as fi:
            reader = DictReader(fi)
            for row in reader:
                for k,v in vals.items(): v.append(row[k])
        a,b,c,d,e,f = tuple([list(x) for x in vals.values()])
        return a,b,c,d,e,f

    pfpb = PyBlock(parse_fitparams,
                   env      = pf_env,
                   outnames = fitcols)

    pop_fitparams =                                                             \
        Gen(name    = 'pop_fitparams',
            desc    = 'Looks in FITPATH for fitting specifications',
            funcs   = [pfpb],
            actions = [Fitparams(insert = True,
                           **{x:pfpb[x] for x in fitcols})])

    ########################################################################
    fcols = ['name','a1','pth','timestamp','decay','nsteps','steps','result','runtime']
    fpb = PyBlock(parse_fit,
                  args     = [Const(fitpth)],
                  outnames = fcols + ['pw','econv','data'] + fitcols + a1msb)

    pop_fit =                                                                   \
        Gen(name    = 'pop_fit',
            desc    = 'looks for completed fitting jobs - grabs identifying info',
            funcs   = [fpb],
            tags    = ['io'],
            actions = [Fit(insert = True,
                           msb    = fpb['msb'],
                           **{x:fpb[x] for x in fcols},
                           calc   = Calc(insert = True,
                                        **{x:fpb[x] for x in ['pw','econv']+a1msb},
                                        functional = Functional(insert = True,
                                                                data   = fpb['data'])),
                           fitparams = Fitparams(insert=True,
                                                 **{x:fpb[x] for x in fitcols}))])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [getlogs,datavasp,kelddata,elemzinfo,calc,contribs,incar,potcar,
            fxbeef,iogen,pop_constraint,pop_fitparams,pop_fit,
            ]#pop_fitstep]

    mod.add(gens)


###deprecated
    ########################################################################
    # def parse_fitstep(pth:str)->T[int,L[int],L[float],L[float]]:
    #     with open(pth,'r') as f: d = load(f)
    #     ntot = len(d)
    #     import pdb;pdb.set_trace()
    #     steps,costs,viols = map(list,zip(*d))
    #     import pdb;pdb.set_trace()
    #     return ntot,steps,costs,viols # type: ignore
    #
    # pfsq    = Query(exprs = {'f':Fit.id(),'o':CONCAT(Fit['pth'](),Lit('/output.json'))})
    # pfsns   = ['n', 'niter', 'cost', 'c_viol']
    # pfspb   = PyBlock(parse_fitstep,
    #                   args     = [pfsq['o']],
    #                   outnames = pfsns)
    #
    # pop_fitstep =                                                      \
    #     Gen(name    = 'pop_fitstep',
    #         desc    = 'Analyzes fitting result log, populates Fitstep',
    #         query   = pfsq,
    #         funcs   = [pfspb],
    #         #tags    = ['fit', 'parallel'],
    #         actions = [Fit(fit      = pfsq['f'],
    #                        nsteps   = pfspb['n']),
    #                    Fit_step(insert = True,
    #                             fit = pfsq['f'],
    #                             **{x:pfspb[x] for x in pfsns[1:]})])
    ########################################################################
    # nlconstcols = ['nlconst_name','description','f','df','hess','lb','ub']
    #
    # ncpb = PyBlock(nlconst,
    #                env      = csv_env,
    #                outnames = nlconstcols)
    #
    # pop_nlconstraint =                                                          \
    #     Gen(name    = 'pop_nlconstraint',
    #         desc    = 'populate nonlinear constraints',
    #         funcs   = [ncpb],
    #         tags    = ['fit'],
    #         actions = [Nl_const(insert = True,
    #                             **{x:ncpb[x] for x in nlconstcols})])
