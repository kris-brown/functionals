# External Modules
from json   import loads

# Internal Modules
from dbgen2 import (Model, Gen, PyBlock, Ref, Query, Expr, AND, GROUP_CONCAT,
                   SUM, COUNT, Literal, NOT, NULL, CONCAT, ABS, MIN, EQ, NE,
                    defaultEnv, Env,Import)

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

def analysis(mod : Model) -> None:
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
    BJ  = Bulk_job.r('job')
    BJE = Bulk_job.r('Expt')
    ########################################################################
    beq = Query(exprs = {'n':Struct['n_atoms']},
                basis   = ['Bulk_job'],
                constr  = EQ(Calc['xc'], 'mBEEF'))
    bulkexpt =                                                                  \
        Gen(name = 'bulkexpt',
            desc = 'All pairs of (bulk) species + calc., if mBEEF calc'\
                    'Links the bulk_job table to Expt table, too',
            actions = [Bulk_job(bulk_job=Ref('bulk_job'),
                                expt=Expt(insert=True,
                                          n_atoms = beq['n'],
                                          species=Ref('species'),
                                          calc = Ref('calc')))],

            query   = beq)
    ########################################################################
    aggs = ['contribs','energies','volumes']

    agg_dict = dict(contribs = Job['contribs'](BJE),
                    energies = Job['energy'](BJE),
                    volumes  = Cell['volume'](BJE))

    paq = Query(exprs = {x : CONCAT('[',GROUP_CONCAT(y),']')
                        for x,y in agg_dict.items()},
                links   = [BJE],
                aggcols = [Expt.id()],
                basis   = ['expt'],
                constr  = Bulk_job['near_min'])

    pop_aggs =                                                                  \
        Gen(name    = 'pop_aggs',
            desc    = 'Store information about closest 5 singlepoints near optimum',
            actions = [Expt(expt=Ref('expt'),**{x:paq[x] for x in aggs})] ,
            query     = paq)

    ########################################################################
    avq = Query(exprs   = {'av':CONCAT('[',GROUP_CONCAT(Bulk_job['dv']),']')},
                links   = [BJE],
                aggcols = [Expt.id()] )
    all_vols =                                                                  \
        Gen(name    = 'all_vols',
            desc    = 'Store all volumes',
            actions = [Expt(expt=Ref('expt'),all_vols=avq['av'])],
            query   = avq)
    ########################################################################
    nmq  = Query(exprs = {'dv':Bulk_job['dv'],'av':Expt['all_vols']})
    nmpb = PyBlock(get_near_min,args=[nmq['dv'],nmq['av']],outnames=['nm'])
    near_min =                                                                  \
        Gen(name = 'near_min',
            desc = '?',
            actions = [Bulk_job(bulk_job=Ref('bulk_job'),near_min=nmpb['nm'])],
            query   = nmq,
            funcs   = [nmpb])

    ########################################################################
    Num = Species_comp['num']
    Eng  = Reference['energy'](multi=((Reference.r('Calc'),   ),
                                      (Reference.r('Element'),)))

    eq = Query(exprs = {'e':Expt['energy_pa'] - SUM(Eng* Num)/ SUM(Num)},
                links   = [Species_comp.r('Species'),
                             Reference.r('Calc'),
                             Reference.r('Element')],
                aggcols = [Expt.id()],
                basis   = ['expt'],
                constr  =  NOT(NULL(Job['contribs'](Reference.r('Job')))), #Edge case - what if job has an energy but x_contribs failed?
                aconstr = EQ(COUNT(Literal(1)), Species['n_elems']))


    eform =                                                                     \
        Gen(name    = 'eform',
            desc    = 'Difference between a relaxed energy (per atom) and '     \
                      'the sum of reference energies (divided by # of atoms)',
            actions = [Expt(expt=Ref('expt'),eform=eq['e'])],
            query   = eq)

    ########################################################################

    ks    = ['kx','ky','kz']
    kds   = ['kptden_x','kptden_y','kptden_z']

    k_env = defaultEnv + Env(Import('ase','Atoms'))
    p_env = defaultEnv + Env(Import('dbgen.core.parsing','parse_line'))

    kq = Query(exprs = {'log':Job['log'],**{x:Cell[x] for x in Cell.ids()}})

    kpb0 =  PyBlock(get_kpts_gpaw,env = p_env,
                    args=[kq['log']], outnames=ks)
    kpb1 = PyBlock(get_kptden, env = k_env,
                   args = [kpb0[x] for x in ks]+[kq[x] for x in Cell.ids()],
                   outnames = kds)
    kpts =                                                                      \
        Gen(name    = 'kpts',
            desc    = 'Kpoints and kpoint density for each job',
            actions = [Job(job  = Ref('job'),
                           **{x:kpb0[x] for x in ks},
                           **{x:kpb1[x] for x in kds})],
            query   = kq,
            funcs   = [kpb1,kpb0])

    ########################################################################

    eos_cols = ['energy_pa', 'bulkmod', 'volume_pa', 'img']
    eos_env  = Env(Import('string','ascii_lowercase'),
                    Import('random','choices'),
                    Import('os','environ','remove'),
                    Import('os.path','join'),
                    Import('base64','b64encode'),
                    Import('ase.eos','EquationOfState'),
                    Import('numpy','polyfit','poly1d','mean','std','array')
                    ) + defaultEnv

    eosq = Query(exprs={'v':GROUP_CONCAT(Cell['volume']),
                        'e':GROUP_CONCAT(Job['energy']),
                        'n':Expt['n_atoms'],
                        'ndata':COUNT(Literal(1))},
                basis   = ['Bulk_job'],
                aggcols = [Expt.id()])

    eospb = PyBlock(eos_func, env = eos_env,
                   args  = [eosq[x] for x in 'ven'], outnames = eos_cols)
    eos =                                                                       \
        Gen(name    = 'eos',
            desc    = 'Uses ASE equation of state over aggregated volumes '
                      'and energies for an experiment',
            actions = [Expt(expt   = Ref('expt'),
                            n_data = eosq['ndata'],
                            **{x:eospb[x] for x in eos_cols})],
            query   = eosq,
            funcs   = [eospb])

    ########################################################################
    lq = Query(exprs={'s':Species['symmetry'],'v':Expt['volume_pa']})
    lpb = PyBlock(conventional_lattice,args=[lq['s'],lq['v']],outnames=['l'])
    lattice =                                                                   \
        Gen(name    = 'lattice',
            desc    = 'Computes lattice constant of conventional unit cell',
            actions = [Expt(expt=Ref('expt'),lattice = lpb['l'])],
            query   = lq,
            funcs   = [lpb])

    ########################################################################

    dcols = ['coefs','name','composition',
             'energy_vector','volume_vector','contrib_vector']
    cq = Query(exprs={'coefs':Calc['coefs'],
                      'name':Species['nickname'],
                      'comp':Species['composition'],
                      'ev' : Expt['energies'],
                      'vv' : Expt['volumes'],
                      'cv' : Expt['contribs']},
            basis   = ['expt'],
            constr  = EQ(Calc['xc'],'mBEEF')
                        |AND| NOT(NULL(Expt['eform']))
                        |AND| (Expt['n_data'] > 5))
    cohesive =                                                                  \
        Gen(name    = 'cohesive',
            desc    = 'Populates DFT data',
            actions = [D(insert = True,
                         coefs  = cq['coefs'],
                         name   = cq['name'],
                         composition = cq['comp'],
                         energy_vector = cq['ev'],
                         volume_vector = cq['vv'],
                         contrib_vector = cq['cv'],
                         expt  = Ref('expt'))],
            query   = cq)



    ########################################################################

    sde_dict = {'cohesive energy' : 'expt_cohesive_energy',
                'bulk modulus'    : 'expt_bm'}

    twogens = []
    for x,y in sde_dict.items():
        q = Query(exprs={'v':SDE['value']},
                  links = [SDE.r('species')],
                  basis = ['dft_data'],
                  constr= EQ(SDE['property'], x))

        twogens.append(Gen(name=x.replace(' ','_'),
                           desc= 'Copying %s info from Kelddata into dataset'%x,
                           actions=[D(dft_data=Ref('dft_data'),**{y:q['v']})],
                           query  = q))

    cohesive_target,bm_target = twogens # unpack


    ########################################################################
    vq = Query({'v':SDE['value'], 's':Species['symmetry']},
                links   = [SDE.r('Species')],
                basis   = ['dft_data'],
                constr  = EQ(SDE['property'], 'lattice parameter'))

    vpb = PyBlock(make_volume,
                   args = [vq['v'],vq['s']],
                   outnames = ['v'])
    vol_target =                                                                \
        Gen(name    = 'vol_target',
            desc    = '?',
            actions = [D(dft_data=Ref('dft_data'),expt_volume=vpb['v'])],
            query   = vq,
            funcs   = [vpb])

    ########################################################################

    def make_dict(attr : str) -> Expr:
        '''
        Makes a serialized python dictionary mapping element_ids to Job
        properties of their reference calculation
        '''
        elemattr = Element['atomic_number'](Species_comp.r('Element'))
        jobattr  = Job[attr](multi=((Reference.r('Element'),),
                                    (Reference.r('Calc'),)))
        group    = GROUP_CONCAT(CONCAT(elemattr,':',jobattr))
        return CONCAT('{',group,'}')

    ceq = Query({'a_c':make_dict('contribs'),'a_e':make_dict('energy')},
                links = [Species_comp.r('Species'),
                             Reference.r('Element'),
                             Reference.r('Calc')],
                basis   = ['dft_data'],
                aggcols = [D.id()])

    cohesive_elems =                                                            \
        Gen(name    = 'cohesive_elems',
            desc    = '?',
            actions = [D(dft_data=Ref('dft_data'),
                        atomic_contribs=ceq['a_c'],
                        atomic_energies=ceq['a_e'])],
            query   = ceq)
    ########################################################################

    vol_pa = Cell['volume'] / Struct['n_atoms']
    gapq = Query({'gap':ABS(vol_pa - Expt['volume_pa']),
                   'dv':vol_pa - Expt['volume_pa']},
                  basis = ['bulk_job'])
    gap =                                                                       \
        Gen(name    = 'gap',
            desc    = 'Computes deviation from optimum for each Bulk_job '
                      'singlepoint',
            actions = [Bulk_job(bulk_job=Ref('bulk_job'),
                                dv=gapq['dv'],gap=gapq['gap'])],
            query   = gapq)

    ########################################################################
    minq = Query({'m':MIN(Bulk_job['gap'])},
                 links   = [BJE],
                 aggcols     = [Expt.id()])
    mingap =                                                                    \
        Gen(name    = 'mingap',
            desc    = 'Stores a value under "Expt" unique to the best corresponding Bulk_job' ,
            actions = [Expt(expt=Ref('expt'),min_gap=minq['m'])] ,
            query   = minq)


    ########################################################################

    ratio = Struct['n_atoms'](BJ) / Species['n_atoms']
    co = Query({'be': Job['energy'](BJ),'bc': Job['contribs'](BJ),'br':ratio},
                links   = [BJE],
                basis   = ['dft_data'],
                constr  = EQ(Bulk_job['gap'], Expt['min_gap']))

    cohesive_optjob =                                                           \
        Gen(name    = 'cohesive_optjob',
            desc    = 'Populates FK using the data populated by generators '
                       '"gap" and "mingap"',
            actions = [D(dft_data      = Ref('dft_data'),
                         best_job      = Ref('job'),
                         bulk_energy   = co['be'],
                         bulk_contribs = co['bc'],
                         bulk_ratio    = co['br'])],
            query   = co)


    gens = [bulkexpt, pop_aggs, all_vols, near_min, eform, kpts, eos, lattice,
            cohesive, cohesive_target, bm_target, vol_target, cohesive_elems,
            gap, mingap, cohesive_optjob,]

    mod.add(gens)
