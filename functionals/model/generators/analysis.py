# External Modules
from json   import loads

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Query, Expr, AND, GROUP_CONCAT,
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
        min_n = (n-1) // 2
    else:
        mid   = n // 2
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
            'species_comp','species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Species, Expt, Bulk_job,\
    Reference, Species_dataset,Species_comp,Species_dataset_element,\
      = map(mod.get, tabs) # type: ignore

    # Abbreviations
    #--------------
    SDE = Species_dataset_element
    BJ  = Bulk_job.r('job')
    BJE = Bulk_job.r('Expt')
    ########################################################################
    beq = Query(exprs = {'b' : Bulk_job.id,
                         'c' : Calc.id,
                         's' : Species.id(Struct.r('species')),
                         'n' : Struct['n_atoms']},
                basis   = ['Bulk_job'],
                constr  = EQ(Calc['xc'], 'mBEEF'))

    iexpt = Expt(insert  = True,
                    n_atoms = beq['n'],
                    species = beq['s'], # otherwise tries to link through FK that we're populating!
                    calc    = beq['c'])
    bulkexpt =                                                                  \
        Gen(name    = 'bulkexpt',
            desc    = 'All pairs of (bulk) species + calc., if mBEEF calc'\
                      'Links the bulk_job table to Expt table, too',
            query   = beq,
            actions = [Bulk_job(bulk_job = beq['b'],
                                expt     = iexpt)])
    ############################################################################

    aggs = ['contribs','energies','volumes']

    agg_dict = dict(contribs = Job['contribs'](BJE),
                    energies = Job['energy'](BJE),
                    volumes  = Cell['volume'](BJE))

    paq = Query(exprs = {'e':Expt.id,
                         **{x : CONCAT('[', GROUP_CONCAT(y), ']')
                             for x,y in agg_dict.items()}},
                links   = [BJE],
                aggcols = [Expt.id],
                basis   = ['expt'],
                constr  = Bulk_job['near_min'])

    pop_aggs =                                                                  \
        Gen(name    = 'pop_aggs',
            desc    = 'Store information about closest 5 singlepoints near optimum',
            query   = paq,
            actions = [Expt(expt=paq['e'],**{x:paq[x] for x in aggs})])

    ########################################################################
    avq = Query(exprs   = {'e':Expt.id,
                           'av':CONCAT('[',GROUP_CONCAT(Bulk_job['dv']),']')},
                basis   = ['expt'],
                links   = [BJE],
                aggcols = [Expt.id] )
    all_vols =                                                                  \
        Gen(name    = 'all_vols',
            desc    = 'Store all volumes',
            query   = avq,
            actions = [Expt(expt = avq['e'], all_vols = avq['av'])])
    ########################################################################
    nmq  = Query(exprs = {'b':Bulk_job.id,
                          'dv':Bulk_job['dv'],
                          'av':Expt['all_vols']},
                 basis = ['bulk_job'])

    nmpb = PyBlock(get_near_min,
                   args=[nmq['dv'],nmq['av']])
    near_min =                                                                  \
        Gen(name    = 'near_min',
            desc    = '?',
            query   = nmq,
            funcs   = [nmpb],
            actions = [Bulk_job(bulk_job = nmq['b'],
                                near_min = nmpb['out'])])

    ########################################################################
    Num = Species_comp['num']
    Eng  = Reference['energy'](multi=((Reference.r('Calc'),   ),
                                      (Reference.r('Element'),)))

    eq = Query(exprs = {'ex':Expt.id,'e':Expt['energy_pa'] - SUM(Eng* Num)/ SUM(Num)},
                links   = [Species_comp.r('Species'),
                           Reference.r('Calc'),
                           Reference.r('Element')],
                aggcols = [Expt.id],
                basis   = ['expt'],
                constr  = NOT(NULL(Job['contribs'](Reference.r('Job')))), #Edge case - what if job has an energy but x_contribs failed?
                aconstr = EQ(COUNT(Literal(1)), Species['n_elems']))


    eform =                                                                     \
        Gen(name    = 'eform',
            desc    = 'Difference between a relaxed energy (per atom) and '     \
                      'the sum of reference energies (divided by # of atoms)',
            actions = [Expt(expt  = eq['ex'],
                            eform = eq['e'])],
            query   = eq)

    ########################################################################

    ks    = ['kx','ky','kz']
    kds   = ['kptden_x','kptden_y','kptden_z']

    k_env = defaultEnv + Env(Import('ase','Atoms'))
    p_env = defaultEnv + Env(Import('dbgen.utils.parsing','parse_line'))

    kq = Query(exprs = {'j'   : Job.id,
                        'log' : Job['log'],
                        **{x:Cell[x] for x in Cell.ids()}},
               basis = ['job'])

    kpb0 =  PyBlock(get_kpts_gpaw,
                    env      = p_env,
                    args     = [kq['log']],
                    outnames = ks)

    kpb1 = PyBlock(get_kptden,
                   env      = k_env,
                   args     = [kpb0[x] for x in ks]+[kq[x] for x in Cell.ids()],
                   outnames = kds)
    kpts =                                                                      \
        Gen(name    = 'kpts',
            desc    = 'Kpoints and kpoint density for each job',
            query   = kq,
            funcs   = [kpb1,kpb0],
            actions = [Job(job = kq['j'],
                           **{x:kpb0[x] for x in ks},
                           **{x:kpb1[x] for x in kds})])

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

    eosq = Query(exprs={'ex':Expt.id,
                        'v':GROUP_CONCAT(Cell['volume']),
                        'e':GROUP_CONCAT(Job['energy']),
                        'n':Expt['n_atoms'],
                        'ndata':COUNT(Literal(1))},
                basis   = ['bulk_job'],
                aggcols = [Expt.id])

    eospb = PyBlock(eos_func,
                    env      = eos_env,
                    args     = [eosq[x] for x in 'ven'],
                    outnames = eos_cols)
    eos =                                                                       \
        Gen(name    = 'eos',
            desc    = 'Uses ASE equation of state over aggregated volumes '
                      'and energies for an experiment',
            query   = eosq,
            funcs   = [eospb],
            actions = [Expt(expt   = eosq['ex'],
                            n_data = eosq['ndata'],
                            **{x:eospb[x] for x in eos_cols})])

    ########################################################################
    lq = Query(exprs={'e' : Expt.id,
                      's' : Species['symmetry'],
                      'v' : Expt['volume_pa']},
               basis = ['expt'])

    lpb = PyBlock(conventional_lattice,
                  args=[lq['s'],lq['v']])
    lattice =                                                                   \
        Gen(name    = 'lattice',
            desc    = 'Computes lattice constant of conventional unit cell',
            actions = [Expt(expt=lq['e'], lattice = lpb['out'])],
            query   = lq,
            funcs   = [lpb])

    ########################################################################

    dcols = ['coefs','name','composition',
             'energy_vector','volume_vector','contrib_vector']

    cq = Query(exprs={'e'     : Expt.id,
                      'coefs' : Calc['coefs'],
                      'name'  : Species['nickname'],
                      'comp'  : Species['composition'],
                      'ev'    : Expt['energies'],
                      'vv'    : Expt['volumes'],
                      'cv'    : Expt['contribs']},
               basis   = ['expt'],
               constr  = EQ(Calc['xc'],'mBEEF')
                            |AND| NOT(NULL(Expt['eform']))
                            |AND| (Expt['n_data'] > 5))
    cohesive =                                                                  \
        Gen(name    = 'cohesive',
            desc    = 'Populates more things in experiment',
            query   = cq,
            actions = [Expt(expt           = cq['e'],
                         coefs          = cq['coefs'],
                         name           = cq['name'],
                         composition    = cq['comp'],
                         energy_vector  = cq['ev'],
                         volume_vector  = cq['vv'],
                         contrib_vector = cq['cv'])])



    ########################################################################

    sde_dict = {'cohesive energy' : 'expt_cohesive_energy',
                'bulk modulus'    : 'expt_bm'}

    twogens = []
    for x,y in sde_dict.items():
        q = Query(exprs={'d':Expt.id,'v':SDE['value']},
                  links = [SDE.r('species')],
                  basis = ['expt'],
                  constr= EQ(SDE['property'], x))

        twogens.append(Gen(name    = x.replace(' ','_'),
                           desc    = 'Copying %s info from Kelddata into dataset'%x,
                           query   = q,
                           actions = [Expt(expt = q['d'], **{y:q['v']})]))

    cohesive_target,bm_target = twogens # unpack


    ########################################################################
    vq = Query(exprs    = {'d': Expt.id,
                           'v':SDE['value'],
                           's':Species['symmetry']},
                links   = [SDE.r('Species')],
                basis   = ['expt'],
                constr  = EQ(SDE['property'], 'lattice parameter'))

    vpb = PyBlock(make_volume,
                   args = [vq['v'],vq['s']])
    vol_target =                                                                \
        Gen(name    = 'vol_target',
            desc    = '?',
            query   = vq,
            funcs   = [vpb],
            actions = [Expt(expt    = vq['d'],
                         expt_volume = vpb['out'])])

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

    ceq = Query({'d'   : Expt.id,
                 'a_c' : make_dict('contribs'),
                 'a_e' : make_dict('energy')},
                links = [Species_comp.r('Species'),
                             Reference.r('Element'),
                             Reference.r('Calc')],
                basis   = ['expt'],
                aggcols = [Expt.id])

    cohesive_elems =                                                            \
        Gen(name    = 'cohesive_elems',
            desc    = '?',
            query   = ceq,
            actions = [Expt(expt        = ceq['d'],
                        atomic_contribs = ceq['a_c'],
                        atomic_energies = ceq['a_e'])])
    ########################################################################

    vol_pa = Cell['volume'] / Struct['n_atoms']
    gapq = Query({'b'   : Bulk_job.id,
                  'gap' : ABS(vol_pa - Expt['volume_pa']),
                  'dv'  : vol_pa - Expt['volume_pa']},
                  basis = ['bulk_job'])
    gap =                                                                       \
        Gen(name    = 'gap',
            desc    = 'Computes deviation from optimum for each Bulk_job singlepoint',
            actions = [Bulk_job(bulk_job = gapq['b'],
                                dv       = gapq['dv'],
                                gap      = gapq['gap'])],
            query   = gapq)

    ########################################################################
    minq = Query({'e':Expt.id,
                  'm':MIN(Bulk_job['gap'])},
                 basis    = ['expt'],
                 links   = [BJE],
                 aggcols = [Expt.id])
    mingap =                                                                    \
        Gen(name    = 'mingap',
            desc    = 'Stores a value under "Expt" unique to the best corresponding Bulk_job' ,
            query   = minq,
            actions = [Expt(expt    = minq['e'],
                            min_gap = minq['m'])])


    ########################################################################

    ratio = Struct['n_atoms'](BJ) / Species['n_atoms']
    co = Query(exprs   = {'d'  : Expt.id,
                          'j'  : Bulk_job.id,
                          'be' : Job['energy'](BJ),
                          'bc' : Job['contribs'](BJ),
                          'br' : ratio},
                links  = [BJE],
                basis  = ['expt'],
                constr = EQ(Bulk_job['gap'], Expt['min_gap']))

    cohesive_optjob =                                                           \
        Gen(name    = 'cohesive_optjob',
            desc    = '???',
            query   = co,
            actions = [Expt(expt      = co['d'],
                             best_job      = co['j'],
                             bulk_energy   = co['be'],
                             bulk_contribs = co['bc'],
                             bulk_ratio    = co['br'])])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [bulkexpt, pop_aggs, all_vols, near_min, eform, kpts, eos, lattice,
            cohesive, cohesive_target, bm_target, vol_target, cohesive_elems,
            gap, mingap, cohesive_optjob,]

    mod.add(gens)
