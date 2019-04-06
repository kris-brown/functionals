# External Modules
from json   import loads

# Internal Modules
from dbgen import (Model, Gen, PyBlock, Query, Expr, AND, GROUP_CONCAT,
                   SUM, COUNT, Literal as Lit, NOT, NULL, CONCAT, ABS, MIN, EQ, NE, GT,
                    defaultEnv, Env,Import, Constraint, JPath, R2)

from functionals.scripts.load.get_kptden             import get_kptden
from functionals.scripts.load.get_kpts_vasp          import get_kpts_vasp
from functionals.scripts.misc.eos_func               import eos_func
from functionals.scripts.misc.conventional_lattice   import conventional_lattice

################################################################################
################################################################################
################################################################################

def make_volume(val : str,sym : str) -> float:
    sg = int(sym.split('_')[-1])
    if sg < 168: raise ValueError
    elif sg < 195:
        return float(val)**2 * 1.633
    else:
        return float(val)**3

def get_near_min(dv : float, all_vols : str) -> bool:
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

    return m in range(min_n - 2, min_n + 3)

################################################################################
################################################################################
################################################################################

def analysis(mod : Model) -> None:
    # Extract tables
    #---------------
    tabs = ['job','atom','element','struct','calc','cell',
            'species','expt','bulk_job','reference','species_dataset',
            'species_comp','species_dataset_element','functional','beef',
            'expt_refs']

    Job, Atom, Element, Struct, Calc, Cell, Species, Expt, Bulk_job,\
    Reference, Species_dataset,Species_comp,Species_dataset_element,Functional,\
    Beef,Expt_refs = map(mod.get, tabs) # type: ignore

    job__calc, bulk_job__job,struct__species, job__struct,bulk_job__expt,   \
        struct__cell,reference__calc, expt__calc,reference__element,        \
        species_comp__element, species_comp__species, expt__species,        \
        reference__job,species_dataset_element__species,calc__functional,   \
        beef__functional,expt_refs__reference,expt_refs__expt =             \
            map(mod.get_rel,[
                Job.r('calc'),Bulk_job.r('job'),Struct.r('species'), Job.r('struct'),
                Bulk_job.r('Expt'),Struct.r('cell'),Reference.r('calc'), Expt.r('calc'),
                Reference.r('element'), Species_comp.r('element'), Species_comp.r('species'),
                Expt.r('species'),Reference.r('job'),Species_dataset_element.r('species'),
                Calc.r('functional'),Beef.r('functional'),Expt_refs.r('reference'),
                Expt_refs.r('expt')])

    # Abbreviations
    #--------------
    SDE = Species_dataset_element
    BJE = bulk_job__expt
    BJ  = bulk_job__job

    xsp = JPath("species",[expt__species])
    xc  = JPath("functional",[calc__functional,expt__calc])
    bjj = JPath('job',[bulk_job__job])

    ############################################################################

    #c = Constraint(Calc)
    #c.find(mod,['bulk_job'],quit=False)

    bjc = JPath("calc", [job__calc, bulk_job__job])
    bfxc = JPath("beef",[beef__functional,calc__functional,job__calc,bulk_job__job])

    #c = Constraint(Species,[Struct.r('species')])
    #c.find(mod,['bulk_job'])

    bjsp = JPath("species", [struct__species, job__struct, bulk_job__job])
    bjst = JPath("struct", [job__struct, bulk_job__job])
    beq  = Query(exprs = dict(b = Bulk_job.id(),
                              c = Calc.id(bjc),
                              s = Species.id(bjsp),
                              n = Struct['n_atoms'](bjst),
                              ),#x = Beef.id(bfxc)),
                 basis   = ['bulk_job'],
                 )#ption  = [Beef.r('functional')])

    iexpt = Expt(insert  = True,
                 n_atoms = beq['n'],
                 species = beq['s'],
                 calc    = beq['c'])
    bulkexpt =                                                                  \
        Gen(name    = 'bulkexpt',
            desc    = 'All pairs of (bulk) species + calc., if mBEEF calc'\
                      'Links the bulk_job table to Expt table, too',
            query   = beq,
            actions = [Bulk_job(bulk_job = beq['b'],
                                expt     = iexpt)])
    ############################################################################

    refpth = JPath("reference", [[[reference__calc, expt__calc],
                                  [reference__element, species_comp__element, species_comp__species, expt__species]]])

    erq = Query(exprs = dict(e = Expt.id(),r = Reference.id(refpth)),
                basis = [Expt])
    exptref =                                                                   \
        Gen(name    = 'expt_refs',
            desc    = 'Populates expt_refs mapping table',
            query   = erq,
            actions = [Expt_refs(insert    = True,
                                 expt      = erq['e'],
                                 reference = erq['r'])])

    ############################################################################

    #c = Constraint(Job,[BJE])
    #c.find(mod,['expt'],links   = [BJE],quit=False)

    jx = JPath("job", [bulk_job__job, bulk_job__expt])
    bjx = JPath("bulk_job", [bulk_job__expt])

    #c = Constraint(Cell,[BJE])
    #c.find(mod,['expt'],links   = [BJE])

    cx = JPath("cell", [struct__cell, job__struct, bulk_job__job, bulk_job__expt])

    aggs = ['contribs','energies','volumes']

    agg_dict = dict(contribs = Job['contribs'](jx),
                    energies = Job['energy'](jx),
                    volumes  = Cell['volume'](cx))

    lbra,rbra = [Lit(x) for x in '[]']

    paq = Query(exprs = {'e':Expt.id(),
                         **{x : CONCAT(lbra, GROUP_CONCAT(y), rbra)
                             for x,y in agg_dict.items()}},
                aggcols = [Expt.id()],
                basis   = ['expt'],
                constr  = Bulk_job['near_min'](bjx),
                aconstr = COUNT(Lit(1)) |EQ| Lit(5))

    pop_aggs =                                                                  \
        Gen(name    = 'pop_aggs',
            desc    = 'Store information about closest 5 singlepoints near optimum',
            query   = paq,
            actions = [Expt(expt=paq['e'],**{x:paq[x] for x in aggs})])

    ########################################################################

    allvolcol = CONCAT(lbra, GROUP_CONCAT(Bulk_job['dv'](bjx)), rbra)

    avq = Query(exprs   = {'e' : Expt.id(), 'av' : allvolcol},
                basis   = ['expt'],
                aggcols = [Expt.id()] )
    all_vols =                                                                  \
        Gen(name    = 'all_vols',
            desc    = 'Store all volumes',
            query   = avq,
            actions = [Expt(expt = avq['e'], all_vols = avq['av'])])

    ############################################################################

    bjexpt = JPath("expt",[bulk_job__expt])

    nmq  = Query(exprs = {'b':Bulk_job.id(),
                          'dv':Bulk_job['dv'](),
                          'av':Expt['all_vols'](bjexpt)},
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

    ################################################################################

    # xbq = Query(exprs  = dict(e = Expt.id(bjexpt),
    #                           x = GROUP_CONCAT(Job['contribs'](bjj))),
    #             basis  = [Bulk_job],
    #             constr = Bulk_job['near_min'](),
    #             aggcols = [Expt.id(bjexpt)] )
    # popxbulk =                                                                  \
    #     Gen(name='xbulk',
    #         desc='Populates exchange contributions from 5 bulk jobs near minimum',
    #         query = xbq,
    #         actions = [Expt(expt=xbq['e'],xbulk=xbq['x'])])
    #
    # popdxba =                                                                   \
    #     Gen(name = 'dxba',
    #         desc = 'Populates dx_bulk_atom')
    #         #Serialized vector of the difference in exchange energy contributions between (optimum) bulk and (stoichiometrically-weighted) atoms

    ############################################################################

    # branch = [Constraint(Reference,[Reference.r('Calc')]),
    #                                   Constraint(Reference,[Reference.r('Element')]
    # findargs = dict(m = mod, basis = ['expt'],
    #                 links = [Species_comp.r('Species'), Reference.r('Calc'), Reference.r('Element')])
    # c=Constraint(Reference, branch = branch)])
    # c=Constraint(Job, [Reference.r('job')], branch = branch)
    # c.find(**findargs)

    # Paths
    refpth = JPath("reference", [[[reference__calc, expt__calc], [reference__element, species_comp__element, species_comp__species, expt__species]]])
    jref   = JPath("job", [reference__job, [[reference__calc, expt__calc], [reference__element, species_comp__element, species_comp__species, expt__species]]])
    spcref =  JPath('species_comp',[species_comp__species, expt__species])

    # Attributes
    Num    = Species_comp['num'](spcref)
    Eng    = Reference['energy'](refpth)
    Norm   = SUM(Eng * Num) / SUM(Num)

    eq = Query(exprs    = {'ex':Expt.id(),
                           'e':Expt['energy_pa']() - Norm},
                aggcols = [Expt.id()],
                basis   = ['expt'])
                #constr  = NOT(NULL(Job['contribs'](jref))), #Edge case - what if job has an energy but x_contribs failed?
                #)#aconstr = COUNT(Lit(1)) |EQ| Species['n_elems'](xsp))


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
    jcell = JPath('cell',[struct__cell,job__struct])
    kq = Query(exprs = {'j'       : Job.id(),
                        'stordir' : Job['stordir'](),
                        **{x:Cell[x](jcell) for x in Cell.ids()}},
               basis = ['job'])

    kpb0 =  PyBlock(get_kpts_vasp,
                    env      = p_env,
                    args     = [kq['stordir']],
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
    bjcell   = JPath('cell',[struct__cell, job__struct, bulk_job__job])
    eos_cols = ['volume_pa','energy_pa', 'bulkmod',  'img','complete']
    eos_env  = Env(Import('string','ascii_lowercase'),
                    Import('numpy','polyfit','poly1d','mean','std','array'),
                    Import('random','choices'),
                    Import('os','environ','remove'),
                    Import('os.path','join'),
                    Import('base64','b64encode'),
                    Import('ase.eos','EquationOfState'),
                    Import('ase.units','kJ'),
                    ) + defaultEnv

    eosq = Query(exprs={'ex' : Expt.id(bjexpt),
                        'v'  : GROUP_CONCAT(Cell['volume'](bjcell)),
                        'e'  : GROUP_CONCAT(Job['energy'](bjj)),
                        'n'  : Expt['n_atoms'](bjexpt),
                        'ndata' : COUNT(Lit(1))},
                basis   = ['bulk_job'],
                aggcols = [Expt.id(bjexpt)])

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
    lq = Query(exprs={'e' : Expt.id(),
                      's' : Species['symmetry'](xsp),
                      'v' : Expt['volume_pa']()},
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

    bxc  = xc + JPath("beef",[beef__functional])

    cq = Query(exprs={'e'     : Expt.id(),
                      'coefs' : Functional['data'](xc),
                      'name'  : Species['nickname'](xsp),
                      'comp'  : Species['composition'](xsp),
                      'ev'    : Expt['energies'](),
                      'vv'    : Expt['volumes'](),
                      'cv'    : Expt['contribs'](),
                      'x'     : Beef.id(bxc)},
               basis   = ['expt'],
               constr  = NOT(NULL(Expt['eform']()))
                            |AND| (Expt['n_data']() |GT| Lit(5)))

    cohesive =                                                                  \
        Gen(name    = 'cohesive',
            desc    = 'Populates more things in experiment',
            query   = cq,
            actions = [Expt(expt        = cq['e'],
                         coefs          = cq['coefs'],
                         name           = cq['name'],
                         composition    = cq['comp'],
                         energy_vector  = cq['ev'],
                         volume_vector  = cq['vv'],
                         contrib_vector = cq['cv'])])

    ########################################################################

    #c = Constraint(SDE)
    #c.find(mod, basis=['expt'], links = [SDE.r('species')])
    xsde = JPath("species_dataset_element", [species_dataset_element__species, expt__species])

    sde_dict = {'cohesive energy' : 'expt_cohesive_energy',
                'bulk modulus'    : 'expt_bm'}

    twogens = []
    for x,y in sde_dict.items():
        q = Query(exprs = dict(d = Expt.id(), v = SDE['value'](xsde)),
                  basis = ['expt'],
                  constr= EQ(SDE['property'](xsde), Lit(x)))

        twogens.append(Gen(name    = x.replace(' ','_'),
                           desc    = 'Copying %s info from Kelddata into dataset'%x,
                           query   = q,
                           actions = [Expt(expt = q['d'], **{y:q['v']})]))

    cohesive_target,bm_target = twogens # unpack

    ########################################################################
    vq = Query(exprs    = {'d' : Expt.id(),
                           'v' : SDE['value'](xsde),
                           's' : Species['symmetry'](xsp)},
                basis   = ['expt'],
                constr  = SDE['property'](xsde) |EQ| Lit('lattice parameter'))

    vpb = PyBlock(make_volume, args = [vq['v'], vq['s']],tests = [(('3','x_225'),27)])
    vol_target =                                                                \
        Gen(name    = 'vol_target',
            desc    = 'Literature value of volume for this experiement',
            query   = vq,
            funcs   = [vpb],
            actions = [Expt(expt    = vq['d'],
                         expt_volume = vpb['out'])])

    ########################################################################

    #links = [Species_comp.r('Species'), Reference.r('Element'), Reference.r('Calc')]
    #c = Constraint(Element,[Species_comp.r('Element')])
    #c.find(mod,['expt'],links = links)

    err  = JPath("reference", [expt_refs__reference])
    erre = err + JPath("element",[reference__element])
    ere  = JPath("expt", [expt_refs__expt])
    def make_dict(attr : str) -> Expr:
        '''
        Makes a serialized python dictionary mapping element_ids to Job
        properties of their reference calculation
        '''
        elemattr = Element['atomic_number'](erre)
        jobattr  = Reference[attr](err)
        group    = GROUP_CONCAT(CONCAT(elemattr,Lit(':'),jobattr))

        return CONCAT(Lit('{'),group,Lit('}'))

    ceq = Query({'d'   : Expt.id(ere),
                 'a_c' : make_dict('contribs'),
                 'a_e' : make_dict('energy')},
                basis   = [Expt_refs],
                aggcols = [Expt.id(ere)])

    cohesive_elems =                                                            \
        Gen(name    = 'cohesive_elems',
            desc    = '?',
            query   = ceq,
            actions = [Expt(expt            = ceq['d'],
                            atomic_contribs = ceq['a_c'],
                            atomic_energies = ceq['a_e'])])
    ########################################################################

    vol_pa = Cell['volume'](bjcell) / Struct['n_atoms'](bjst)
    gapq = Query({'b'   : Bulk_job.id(),
                  'gap' : ABS(vol_pa - Expt['volume_pa'](bjexpt)),
                  'dv'  : vol_pa - Expt['volume_pa'](bjexpt)},
                  basis = ['bulk_job'])
    gap =                                                                       \
        Gen(name    = 'gap',
            desc    = 'Computes deviation from optimum for each Bulk_job singlepoint',
            actions = [Bulk_job(bulk_job = gapq['b'],
                                dv       = gapq['dv'],
                                gap      = gapq['gap'])],
            query   = gapq)

    ########################################################################
    #links   = [BJE],
    minq = Query({'e':Expt.id(), 'm':MIN(Bulk_job['gap'](bjx))},
                 basis    = ['expt'],
                 aggcols = [Expt.id()])
    mingap =                                                                    \
        Gen(name    = 'mingap',
            desc    = 'Stores a value under "Expt" unique to the best corresponding Bulk_job' ,
            query   = minq,
            actions = [Expt(expt    = minq['e'],
                            min_gap = minq['m'])])


    ########################################################################
    exptstructs = JPath("struct", [job__struct, bulk_job__job,bulk_job__expt])
    exptsjobs   = JPath("job",[bulk_job__job,bulk_job__expt])

    ratio = Struct['n_atoms'](exptstructs) / Species['n_atoms'](xsp)
    co = Query(exprs   = {'d'  : Expt.id(),
                          'j'  : Bulk_job.id(bjx),
                          'be' : Job['energy'](exptsjobs),
                          'bc' : Job['contribs'](exptsjobs),
                          'br' : ratio},
                basis    = ['expt'],
                constr   = EQ(Bulk_job['gap'](bjx), Expt['min_gap']()),
                opt_attr = [Job['contribs'](exptsjobs)])

    cohesive_optjob =                                                           \
        Gen(name    = 'cohesive_optjob',
            desc    = '???',
            query   = co,
            actions = [Expt(expt          = co['d'],
                            best_job      = co['j'],
                            bulk_energy   = co['be'],
                            bulk_contribs = co['bc'],
                            bulk_ratio    = co['br'])])

    ########################################################################
    fpth  = JPath('calc',[expt__calc])
    qr2bm = Query(exprs = dict(c = Calc.id(fpth),
                               xb = GROUP_CONCAT(Expt['expt_bm']()),
                               b  = GROUP_CONCAT(Expt['bulkmod']()),
                               xl = GROUP_CONCAT(Expt['expt_lattice']())),
                  basis = [Expt],
                  aggcols   = [Calc.id(fpth)])

    r2bm =                                                                      \
        Gen(name    = 'r2bm',
            desc    = 'Calculates R2 for bulk modulus for functionals',
            query   = qr2bm,
            actions = [Calc(calc=qr2bm['c'],r2_bm=qr2bm['r'])])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [bulkexpt, pop_aggs, all_vols, near_min, eform, kpts, eos, lattice,
            cohesive, cohesive_target, bm_target, vol_target, cohesive_elems,
            gap, mingap, cohesive_optjob, exptref, r2bm]

    mod.add(gens)
