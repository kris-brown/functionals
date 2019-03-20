# External Modules
from typing     import Tuple as T, List as L
from re         import search
from ase.data   import chemical_symbols # type: ignore
from ast        import literal_eval

# Internal Modules
from dbgen import (Model, Gen, Query, PyBlock, AND, Env, defaultEnv, Import,
                     Literal as Lit,  EQ, JPath, Constraint)

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
from functionals.scripts.atoms.traj_to_json     import traj_to_json

'''
Pure functions to extract information loaded by IO actions
'''

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

eng_env = Env(Import('re','search'))

def eng_gpaw(s:str)->float:
    '''Parse the result free energy from a GPAW logfile'''
    pat = r'Free energy:\s+([-+]?\d+\.\d+)'
    from re import findall
    match = findall(pat, s); assert match
    return float(match[-1])

def eng_vasp(s:str)->float:
    '''Parse the result free energy from a GPAW logfile'''
    from re import findall
    pat = r'TOTEN\s+=\s+([-+]?\d+\.\d+)'
    match = findall(pat, s); assert match
    out = float(match[-1])
    return out
##############################################################################
##############################################################################
##############################################################################

def load(mod : Model) -> None:
    # Extract tables
    tabs = ['job','atom','element','struct','calc','cell','pure_struct',
            'species','bulk_job','reference','species_comp','species_dataset',
            'species_dataset_element']

    Job, Atom, Element, Struct, Calc, Cell, Pure_struct, Species, Bulk_job,\
    Reference, Species_comp, Species_dataset,         \
    Species_dataset_element = map(mod.get, tabs) # type: ignore

    struct__pure_struct,job__calc,atom__struct, job__struct = map(mod.get_rel,[
        Struct.r('pure_struct'),Job.r('calc'),Atom.r('struct'),Job.r('struct')])

    sdq   = Query(exprs={'j':Job.id(),    'sd'  : Job['stordir']()})
    logq  = Query(exprs={'j':Job.id(),    'log' : Job['log']()})
    rawq  = Query(exprs={'s':Struct.id(), 'raw' : Struct['raw']()})
    spcq  = Query(exprs={'s':Species.id(),'c'   : Species['composition']()})

    ########################################################################

    psq = Query(exprs = {'p'  : Pure_struct.id(),
                         'pt' : Pure_struct['prototype']()})

    pspb = PyBlock((lambda x: {'AB_1_a_b_225'       : 'rocksalt',
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
                    args=[psq['pt']])
    psnick =                                                                    \
        Gen(name    = 'psnick',
            desc    = 'Use a dict to map nicknames to some Pure structs',
            actions = [Pure_struct(pure_struct = psq['p'],
                                   nickname    = pspb['out'])],
            query   = psq,
            funcs   = [pspb],
            tags    = ['species'])

    ########################################################################
    anytraj_env = Env(Import('ase.io','read'),
                       Import('glob','glob'),
                       Import('os.path','getsize'))

    ttjenv = Env(Import('ase.constraints','FixAtoms'),
                 Import('ase','Atoms')) + defaultEnv

    spb1 = PyBlock(anytraj,
                   env  = anytraj_env,
                   args = [sdq['sd']])
    spb2 = PyBlock(traj_to_json,
                    env  =  ttjenv,
                    args = [spb1['out']])

    struct =                                                                    \
        Gen(name    = 'struct',
            desc    = "Assumes there's only one traj in the directory",
            actions = [Job(job    = sdq['j'],
                           struct = Struct(insert = True,
                                           raw    = spb2['out']))],
            query   = sdq,
            funcs   = [spb1,spb2])

    ########################################################################

    jepb = PyBlock(eng_vasp,
                   env  = eng_env,
                   args = [logq['log']])

    jobeng = Gen(name    = 'jobeng',
                 desc    = '?',
                 actions = [Job(job    = logq['j'],
                                energy = jepb['out'])],
                 query   = logq,
                 funcs   = [jepb])

    ########################################################################

    md_env = Env(Import('os', 'stat'),
                 Import('pwd', 'getpwuid'),
                 Import('os.path', 'getmtime')) + defaultEnv

    mdpb = PyBlock(metadata,
                   env      = md_env,
                   args     = [sdq['sd']],
                   outnames = ['u','t'])

    jobmetadata =                                                               \
        Gen(name    = 'jobmetadata',
            desc    = 'Scrape metadata about job',
            actions = [Job(job       = sdq['j'],
                           user      = mdpb['u'],
                           timestamp = mdpb['t'])],
            query   = sdq,
            funcs   = [mdpb])

    ########################################################################

    atomdetails = ['num','x','y','z','const','mag','tag','ind']
    apb = PyBlock(get_atoms, args= [rawq['raw']],
                  outnames = atomdetails)

    iatom = Atom(insert = True,
                 element= Element(atomic_number = apb['num']),
                 struct = rawq['s'],magmom=apb['mag'], ind = apb['ind'],
                 number = apb['num'], constrained=apb['const'],
                 x=apb['x'],y=apb['y'],z=apb['z'],tag=apb['tag'])
    atom =                                                                      \
        Gen(name    = 'atom',
            desc    = 'Initializes atoms from knowing n_atoms',
            actions = [iatom],
            query   = rawq,
            funcs   = [apb])

    ########################################################################
    epb = PyBlock(lambda:list(range(1,100)),outnames=['z'])

    elems =                                                                     \
        Gen(name    = 'elems',
            desc    = 'Seeds the Element table by providing atomic numbers',
            actions = [Element(insert=True,
                               atomic_number=epb['z'])],
            funcs   = [epb])


    ########################################################################
    cpb = PyBlock(get_cell, args=[rawq['raw']], outnames = Cell.ids())

    cell =                                                                      \
        Gen(name    = 'cell',
            desc    = 'Populate cells from Structs',
            actions = [Struct(struct = rawq['s'],
                              cell   = Cell(insert = True,
                                        **{x:cpb[x] for x in Cell.ids()}))],
            query   = rawq,
            funcs   = [cpb])

    ########################################################################
    cellq   = Query(exprs = {'c':Cell.id(), **{x : Cell[x]() for x in Cell.ids()}})
    cpbout  = ['surface_area', 'volume', 'a', 'b', 'c']
    cellpb  = PyBlock(cell_info,
                      args     = [cellq[x] for x in Cell.ids()],
                      outnames = cpbout)
    celinfo =                                                               \
        Gen(name    = 'celinfo',
            desc    = 'Basic geometric formulas applied to 3x3 cell representation',
            actions = [Cell(cell=cellq['c'],
                            **{x:cellpb[x] for x in cpbout})],
            query   = cellq,
            funcs   = [cellpb])

    ########################################################################
    jtt_env = Env(Import('ase.constraints','FixAtoms'),
                  Import('ase','Atoms')) + defaultEnv

    jtt = PyBlock(json_to_traj,env=jtt_env,args=[rawq['raw']],outnames=['t'])
    gst = PyBlock(get_system_type,args=[jtt['t']],outnames=['st'])
    systype =                                                                   \
        Gen(name    = 'systype',
            desc    = 'Apply function to identify whether system is bulk/mol/surf',
            actions = [Struct(struct      = rawq['s'],
                              system_type = gst['st'])],
            query   =  rawq,
            funcs   =  [jtt,gst])

    ########################################################################
    gb_env = Env(Import('bulk_enumerator.bulk', 'BULK'),
                  Import('ase.io','write'),
                  Import('os','remove','environ'),
                  Import('os.path','join'),
                  Import('random','choices'),
                  Import('string', 'ascii_lowercase'))

    psq = Query(exprs  = {'s': Struct.id(), 'raw' : Struct['raw']()},
                constr = EQ(Struct['system_type'](), Lit('bulk')))

    psb1 = PyBlock(json_to_traj,env=jtt_env,args=[psq['raw']],outnames=['traj'])
    psb2 = PyBlock(get_bulk,env=gb_env,args=[psb1['traj']],outnames=['b'])
    psb3 = PyBlock(get_pure_struct,args=[psb2['b']],outnames=['ptype'])

    ips  = Pure_struct(insert    = True,
                       prototype = psb3['ptype'])
    ps =                                                                        \
        Gen(name    = 'ps',
            desc    = 'Determine the "pure structure" code for each Bulk',
            actions = [Struct(struct      = psq['s'],
                              pure_struct = ips)],
            query   = psq,
            funcs   = [psb1,psb2,psb3],
            tags    = ['long','parallel'])
    ########################################################################
    le_env = defaultEnv + Env(Import('ast', 'literal_eval'))
    cs_env = Env(Import('ase.data','chemical_symbols'))

    spnq  = Query(exprs={'sp' : Species.id(),
                         'c'  : Species['composition'](),
                         's'  : Species['symmetry']()})

    spnpb =  PyBlock(species_nick,
                     env      = le_env + cs_env,
                     args     = [spnq['c'],spnq['s']],
                     outnames = ['nn'])
    spnick =                                                                    \
        Gen(name    = 'spnick',
            desc    = "Determine species' nicknames by analyzing their "        \
                       "composition and symmetry" ,
            actions = [Species(species  = spnq['sp'],
                               nickname = spnpb['nn'])],
            query   = spnq,
            funcs   = [spnpb])
    ########################################################################
    def getNum(d:str) -> T[L[int], L[int]]:
        elems, counts = zip(*literal_eval(d).items())
        return list(elems), list(counts)

    spf = PyBlock(getNum, env=le_env, args = [spcq['c']], outnames = ['z','num'])
    sp_comp =                                                                  \
        Gen(name    = 'sp_comp',
            desc    = 'Populates species-composition mapping table',
            actions = [Species_comp(insert  = True,
                                    num     = spf['num'],
                                    species = spcq['s'],
                                    element = Element(atomic_number=spf['z']))],
            query   = spcq,
            funcs   = [spf])

    ########################################################################

    sps = JPath(Pure_struct,[struct__pure_struct])

    pops_q = Query(exprs={'str' : Struct.id(),
                          'c'   : Struct['composition_norm'](),
                          's'   : Pure_struct['prototype'](sps),
                          'n'   : Struct['n_elems']()},
                    basis = ['struct'])

    ispecies = Species(insert      = True,
                       symmetry    = pops_q['s'],
                       composition = pops_q['c'],
                       n_elems     = pops_q['n'])
    pop_species =                                                               \
        Gen(name    = 'pop_species',
            desc    = 'Populate species from Struct and Purestruct tables',
            actions = [Struct(struct  = pops_q['str'],
                              species = ispecies)],
            query   = pops_q)

    ############################################################################

    def sum_values(x:str)->int:
        return sum(dict(literal_eval(x)).values())

    snapb = PyBlock(sum_values, env = le_env,args=[spcq['c']],outnames=['n'])

    species_natoms =                                                            \
        Gen(name    = 'species_natoms',
            desc    = 'Total number of atoms in the most reduced stoichiometric ratios',
            actions = [Species(species=spcq['s'],n_atoms=snapb['n'])],
            query   = spcq,
            funcs   = [snapb])
    ############################################################################

    #c = Constraint(Atom)
    #c.find(mod,basis=['job'],links=[Atom.r('struct')],quit=False)

    jc    = JPath(Calc,[job__calc])
    jatom = JPath("atom", [atom__struct, job__struct])
    sj    = JPath("struct", [job__struct])

    refq = Query(exprs  = dict(j = Job.id(),
                               c = Calc.id(jc),
                               z = Atom['number'](jatom),
                               e = Job['energy'](),
                               t = Job['contribs']()),
                 basis  = ['job'],
                 constr = EQ(Struct['system_type'](sj), Lit('molecule'))
                              |AND| EQ(Struct['n_atoms'](sj),Lit(1)))
    refr =                                                                      \
        Gen(name    = 'refr',
            desc    = 'Collects calculations of isolated atoms',
            query   = refq,
            actions = [Reference(insert   = True,
                                 job      = refq['j'],
                                 calc     = refq['c'],
                                 element  = Element(atomic_number = refq['z']),
                                 energy   = refq['e'],
                                 contribs = refq['t'])])

    ########################################################################

    bjq = Query(exprs  = {'j':Job.id()},
                basis  = ['job'],
                constr = EQ(Struct['system_type'](sj), Lit('bulk')))
    bulkjob =                                                                   \
        Gen(name    = 'bulkjob',
            desc    = "Subset of jobs that are bulk calculations",
            actions = [Bulk_job(insert = True, job=bjq['j'])],
            query   =  bjq)

    ###########################################################################

    spinb = PyBlock(lambda x: int('Spin-polarized calculation' in x),
                    args    =[logq['log']])

    spinpol =                                                               \
        Gen(name    = 'spinpol',
            desc    = 'Populates Job.spinpol',
            actions = [Job(job     = logq['j'],
                           spinpol = spinb['out'])],
            query   = logq,
            funcs   = [spinb])

    ########################################################################

    compcols = ['n_atoms','n_elems','composition',
                     'composition_norm','metal_comp','str_symbols',
                     'str_constraints']
    ca_env = Env(Import('dbgen.utils.lists','normalize_list','nub'))+defaultEnv

    capb = PyBlock(countatm,
                   env      = ca_env,
                   args     =[rawq['raw']],
                   outnames = compcols)

    countatoms =                                                            \
        Gen(name    = 'countatoms',
            desc    = 'Get stoichiometric properties of Structure from raw',
            actions = [Struct(struct = rawq['s'],
                              **{k:capb[k] for k in compcols})],
            query   =  rawq,
            funcs   = [capb])

    ########################################################################
    ########################################################################
    ########################################################################

    gens = [psnick,struct,jobeng,jobmetadata,atom,elems,cell,celinfo,
            systype,ps,spnick,sp_comp,pop_species,species_natoms,refr,bulkjob,
            spinpol, countatoms]

    mod.add(gens)
