from typing     import Any, List as L, Dict as D, Tuple as T
from argparse   import ArgumentParser
from os.path    import join, exists
from subprocess  import check_output, call
from itertools  import product as prod
from os         import environ, listdir
from shutil     import copyfile, rmtree
from pathlib    import Path
from random     import choice
from ast        import literal_eval
from json       import load,dump, loads
from numpy      import array,inf,sum,vstack,ones # type: ignore
from numpy.linalg import inv   # type: ignore
from psycopg2 import connect # type: ignore
from functionals.fit.fit import Fit,sqlselect

'''
Submit fitting jobs, provided a DB connection to the DB with DFT jobs
'''
###############################################################################
# Constants
###########
home = '/'+join(*__file__.split('/')[:-3])

def str2bool(v:str)->bool:  return v.lower() in ("yes", "true", "t", "1")

def mkconn(db:str)->Any:
    with open(db,'r') as f:
        kwargs = load(f)
        kwargs['dbname']=kwargs.pop('db'); kwargs['password']=kwargs.pop('passwd')
        return connect(**kwargs)

# increase weight on magnetic lattice constants

spec = dict(
    consts    = ['lda hnorm pos liebox scan11'],
    reg       = [0.,0.8,1.],
    ce_scale  = [0.01], #[0,0.01],
    bm_scale  = [0.5],
    lc_scale  = [0.1], #[0,0.1,0.5]
    mag_scale = [1,10],
) # type: D[str,L]

# Queries
#------------
q2 = '''SELECT calc_id FROM calc
         JOIN functional ON functional = functional_id WHERE beef  --and done'''

def main(db_ : str, submit: bool, retry:bool) -> None:
    '''
    For every completed "expt" and every single set of fitting parameters,
    set up 5 fitting jobs
    '''
    pth  = join(home, 'data/fit')
    conn = mkconn(db_)

    calcs = sqlselect(conn,q2)
    params =  list(prod(*[set(x) for x in spec.values()]))
    for (constr,reg,ce,bm,lc,mc),(calc,) in sorted(list(prod(params,calcs))):
        fit = Fit.from_db(db=db_,calc_id=calc,constr=constr,reg=reg,ce=ce,bm=bm,lc=lc,mc=mc)
        root = join(pth,fit.metadata()['uid'][:10])
        Path(root).mkdir(parents=True, exist_ok=True)
        fit.write(root)
        fit2 = fit.from_json(root); print(fit==fit2) # CHECK SERIALIZATION
        if not fit==fit2: import pdb;pdb.set_trace(); assert False

    if submit:
        sub(1,False)

def sub(time:int, retry:bool) -> None:
    args = ['functionals.CLI.fit', time, 'T' if retry else 'F']
    cmds = ';'.join(['cds','cd functionals',
                     'source /nfs/slac/g/suncatfs/sw/py3.7.4/env.csh', 'act',
                     'python -m {} --sub=T --time={} --retry={} --local=F'.format(*args)])
    import pdb;pdb.set_trace()
    call('rsync -avz --exclude "*.env"  /Users/ksb/functionals ksb@suncatls1.slac.stanford.edu:/nfs/slac/g/suncatfs/ksb/', shell=True)
    call('ssh ksb@suncatls1.slac.stanford.edu "{}"'.format(cmds), shell=True)

def sub_suncat(time : int, retry : bool) -> None:
    pth     = join(home, 'data/fit')
    dirs    = listdir(pth)
    currstr = 'bjobs -o "exec_cwd jobid" -noheader'
    output  = check_output(currstr, encoding='UTF-8',shell=True).split('\n')
    working  = dict(map(str.split, filter(None,output)))  # type: ignore

    for d in dirs:
        dir    = join(pth,d)
        root   = '/nfs/slac/g/suncatfs/ksb/beefjobs/fit/'
        args   = [root+d,time,join(root+d,'subfit.sh')]
        inprog = root+d in working
        done   = exists(join(pth,dir,'x9.json'))
        if retry or (not (done or inprog)):
            if retry:
                if inprog: call('bkill '+working[d], shell=True)
                call('rm -rf '+root+d)
            call('cp -r {} {}'.format(dir, root), shell=True)
            act = 'cd {}; bsub -n 1 -W0{}:30 -q suncat {}'.format(*args)
            call(act, shell=True)


def sub_(local:str,time:int,retry:bool) -> None:
    assert time < 10 and time > 0
    if   local=='T': sub(time,retry)
    elif local=='F': sub_suncat(time,retry)
    else: raise ValueError(local)
# Parser
########
parser = ArgumentParser(description  = 'Submit some fitting jobs', allow_abbrev = True)
parser.add_argument('--db', type = str, help='JSON with DB connection info',
                    default = join(home, 'data/functionals.json'),)

parser.add_argument('--sub',default=False,
                    help='Add anything to do submit command instead')
parser.add_argument('--local',default='T',
                    help='Add anything to do submit command instead')

parser.add_argument('--create',default=False,
                    help='Add anything to do create command instead')

parser.add_argument('--time', type    = int, default = 1,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--retry', type    = str2bool, default = False,
                    help    = 'Ignore existence of result.json')

parser.add_argument('--clean', type    = bool, default = False,
                    help    = 'Remove things not in DB')

def clean(db:str) -> None:
    conn = mkconn(db)
    names  = set([x[0] for x in sqlselect(conn,'SELECT name FROM fit')])
    import pdb;pdb.set_trace()
    for root in [join(home,'data','fit'), '/Users/ksb/scp_tmp/fitplt']:
        for d in listdir(root):
            if d not in names:
                print('removing ' + join(root,d))
                rmtree(join(root,d))

    raise NotImplementedError('DELETE DIRS NOT IN NAMES')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.create:
        main(args.db,args.sub,args.retry)
    elif args.sub:
        sub_(args.local,args.time,args.retry)
    elif args.clean:
        clean(args.db)
    else:
        raise ValueError("Run with CREATE SUB or MOVE")
