from typing     import Any, List as L, Dict as D, Tuple as T
from argparse   import ArgumentParser
from os.path    import join, exists
from subprocess  import check_output
from itertools  import product as prod
from os         import environ, system, listdir
from shutil     import copyfile
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
home = __file__.split('functionals')[0]+'functionals/'

# Queries
#------------
q1 = 'SELECT fitparams_id  FROM fitparams'
q2 = '''SELECT calc_id FROM calc
         JOIN functional ON functional = functional_id WHERE beef  --and done'''

def main(db_ : str, sub: bool, retry:bool) -> None:
    '''
    For every completed "expt" and every single set of fitting parameters,
    set up 5 fitting jobs
    '''
    pth  = join(home, 'data/fit')

    with open(db_,'r') as f:
        kwargs = load(f)
        kwargs['dbname']=kwargs.pop('db'); kwargs['password']=kwargs.pop('passwd')
        conn = connect(**kwargs)

    params,calcs = [sqlselect(conn,x) for x in [q1,q2]]
    for (fp,),(calc,) in sorted(list(prod(params,calcs))):
        print(fp,calc)
        fit = Fit.from_db(db=db_,fp_id=fp,calc_id=calc)
        root = join(pth,fit.metadata()['uid'][:10])
        Path(root).mkdir(parents=True, exist_ok=True)
        fit.write(root)
        #fit2 = fit.from_json(root); print(fit==fit2)
        if sub and (retry or not exists(join(root,'result.json'))):
            act   = 'python runfit.py'
            cmd   = 'cd {}; '.format(root) + act
            system(cmd)

def sub(time : int, retry : bool) -> None:
    pth     = join(home, 'data/fit')
    dirs    = listdir(pth)
    currstr = 'bjobs -o exec_cwd -noheader'
    working = set([x.split('/')[-1] for x in check_output(currstr, encoding='UTF-8',shell=True).split()])

    for d in dirs:
        dir    = join(pth,d)
        inprog = d in working
        done   = exists(join(pth,dir,'result.json'))
        if retry or (not (done or inprog)):
            act   = 'bsub -n 1 -W{}:30 -q suncat3 subfit.sh'.format(time)
            cmd   = 'cd {}; '.format(dir) + act
            system(cmd)


# Parser
########
parser = ArgumentParser(description  = 'Submit some fitting jobs',
                        allow_abbrev = True)

parser.add_argument('--db', type = str,
                    default = join(home, 'data/functionals.json'),
                    help    = 'Path to JSON with DB connection info')

parser.add_argument('--sub',default=False,
                    help='Add anything to do submit command instead')

parser.add_argument('--create',default=False,
                    help='Add anything to do create command instead')

parser.add_argument('--time', type    = int, default = 0,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--retry', type    = bool, default = False,
                    help    = 'Ignore existence of result.json')

parser.add_argument('--move', type    = bool, default = False,
                    help    = 'Move jobs from scp_tmp/vauto/fit to /functionals/data/fit')

def move() -> None:
    s='/ksb/functionals/data/fit/'
    cmd = 'rsync -avz ksb@suncatls1.slac.stanford.edu:/nfs/slac/staas/fs1/g/suncat{0} /Users{0}'.format(s)
    system(cmd)

if __name__ == '__main__':
    args = parser.parse_args()
    if args.create:
        main(args.db,args.sub,args.retry)
    elif args.sub:
        sub(args.time,args.retry)
    elif args.move:
        move()
    else:
        raise ValueError("Run with CREATE SUB or MOVE")
