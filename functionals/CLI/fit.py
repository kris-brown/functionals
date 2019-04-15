from typing     import Any, List as L, Dict as D, Tuple as T
from argparse   import ArgumentParser
from os.path    import join, exists
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
home = environ['FUNCTIONALS_ROOT']
# Queries
#------------
q1 = 'SELECT fitparams_id  FROM fitparams'
q2 = '''SELECT calc_id FROM calc
         JOIN functional ON functional = functional_id WHERE beef''' # also need to add WHERE calc.done once we have all the data


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

    for (fp,),(calc,) in list(prod(params,calcs)):
        print('\n',fp,calc)
        if fp!=7: continue
        fit = Fit.from_db(db=db_,fp_id=fp,calc_id=calc)
        root = join(pth,fit.metadata()['uid'][:10])
        Path(root).mkdir(parents=True, exist_ok=True)
        fit.write(root)
        if sub and (retry or not exists(join(root,'result.json'))):
            act   = 'python runfit.py'
            cmd   = 'cd {}; '.format(root) + act
            system(cmd)

def sub(time : int, retry : bool, local : bool) -> None:
    pth  = join(home, 'data/fit')
    dirs = listdir(pth)
    for d in dirs:
        dir   = join(pth,d)
        if retry or not exists(join(dir,'result.json')):
            sunc  = choice(['','','','','2','2','3'])
            act   = 'python runfit.py' if local else 'bsub -n 1 -W{}:30 -q suncat{} subfit.sh'.format(time,sunc)
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
                    help='Add anything to do submit command instead')

parser.add_argument('--time', type    = int,
                    default = 0,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--retry', type    = bool,
                    default = False,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--local', type    = bool,
                    default = False,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--move', type    = bool,
                    default = False,
                    help    = 'Move jobs from scp_tmp/vauto/fit to /functionals/data/fit')

def move()->None:
    r= '/Users/ksb/'
    r1,r2=[r+x for x in ['scp_tmp/vauto/fit','functionals/data/fit']]
    dirs = set( listdir(r2))
    for dir in listdir(r1):
        res = join(r1,dir,'result.json')
        if dir in dirs and exists(res):
            copyfile(res,join(r2,dir,'result.json'))

if __name__ == '__main__':
    args = parser.parse_args()
    if args.create:
        main(args.db,args.sub,args.retry)
    elif args.sub:
        sub(args.time,args.retry,args.local)
    elif args.move:
        move()
