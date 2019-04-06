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
q2 = '''SELECT DISTINCT E.calc FROM expt E JOIN calc C on E.calc=C.calc_id
         JOIN functional F ON C.functional = F.functional_id WHERE F.beef'''

def main(db_ : str, pth : str) -> None:
    '''
    For every completed "expt" and every single set of fitting parameters,
    set up 5 fitting jobs
    '''


    with open(db_,'r') as f:
        kwargs = load(f)
        kwargs['dbname']=kwargs.pop('db'); kwargs['password']=kwargs.pop('passwd')
        conn = connect(**kwargs)

    params,calcs = [sqlselect(conn,x) for x in [q1,q2]]

    for (fp,),(calc,) in prod(params,calcs):
        for decay in range(5):
            fit = Fit.from_db(db=db_,fp_id=fp,calc_id=calc,decay=decay)
            root = join(pth,fit.uid()[:10],str(decay))
            Path(root).mkdir(parents=True, exist_ok=True)
            fit.write(root)

def sub(pth : str, time : int, retry : bool, local : bool) -> None:
    dirs = listdir(pth)
    for d in dirs:
        for dd in listdir(join(pth,d)):
            dir   = join(pth,d,dd)
            if retry or not exists(join(dir,'result.json')):
                sunc  = choice(['','','','2','2','3'])
                act   = 'python runfit.py' if local else 'bsub -n 1 -W{}:09 -q suncat{} subfit.sh'.format(time,sunc)
                cmd   = 'cd {}; '.format(dir) + act
                system(cmd)


# Parser
########
parser = ArgumentParser(description  = 'Submit some fitting jobs',
                        allow_abbrev = True)

parser.add_argument('--db', type = str,
                    default = join(home, 'data/functionals.json'),
                    help    = 'Path to JSON with DB connection info')

parser.add_argument('--pth', type = str,
                    default = join(home, 'data/fit'),
                    help    = 'Path to where fitting jobs will be performed')

parser.add_argument('--sub',
                    help='Add anything to do submit command instead')

parser.add_argument('--time', type    = int,
                    default = 1,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--retry', type    = bool,
                    default = False,
                    help    = 'Walltime for batch jobs')

parser.add_argument('--local', type    = bool,
                    default = False,
                    help    = 'Walltime for batch jobs')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.sub:
        sub(args.pth,args.time,args.retry,args.local)
    else:
        main(args.db,args.pth)