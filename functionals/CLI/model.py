# External
from os import environ
from os.path import join
# Internal Modules
from dbgen             import parser, ConnectInfo
from functionals.model  import make_model
"""
Run the model defined in /functionals/model
"""
################################################################################

root = join(environ['FUNCTIONALS_ROOT'],'data/')

def main(args:dict)->None:
    """
    Run the model with no extensions from command line.
    """

    m    = make_model()
    db   = ConnectInfo.from_file(root+'functionals.json')
    mdb  = ConnectInfo.from_file(root+'functionals_log.json')

    m.run(db, mdb, **args)

if __name__=='__main__':
    args = parser.parse_args()
    main(vars(args))
