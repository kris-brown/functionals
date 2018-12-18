# Internal Modules
from dbgen2             import parser, ConnectInfo
from functionals.model  import make_model
"""
Run the model defined in /functionals/model
"""
################################################################################

root = '/Users/ksb/Documents/JSON/'

def main(args:dict)->None:
    """
    Run the model with no extensions from command line.
    """

    m    = make_model()
    db   = ConnectInfo.from_file(root+'functionals.json')
    mdb  = ConnectInfo.from_file(root+'functionals_log.json')

    for x in ['only','xclude']:
        args[x] = set(args[x].split())

    m.run(db, mdb, **args)

if __name__=='__main__':
    args = parser.parse_args()
    main(vars(args))
