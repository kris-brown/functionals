# External Modules
from os     import environ

from dbgen.core.parsing           import parser
from dbgen.support.misc           import ConnectInfo
from functionals.model            import make_model
"""
Work in progress

Preliminary sketches for a DB to be used in association with functional
development

"""
################################################################################



def main()->None:
    """
    Run the model with no extensions from command line.
    """
    args = parser.parse_args()
    m    = make_model()
    db   = ConnectInfo.from_file('/Users/ksb/Documents/JSON/functionals.json')
    mdb  = ConnectInfo.from_file('/Users/ksb/Documents/JSON/functionals_log.json')

    only,xclude = [set(x.split()) for x in [args.only,args.xclude]]

    m._run(db, mdb, nuke=args.nuke, add=args.add, retry=args.retry, only=only,
          xclude=xclude, start=args.start, until=args.until)

if __name__=='__main__':
    main()
