# External Modules
from typing import Any,Type,TYPE_CHECKING
from os     import environ

# Internal Modules
from dbgen.support.model          import Model,new_model
from dbgen.core.parsing           import parser
from dbgen.support.misc           import ConnectInfo
from functionals.model.objects    import add_objects
from functionals.model.generators import add_generators

"""
Work in progress

Preliminary sketches for a DB to be used in association with functional
development

"""
################################################################################

def make_model()->Type[Model]:
    m = new_model('johannes') # type: Type['Model']
    add_objects(m)
    add_generators(m)
    return m

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
