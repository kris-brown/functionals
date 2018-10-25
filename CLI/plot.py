# Internal Modules
from dbgen.support.misc import ConnectInfo
from dbplot.main        import main as plt,parser
################################################################################
def main(args:dict)->None:
    plt(args)

if __name__=='__main__':
    args = vars(parser.parse_args())
    args['db'] = '/Users/ksb/Documents/JSON/functionals.json'
    main(args)
