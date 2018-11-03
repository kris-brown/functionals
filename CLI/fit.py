from os                 import environ
from argparse           import ArgumentParser
from distutils.util     import strtobool
# Internal
from dbgen.support.misc import ConnectInfo as Conn
from functionals.fit    import Fit
##############################################################################
# Utility functions
#------------------

def main(args:dict)->None:
    db  = Conn.from_file('/Users/ksb/Documents/JSON/functionals.json')
    fit = Fit(db, size=args['size'], bound = args['bound'], norm = args['norm'],
              initfit = args['initfit'], gridden = args['gridden'],
              maxiter = args['maxiter'])
    fit.submit(environ["FITPATH"])

if __name__ == '__main__':
    parser = ArgumentParser(description  = 'Submit fitting jobs',
                            allow_abbrev = True)
    parser.add_argument('--size', default = 5, type = int,
                        help = 'NxN size of fitted functional')
    parser.add_argument('--maxiter', default = 1000, type = int,
                        help = 'Max number of steps')
    parser.add_argument('--norm', default = 0.1, type = float,
                        help = 'Regularization parameter ')
    parser.add_argument('--bound', default = 0.1, type = int,
                        help = 'Max |value| of fitted coefs (except <1,1>)')
    parser.add_argument('--initfit', default = True, type  = strtobool,
                        help  = 'Initialize with unconstrained least squares fit')
    parser.add_argument('--gridden', default = 5, type = int,
                        help = 'Density of points for linear constraints')


    main(vars(parser.parse_args()))
