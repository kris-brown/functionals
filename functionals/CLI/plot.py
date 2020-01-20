# External
from os.path import join
# Internal Modules
from dbplot.main import main as plt, parser

'''
Visualize results from the DB

e.g. >> python CLI/plot.py --pltpth=plots/example.json
'''
##############################################################################

if __name__ == '__main__':
    args = vars(parser.parse_args())
    args['db'] = '/'+join(*__file__.split('/')[:-3], 'data/functionals.json')
    args['open'] = True
    plt(args)
