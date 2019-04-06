# External
from os import environ
from os.path import join
# Internal Modules
from dbplot.main  import main as plt, parser

'''
Visualize results from the DB

e.g. >> python CLI/plot.py --pltpth=plots/example.json
'''
################################################################################
def main(args:dict)->None:
    plt(args)

if __name__=='__main__':
    args = vars(parser.parse_args())
    args['db'] = join(environ['FUNCTIONALS_ROOT'],'data/functionals.json')
    main(args)