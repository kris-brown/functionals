from argparse       import ArgumentParser

###########################
# Command line parsing
########################
parser = ArgumentParser(description  = 'Submit some jobs'
                       ,allow_abbrev = True)

parser.add_argument('--time'
                   ,default = 5
                   ,type    = int
                   ,help    = 'Walltime for batch jobs')

parser.add_argument('--sigma'
                   ,default = 0.01
                   ,type    = float
                   ,help    = 'Walltime for batch jobs')

parser.add_argument('--econv'
                   ,default = 0.001
                   ,type    = float
                   ,help    = 'Walltime for batch jobs')

parser.add_argument('--dconv'
                   ,default = 0.001
                   ,type    = float
                   ,help    = 'Walltime for batch jobs')

parser.add_argument('--src'
                   ,default = ''
                   ,type    = str
                   ,help    = 'Path to bulk .traj files')

parser.add_argument('--target'
                   ,default = ''
                   ,type    = str
                   ,help    = 'Path to where jobs will be submitted from')

parser.add_argument('--elems'
                   ,default = ''
                   ,type    = lambda x: list(map(int,x.split()))
                   ,help    = 'Path to bulk .traj files')
