# External
from unittest import main
from sys      import exit
################################################################################
pre  = 'functionals.test.tests.'
mods = ['UnitSimple']
args = dict(failfast = True, exit = False, warnings = 'ignore')

if __name__ == '__main__':

    for m in mods:
        print('\nTesting ',m)
        res = main(module=pre+m, **args) # type: ignore
        if res.result.errors or res.result.failures:
            exit()
