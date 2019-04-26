# External
from json import dump
from os   import getcwd
# Internal
from functionals.fit.fit import Fit
################################################################################

def main() -> None:
    fit = Fit.from_json(getcwd())
    x= fit.fit()
    with open('result.json','w') as fi:
        dump(x, fi)

if __name__=='__main__':
    main()
