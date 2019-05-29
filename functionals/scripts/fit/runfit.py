# External
from json import dump
from os   import getcwd
# Internal
from functionals.fit.fit import Fit
################################################################################

def main() -> None:
    fit = Fit.from_json(getcwd())
    x   = [fit.fit(i).to_dict() for i in range(5)]
    with open('result.json','w') as fi:
        dump(x, fi)

if __name__=='__main__':
    main()
