# External
from json import dump
from os   import getcwd
# Internal
from functionals.fit.fit import Fit
################################################################################

def main() -> None:
    fit = Fit.from_json(getcwd())
    fit.fit(getcwd())
if __name__=='__main__':
    main()
