# Internal Modules
from dbgen import Model

from functionals.model.generators.io import io
from functionals.model.generators.fit import fit
from functionals.model.generators.load import load
from functionals.model.generators.analysis import analysis
#############################################################
def add_generators(mod : Model) -> None:
    io(mod)
    fit(mod)
    load(mod)
    analysis(mod)
