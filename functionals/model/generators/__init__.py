# Internal Modules
from dbgen import Model

from functionals.model.generators.io import io
from functionals.model.generators.fit import fit
from functionals.model.generators.analysis import analysis
from functionals.model.generators.bulk_analysis import bulk_analysis
#############################################################


def add_generators(mod: Model) -> None:
    io(mod)
    fit(mod)
    analysis(mod)
    bulk_analysis(mod)
