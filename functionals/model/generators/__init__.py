# External Modules
from typing import Any, Type, TYPE_CHECKING

# Internal Modules
if TYPE_CHECKING:
    from dbgen.support.model     import Model

from functionals.model.generators.io import io
from functionals.model.generators.load import load
from functionals.model.generators.analysis import analysis
#############################################################
def add_generators(mod:Type['Model'])->None:
    io(mod)
    load(mod)
    analysis(mod)
