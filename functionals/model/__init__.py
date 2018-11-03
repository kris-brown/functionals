from typing import Any,Type

# Internal Modules
from dbgen.support.model          import Model,new_model
from functionals.model.objects    import add_objects
from functionals.model.generators import add_generators

def make_model()->Type[Model]:
    m = new_model('johannes') # type: Type['Model']
    add_objects(m)
    add_generators(m)
    return m
