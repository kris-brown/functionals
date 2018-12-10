# Internal Modules
from dbgen                          import Model
from functionals.model.objects    import add_objects
from functionals.model.generators import add_generators

def make_model()->Model:
    m = Model('johannes')
    add_objects(m)
    add_generators(m)
    return m
