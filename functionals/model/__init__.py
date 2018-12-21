# Internal Modules
from dbgen                        import Model
from functionals.model.objects    import new_model
from functionals.model.generators import add_generators

def make_model() -> Model:
    m = new_model()
    add_generators(m)
    return m
