# External Modules
from typing import Any,Type,TYPE_CHECKING
from os import environ
# Internal Modules
from dbgen.support.model     import Model,new_model
from dbgen.support.object    import DEFAULT,DESC
from dbgen.support.sqltypes  import Int,Varchar,Text,Decimal,Date
from dbgen.core.parsing                 import parser
from dbgen.support.misc                 import ConnectInfo

"""
Work in progress

Preliminary sketches for a DB to be used in association with functional
development

"""
################################################################################

# Add objects and relations
#-------------------------
def add_objects(mod:Type['Model'])->None:

    model = mod # type: Any                             ### WHY, MYPY, WHY ###

    ###############
    # Job related #
    ###############
    class Calc(model):
        beef = Int('tiny'), DESC('Whether or not functional is BEEF-style (i.e. '
                                 ' whether or not the coefs attr is NULL)')
        coefs = Text(), DESC('JSON dumped 8x8 matrix of exchange coefficients '
                             ' which defines a functional as a sum of bases')
        xc   = Varchar(), DESC("Name of functional")
        psp  = Varchar(), DESC("Path to PSPs used")
        kx   = Int(), DESC('Kpoints in x direction')
        ky   = Int(), DESC('Kpoints in y direction')
        kz   = Int(), DESC('Kpoints in z direction')
        spinpol = Int('tiny'), DESC('Whether calculation was spin polarized')
        econv   = Decimal(), DESC('Energy convergance criterion')

    class Job(object):
        logfile    = Varchar(),   DESC("Path to primary log file")
        log        = Text('long'),DESC('Content of primary log file')
        timestamp  = Date(),      DESC('Unix timestamp on log file')

        _init = logfile
        _components = Calc

    class Cell(model):
        """
        Periodic cells defined by three vectors
        """
        a0,a1,a2      = Decimal(),Decimal(),Decimal()
        b0,b1,b2      = Decimal(),Decimal(),Decimal()
        c0,c1,c2      = Decimal(),Decimal(),Decimal()
        a,b,c         = Decimal(),Decimal(),Decimal()
        surface_area  = Decimal()
        volume        = Decimal()

        _init         = a0,a1,a2,b0,b1,b2,c0,c1,c2

    class Species(model):
        """
        Abstraction of a struct which throws away position info.

        Like Pure_struct, but contains stoich information and special considerations
        for molecules which are defined by pointgroup rather than spacegroup

        Composition is a python dictionary mapping atomic number to NORMALIZED stoich
        """
        composition = Varchar(), DESC('Stringified python dict (ordered, normalized)')
        symmetry    = Varchar(), DESC('For molecules, pointgroup; '
                                      'for bulks, Prototype name;'
                                      'for surfaces, underlying prototype+facet')
        nickname    = Varchar(), DESC('Human-readable name, ought (but need not) be unique')

        _init = composition,symmetry

    class Element(model):
        """
        Chemical element
        """
        # Subset of Mendeleev data
        #-------------------------
        atomic_number           = Int()
        symbol                  = Varchar()
        atomic_weight           = Decimal()
        name                    = Varchar()
        atomic_radius           = Int(),     DESC("Angstrom")
        phase                   = Varchar()
        group_id                = Int(),     DESC("Column in periodic table")
        period                  = Int(),     DESC("Row in periodic table")
        evaporation_heat        = Decimal(), DESC('kJ/mol')
        melting_point           = Decimal(), DESC('K')
        is_radioactive          = Int(),     DESC('eV')
        lattice_struct          = Varchar()
        econf                   = Varchar(), DESC('electron configuration')
        heat_of_formation       = Decimal(), DESC('kJ/mol')
        electron_affinity       = Decimal(), DESC('eV')
        boiling_point           = Decimal(), DESC('K')
        proton_affinity         = Decimal(), DESC('kJ/mol')
        en_pauling              = Decimal()

        _init = atomic_number

    class Struct(model):
        """
        Chemical structure defined in periodic cell
        """
        raw              = Text(),      DESC('JSON encoding of ASE atoms object')
        system_type      = Varchar(),   DESC('One of: bulk, molecule, surface')

        _init            = raw
        _components      = Cell,Species

    class Atom(model):
        """
        An atom, considered within a specific chemical structure
        """
        ind          = Int(), DESC('ASE atom index')
        number       = Int(), DESC('Atomic number')
        x,y,z        = Decimal(),Decimal(),Decimal()
        constrained  = Int(), DESC('Whether or not there was a FixAtoms constraint')
        magmom       = Decimal(), DESC('Units: Bohr')
        tag          = Int(), DESC('ASE atom tag')
        layer        = Int(), DESC('Autodetected layer: NULL if not a surface atom')
        ads          = Int(), DESC('Which adsorbate, if any (0 or NULL = Not an adsorbate atom)')

        _init       = ind
        _parents    = Struct,
        _components = Element,

def add_generators(mod:Type['Model'])->None:
    pass

def make_model()->Type[Model]:
    m = new_model('johannes') # type: Type['Model']
    add_objects(m)
    add_generators(m)
    return m

def main()->None:
    """
    Run the model with no extensions from command line.
    """
    args = parser.parse_args()
    m    = new_model('catalysis') # type: Type['Model']
    add_objects(m)
    add_generators(m)

    # Assume connection info in Environment Variables
    #-----------------------------------------------
    db   = ConnectInfo.from_file(environ['DEV'])
    mdb  = ConnectInfo.from_file(environ['LOG'])

    only,xclude = [set(x.split()) for x in [args.only,args.xclude]]

    m._run(db,mdb,nuke=args.nuke,add=args.add,retry=args.retry,only=only,
          xclude=xclude,start=args.start,until=args.until)

if __name__=='__main__':
    main()
