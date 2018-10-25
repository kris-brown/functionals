# External Modules
from typing import Any, Type, TYPE_CHECKING

from dbgen.support.sqltypes  import Int, Varchar, Text, Decimal, Date

from dbgen.support.object    import DEFAULT, DESC

if TYPE_CHECKING:
    from dbgen.support.model     import Model
################################################################################
# Add objects and relations
#-------------------------
def add_objects(mod:Type['Model'])->None:

    model = mod # type: Any                             ### WHY, MYPY, WHY ###

    ###############
    # Job related #
    ###############

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

    class Setup_family(model):
        '''Class of Setups that are the same 'type' '''
        name = Varchar()
        kind = Varchar(),
        xc   = Varchar(), DESC("XC Functional used to generate")

    class Setup(model):
        """
        Pseudopotential
        """
        checksum = Varchar(), DESC("MD5 hash of file")
        val = Int(), DESC("number of valence electrons")

        _init       = checksum
        _components = Element, Setup_family

    class Calc(model):
        beef    = Int('tiny'), DESC('Whether or not functional is BEEF-style '
                                    '(i.e. whether or not the coefs attr is meaningful)')
        coefs   = Text(),      DESC('JSON dumped 8x8 matrix of exchange coefficients '
                                    ' which defines a functional as a sum of bases')
        xc      = Varchar(),   DESC("Name of functional")
        pw      = Int(),       DESC("Planewave cutoff, eV")
        econv   = Decimal(),   DESC('Energy convergance criterion')

        _init = beef, coefs, xc, econv # kx,ky,kz need to change w/ cell size, atomic calculations are spinpol while the bulk jobs we want to compare to are not

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

        _init         = a0, a1, a2, b0, b1, b2, c0, c1, c2

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
        n_elems     = Int(),     DESC("# of distinct chemical elements in species")
        _init = composition, symmetry

    class Species_dataset(model):
        """

        """
        dataset_name = Varchar()

    class Species_dataset_element(model):
        """

        """
        property = Varchar()
        value    = Text()
        datatype = Varchar()

        _init    = property
        _parents = Species, Species_dataset

    class Species_comp(model):
        """
        Mapping table between species and element to show composition
        """
        num = Int(), DESC('Number of this element in lowest integer terms')
        _parents = Species, Element

    class Pure_struct(model):
        """
        Structure abstraction based on AFLOW prototypes refined by Ankit Jain
        """
        prototype  = Varchar()
        spacegroup = Int()
        free       = Int(),    DESC('Number of free parameters')
        nickname   = Varchar()

        _init = prototype

    class Struct(model):
        """
        Chemical structure defined in periodic cell
        """
        raw              = Text(),      DESC('JSON encoding of ASE atoms object')
        system_type      = Varchar(),   DESC('One of: bulk, molecule, surface')
        sg               = Int(),       DESC("spacegroup")
        composition_norm = Text()
        n_atoms          = Int()
        n_elems          = Int()
        composition      = Text()
        metal_comp       = Text()
        str_symbols      = Text()
        str_constraints  = Text()

        _init            = raw
        _components      = Cell, Species, Pure_struct


    class Job(model):
        logfile    = Varchar(),    DESC("Path to primary log file")
        stordir    = Varchar()
        log        = Text('long'), DESC('Content of primary log file')
        timestamp  = Date(),       DESC('Unix timestamp on log file')
        user       = Varchar()
        kx         = Int(), DESC('Kpoints in x direction')
        ky         = Int(), DESC('Kpoints in y direction')
        kz         = Int(), DESC('Kpoints in z direction')
        kptden_x   = Decimal()
        kptden_y   = Decimal()
        kptden_z   = Decimal()
        spinpol    = Int('tiny'), DESC('Whether calculation was spin polarized')

        energy     = Decimal()

        _init       = logfile
        _components = Calc, Struct

    class Job_setup(model):
        """mapping table between Setup and Job"""
        _parents = Job, Setup

    class Atom(model):
        """
        An atom, considered within a specific chemical structure
        """
        ind          = Int(), DESC('ASE atom index')
        number       = Int(), DESC('Atomic number')
        x,y,z        = Decimal(), Decimal(), Decimal()
        constrained  = Int(),     DESC('Whether or not there was a FixAtoms constraint')
        magmom       = Decimal(), DESC('Units: Bohr')
        tag          = Int(),     DESC('ASE atom tag')

        _init       = ind
        _parents    = Struct,
        _components = Element,


    class Expt(model):
        """
        Set of single points on a particular material (with a particular calc)

        This should be a bulk calculation
        """
        energy  = Decimal(), DESC("Per atom, eV")
        bulkmod = Decimal()
        volume  = Decimal(), DESC("Per atom, A^3")
        img     = Text('long')
        eform   = Decimal(), DESC("Per atom, eV")

        _init    = ()
        _parents = Species, Calc

    class Bulk_job(model):
        """
        A subset of jobs which have a many-one relationship linking jobs to an experiment
        """

        _components = Expt
        _parents    = Job

    class Reference(model):
        """
        A single calculation that gives the energy of an isolated atom (w/ a calc)
        """
        energy   = Decimal()

        _init       = ()
        _parents    = Job
        _components = Calc, Element
