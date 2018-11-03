# External Modules
from typing import Any, Type, TYPE_CHECKING

# Internal Modules
if TYPE_CHECKING:
    from dbgen.support.model import Model

from dbgen.support.sqltypes  import Int, Varchar, Text, Decimal, Date
from dbgen.support.object    import DEFAULT, DESC

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
        fusion_heat             = Decimal(), DESC('kJ/mol')
        melting_point           = Decimal(), DESC('K')
        is_radioactive          = Int(),     DESC('eV')
        lattice_struct          = Varchar(), DESC("e.g. HEX, SC")
        econf                   = Varchar(), DESC('electron configuration')
        heat_of_formation       = Decimal(), DESC('kJ/mol')
        electron_affinity       = Decimal(), DESC('eV')
        boiling_point           = Decimal(), DESC('K')
        proton_affinity         = Decimal(), DESC('kJ/mol')
        en_pauling              = Decimal(),
        pointgroup              = Varchar()
        spacegroup              = Int()
        metallic_radius         = Decimal()
        vdw_radius              = Decimal()
        density                 = Decimal()
        en_allen                = Decimal()
        en_ghosh                = Decimal()
        covalent_radius_bragg   = Decimal()
        covalent_radius_slater  = Decimal()
        geochemical_class       = Varchar()
        abundance_crust         = Decimal()
        abundance_sea           = Decimal()
        atomic_volume           = Decimal()
        lattice_constant        = Decimal()
        dipole_polarizability   = Decimal()
        thermal_conductivity    = Decimal()
        gas_basicity            = Decimal()

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
        checksum = Varchar(), DESC("MD5 hash of fiwle")
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

        _init = beef, coefs, xc, pw, econv # kx,ky,kz need to change w/ cell size, atomic calculations are spinpol while the bulk jobs we want to compare to are not

    class Cell(model):
        """
        Periodic cells defined by three vectors
        """
        a0,a1,a2     = Decimal(), Decimal(), Decimal()
        b0,b1,b2     = Decimal(), Decimal(), Decimal()
        c0,c1,c2     = Decimal(), Decimal(), Decimal()
        a,b,c        = Decimal(), Decimal(), Decimal()
        surface_area = Decimal()
        volume       = Decimal()

        _init = a0, a1, a2, b0, b1, b2, c0, c1, c2

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
        n_atoms     = Int(),     DESC("Total # of atoms in normalized stoich")

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

        _init       = raw
        _components = Cell, Species, Pure_struct


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
        contribs   = Text(),      DESC('output exchange contribs, if beef calculation')
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
        n_atoms    = Int(),        DESC("Value for all jobs in this experiment")
        energy_pa  = Decimal(),    DESC("Per atom, eV")
        bulkmod    = Decimal(),    DESC("Bulk modulus, GPa")
        volume_pa  = Decimal(),    DESC("Per atom, A^3")
        img        = Text('long'), DESC('base64 encoded image of EOS fit')
        eform      = Decimal(),    DESC("Per atom, eV")
        lattice    = Decimal(),    DESC("Conventional unit cell lattice (optimized)")
        n_data     = Int(),        DESC('number of aggregated data points')
        min_gap    = Decimal(),    DESC('Absolute difference between best singlepoint '
                                        'volume_pa and the fitted optimum')
        _init    = n_atoms # want to distinguish fits with different scales b/c DFT isn't perfect
        _parents = Species, Calc

    class Bulk_job(model):
        """
        A subset of jobs which have a many-one relationship linking jobs to an experiment
        """
        gap = Decimal(), DESC("Absolute difference in volume from 'minimum'")

        _init       = ()
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

    class Cohesive_data(model):
        """
        Data that are relevant to fitting BEEF coefs using cohesive energy
        """
        name            = Varchar(),    DESC("Species nickname")
        coefs           = Text('long'), DESC('Calc Coefs')
        composition     = Varchar(),    DESC("Species composition")
        atomic_contribs = Text('long'), DESC("Serialized dict of relevant reference info")
        atomic_energies = Text(),       DESC("Serialized dict of relevant reference info")
        bulk_contribs   = Text('long'), DESC("Best job xc contribs")
        bulk_energy     = Decimal(),    DESC("Best job energy")
        target          = Decimal(),    DESC("Experimental cohesive energy")
        bulk_ratio      = Int(),        DESC("Ratio of bulk system to size of normalized species")

        _init       = ()
        _parents    = Expt
        _components = Job # best job

    class Fit(model):
        """
        A fit to some subset of cohesive data + lattice data
        """
        name      = Varchar()
        timestamp = Date(),     DESC("Timestamp of fitting")
        runtime   = Decimal(),  DESC("Duration of fitting, s")
        result    = Text()
        basis     = Int(),       DESC('Size of fitted functional')
        norm      = Decimal(),   DESC('Regularization term')
        initfit   = Int('tiny'), DESC('If true: initialize with lstsq fit w/o '
                                        'constraints (else with 0)')
        bound     = Decimal(),  DESC('Range over which to search for coefficients')
        maxiter   = Int(),      DESC('Stop nonlinear fitting after this step')
        constden  = Int(),      DESC("Number of s or alpha points linearly  "
                                     "constrained between logscale(-2,2) ")

        resid      = Decimal(), DESC('Final value of cost function')

        _init = name

    class Fit_step(model):
        '''
        A single iteration in a fitting process
        '''
        niter  = Int(),     DESC("Iteration number")
        cost   = Decimal(), DESC("Objective function cost")
        c_viol = Decimal(), DESC("Constraint cost")

        _init    = niter
        _parents = Fit

    class Fit_data(model):
        '''
        Mapping table specifying which Cohesive_data was used in a given fit
        '''
        _parents = Fit, Cohesive_data

    class Const(model):
        '''
        Put a linear constraint (LT/GT/EQ) to Fx(s,a)
        '''
        const_name  = Varchar(), DESC('Name of constraint')
        description = Text(),    DESC("Description of constraint ")
        val         = Decimal(), DESC('Value of Fx(s,alpha)')
        kind        = Varchar(), DESC("GT/LT/EQ")
        s           = Decimal(), DESC("If valid only for particular s, else None")
        alpha       = Decimal(), DESC("If valid only for particular alpha, else None")

        _init = const_name

    class Fit_const(model):
        '''
        Mapping table to denote which constraints were used in a given fit
        '''
        _parents = Fit, Const
