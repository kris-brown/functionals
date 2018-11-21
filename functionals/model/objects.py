# External Modules
from typing import Any, Type

# Internal Modules
from dbgen import Model, Int, Varchar, Text, Decimal, Date, DEFAULT

################################################################################
# Add objects and relations
#-------------------------
def add_objects(mod:Type[Model])->None:

    model : Any = mod # WHY, MYPY, WHY ###

    ###############
    # Job related #
    ###############
    class Globals(model):
        '''Properties of entire database. There will be just one row'''
        all_data            = Text('long') # All cohesive data
        all_constraints     = Text('long') # All linear constraints
        all_nlconstraints   = Text('long') # All nonlinear constraints

        _init = ()

    class Element(model):
        """
        Chemical element
        """
        # Subset of Mendeleev data
        #-------------------------
        atomic_number           = Int()      # Number of protons
        symbol                  = Varchar()  # E.g. He, K, Li
        atomic_weight           = Decimal()  # Atomic units
        name                    = Varchar()
        atomic_radius           = Int(),     # Angstrom
        phase                   = Varchar()
        group_id                = Int(),     # Column in periodic table
        period                  = Int(),     # Row in periodic table
        evaporation_heat        = Decimal(), # kJ/mol
        fusion_heat             = Decimal(), # kJ/mol
        melting_point           = Decimal(), # K
        is_radioactive          = Int(),     # eV
        lattice_struct          = Varchar(), # e.g. HEX, SC
        econf                   = Varchar(), # electron configuration
        heat_of_formation       = Decimal(), # kJ/mol
        electron_affinity       = Decimal(), # eV
        boiling_point           = Decimal(), # K
        proton_affinity         = Decimal(), # kJ/mol
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
        kind = Varchar()
        xc   = Varchar() #XC Functional used to generate

    class Setup(model):
        """
        Pseudopotential
        """
        checksum = Varchar() # MD5 hash of file
        val      = Int()     # number of valence electrons

        _init       = checksum
        _components = Element, Setup_family

    class Calc(model):
        ''' Calculator details '''
        beef    = Int('tiny') # Whether or not functional is BEEF-style
                              #  (i.e. whether or not the coefs attr is meaningful)

        coefs   = Text()      # JSON dumped 8x8 matrix of exchange coefficients
                              #  which defines a functional as a sum of bases

        xc      = Varchar()   # Name of functional
        pw      = Int()       # Planewave cutoff, eV
        econv   = Decimal()   # Energy convergance criterion

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
        composition = Varchar() # Stringified python dict (ordered, normalized)
        symmetry    = Varchar() # For molecules, pointgroup;
                                 # for bulks, Prototype name;
                                 # for surfaces, underlying prototype+facet
        nickname    = Varchar() #Human-readable name, ought (but need not) be unique
        n_elems     = Int() ## of distinct chemical elements in species
        n_atoms     = Int() #Total # of atoms in normalized stoich

        _init = composition, symmetry

    class Species_dataset(model):
        """
        Datasets that contain information about chemical species
        """
        dataset_name = Varchar()

    class Species_dataset_element(model):
        """
        Datum about a particular species
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
        num = Int() # Number of this element in lowest integer terms

        _parents = Species, Element

    class Pure_struct(model):
        """
        Structure abstraction based on AFLOW prototypes refined by Ankit Jain
        """
        prototype  = Varchar()
        spacegroup = Int()
        free       = Int() #Number of free parameters
        nickname   = Varchar()

        _init = prototype

    class Struct(model):
        """
        Chemical structure defined in periodic cell
        """
        raw              = Text()    # JSON encoding of ASE atoms object
        system_type      = Varchar() # One of: bulk, molecule, surface
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
        logfile    = Varchar()      # Path to primary log file
        stordir    = Varchar()
        log        = Text('long')   # Content of primary log file
        timestamp  = Date()         # Unix timestamp on log file
        user       = Varchar()
        kx,ky,kz   = Int(),Int(),Int() # K-points
        kptden_x   = Decimal()
        kptden_y   = Decimal()
        kptden_z   = Decimal()
        spinpol    = Int('tiny')   # Whether calculation was spin polarized
        has_contribs = Int('tiny') # Whether xc_contribs.txt exists
        contribs   = Text()        # output exchange contribs, if beef calculation
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
        ind          = Int() # ASE atom index
        number       = Int() # Atomic number
        x,y,z        = Decimal(), Decimal(), Decimal() # position
        constrained  = Int()     # Whether or not there was a FixAtoms constraint
        magmom       = Decimal() # Units: Bohr
        tag          = Int()     # ASE atom tag

        _init       = ind
        _parents    = Struct,
        _components = Element,

    class Expt(model):
        """
        Set of single points on a particular material (with a particular calc)

        This should be a bulk calculation
        """
        n_atoms    = Int()        # Value for all jobs in this experiment
        energy_pa  = Decimal()    # Per atom, eV
        bulkmod    = Decimal()    # Bulk modulus, GPa
        volume_pa  = Decimal()    # Per atom, A^3
        all_vols   = Text('long') # every volume of the related bulk jobs
        volumes    = Text('long') # volumes of the 5 most optimal jobs
        energies   = Text('long') # energies of the 5 most optimal jobs
        contribs   = Text('long') # exchange contribs of the 5 most optimal jobs
        img        = Text('long') # base64 encoded image of EOS fit
        eform      = Decimal()    # Per atom, eV
        lattice    = Decimal()    # Conventional unit cell lattice (optimized)
        n_data     = Int()        # Number of aggregated data points
        min_gap    = Decimal()    # Absolute difference between best singlepoint
                                   # volume_pa and the fitted optimum

        _init    = n_atoms       # distinguish fits with different scales b/c DFT isn't perfect
        _parents = Species, Calc

    class Bulk_job(model):
        """
        A subset of jobs which have a many-one relationship linking jobs to an experiment
        """
        dv  = Decimal() # Difference in volume from 'minimum'
        gap = Decimal() # Abs(dv)
        near_min = Int('tiny') # Whether or not job is in the bottom 5 data points

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

    class Dft_data(model):
        """
        Data that are relevant to fitting BEEF coefs using cohesive energy
        """
        name                 = Varchar()     # Species nickname
        coefs                = Text('long')  # Calc Coefs
        composition          = Varchar()     # Species composition
        atomic_contribs      = Text('long')  # Serialized dict of relevant reference info
        atomic_energies      = Text()        # Serialized dict of relevant reference info
        bulk_contribs        = Text('long')  # Best job xc contribs
        bulk_energy          = Decimal()     # Best job energy
        expt_cohesive_energy = Decimal()     # Experimental cohesive energy
        expt_bm              = Decimal()     # Experimental bulk modulus
        expt_volume          = Decimal()     # Experimental volume of reference stoichiometry
        energy_vector        = Text('long')  # JSON'd vector of 5 energies
        volume_vector        = Text('long')  # JSON'd vector of 5 volumes
        contrib_vector       = Text('long')  # JSON'd 5x64 matrix with exchange contributions
        bulk_ratio           = Int()         # Ratio of bulk system to size of normalized species

        _init       = ()
        _parents    = Expt
        _components = Job # best job

    class Fit(model):
        """
        A fit to some subset of cohesive data + lattice data
        """
        # Input params
        name         = Varchar()
        bm_weight    = Decimal()   # Relative weight of bulk-modulus data to cohesive energy
        lat_weight   = Decimal()   # Relative weight of lattice data to cohesive energy
        constconst   = Text()      # SQL const on what linear constraints should be included
        nlconstconst = Text()      # SQL const on what linear constraints should be included
        dataconst    = Text()      # SQL const on what data should be included
        basis        = Int()       # Size of fitted functional
        initfit      = Int('tiny') # If true: initialize with lstsq fit w/o
                                    # constraints (else with 0)
        bound     = Decimal()   # Range over which to search for coefficients
        maxiter   = Int()       # Stop nonlinear fitting after this step
        constden  = Int()       # Number of s or alpha points linearly
                                # constrained between logscale(-2,2)
        # Intermediate computation
        raw_data    = Text('long') # Data to be fed to fitting script
        raw_const   = Text('long') # Data to be fed to fitting script
        raw_nlconst = Text('long') # Data to be fed to fitting script
        n_const     = Int(),DEFAULT(0)  # Number of constraints
        n_data      = Int(),DEFAULT(0)  # Number of data points
        n_nlconst   = Int(),DEFAULT(0)  # Number of nonlinear constraints

        # Result params
        timestamp = Date()       # Timestamp of fitting
        runtime   = Decimal()    # Duration of fitting, s
        r2_ce     = Decimal()    # R2 fit of cohesive energies
        r2_bm     = Decimal()    # R2 fit of bulk moduli
        r2_lat    = Decimal()    # R2 fit of lattice constants
        c_viol    = Decimal()
        score     = Decimal()    # Arbitrary combination of R2's and c_viol
        result    = Text()       # Flattened NxN fitted coefficients
        log       = Text('long') # Output from fitting
        err       = Text()       # Error during scipy.minimize()
        beefdist  = Decimal()    # Score for cartesian distance btw output and BEEF

        _init = name

    class Fit_step(model):
        '''
        A single iteration in a fitting process
        '''
        niter  = Int()     # Iteration number
        cost   = Decimal() # Objective function cost
        c_viol = Decimal() # Constraint cost

        _init    = niter
        _parents = Fit

    class Fit_data(model):
        '''
        Mapping table specifying which Cohesive_data was used in a given fit
        '''
        _parents = Fit, Dft_data


    class Const(model):
        '''
        Put a linear constraint (LT/GT/EQ) to Fx(s,a)
        '''
        const_name  = Varchar() # Name of constraint
        description = Text()    # Description of constraint
        val         = Decimal() # Value of Fx(s,alpha)
        kind        = Varchar() # GT/LT/EQ
        s           = Decimal() # If valid only for particular s, else None
        alpha       = Decimal() # If valid only for particular alpha, else None
        vec         = Text()    # Raw 64-float JSON, ignore s/a if defined

        _init = const_name

    class Fit_const(model):
        '''
        Mapping table to denote which constraints were used in a given fit
        '''
        _parents = Fit, Const

    class Nonlin_const(model):
        '''
        Nonlinear constraints
        '''
        nlconst_name = Varchar() # Name of nonlinear constraint
        description  = Text()    # Description of nonlinear constraint
        f            = Text()    # Source code for nonlin function to be minimized
                                  # expecting an 'x' input vector
        df           = Text()    # Source code for derivative of f
                                  # expecting an 'x' input vector
        hess         = Text()    # Source code for hessian of  f
                                  # expecting an 'x' input vector
        lb           = Decimal() # Lower bound
        ub           = Decimal() # Upper bound

        _init = nlconst_name

    class Fit_nonlin_const(model):
        '''
        Mapping table to denote which constraints were used in a given fit
        '''
        _parents = Fit, Nonlin_const
