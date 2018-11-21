from typing     import List,Tuple
from json       import load

def parse_mendeleev(i:int,
                    pth:str)->Tuple[List[str],List[float],List[str]
                                 ,List[int],List[str],  List[str]]:
    """
    Extracts information of elements
    """
    # Get data from public JSON file
    #-------------------------------
    with open(pth,'r') as f:
        data = load(f)[i - 1] # access the i'th element of the list
    ############################
    # Constants
    #----------
    #cols = ['symbol','atomic_weight','lattice_structure','econf'] # etc

    allcols = [ 'symbol', 'name', 'atomic_weight','atomic_radius'
              , 'phase','evaporation_heat', 'pointgroup','spacegroup'
              , 'melting_point', 'metallic_radius', 'vdw_radius'
              , 'density', 'en_allen' , 'is_radioactive'
              , 'lattice_structure' , 'fusion_heat'
              , 'econf', 'period', 'covalent_radius_bragg'
              , 'geochemical_class', 'abundance_crust', 'heat_of_formation'
              , 'electron_affinity', 'atomic_volume',  'boiling_point'
              , 'proton_affinity', 'covalent_radius_slater'
              , 'lattice_constant', 'dipole_polarizability'
              , 'en_ghosh', 'thermal_conductivity', 'group_id', 'en_pauling'
              , 'gas_basicity'
              ,'abundance_sea']

    # Initialize Variables
    #----------------------
    output = []
    for k in allcols:
        v = data.get(k)
        if isinstance(v,float):
            output.append(round(v,3))
        else:
            output.append(v)

    return tuple(output) # type: ignore
