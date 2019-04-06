from typing     import List as L,Tuple as T
from ase.data   import chemical_symbols # type: ignore
from re         import findall
from json       import load
################################################################################

def parse_keld(name:str,keldpth:str) -> T[float,float,float,float]:
    """
    Extracts information of materials into triples for struct_dataset_element
    If the species doesn't exist in the database, we need to add it.
    """

    # Get data from public JSON file
    #--------------------------------
    with open(keldpth,'r') as f: data = load(f)
    keld = data[name.split('_')[0]]
    return (keld.get('cohesive energy'),
            keld.get('bulk modulus'),
            keld.get('lattice parameter'),
            keld.get('magmom',0.))
