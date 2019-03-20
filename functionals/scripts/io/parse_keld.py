from typing     import List as L,Tuple as T
from ase.data   import chemical_symbols # type: ignore
from re         import findall
from json       import load
################################################################################

def parse_keld(pth:str) -> T[str,L[str],L[str],L[int],int,L[str],L[str],str]:
    """
    Extracts information of materials into triples for struct_dataset_element
    If the species doesn't exist in the database, we need to add it.
    """
    # Constants
    #----------
    struct_dict = {'fcc'    :'A_1_a_225'
                  ,'bcc'    :'A_1_a_229'
                  ,'hcp'    :'A_2_c_194'
                  ,'diamond':'A_2_a_227'
                  ,'b1'     :'AB_1_a_b_225'  # rocksalt
                  ,'b2'     :'AB_1_a_b_221'  # cesium chloride
                  ,'b3'     :'AB_1_a_c_216'} # zincblende

    # Initialize Variables
    #----------------------
    output = [] # type: L[T[str,str,int,int,str,str]]

    # Get data from public JSON file
    #--------------------------------
    with open(pth,'r') as f:
        keld = load(f)

    cols      = ['lattice parameter',
                 'cohesive energy',
                 'bulk modulus',
                 'magmom',]

    # Extract data
    #------------
    for _,v in keld.items():
        dataset     = [] # type: L[T[str,str]]
        composition = [] # type: L[T[str,str,int,int]]

        k = v['name']

        bm = 'bulk modulus'
        if bm not in v and bm + ' kittel' in v:
            v[bm] = v[bm + ' kittel']

        elemstr,crystal = k.split('-')                  # 'LiH','b1'
        sym     = struct_dict[crystal]             # 'b1' -> 'AB_1_a_b_225'
        elems_ = findall('[A-Z][^A-Z]*', elemstr)        # ['Li','H']
        elems  = list(sorted([chemical_symbols.index(x) for x in elems_]))

        compstr     = str({e:1 for e in elems})
        composition = [(compstr,sym,num,len(elems)) for num in elems]
        dataset     = [(c,v.get(c)) for c in cols]

        # need to combine the composition list of tuples and the dataset list of tuples
        output.extend([(*c,*dataset[0]) for c in composition]
                     +[(*composition[0],*d) for d in dataset])

    # Organize data for output
    #-------------------------
    compositions,symmetries,atomic_numbers,n_elems,properties,values  = map(list,zip(*output))

    return ('keld',compositions,symmetries,atomic_numbers,1,n_elems,properties,values,'float') # type: ignore

if __name__=='__main__':
    from os import getenv
    print(parse_keld(getenv('KELDPATH','')))
