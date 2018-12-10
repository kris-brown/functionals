from typing  import Tuple as T,List as L, Any
from gzip    import open as gopen
from hashlib import md5
from os      import listdir
from os.path import isdir,join
from xml.etree.ElementTree import fromstring
from ase.data import chemical_symbols # type: ignore

################################################################################
def parse_setup(root:str)->T[L[str],L[str],L[str],
                             L[str],L[int],L[int]]:
    '''Applies parse to a directory'''

    def get(tree:Any,key:str,attr:str)->str:
        e = tree.find(key)
        assert e is not None
        return e.get(attr).strip()

    def parse_gz(tree:Any)->T[str,str,str,int,int]:
        """
        Reads a gzipped XML GPAW setup file
        """
        xc   = get(tree,'xc_functional','type')
        name = get(tree,'generator','name')
        z    = int(get(tree,'atom','Z'))
        val  = int(get(tree,'atom','valence'))
        kind = 'GPAW setup'
        return xc,kind,name,z,val

    def parse_upf(tree:Any)->T[str,str,str,int,int]:

        xc    = get(tree,'PP_HEADER','functional')
        z     = chemical_symbols.index(get(tree,'PP_HEADER','element'))
        val   = int(float(get(tree,'PP_HEADER','z_valence')))
        gen   = get(tree,'PP_HEADER','generated')
        sg15  = 'ONCVPSP' in gen
        date  = get(tree,'PP_HEADER','date')
        assert sg15, "New class of UPF pseudopotentials?"+gen
        name = 'sg15_'+date
        kind = 'sg15'
        return xc,kind,name,z,val

    def parse_general(pth:str)->T[str,str,str,str,int,int]:
        full = join(root,pth)
        gz   = pth[-3:]=='.gz'
        upf  = pth[-4:]=='.upf'
        assert gz != upf, pth

        # Get specific file opener and parser
        parser = parse_gz if gz else parse_upf
        opener = gopen if gz else open

        # Common parsing actions
        with opener(full,'rb') as f: # type: ignore
            content = f.read()

        tree = fromstring(content)
        hash = md5(bytes(content)).hexdigest()

        # return result of specific parser
        return (hash,*parser(tree))

    filt   = [x for x in listdir(root) if x[-3:] in ['upf','.gz']]
    parsed = [parse_general(x) for x in filt if 'basis' not in x]

    hash,xc,kind,name,z,val = list(map(list,zip(*parsed)))

    return hash,xc,kind,name,z,val # type: ignore

if __name__=='__main__':
    from sys     import argv
    root = argv[1]
    assert isdir(root), 'parse_setup operates on a DIRECTORY'
    print(parse_setup(root))
