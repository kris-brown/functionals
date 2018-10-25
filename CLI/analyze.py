# External modules
from typing  import Tuple,Any,Optional as O
from sys     import argv
from os      import listdir
from os.path import join,isdir,exists
from numpy   import polyfit,poly1d,mean,std,array # type: ignore
from ase.io  import read # type: ignore
from ase.eos import EquationOfState # type: ignore

# Internal Modules
from dbgen.support.misc import ConnectInfo
from dbgen.core.sql     import sqlselect

"""
A CLI interface for viewing the results of jobs
"""
################################################################################
default = ConnectInfo.from_file('/Users/ksb/Documents/JSON/functionals.json')

def decode_base64(data:bytes)->bytes:
    """Decode base64, padding being optional."""
    import base64
    missing_padding = len(data) % 4
    if missing_padding != 0:
        data += b'='* (4 - missing_padding)
    return base64.decodestring(data)

def viz_(img:bytes)->None:
    import base64
    import io
    from matplotlib import pyplot as plt # type: ignore
    import matplotlib.image as mpimg # type: ignore
    i = io.BytesIO(img)
    i2 = mpimg.imread(i, format='png')

    plt.imshow(i2, interpolation='nearest')
    plt.show()

def viz(material:str,cxn:ConnectInfo=default)->None:

    q = '''SELECT img
            FROM functionals.expt
                JOIN species ON expt__species = species_id
            WHERE species.nickname=%s
            LIMIT 1
            '''

    from os import environ
    conn = cxn.connect()
    image = sqlselect(conn,q,[material])[0][0]

    i1 = decode_base64(bytes(image,encoding='UTF-8')) # type: ignore

    viz_(i1)

class Bulk(object):
    """
    Result of a Bulk job, initialized with path to folder containing subfolders
    """
    def __init__(self, pth : str) -> None:
        self.pth=pth

        eng_,vol_ = [],[]
        for strain in listdir(self.pth):
            efile = join(self.pth,strain,'energy.txt')
            if exists(efile): # the calculation has finished
                with open(efile,'r') as f:
                    eng_.append(float(f.read()))
                atoms = read(join(self.pth,strain,'init.traj'))
                vol_.append(atoms.get_volume())

        eng,vol = array(eng_),array(vol_)

        # Determine if there are any outliers
        p      = poly1d(polyfit(vol,eng,2))
        resid  = [abs(p(v)-e) for v,e in zip(vol,eng)]
        maxres = mean(resid)+std(resid)
        inds   = [i for i,r in enumerate(resid) if r < maxres]

        self.eng = eng[inds]
        self.vol = vol[inds]

    def fit(self) -> Tuple[float,float,float]:
        eos = EquationOfState(self.vol,self.eng)
        return eos.fit()

    def plot(self, pth : O[str] = None) -> None:
        import matplotlib   # type: ignore
        matplotlib.use('Qt5Agg')
        eos = EquationOfState(self.vol,self.eng)
        print('Volume %f, Energy %f, Bulkmod %f'%self.fit())
        eos.plot(pth)


def main() -> None:
    import sys
    viz(sys.argv[1])

if __name__=='__main__':
    main()
