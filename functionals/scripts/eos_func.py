from typing import Tuple
from string import ascii_lowercase
from random import choices
from os     import environ,remove
from os.path import join
from base64  import b64encode
from ase.eos  import EquationOfState # type: ignore
import matplotlib.pyplot as plt # type: ignore
#####################################################
def eos_func(vols:str,engs:str,n_atoms:int)->Tuple[float,float,float,str]:
    # Constants
    #----------
    randroot = environ['HOME']
    suffix   = 'tmp_'+''.join(choices(ascii_lowercase,k=8))+'.png'
    pth      = join(randroot,suffix)
    floats   = lambda x: list(map(float,x.split(',')))
    eos      = EquationOfState(floats(vols),floats(engs)) # type: ignore
    vol,eng,bmod = eos.fit() # type: Tuple[float,float,float]
    eos.plot(pth)
    encoded = b64encode(open(pth, "rb").read()).decode("utf-8")
    remove(pth)
    plt.clf()
    return vol/n_atoms,eng/n_atoms,bmod,encoded
