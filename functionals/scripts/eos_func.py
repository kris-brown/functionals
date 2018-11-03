from typing   import Tuple
from string   import ascii_lowercase
from random   import choices
from os       import environ,remove
from os.path  import join
from base64   import b64encode
from ase.eos  import EquationOfState # type: ignore
from numpy   import polyfit,poly1d,mean,std,array # type: ignore

#####################################################
def eos_func(vols_:str,engs_:str,n_atoms:int)->Tuple[float,float,float,str]:
    '''
    Use fancy equation of state fitter to series of vol,eng pairs
    '''
    import matplotlib.pyplot as plt # type: ignore

    # Constants
    #----------
    randroot = environ['HOME']
    suffix   = 'tmp_'+''.join(choices(ascii_lowercase,k=8))+'.png'
    pth      = join(randroot,suffix)
    floats   = lambda x: array(list(map(float,x.split(','))))
    vol,eng  = map(floats,[vols_,engs_])

    # Remove any outliers
    p      = poly1d(polyfit(vol,eng,2))
    resid  = [abs(p(v)-e) for v,e in zip(vol,eng)]
    maxres = mean(resid)+std(resid)
    inds   = [i for i,r in enumerate(resid) if r < maxres]

    vols,engs = vol[inds],eng[inds]

    eos = EquationOfState(vols,engs) # type: ignore
    bestvol,besteng,bmod = eos.fit() # type: Tuple[float,float,float]
    eos.plot(pth)
    encoded = b64encode(open(pth, "rb").read()).decode("utf-8")
    remove(pth)
    plt.clf()
    return bestvol/n_atoms, besteng/n_atoms, bmod, encoded

if __name__=='__main__':
    e = '-242.913,-242.323,-242.568,-242.765,-241.563,-242.995,-240.974,-242.039,-243.003,-243.033,-241.717,-242.885,-242.039,-242.417,-242.695'
    v = '54.620,59.595,57.912,56.272,42.752,49.949,41.405,61.311,53.032,51.475,63.059,48.454,44.129,45.563,46.988'
    print(eos_func(v, e, 4))
