from typing    import Tuple as T, Optional as O
from string    import ascii_lowercase
from random    import choices
from os        import environ,remove
from os.path   import join
from base64    import b64encode
from ase.eos   import EquationOfState   # type: ignore
from ase.units import kJ                # type: ignore
from numpy     import polyfit,poly1d,mean,std,array # type: ignore

#####################################################
def eos_func(vols_:str,engs_:str,n_atoms:int) -> T[O[float],O[float],O[float],O[str],bool]:
    '''
    Use fancy equation of state fitter to series of vol,eng pairs
    '''
    import warnings; warnings.filterwarnings("ignore")
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

    try:
        bestvol,besteng,bmod = eos.fit() # type: T[float,float,float]
    except:
        return None,None,None,None,False

    eos.plot(pth)
    encoded = b64encode(open(pth, "rb").read()).decode("utf-8")
    remove(pth)
    plt.clf()

    return bestvol/n_atoms, besteng/n_atoms, bmod/ kJ * 1.0e24, encoded, True

if __name__=='__main__':
    #e = '-242.913,-242.323,-242.568,-242.765,-241.563,-242.995,-240.974,-242.039,-243.003,-243.033,-241.717,-242.885,-242.039,-242.417,-242.695'
    #v = '54.620,59.595,57.912,56.272,42.752,49.949,41.405,61.311,53.032,51.475,63.059,48.454,44.129,45.563,46.988'
    #print(eos_func(v, e, 4))

    v = '42.50854900000001,43.800328125,42.47201916899999,41.242421375,42.545099771000004,32.036235776000005,42.36255525600001,38.78609162499999,42.399022303,33.107082930999994,42.289683904,34.201530936,42.289683904,45.118016000000004,42.545099771000004,36.462258496,42.47201916899999,40.001687999999994,42.618264157,37.595375000000004,42.50854900000001,30.988732221,42.50854900000001,46.461869875000005,42.50854900000001,47.832147,42.47201916899999,52.062250904,42.545099771000004,49.229104625,42.581671488,53.540005609000005,42.581671488,55.04546246399999,42.545099771000004,56.578878718999995,42.545099771000004'
    e = '-3.601164,-3.598879,-3.601192,-3.599922,-3.601143,-3.446590,-3.600847,-3.587052,-3.601095,-3.482425,-3.600894,-3.513082,-3.600938,-3.593725,-3.601192,-3.559190,-3.600998,-3.595246,-3.601080,-3.574979,-3.601017,-3.404968,-3.600774,-3.586008,-3.601212,-3.575962,-3.601226,-3.531973,-3.601243,-3.563194,-3.601218,-3.513245,-3.600965,-3.492687,-3.600981,-3.470342,-3.600994'
    print(eos_func(v, e, 1)[:-1])
