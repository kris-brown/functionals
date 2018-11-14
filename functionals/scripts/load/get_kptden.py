from typing import Tuple
import json
import ase # type: ignore


def get_kptden(kx : int
              ,ky : int
              ,kz : int
              ,ax : float
              ,ay : float
              ,az : float
              ,bx : float
              ,by : float
              ,bz : float
              ,cx : float
              ,cy : float
              ,cz : float
              ) -> Tuple[float,float,float]:
    """
    Returns kpt density:
        https://wiki.physics.udel.edu/phys824/About_k-point_sampling
    """

    atoms = ase.Atoms(numbers   = [1] # fake Atoms object with correct cell
                     ,positions = [[0,0,0]]
                     ,cell      = [[ax,ay,az],[bx,by,bz],[cx,cy,cz]])

    recipcell = atoms.get_reciprocal_cell()

    denoms = [2. * 3.1415 * (recipcell[i]**2).sum()**(0.5) for i in range(3)]

    kptdenx = round(kx/denoms[0],2)
    kptdeny = round(ky/denoms[1],2)
    kptdenz = round(kz/denoms[2],2)

    return kptdenx,kptdeny,kptdenz
