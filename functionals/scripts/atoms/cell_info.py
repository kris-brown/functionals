from typing import Tuple as T
######################################################################################################
def cell_info(a0:float,a1:float,a2:float,
              b0:float,b1:float,b2:float,
              c0:float,c1:float,c2:float
              ) -> T[float,float,float,float,float]:
    """ Basic geometry (inputs are 3 <x,y,z> vectors) """
    surface_area = float((a0*b2-a2*b1) * (a1*b2 - a2*b1)
                     +(a0*b2-a2*b0) * (a0*b2 - a2*b0)
                     +(a0*b1-a1*b0) * (a0*b1 - a1*b0))**0.5
    volume = a0*(b1*c2-b2*c1)\
            -a1*(b0*c2-b2*c0)\
            +a2*(b0*c1-b1*c0)
    a = float(a0*a0+a1*a1+a2*a2)**0.5
    b = float(b0*b0+b1*b1+b2*b2)**0.5
    c = float(c0*c0+c1*c1+c2*c2)**0.5
    return surface_area,volume,a,b,c
