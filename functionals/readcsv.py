import numpy as np # type: ignore
import csv
from ast import literal_eval

################################################################################
def readcsv(fi:str)->tuple:
    names = []
    target = []
    nonxc = [] #everything in total energy but xc (Hartree, electron-ion); not fitted
    x = []     #exchange
    cf = open(fi, 'r')
    reader = csv.DictReader(cf, delimiter=',', quotechar='"')
    for row in reader:
        names.append(row['Element'])
        target.append(float(row['keldAE']))
        nonxc.append(float(row['DFT atomE']))
        x.append(np.array(literal_eval(row['xc_contribs'])))
    cf.close()
    return names,np.array(target),np.array(nonxc),np.array(x)
