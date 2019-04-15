# External Modules
from typing import Any, Type

# Internal Modules
from dbgen import Model, Obj, Attr, Rel, Int, Varchar, Text, Decimal, Boolean, Date, PathEQ, Path

####################
# VASP INPUT FILES # --- curently no real need to parse these?
####################

###############################################################################
potcar = Obj(
    name  = 'potcar',
    desc  = 'POTCAR pseudopotential file',
    attrs = [Attr('titel', Varchar(),   desc='title? [sic]', id = True),
              Attr('lultra',Boolean(),  desc='use ultrasoft pp'),
              Attr('iunscr',Int(),      desc='unscreen: 0-lin, 1-nonlin, 2-no'),
              Attr('rpacor',Decimal(),  desc='partial core radius'),
              Attr('pomass',Decimal(),  desc='mass'),
              Attr('zval',  Decimal(),  desc='valence'),
              Attr('rcore', Decimal(),  desc='outmost cutoff radius'),
              Attr('rwigs', Decimal(),  desc='wigner-seitz radius (au A)'),
              Attr('enmax', Decimal()),
              Attr('enmin', Decimal()),
              Attr('lcor',  Boolean(),  desc='correct aug charges'),
              Attr('lpaw',  Boolean(),  desc='paw pp'),
              Attr('eaug',  Decimal()),
              Attr('rmax',  Decimal(),  desc='core radius for proj-oper'),
              Attr('raug',  Decimal(),  desc='factor for augmentation sphere'),
              Attr('rdep',  Decimal(),  desc='radius for radial grids'),
              Attr('rdept', Decimal(),  desc='core radius for aug-charge')])
###############################################################################

ams = ['a1'+x for x in '12345']+ ['msb']

idecs = ['ediff','encut','sigma','magmom'] + ams
iints = ['ismear','npar','nelm','ispin','ibrion']
ibool = ['lbeefens','addgrid','lasph','lwave']
ivars = ['metagga','gga','prec','algo']
idict = {Decimal() : idecs, Int() : iints, Boolean() : ibool, Varchar() : ivars}

incar = Obj(
    name = 'incar',
    desc  = 'VASP input file',
    attrs = [Attr(x,dt) for dt,xs in idict.items() for x in xs])
###############################################################################
