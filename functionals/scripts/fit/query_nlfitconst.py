# External
from typing import List as L
# Internal
from dbgen.core.misc import ConnectInfo as Conn
from dbgen.core.sql     import sqlselect

################################################################################
def query_nlfitconst(cxnpth:str,where:str)->L[str]:
    db   = Conn.from_file(cxnpth)
    conn = db.connect()
    q    = 'SELECT nlconst_name FROM nonlin_const WHERE '+(where or '1')
    return list([x[0] for x in sqlselect(conn,q)])

if __name__ == '__main__':
    print(query_nlfitconst('/Users/ksb/Documents/JSON/functionals.json',''))
