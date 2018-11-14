# External
from typing import List as L, Tuple as T
# Internal
from dbgen.support.misc import ConnectInfo as Conn
from dbgen.core.sql     import sqlselect


BigTuple = T[L[str],L[str],L[str],L[str],L[str],L[float],L[float],L[int]]
################################################################################
def query_fitdata(cxnpth : str, where : str) -> BigTuple:
    db   = Conn.from_file(cxnpth)
    conn = db.connect()
    q    = '''SELECT species.composition,
                     symmetry,
                     beef,
                     calc.coefs,
                     xc,
                     pw,
                     econv,
                     expt.n_atoms
                FROM dft_data
                    JOIN expt ON dft_data__expt = expt_id
                    JOIN species ON expt__species = species_id
                    JOIN calc ON expt__calc = calc_id
                 WHERE '''+(where or '1')
    output = sqlselect(conn,q)

    if output:
        return tuple(map(list, zip(*output))) # type: ignore
    else:
        return ([],[],[],[],[],[],[],[])

if __name__=='__main__':
    print(query_fitdata('/Users/ksb/Documents/JSON/functionals.json',''))