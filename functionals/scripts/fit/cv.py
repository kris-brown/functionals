from typing import List
import numpy as np
import functionals.fit.functional as fx
from functionals.scripts.fit.count_bumps import count_bumps, all_bumps
from functionals.scripts.fit.fitfun import fitfun
from functionals.CLI.analyze import default, mk_dataset
import dbgen


def cv() -> None:

    n_split = 5

    # Assemble data
    ms2_id = dbgen.sqlselect(default.connect(),
                             "SELECT calc_id FROM calc WHERE name='ms2'")

    data = mk_dataset(ms2_id[0][0])
    datas = data.split(n_split)

    # Initialize functionals
    x_init = fx.FromMatrix.frompath('ms2').x

    # Digression, what is mCAML data error based on ms2 wavefuns
    ms2err_ce = data.mae(x_init, 'ce', rel=False)
    ms2err_bm = data.mae(x_init, 'bm')
    ms2err_lat = data.mae(x_init, 'lc')
    mcaml_x = fx.FromMatrix.frompath('msurf').x
    mcamlerr_ce = data.mae(mcaml_x, 'ce', rel=False)
    mcamlerr_bm = data.mae(mcaml_x, 'bm')
    mcamlerr_lat = data.mae(mcaml_x, 'lc')
    print(ms2err_ce, ms2err_bm, ms2err_lat)
    print(mcamlerr_ce, mcamlerr_bm, mcamlerr_lat)
    data.resid('ce', mcaml_x, False)
    breakpoint()
    # END digression
    initial_bumps = all_bumps(x_init)
    xs = [x_init for _ in range(n_split)]

    # initialize result containers
    bumps: List[List[int]] = []
    trainerrs: List[List[float]] = []
    testerrs: List[List[float]] = []

    for i in range(n_split):
        train, test = datas[i]
        trainA, trainB = train.xy()

        # Reference overfit solution
        res = np.linalg.lstsq(trainA, trainB, rcond=-1)[0]

        # explore tradeoff off regularization vs error on subset of data
        for bump in [5, 10, 15, 100, 10000, 10000000]:

            xs[i] = fitfun(xs[i], bump, trainA, trainB, initial_bumps)
            trainerrs[i].append(train.mae(xs[i], 'ce'))
            testerrs[i].append(test.mae(xs[i], 'ce'))
            bumps[i].append(count_bumps(xs[i], initial_bumps))


if __name__ == '__main__':
    cv()
