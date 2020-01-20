# External
from typing import (Any,
                    Dict as D,
                    List as L,
                    Union as U,
                    Tuple as T)

import os
import numpy as np
import shutil
import json
import hashlib
import tempfile
import scipy as sp
import psycopg2
import plotly
import plotly.graph_objs as go

# Internal
import functionals.fit.constraint as constr
import functionals.fit.data as fdata
import functionals.fit.functional as fx

# Type synonyms
Connection = Any
Arrs = L[np.ndarray]
Metric = U[str, D[str, float]]
###############################################################################
# CONSTANTS
###########
errtypes = ['ce', 'bm', 'lc']

# unequally weigh matrix elements in imporance
bias = np.array([1.2**(x % 8 + x//8) for x in range(64)])

root = '/'+os.path.join(*__file__.split('/')[:-3])

q3 = '''SELECT a_ce,a_bm,a_l,b_ce,b_bm,b_l,np.expt_ce,np.expt_bm,np.expt_vol,name,ce
        FROM bulks JOIN job ON job=job_id
        WHERE calc = %s'''


q4 = '''SELECT pw,data,fitdata
        FROM calc JOIN functional on functional=functional_id
        WHERE calc_id = %s'''

###############################################################################
# Helper functions
# ----------------


def pltstr(x: go.Figure, show: bool = True) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        plotly.offline.plot(x, filename=tmp+'/temp0.html', auto_open=show)
        with open(tmp+'/temp0.html', 'r') as f:
            out = f.read()
    return out


def sqlselect(conn: Connection, q: str,
              binds: L[Any] = []) -> L[T[Any, ...]]:
    with conn.cursor() as cxn:
        cxn.execute(q, vars=binds)
        return list(cxn.fetchall())


def hash_(x: Any) -> str:
    return hashlib.sha512(str(x).encode()).hexdigest()[:10]

###############################################################################
# Classes
#########

    #################
    # VISUALIZATION #
    #################

    def fxviz(self, name: str = 'fit', show: bool = True) -> str:
        steps, data = [], []
        allerr = zip(*[self.test_err[x] for x in errtypes])
        cviols = self.cviols().tolist()
        mincv, maxcv = min(cviols), max(cviols)
        it = zip(self.a0, self.a1, allerr, cviols, self.n_sat)
        minmax = [[self.test_err[x][i] for i in [0, -1]] for x in errtypes]
        heights = [1.9, 2., 2.1, 2.2]
        names = ['(α=0)', '(α=1)']
        # import pdb; pdb.set_trace()
        for i, (a0, a1, err, cv, n_sat) in enumerate(it):
            As, stys = [a0, a1], ['solid', 'dot']
            zip1 = zip(As, stys, names)
            zip2 = zip(errtypes, err, heights, minmax)
            plts = [go.Scatter(x=np.linspace(0, 5, 20).tolist(), y=a, name=n,
                               line=dict(dash=sty, shape='spline'))
                    for a, sty, n in zip1]
            errs = [go.Scatter(x=[0, 5*(e-mn)/((mx-mn) or 1)], y=[h, h],
                               mode='lines', name='%s (∆=%.2E)' % (x, mx-mn))
                    for x, e, h, (mn, mx) in zip2]
            x = 5*(cv-mincv)/((maxcv-mincv) or 1)
            cv = [go.Scatter(x=[0, x], y=[2.2, 2.2], mode='lines',
                             name='Const viol: %f -> %f' % (maxcv, mincv))]
            data.extend(plts+errs+cv)
            titles = [' %s:%f ' % (x, e) for x, e in zip(errtypes, err)]
            sat = "Sat: " + ', '.join(self.cons[:n_sat])
            title = sat+'   '+' '.join(titles)
            args = [dict(visible=[False] * 6*len(self)), {'title.text': title}]
            trues = [True, True, True, True, True, True]
            args[0]['visible'][6*i:6*i+6] = trues
            steps.append(dict(args=args, method='update'))

        layout = go.Layout(hovermode='closest', legend=dict(x=1, y=1),
                           xaxis=dict(title='s', rangemode="tozero"),
                           yaxis=dict(title='Fx', rangemode="tozero"),
                           sliders=[dict(steps=steps)],)
        fig = go.Figure(data=data, layout=layout)
        return pltstr(fig, show=show)

    def plot(self, met: Metric, cons: L[str] = None, show: bool = True) -> str:
        '''Plots a metric over the trajectory.'''
        xs = self.cviols(cons)

        trainys = self.weighted_err(met, train=True)
        data = [go.Scatter(x=xs, y=trainys, name='train', mode='markers',
                           text=[str(x) for x in range(len(xs))],
                           hoverinfo='text')]

        if self.test != self.train:
            testys = self.weighted_err(met, train=False)
            data.append(go.Scatter(x=xs, y=testys, name='test', mode='markers',
                                   text=[str(x) for x in range(len(xs))],
                                   hoverinfo='text'))
        annotations = []  # type: L[Any]

        layout = go.Layout(title='Trajectory of constrained optimization',
                           hovermode='closest',
                           xaxis=dict(title='Constraint Violation', ticklen=5,
                                      zeroline=False, gridwidth=2),
                           yaxis=dict(title=str(met), ticklen=5, gridwidth=2),
                           annotations=annotations)

        fig = go.Figure(data=data, layout=layout)

        return pltstr(fig, show=show)

    def resid(self, step: int = -1, train: bool = False, show: bool = True
              ) -> str:
        '''Plots residual BarPlot across all materials.'''
        data = self.train if train else self.test
        datalist = []
        x = np.array(self.xs[step])
        for m in errtypes:
            mats, ys = [], []
            for d in getattr(data, m):
                mats.append(d.mat)
                ys.append(d.err(x))
            datalist.append(go.Bar(x=mats, y=ys, name=m))

        layout = go.Layout(title='Residuals ', hovermode='closest',
                           xaxis=dict(title='Material'),
                           yaxis=dict(title='Error', ticklen=5, gridwidth=2))
        fig = go.Figure(data=datalist, layout=layout)
        return pltstr(fig, show=show)





def cost(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ x - B


def dcost(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A


def cost_scalar(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    res = A @ x - B
    return res @ res / 10000000


def dcost_scalar(x: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    res = A @ x - B
    return 2 * res.T @ A / 10000000


class Fit(object):
    """
    Fit a 8 x 8 matrix of coefficients corresponding to Legendre polynomial
    basis functions to DFT data + constraints
    """
    def __init__(self, calc: D[str, Any], data: fdata.Data, cons: L[str],
                 cescale: float, bmscale: float, lcscale: float) -> None:
        self.data = data
        self.calc = calc
        self.cons = cons
        self.cescale = cescale
        self.bmscale = bmscale
        self.lcscale = lcscale

        # Constant parameters for sp.optimize.least_squares
        self.kwargs = dict(
            method='trf', ftol=1e-15, xtol=1e-15, gtol=1e-15,
            # bounds=([-float('inf')]*64, [float('inf')]*64),
            bounds=([-2.]*64, [2.]*64),
            x_scale='jac', loss='linear', max_nfev=50000,
            verbose=0)

    def __eq__(self, other: object) -> bool:
        fields = ['calc', 'cons', 'cescale', 'bmscale', 'lcscale', 'data']
        return False if not isinstance(other, Fit) else \
            all([getattr(self, x) == getattr(other, x) for x in fields])

    ###############
    # MAIN METHOD #
    ###############
    def fit(self, pth: str, n: int = 30) -> None:
        '''Cross validate on 30 random 50/50 splits of the data.'''
        tpth = os.path.join(pth, 'fit.json')
        if not os.path.exists(tpth):
            t = self.fittraj2().to_dict()
            with open(tpth, 'w') as f:
                json.dump(t, f)
            print('finished main fit')

        # Do cross validation
        for i in range(3):
            tpth = os.path.join(pth, ' x%d.json' % i)
            trn, tst = self.data.split(i)
            if not os.path.exists(tpth):
                t = self.fittraj2(train=trn, test=tst).to_dict()
                with open(tpth, 'w') as f:
                    json.dump(t, f)
                print('finished CV fit %d' % i)

    def fittraj2(self, train: fdata.Data = None, test: fdata.Data = None
                 ) -> Traj:
        if train is None and test is None:
            train, test = (self.data, self.data)
        assert train is not None and test is not None
        t = Traj(train, test, cons=self.cons)
        cons = [constr.consts[c] for c in self.cons]  # ORDER IMPORTANT
        cons_eq = [c for c in cons if c.eq]  # ORDER IMPORTANT
        cons_lt = [c for c in cons if not c.eq]  # ORDER IMPORTANT
        eq_ab, lt_ab = [[c.ab() for c in cs] for cs in [cons_eq, cons_lt]]
        A_expt, B_expt = train.xy(self.cescale, self.bmscale, self.lcscale)

        A_const = np.vstack([a for (a, _) in eq_ab+lt_ab])
        B_const = np.concatenate([b for (_, b) in eq_ab] +
                                 [c.default*np.ones(len(c)) for c in cons_lt])

        x0 = sp.optimize.least_squares(
            cost, fx.BEEF().x, jac=dcost,
            kwargs={'A': A_const, 'B': B_const}, **self.kwargs).x
        # print('Const err: ', np.linalg.norm(A_const@x0-B_const))
        constraints = constr.Constraint.lc(cons)
        x = sp.optimize.minimize(
            method='slsqp', x0=x0, args=(A_expt, B_expt),
            fun=cost_scalar, jac=dcost_scalar,
            constraints=constraints,  # bounds=[(-2., 2.)]*64,
            options=dict(disp=True, maxiter=5000)).x
        # for k in errtypes: print(k, train.mse(x, k))
        # pltstr(fx.plots([fx.FromMatrix(x)]))
        # pltstr(fx.FromMatrix(x).plot3d())
        return t

    ################
    # CONSTRUCTORS #
    ################
    @classmethod
    def from_db(cls, db: str, calc_id: int, constr: str,
                ce: float, bm: float, lc: float) -> 'Fit':
        with open(db, 'r') as fi:
            kwargs = json.load(fi)
            kwargs['dbname'] = kwargs.pop('db')
            kwargs['password'] = kwargs.pop('passwd')
            conn = psycopg2.connect(**kwargs)

        # Extract fitting parameters
        constraints = constr.split()
        ce_scale, bm_scale, lc_scale = map(float, [ce, bm, lc])

        # Extract calculation parameters
        pw_, fx_, data_ = sqlselect(conn, q4, [calc_id])[0]
        pw = float(pw_)
        calc = dict(pw=pw, fx=json.loads(fx_))

        # Extract fitted data
        data = fdata.Data.from_list(json.loads(data_), full=True)

        return cls(cons=constraints, calc=calc, cescale=ce_scale,
                   bmscale=bm_scale, lcscale=lc_scale, data=data)

    @classmethod
    def from_json(cls, pth: str) -> 'Fit':
        with open(os.path.join(pth, 'metadata.json'), 'r') as fi:
            md = json.load(fi)
        with open(os.path.join(pth, 'data.json'), 'r') as fi:
            data = json.load(fi)
        p = md['params']
        return cls(cons=p['c'], calc=md['calc'], lcscale=p['lc_scale'],
                   data=fdata.Data.from_list(data),
                   cescale=p['ce_scale'], bmscale=p['bm_scale'])

    ##########
    # HELPER #
    ##########

    def metadata(self) -> D[str, Any]:
        '''Summary of parameters used to determine fit inputs.'''
        uid = hash_([self.calc[x] for x in ['pw', 'fx']] +
                    [self.cescale, self.bmscale, self.lcscale, self.cons])
        md = dict(calc=self.calc, uid=uid,
                  params=dict(ce_scale=self.cescale, bm_scale=self.bmscale,
                              lc_scale=self.lcscale, c=self.cons))
        return md

    def write(self, pth: str) -> None:
        rootpth = os.path.join(root, 'functionals/scripts/fit/')

        def write(fi: str, x: Any) -> None:
            with open(os.path.join(pth, fi)+'.json', 'w') as file:
                json.dump(x, file)

        # Write to directory
        # ------------------
        write('metadata', self.metadata())
        write('data', self.data.to_list())
        shutil.copyfile(rootpth+'runfit.py', os.path.join(pth, 'runfit.py'))
        shutil.copyfile(rootpth+'subfit.sh', os.path.join(pth, 'subfit.sh'))
        os.system('chmod 755 '+os.path.join(pth, 'subfit.sh'))


if __name__ == '__main__':
    fr = FitResult.from_pth('/Users/ksb/scp_tmp/vauto/fit/?')
    fr.full.fxviz()



# class ResX(object):
#     '''Pair of results, one constrained, other unconstrained.'''

#     def __init__(self, train: fdata.Data, test: fdata.Data,
#                  cons: L[str], opt: Optional[int] = None) -> None:
#         self.train = train
#         self.test = test
#         self.cons = cons
#         # Initialize state
#         self.xs = []  # type: L[np.ndarray]
#         self.train_err = collections.defaultdict(list)  # type: D[str,L[float]]
#         self.test_err = collections.defaultdict(list)  # type: D[str,L[float]]
#         self.cviol = collections.defaultdict(list)  # type: D[str,L[float]]
#         self.a0 = []  # type: L[L[float]]
#         self.a1 = []  # type: L[L[float]]
#         self._opt = opt
#         self.n_sat = []  # type: L[int]

#     def __len__(self) -> int: return len(self.xs)

#     def __eq__(self, other: object) -> bool:
#         return False if not isinstance(other, Traj) else \
#                 (vars(self)) != (vars(other))

#     # Serialization/Deserialization
#     def to_dict(self) -> D[str, Any]:
#         keys = ['train_err', 'test_err', 'cviol', 'a0', 'a1',
#                 'n_sat', 'cons']
#         return dict(train=self.train.to_list(), test=self.test.to_list(),
#                     xs=[x.tolist() for x in self.xs],
#                     **{x: getattr(self, x) for x in keys})

#     @classmethod
#     def from_dict(cls, d: D[str, Any]) -> 'Traj':
#         '''Inverse of to_dict.'''
#         t = cls(train=fdata.Data.from_list(d['train']),
#                 test=fdata.Data.from_list(d['test']),
#                 opt=d.get('opt'), cons=d['cons'])
#         t.xs = [np.array(x) for x in d['xs']]
#         for x in ['train_err', 'test_err', 'cviol', 'n_sat', 'a0', 'a1']:
#             setattr(t, x, d[x])
#         return t

#     def add(self, x: np.ndarray, cv: D[str, float], n_sat: int) -> None:
#         '''Add a step to the trajectory. Compute metrics upon insert.'''
#         self.xs.append(x)
#         self.n_sat.append(n_sat)
#         for k in errtypes:
#             self.train_err[k].append(self.train.mse(x, k))
#             self.test_err[k].append(self.test.mse(x, k))
#             # print('err '+k, self.train_err[k][-1])

#         for c, err in cv.items():
#             self.cviol[c].append(err)

#         fxl = fx.FromMatrix(x)
#         srange = np.linspace(0, 5)
#         for a in [0, 1]:
#             getattr(self, 'a%d' % a).append([fxl.apply(s, a) for s in srange])

#     def gap(self, met: Metric) -> np.ndarray:
#         '''Gap between train and test.'''
#         return np.maximum(self.weighted_err(met, train=False)
#                           - self.weighted_err(met, train=True), 0.0)

#     def weighted_err(self, met: Metric, train: bool) -> np.ndarray:
#         metdic = {met: 1.} if isinstance(met, str) else met
#         it = self.train_err if train else self.test_err
#         errs = [metdic.get(k, 0)*np.array(v) for k, v in it.items()]
#         return np.sum(errs, axis=0)

#     def cviols(self, cons: L[str] = None) -> np.ndarray:
#         '''np.sum all constraint violations found in cons'''
#         if not self.cviol:
#             return np.zeros(len(self))
#         return np.sum([np.array(x) for k, x in self.cviol.items()
#                        if (cons is None) or (k in cons)], axis=0)

#     def opt_err(self) -> T[Optional[float], Optional[float], Optional[float]]:
#         opt = self.opt
#         if opt is None:
#             return None, None, None
#         else:
#             a, b, c = [self.test_err[x][opt] for x in errtypes]
#             return a, b, c

#     @property
#     def opt(self) -> Optional[int]:
#         return -1 if self._opt is None else self._opt

    # def fittraj(self, train: fdata.Data = None, test: fdata.Data = None
    #             ) -> Traj:
    #     '''Trajectory of constrained solutions to Ax=b until convergence.'''
    #     # Initialize fitting
    #     #####################
    #     if train is None and test is None:
    #         train, test = (self.data, self.data)
    #     assert train is not None and test is not None
    #     print('FORCE CONS AND SCALE')
    #     self.cons = ['lda', 'zerocurv', 'pos', 'liebox',
    #                  'curvpos', 'curvneg', 'hnorm', 'scan11', 'decay']
    #     self.cescale, self.bmscale, self.lcscale = 0.01, .5, .1
    #     t = Traj(train, test, cons=self.cons)
    #     cons = [constr.consts[c] for c in self.cons]  # ORDER IMPORTANT
    #     cons_eq = [c for c in cons if c.eq]  # ORDER IMPORTANT
    #     cons_lt = [c for c in cons if not c.eq]  # ORDER IMPORTANT
    #     eq_ab, lt_ab = [[c.ab() for c in cs] for cs in [cons_eq, cons_lt]]
    #     weights = {c: 1e-15 for c in self.cons}
    #     satisfied_constraints, no_change, changed = 0, 0, 0
    #     A_expt, B_expt = train.xy(self.cescale, self.bmscale, self.lcscale)
    #     snapshot = np.ones(64) * 0.001
    #     #  x = np.linalg.lstsq(A_expt, B_expt, rcond=None)[0]  # init
    #     x = sp.optimize.least_squares(
    #             cost, fx.BEEF().x, jac=dcost,
    #             kwargs={'A': A_expt, 'B': B_expt}, **self.kwargs).x
    #     print('INIT ERR', *[(k, train.mse(x, k)) for k in errtypes])
    #     # Constrain
    #     ###########
    #     for counter in itertools.count():
    #         # Equality constraints
    #         A_const_eq = np.vstack(
    #             [a * weights[c.name] for c, (a, _) in zip(cons_eq, eq_ab)]
    #             or [np.zeros((0, 64))])
    #         B_const_eq = np.concatenate(
    #             [b * weights[c.name] for c, (_, b) in zip(cons_eq, eq_ab)]
    #             or [[]])
    #         # Treat inequality constraints
    #         AB_inds = [a@x-b > 0 for a, b in lt_ab]  # ONLY want these rows
    #         A_const_lt = np.vstack(
    #             [a[i, :] * weights[c.name]
    #              for c, i, (a, _) in zip(cons_lt, AB_inds, lt_ab)])
    #         B_const_lt = np.concatenate(
    #             [b[i] * weights[c.name]
    #              for c, i, (_, b) in zip(cons_lt, AB_inds, lt_ab)])

    #         # Assemble A, B matrices
    #         A = np.vstack((A_expt, A_const_eq, A_const_lt))
    #         B = np.concatenate((B_expt, B_const_eq, B_const_lt))
    #         # Do fitting
    #         kwargs = dict(A=A, B=B)
    #         res = sp.optimize.least_squares(
    #             cost, x, jac=dcost, kwargs=kwargs, **self.kwargs)

    #         # Double constraint weights for next unsatisifed constraint
    #         const_eq_err = {c.name: np.square(a@x - b)
    #                         for c, (a, b) in zip(cons_eq, eq_ab)}
    #         const_lt_err = {
    #             c.name: np.square(a[i, :]@x-b[i]) for c, i, (a, b)
    #             in zip(cons_lt, AB_inds, lt_ab)}
    #         const_viol = {**const_eq_err, **const_lt_err}
    #         cv = {k: np.sum(v) for k, v in const_viol.items()}
    #         for con in cons[:satisfied_constraints+1]:
    #             last = con.name == self.cons[satisfied_constraints]
    #             # the 'initial' is needed in case of empty list
    #             # (i.e. a LT constraint where no point is active)
    #             if np.max(const_viol[con.name], initial=0.) > con.tol:
    #                 weights[con.name] *= 1+con.dw if last else 1+5*con.dw
    #             elif last:
    #                 satisfied_constraints += 1
    #                print('Sat: '+' '.join(self.cons[:satisfied_constraints]))

    #                 if satisfied_constraints == len(cons):
    #                     t.add(x, cv, satisfied_constraints)
    #                     return t
    #                 else:
    #                     with open('tmp.json', 'w') as f:
    #                         json.dump(t.to_dict(), f)
    #                 snapshot = np.ones(64)*0.001  # make sure this is added
    #         x = res.x

    #         # how much has solution changed
    #         dx = np.mean(np.abs(np.divide(snapshot - x, snapshot)))
    #         print('dx %.2E (%d) %s viol: %.2E, dw=%.2E' % (dx, len(t),
    #               self.cons[satisfied_constraints],
    #               np.max(const_viol[self.cons[satisfied_constraints]]),
    #               cons[satisfied_constraints].dw))
    #         if dx > 0.2:
    #             t.add(x, cv, satisfied_constraints)  # update traj
    #             snapshot = x
    #             changed = counter
    #             # Check if we have been changing TOO frequently
    #             if counter - no_change > 1:
    #                 no_change = counter
    #                 print('2 consecutive steps of big dx change')
    #                 cons[satisfied_constraints].dw /= 1.1
    #         else:
    #             no_change = counter  # record last step of no change
    #             # Check if we haven't been changing enough
    #             if counter - changed > 100:
    #                 changed = counter
    #                 print('100 consecutive steps without dx change')
    #                 cons[satisfied_constraints].dw *= 1.1

    #     raise ValueError()
# class FitResult(object):
#     def __init__(self, full: Traj, cv: L[Traj]) -> None:
#         self.full = full
#         self.cv = cv

#     @classmethod
#     def from_pth(cls, pth: str) -> 'FitResult':
#         cv = []
#         with open(os.path.join(pth, 'fit.json'), 'r') as f:
#             data = Traj.from_dict(json.load(f))
#         for x in filter(lambda x: x[0] == 'x', os.listdir(pth)):
#             with open(os.path.join(pth, x), 'r') as f:
#                 try:
#                     cv.append(Traj.from_dict(json.load(f)))
#                 except OSError:
#                     print('error loading x ', pth, x)
#         return cls(data, cv)

#     def plot_cv(self, met: Metric, show: bool = True) -> str:
#         '''Plots a metric over the trajectory.'''
#         kwargs = dict(mode='markers', hoverinfo='text')
#         data = []
#         for i, t in enumerate(self.cv):
#             x = t.cviols()
#             y = t.weighted_err(met, train=True)
#             z = t.weighted_err(met, train=False)
#             text = [str(x) for x in range(len(x))]
#             data.extend([go.Scatter(x=x, y=y, name='train_%d' % i,
#                                     text=text, **kwargs),
#                          go.Scatter(x=x, y=z, name='test_%d' % i, **kwargs)])

#         layout = go.Layout(
#             title='Trajectory of constrained opt', hovermode='closest',
#             xaxis=dict(title='∆Train/Test', ticklen=5, zeroline=False,
#                        gridwidth=2),
#             yaxis=dict(title='Test error', ticklen=5, gridwidth=2,))

#         fig = go.Figure(data=data, layout=layout)
#         return pltstr(fig, show=show)

#     def plot_transfer(self, met: Metric, cons: L[str] = None,
#                       show: bool = True) -> str:
#         fullx = self.full.cviols()
#         fully = self.full.weighted_err(met, train=False)
#         tx, ty, tx_, ty_ = self.transferability(met, cons)

#         kwargs = dict(mode='markers', text=[str(x) for x in range(len(fullx))],
#                       hoverinfo='text')
#         data = [go.Scatter(name='full', x=fullx, y=fully, **kwargs),
#                 go.Scatter(name='∆Err', x=tx, y=ty, yaxis='y2', **kwargs),
#                 go.Scatter(name='∆Err_smooth', x=tx_, y=ty_, yaxis='y2',
#                            **kwargs), ]
#         layout = go.Layout(
#             title='Transferability', hovermode='closest',
#             xaxis=dict(title='Constraint Violation', ticklen=5, zeroline=False,
#                        gridwidth=2),
#             yaxis=dict(title=str(met), ticklen=5, gridwidth=2,),
#             yaxis2=dict(title='∆Err', overlaying='y', side='right'))

#         fig = go.Figure(data=data, layout=layout)
#         return pltstr(fig, show=show)

#     def transferability(self, met: Metric, cons: L[str] = None
#                         ) -> T[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
#         '''For each data point, compute average ErrTest / ErrTrain ratio'''
#         xs, tsterr, trnerr = [], [], []  # type: ignore
#         for traj in self.cv:
#             xs.extend(traj.cviols(cons))
#             tsterr.extend(traj.weighted_err(met, train=False))
#             trnerr.extend(traj.weighted_err(met, train=True))

#         # KNN regression
#         knn = neighbors.KNeighborsRegressor(20)
#         ys = np.maximum(np.array(tsterr)-np.array(trnerr), 0.0)
#         x_ = np.array(xs).reshape((-1, 1))

#         outx = self.full.cviols()

#         outy = knn.fit(X=x_, y=ys).predict(outx.reshape((-1, 1)))

#         return xs, ys, outx, outy