import numpy as np
import pandas as pd
import random
from scipy.optimize import minimize
from rich import print


def read_edges(filename: str):
    print("Reading edges from", filename)
    df = pd.read_csv(filename)
    # keep only the relevant columns
    df = df[
        [
            "ATOM_SERIAL_NUMBER_I",
            "ATOM_SERIAL_NUMBER_J",
            "LOWER_BOUND_IJ",
            "UPPER_BOUND_IJ",
        ]
    ]
    # change column names
    df.columns = ["i", "j", "l", "u"]
    # V is the union of E['i'] and E['j']
    V = np.unique(df[["i", "j"]].values)
    # renumber the atoms
    V = {i: idx for idx, i in enumerate(V)}
    df["i"] = df["i"].map(V)
    df["j"] = df["j"].map(V)
    # convert to dict
    E = {}
    for i, j, l, u in df.itertuples(index=False):
        E[i, j] = (l, u)
    return V, E


def phi(lamb: float, tau: float, y: float) -> tuple[float, float]:
    kappa = np.sqrt(lamb**2 * y**2 + tau**2)
    f = lamb * y + kappa
    df = lamb + (lamb**2 * y) / kappa
    return f, df


def fobj(
    x: np.ndarray,
    V: dict[int, int],
    E: dict[tuple[int, int], tuple[float, float]],
    lamb: float,
    tau: float,
) -> tuple[float, np.ndarray]:
    N = len(V)
    f = 0.0
    x = x.reshape(N, 3)
    g = np.zeros_like(x)
    for (i, j), (lij, uij) in E.items():
        xij = x[i] - x[j]
        tij = np.sqrt(tau**2 + np.linalg.norm(xij) ** 2)
        phi_lij, dphi_lij = phi(lamb, tau, lij - tij)
        phi_uij, dphi_uij = phi(lamb, tau, tij - uij)
        f += phi_lij + phi_uij
        gij = ((dphi_uij - dphi_lij) / tij) * xij
        g[i] += gij
        g[j] += -gij
    return f, g.ravel()


def numerical_gradient(func, x0, eps=1e-6):
    """
    Approximate gradient of func at x0 via finite differences.
    func: callable that takes a 1D numpy array and returns a float.
    x0: 1D numpy array (point of evaluation).
    eps: step size for finite differences.
    """
    grad_approx = np.zeros_like(x0)
    for i in range(len(x0)):
        old_val = x0[i]
        # f(x + eps)
        x0[i] = old_val + eps
        f_plus = func(x0)
        # f(x - eps)
        x0[i] = old_val - eps
        f_minus = func(x0)
        # restore
        x0[i] = old_val
        grad_approx[i] = (f_plus - f_minus) / (2 * eps)
    return grad_approx


def test_gradient():
    np.random.seed(0)

    # Example data
    N = 5  # number of vertices
    lamb = 2.0
    tau = 0.5

    # Random geometry (N points in 3D)
    x = np.random.randn(N, 3).ravel()

    # Example edges: (i, j, lij, uij)
    # We'll just define some random i, j and random bounds
    V = {i: i for i in range(N)}
    E = dict()
    for _ in range(4 * N):  # just pick 4N random edges
        i = np.random.randint(0, N)
        j = np.random.randint(0, N)
        if i == j:  # skip self-edges
            continue
        lij = np.random.rand()  # random lower bound
        uij = lij + np.random.rand()  # random upper bound
        E[i, j] = (lij, uij)

    # Get the analytic gradient
    _, g_analytic = fobj(x, V, E, lamb, tau)

    # Get the numerical gradient
    g_numerical = numerical_gradient(lambda z: fobj(z, V, E, lamb, tau)[0], x)

    # Compare them
    diff = np.linalg.norm(g_analytic - g_numerical)
    rel_diff = diff / (1e-14 + np.linalg.norm(g_analytic) + np.linalg.norm(g_numerical))

    print(f"Analytic gradient: {g_analytic}")
    print(f"Numerical gradient: {g_numerical}")
    print(f"Absolute diff: {diff:g}")
    print(f"Relative diff: {rel_diff:g}")


def random_sphere() -> np.ndarray:
    x = np.random.randn(3)
    return x / np.linalg.norm(x)


def init_x(
    V: dict[int, int],
    E: dict[tuple[int, int], tuple[float, float]],
    seed: int = 0,
) -> np.ndarray:
    random.seed(seed)
    N = len(V)
    x = np.zeros((N, 3))
    L = set(i for i in range(N))
    # M is the set of neighbors of i
    M = {i: set() for i in range(N)}
    for i, j in E:
        M[i].add(j)
        M[j].add(i)
    while len(L) > 0:
        # random vertex from L
        i = random.choice(list(L))
        M_i = [j for j in M[i] if j in L]
        xi = x[i]
        for j in M_i:
            lij, uij = E[i, j] if (i, j) in E else E[j, i]
            dij = (lij + uij) / 2
            # xj lies on the sphere of center xi and radius dij
            x[j] = xi + dij * random_sphere()
        L.remove(i)
    return x.ravel()


def sph(V, E, lamb, rho, x, eps=1e-6):
    # tau is the third quartile of the distances
    D = [dij for (i, j), (lij, uij) in E.items() for dij in [lij, uij]]
    tau = np.percentile(D, 75)
    while fobj(x, V, E, lamb, tau)[0] > eps:
        x = minimize(fobj, x, args=(V, E, lamb, tau), jac=True)
        tau = rho * tau
    f = fobj(x, V, E, lamb, tau)[0]
    return x, f


def locmin(V, E, x, eps=1e-6):
    # Naive local minimizer
    lamb = 0.5
    tau = 1e-8
    x = minimize(fobj, x, args=(V, E, lamb, tau), jac=True)
    f = fobj(x, V, E, lamb, tau)[0]
    return x, f


def main():
    V, E = read_edges("instances/1GPV_A_N_100_eps_0.04.txt")
    num_trials = 10
    for seed in range(num_trials):
        x = init_x(V, E, seed=seed)
        x_loc, f_loc = locmin(V, E, x=x.copy())
        x_sph, f_sph = sph(V, E, lamb=0.5, rho=0.99, x=x.copy())
        print(f_sph, f_loc)


if __name__ == "__main__":
    test_gradient()
    main()
