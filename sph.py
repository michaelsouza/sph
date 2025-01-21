import numpy as np
import pandas as pd


def read_edges(filename: str) -> tuple[dict[int, int], pd.DataFrame]:
    print("Reading edges from", filename)
    E = pd.read_csv(filename)
    # keep only the relevant columns
    E = E[
        [
            "ATOM_SERIAL_NUMBER_I",
            "ATOM_SERIAL_NUMBER_J",
            "LOWER_BOUND_IJ",
            "UPPER_BOUND_IJ",
        ]
    ]
    # change column names
    E.columns = ["i", "j", "l", "u"]
    # V is the union of E['i'] and E['j']
    V = np.unique(E[["i", "j"]].values)
    # renumber the atoms
    V = {i: idx for idx, i in enumerate(V)}
    E["i"] = E["i"].map(V)
    E["j"] = E["j"].map(V)
    return V, E


def phi(lamb: float, tau: float, y: float) -> tuple[float, float]:
    kappa = np.sqrt(lamb**2 * y**2 + tau**2)
    f = lamb * y + kappa
    df = lamb + (lamb**2 * y) / kappa
    return f, df


def fobj(
    E: list[tuple[int, int, float, float]], lamb: float, tau: float, x: np.ndarray
) -> tuple[float, np.ndarray]:
    f = 0.0
    g = np.zeros_like(x)
    for i, j, lij, uij in E:
        xij = x[i] - x[j]
        tij = np.sqrt(tau**2 + np.linalg.norm(xij) ** 2)
        phi_lij, dphi_lij = phi(lamb, tau, lij - tij)
        phi_uij, dphi_uij = phi(lamb, tau, tij - uij)
        f += phi_lij + phi_uij
        gij = ((dphi_uij - dphi_lij) / tij) * xij
        g[i] += gij
        g[j] += -gij
    return f, g


def fobj_flat(
    flat_x: np.ndarray, E: pd.DataFrame, lamb: float, tau: float, N: int
) -> float:
    x = flat_x.reshape(N, 3)
    f_val, _ = fobj(E, lamb, tau, x)
    return f_val


def grad_flat(
    flat_x: np.ndarray, E: pd.DataFrame, lamb: float, tau: float, N: int
) -> np.ndarray:
    """
    Returns the analytic gradient as a flattened vector.
    """
    x = flat_x.reshape(N, 3)
    _, g = fobj(E, lamb, tau, x)
    return g.ravel()


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
    x_init = np.random.randn(N, 3)

    # Example edges: (i, j, lij, uij)
    # We'll just define some random i, j and random bounds
    E = []
    for _ in range(4 * N):  # just pick 4N random edges
        i = np.random.randint(0, N)
        j = np.random.randint(0, N)
        if i == j:  # skip self-edges
            continue
        lij = np.random.rand()  # random lower bound
        uij = lij + np.random.rand()  # random upper bound
        E.append((i, j, lij, uij))

    # Flatten x for the numerical derivative
    x_flat = x_init.ravel()

    # Get the analytic gradient
    g_analytic = grad_flat(x_flat, E, lamb, tau, N)

    # Get the numerical gradient
    g_numerical = numerical_gradient(lambda z: fobj_flat(z, E, lamb, tau, N), x_flat)

    # Compare them
    diff = np.linalg.norm(g_analytic - g_numerical)
    rel_diff = diff / (1e-14 + np.linalg.norm(g_analytic) + np.linalg.norm(g_numerical))

    print("Analytic gradient:", g_analytic)
    print("Numerical gradient:", g_numerical)
    print(f"Absolute diff: {diff:g}")
    print(f"Relative diff: {rel_diff:g}")


def main():
    V, E = read_edges("instances/1GPV_A_N_100_eps_0.04.txt")
    print("Number of atoms:", len(V))
    print("Number of edges:", len(E))


if __name__ == "__main__":
    test_gradient()
    main()
