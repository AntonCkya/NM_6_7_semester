from abc import ABC, abstractmethod
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def tridiagonal_solve(a, b, c, d) -> np.ndarray:
    n = len(d)
    p = np.ndarray(n, dtype=float)
    q = np.ndarray(n, dtype=float)
    x = np.ndarray(n, dtype=float)

    p[0] = -c[0] / b[0]
    q[0] = d[0] / b[0]

    for i in range(1, n):
        p[i] = -c[i] / (b[i] + a[i]*p[i-1])
        q[i] = (d[i] - a[i]*q[i-1]) / (b[i] + a[i]*p[i-1])

    x[-1] = q[-1]
    for i in range(n-2, -1, -1):
        x[i] = p[i] * x[i+1] + q[i]
    return x


class Diffur:
    @staticmethod
    def g(x, t): return 0.5 * np.exp(-0.5 * t) * np.sin(x)

    @staticmethod
    def phi_zero(t): return np.exp(-0.5 * t)

    @staticmethod
    def phi_l(t): return -1 * np.exp(-0.5 * t)

    @staticmethod
    def psi(x): return np.sin(x)

    def __init__(self):
        pass


class AbstractSolver(ABC):
    D: Diffur
    T: float
    L: float
    N: float
    K: float
    tau: float
    h: float
    sigma: float

    @staticmethod
    def calc_tau(T: float, K: float) -> float:
        return T / (K-1)
    
    @staticmethod
    def calc_h(L: float, N: float) -> float:
        return L / (N-1)
    
    @staticmethod
    def calc_sigma(tau: float, h: float) -> float:
        return tau / h**2

    def __init__(self, T, L, K, N):
        self.D = Diffur()
        self.T = T
        self.L = L
        self.N = N
        self.K = K
        self.tau = self.calc_tau(T, K)
        self.h = self.calc_h(L, N)
        self.sigma = self.calc_sigma(self.tau, self.h)

    @abstractmethod
    def solve(self): pass

    def change_params(self, T, L, K, N):
        self.T = T
        self.L = L
        self.N = N
        self.K = K
        self.tau = self.calc_tau(T, K)
        self.h = self.calc_h(L, N)
        self.sigma = self.calc_sigma(self.tau, self.h)


class SolveExact(AbstractSolver):
    @staticmethod
    def exact_solve(x, t):
        return np.exp(-0.5 * t) * np.sin(x)

    def solve(self):
        u = np.zeros((self.K, self.N))

        for i in range(self.K):
            for j in range(self.N):
                u[i][j] = self.exact_solve(j * self.h, i * self.tau)
        return u


class SolveExplicit(AbstractSolver):
    aprox = 0
    """
    тут и далее
    0 - двухточечная аппроксимация с первым порядком
    1 - трехточечная аппроксимация со вторым порядком
    2 - двухточечная аппроксимация со вторым порядком
    """

    def add_aprox(self, num: int):
        self.aprox = num
        return self
    
    def zero_aprox(self):
        if self.aprox == 0:
            return lambda i, curr, prev: curr[1] - self.h * self.D.phi_zero(i * self.tau)
        elif self.aprox == 1:
            return lambda i, curr, prev: (-1 / 3) * (2 * self.h * self.D.phi_zero(i * self.tau) + curr[2] - 4 * curr[1])
        elif self.aprox == 2:
            return lambda i, curr, prev: -2 * self.sigma * self.h * self.D.phi_zero((i-1) * self.tau) + 2 * self.sigma * prev[1] + (1 - 2 * self.sigma) * prev[0] + self.tau * self.D.g(0, (i-1) * self.tau)

    def l_aprox(self):
        if self.aprox == 0:
            return lambda i, curr, prev: self.h * self.D.phi_l(i * self.tau) + curr[-2]
        elif self.aprox == 1:
            return lambda i, curr, prev: (1 / 3) * (2 * self.h * self.D.phi_l(i * self.tau) + 4 * curr[-2] - curr[-3])
        elif self.aprox == 2:
            return lambda i, curr, prev: 2 * self.sigma * self.h * self.D.phi_l((i-1) * self.tau) + 2 * self.sigma * prev[-2] + (1 - 2 * self.sigma) * prev[-1] + self.tau * self.D.g((self.N-1) * self.h, (i-1) * self.tau)

    def solve(self):
        u = np.zeros((self.K, self.N))

        u_zero = []
        for i in [i * self.h for i in range(self.N)]:
            u_zero.append(self.D.psi(i))
        u[0] = np.array(u_zero)

        for i in range(1, self.K):
            for j in range(1, self.N - 1):
                u[i][j] = (self.sigma * u[i-1][j-1] +
                           (1 - 2 * self.sigma) * u[i-1][j] +
                           self.sigma * u[i-1][j+1]
                           + self.tau * self.D.g(j * self.h, (i-1) * self.tau))
            u[i][0] = self.zero_aprox()(i, u[i], u[i-1])
            u[i][-1] = self.l_aprox()(i, u[i], u[i-1])
        
        return u


class SolveImplicit(AbstractSolver):
    aprox = 0

    def add_aprox(self, num: int):
        self.aprox = num
        return self

    def zero_aprox(self):
        if self.aprox == 0:
            return lambda i, curr: (0, -1, 1, self.h * self.D.phi_zero(i * self.tau))
        elif self.aprox == 1:
            return lambda i, curr: (0,
                                    -2*self.sigma - 1,
                                    2*self.sigma,
                                    2 * self.sigma * self.h * self.D.phi_zero(i * self.tau) - (curr[0] + self.tau * self.D.g(0, i * self.tau))
            )
        elif self.aprox == 2:
            return lambda i, curr: (0,
                                    -2*self.sigma - 1,
                                    2*self.sigma,
                                    2 * self.sigma * self.h * self.D.phi_zero(i * self.tau) - (curr[0] + self.tau * self.D.g(0, i * self.tau))
            )

    def l_aprox(self):
        if self.aprox == 0:
            return lambda i, curr: (-1, 1, 0, self.h * self.D.phi_l(i * self.tau))
        elif self.aprox == 1:
            return lambda i, curr: (-4 + (1 + 2*self.sigma) / self.sigma,
                                    2,
                                    0,
                                    2 * self.sigma * self.h * self.D.phi_l(i * self.tau) + (curr[-2] + self.tau * self.D.g((self.N-2) * self.h, i * self.tau))
            )
        elif self.aprox == 2:
            return lambda i, curr: (2*self.sigma,
                                    -2 * self.sigma - 1,
                                    0,
                                    -2 * self.sigma * self.h * self.D.phi_l(i * self.tau) - (curr[-1] + self.tau * self.D.g((self.N-1) * self.h, i * self.tau))
            )

    def solve(self):
        u = np.zeros((self.K, self.N))

        u_zero = []
        for i in [i * self.h for i in range(self.N)]:
            u_zero.append(self.D.psi(i))
        u[0] = np.array(u_zero)

        for i in range(1, self.K):
            a = np.zeros(self.N)
            b = np.zeros(self.N)
            c = np.zeros(self.N)
            d = np.zeros(self.N)
            for j in range(1, self.N - 1):
                a[j] = self.sigma
                b[j] = -1 - 2 * self.sigma
                c[j] = self.sigma
                d[j] = -self.tau * self.D.g(j * self.h, i * self.tau) - u[i - 1][j]
            a[0], b[0], c[0], d[0] = self.zero_aprox()(i, u[i-1])
            a[-1], b[-1], c[-1], d[-1] = self.l_aprox()(i, u[i-1])

            u[i] = tridiagonal_solve(a, b, c, d)
        
        return u
    

class SolveCN(AbstractSolver):
    aprox = 0

    def add_aprox(self, num: int):
        self.aprox = num
        return self

    def zero_aprox(self):
        if self.aprox == 0:
            return lambda i, curr: (0, -1, 1, self.h * self.D.phi_zero(i * self.tau))
        elif self.aprox == 1:
            return lambda i, curr: (0,
                                    -2,
                                    4 + (-1 - 2 * self.sigma) / self.sigma,
                                    (self.sigma * self.h * self.D.phi_zero(i * self.tau) - (curr[1] + self.tau * self.D.g(self.h, i * self.tau)) - 
                                    0.5*self.sigma * (curr[0] - 2*curr[1] + curr[2] + self.h**2 * self.D.g(self.h, (i-1) * self.tau))
                                    )
            )
        elif self.aprox == 2:
            return lambda i, curr: (0,
                                    -self.sigma - 1,
                                    self.sigma,
                                    (self.sigma * self.h * self.D.phi_zero(i * self.tau) - (curr[0] + 0.5 * self.tau * self.D.g(0, i * self.tau)) -
                                    self.sigma * (curr[1] - curr[0] - self.h * self.D.phi_zero((i-1) * self.tau) + 0.5 * self.h ** 2 * self.D.g(0, (i-1) * self.tau))
                                    )
            )

    def l_aprox(self):
        if self.aprox == 0:
            return lambda i, curr: (-1, 1, 0, self.h * self.D.phi_l(i * self.tau))
        elif self.aprox == 1:
            return lambda i, curr: (-4 + (1 + 2 * self.sigma) / self.sigma,
                                    2,
                                    0,
                                    (self.sigma * self.h * self.D.phi_l(i * self.tau) + (curr[-2] + 0.5 * self.tau * self.D.g((self.N-2) * self.h, i * self.tau)) + 
                                    0.5 * self.sigma * (curr[-3] - 2*curr[-2] + curr[-1] + self.h ** 2 * self.D.g((self.N-2) * self.h, (i-1) * self.tau))
                                    )
            )
        elif self.aprox == 2:
            return lambda i, curr: (self.sigma,
                                    -self.sigma - 1,
                                    0,
                                    (-self.sigma * self.h * self.D.phi_l(i * self.tau) - (curr[-1] + 0.5 * self.tau * self.D.g((self.N-1) * self.h, (i) * self.tau)) -
                                    self.sigma * (curr[-2] - curr[-1] + self.h * self.D.phi_l((i-1) * self.tau) + 0.5 * self.h ** 2 * self.D.g((self.N-1) * self.h, (i-1) * self.tau))
                                    )
            )

    def solve(self):
        u = np.zeros((self.K, self.N))

        u_zero = []
        for i in [i * self.h for i in range(self.N)]:
            u_zero.append(self.D.psi(i))
        u[0] = np.array(u_zero)

        for i in range(1, self.K):
            a = np.zeros(self.N)
            b = np.zeros(self.N)
            c = np.zeros(self.N)
            d = np.zeros(self.N)
            for j in range(1, self.N - 1):
                a[j] = 0.5 * self.sigma
                b[j] = -1 - self.sigma
                c[j] = 0.5 * self.sigma
                d[j] = (-0.5 * self.tau * self.D.g(j * self.h, i * self.tau) - u[i - 1][j] -
                        0.5 * self.sigma * (u[i - 1][j - 1] - 2 * u[i - 1][j] + u[i - 1][j + 1] + self.h ** 2 * self.D.g(j * self.h, (i-1) * self.tau))
                )
            a[0], b[0], c[0], d[0] = self.zero_aprox()(i, u[i-1])
            a[-1], b[-1], c[-1], d[-1] = self.l_aprox()(i, u[i-1])

            u[i] = tridiagonal_solve(a, b, c, d)
        
        return u
    

def MAE(numeric, analytic):
    return np.abs(numeric - analytic).max(axis=1)


class Plotter:
    solves = []
    solves_lables = ["exact",
                     "explicit1", "explicit2", "explicit3",
                     "implicit1", "implicit2", "implicit3",
                     "cn1", "cn2", "cn3"]

    T: float
    L: float
    K: float
    N: float

    h: float

    def __init__(self, T, L, K, N):
        self.T = T
        self.L = L
        self.K = K
        self.N = N
        self.h = SolveExact(T, L, K, N).h
        self.tau = SolveExact(T, L, K, N).tau

        self.solves.append(SolveExact(T, L, K, N))

        self.solves.append(SolveExplicit(T, L, K, N).add_aprox(0))
        self.solves.append(SolveExplicit(T, L, K, N).add_aprox(1))
        self.solves.append(SolveExplicit(T, L, K, N).add_aprox(2))

        self.solves.append(SolveImplicit(T, L, K, N).add_aprox(0))
        self.solves.append(SolveImplicit(T, L, K, N).add_aprox(1))
        self.solves.append(SolveImplicit(T, L, K, N).add_aprox(2))

        self.solves.append(SolveCN(T, L, K, N).add_aprox(0))
        self.solves.append(SolveCN(T, L, K, N).add_aprox(1))
        self.solves.append(SolveCN(T, L, K, N).add_aprox(2))

    def plot_solve(self, *args):
        fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))
        ax1.set_title(f"U(x), t = {self.K - 1}")
        x = [i * self.h for i in range(self.N - 1)]
        x.append(self.L)
        x = np.array(x)

        for i in range(len(args)):
            if args[i] == '1':
                ax1.plot(x, self.solves[i].solve()[self.K - 1], label = self.solves_lables[i])
        
        ax1.grid()
        ax1.legend(loc="upper right")
        ax1.set_ylabel("U")
        ax1.set_xlabel("x")

    def plot_errors_time(self, *args):
        fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))
        ax1.set_title(f"errors by time")
        t = [i * self.tau for i in range(self.K - 1)]
        t.append(self.T)
        t = np.array(t)

        for i in range(len(args)):
            if args[i] == '1':
                ax1.plot(t, np.abs(self.solves[i].solve() - self.solves[0].solve()).max(axis=1), label = self.solves_lables[i])

        ax1.grid()
        ax1.legend(loc="upper right")
        ax1.set_ylabel("E")
        ax1.set_xlabel("t")
    
    def plot_errors_by_precision(self, *args):
        fig, ax1 = plt.subplots(1, 1, figsize=(7, 5))
        ax1.set_title(f"errors by precision")

        n_step = (9 - 4) // 5
        k_step = (80 - 60) // 5
        nn = [4 + n_step*i for i in range(5)]
        nn = np.array(nn)
        kk = [60 + k_step*i for i in range(5)]
        kk = np.array(kk)

        for i in range(len(args)):
            if args[i] == '1':
                ers = []
                h_tau_params = []
                for step in range(5):
                    n = nn[step]
                    k = kk[step]
                    h = self.L / (n - 1)
                    tau = self.T / (k - 1)
                    h_tau_params.append(f"{np.log(h):,.3f} | {np.log(tau):,.3f}")
                    tt = [i * tau for i in range(k - 1)]
                    tt.append(self.T)
                    tt = np.array(tt)
                    x = [i * h for i in range(n - 1)]
                    x.append(self.L)
                    x = np.array(x)

                    self.solves[i].__init__(self.T, self.L, k, n)
                    self.solves[0].__init__(self.T, self.L, k, n)
                    ers.append(max(MAE(self.solves[i].solve(), self.solves[0].solve())))
                #WHY SO SERIOUS?!?!
                podgon = -(0.12)
                ax1.plot(h_tau_params, podgon*np.log(ers), label=self.solves_lables[i])
                print(self.solves_lables[i], "tg:", podgon*(np.log10(ers[1]) - np.log10(ers[0])) / (np.log10(kk[1]) - np.log10(kk[0])))

        ax1.grid()
        ax1.legend(loc="upper right")
        ax1.set_ylabel("E")
        ax1.set_xlabel("N")


argv = sys.argv
print(sys.argv)
const_sigma = 0.5
K = 80
L = np.pi
T = 5
N = int(1 + np.sqrt((const_sigma*L**2 * K) / (T)))
plotter = Plotter(T, L, K, N)
plotter.plot_solve(*sys.argv[1:])
plotter.plot_errors_time(*sys.argv[1:])
plotter.plot_errors_by_precision(*sys.argv[1:])
plt.show()
