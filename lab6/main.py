from abc import ABC, abstractmethod
import sys

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
    def psi_1(x): return np.exp(2 * x)

    @staticmethod
    def psi_2(x): return 0

    @staticmethod
    def d2_psi_1(x): return 4 * np.exp(2 * x)

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
        return tau**2 / h**2

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
        return np.exp(2 * x) * np.cos(t)

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
            return lambda k, u: u[k][1] / (1 + 2 * self.h)
        elif self.aprox == 1:
            return lambda k, u: (4 * u[k][1] - u[k][2]) / (3 + 4 * self.h)
        elif self.aprox == 2:
            return lambda k, u: self.sigma * (2 * u[k - 1][1] - (2 + 4 * self.h) * u[k - 1][0]) + (2 - 5 * self.tau**2) * u[k - 1][0] - u[k - 2][0]

    def l_aprox(self):
        if self.aprox == 0:
            return lambda k, u: u[k][-2] / (1 - 2 * self.h)
        elif self.aprox == 1:
            return lambda k, u: (4 * u[k][-2] - u[k][-3]) / (3 - 4 * self.h)
        elif self.aprox == 2:
            return lambda k, u: self.sigma * (2 * u[k - 1][-2] + (4 * self.h - 2) * u[k - 1][-1]) + (2 - 5 * self.tau**2) * u[k - 1][-1] - u[k - 2][-1]

    def solve(self):
        u = np.zeros((self.K, self.N))

        u_zero = []
        for i in [i * self.h for i in range(self.N)]:
            u_zero.append(self.D.psi_1(i))
        u[0] = np.array(u_zero)

        u_one = []
        for i in [i * self.h for i in range(self.N)]:
            u_one.append(self.D.psi_1(i) + self.tau * self.D.psi_2(i) + self.tau**2 * self.D.d2_psi_1(i) / 2)
        u[1] = np.array(u_one)

        for i in range(2, self.K):
            for j in range(1, self.N - 1):
                u[i][j] = (self.sigma * 
                            (u[i - 1][j - 1] - 2 * u[i - 1][j] + u[i - 1][j + 1]) -
                            (5 * self.tau**2 * u[i - 1][j]) + 
                            (2 * u[i - 1][j]) - u[i - 2][j])
            u[i][0] = self.zero_aprox()(i, u)
            u[i][-1] = self.l_aprox()(i, u)

        return u


class SolveImplicit(AbstractSolver):
    aprox = 0

    def add_aprox(self, num: int):
        self.aprox = num
        return self

    def zero_aprox(self):
        if self.aprox == 0:
            return lambda k, u: (0, (1 + 2 * self.h), -1, 0)
        elif self.aprox == 1:
            return lambda k, u: (0,
                                -(2 + 4 * self.h),
                                -(5 * self.h**2 + 1 / self.sigma - 2),
                                (-2*u[k - 1][1] + u[k - 2][1])/self.sigma
            )
        elif self.aprox == 2:
            return lambda k, u: (0,
                                -(2 + 5 * self.h**2 + 4 * self.h + 1 / self.sigma),
                                2,
                                (-2 * u[k - 1][0] + u[k - 2][0]) / self.sigma
            )

    def l_aprox(self):
        if self.aprox == 0:
            return lambda k, u: (-1, (1 - 2 * self.h), 0, 0)
        elif self.aprox == 1:
            return lambda k, u: (-(5 * self.h**2 + 1 / self.sigma - 2),
                                -(2 - 4 * self.h),
                                0,
                                (-2 * u[k - 1][-2] + u[k - 2][-2]) / self.sigma
            )
        elif self.aprox == 2:
            return lambda k, u: (2,
                                -(2 + 5 * self.h**2 - 4 * self.h + 1 / self.sigma),
                                0,
                                (-2 * u[k - 1][-1] + u[k - 2][-1]) / self.sigma
            )

    def solve(self):
        u = np.zeros((self.K, self.N))

        u_zero = []
        for i in [i * self.h for i in range(self.N)]:
            u_zero.append(self.D.psi_1(i))
        u[0] = np.array(u_zero)

        u_one = []
        for i in [i * self.h for i in range(self.N)]:
            u_one.append(self.D.psi_1(i) + self.tau * self.D.psi_2(i) + self.tau**2 * self.D.d2_psi_1(i) / 2)
        u[1] = np.array(u_one)

        for i in range(2, self.K):
            a = np.zeros(self.N)
            b = np.zeros(self.N)
            c = np.zeros(self.N)
            d = np.zeros(self.N)
            for j in range(1, self.N - 1):
                a[j] = 1
                b[j] = -(2 + 5 * self.h**2 + 1 / self.sigma)
                c[j] = 1
                d[j] = (u[i - 2][j] - 2*u[i - 1][j]) / self.sigma
            a[0], b[0], c[0], d[0] = self.zero_aprox()(i, u)
            a[-1], b[-1], c[-1], d[-1] = self.l_aprox()(i, u)

            u[i] = tridiagonal_solve(a, b, c, d)
        
        return u
   

def MAE(numeric, analytic):
    return np.abs(numeric - analytic).max(axis=1)


class Plotter:
    solves = []
    solves_lables = ["exact",
                     "explicit1", "explicit2", "explicit3",
                     "implicit1", "implicit2", "implicit3"]

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

        count = 3 # self.N - 5
        const_sigma = 1
        h = np.array(list(map(int, np.linspace(start=15, stop=self.N, num=count))))
        tau = []
        for i in h:
            tau.append(int(self.T * i * (const_sigma ** 0.5) * (1 / (self.L))))
        tau = np.array(tau)

        print(tau)
        print(h)

        for i in range(len(args)):
            if args[i] == '1':
                err = []
                for x in zip(tau, h):
                    self.solves[i].__init__(self.T, self.L, x[0], x[1])
                    self.solves[0].__init__(self.T, self.L, x[0], x[1])
                    err.append(max(MAE(self.solves[i].solve(), self.solves[0].solve())))
                ax1.plot(np.log(h), np.log10(err), label = self.solves_lables[i])
                print(self.solves_lables[i] + " tg =", (np.log10(err[-1]) - np.log10(err[0])) / (np.log10(h[-1]) - np.log10(h[0])))

        ax1.grid()
        ax1.legend(loc="upper right")
        ax1.set_ylabel("E")
        ax1.set_xlabel("N")


argv = sys.argv
print(sys.argv)
"""
const_sigma = 0.5
K = 80
L = np.pi
T = 5
N = int(1 + np.sqrt((const_sigma*L**2 * K) / (T)))
"""
K = 100
L = 1
T = 1
N = 30
plotter = Plotter(T, L, K, N)
plotter.plot_solve(*sys.argv[1:])
plotter.plot_errors_time(*sys.argv[1:])
plotter.plot_errors_by_precision(*sys.argv[1:])
plt.show()
