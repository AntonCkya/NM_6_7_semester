from math import *
from abc import ABC, abstractmethod

import numpy as np
import matplotlib.pyplot as plt

class Diffur:
    a = 2
    b = 2
    c = 4

    alpha_1 = 0
    beta_1 = 1

    alpha_2 = 0
    beta_2 = 1

    alpha_3 = 0
    beta_3 = 1

    alpha_4 = 0
    beta_4 = 1


    @staticmethod
    def f(x, y): return 0

    @staticmethod
    def phi_1(y): return exp(-y) * cos(y)

    @staticmethod
    def phi_2(y): return 0

    @staticmethod
    def phi_3(x): return exp(-x) * cos(x)

    @staticmethod
    def phi_4(x): return 0

    @staticmethod
    def exact(x, y): return exp(-x - y) * cos(x) * cos(y)

    def __init__(self):
        pass


def interpolate_2d(n, m, a):
    up = a[0, 1:-1].any()
    down = a[-1, 1:-1].any()
    left = a[1:-1, 0].any()
    right = a[1:-1, -1].any()

    if up and left:
        for i in range(1, n):
            for j in range(1, m):
                a[i, j] = np.linspace(a[0, j], a[i, 0], i + 1 + j + 1 - 1)[i]
    elif up and right:
        for i in range(1, n):
            for j in range(0, m - 1):
                a[i, j] = np.linspace(a[0, j], a[i, -1], i + 1 + (m - j) - 1)[i]
    elif down and right:
        for i in range(0, n - 1):
            for j in range(0, m - 1):
                a[i, j] = np.linspace(a[-1, j], a[i, -1], (n - i) + (m - j) - 1)[n - i - 1]
    elif down and left:
        for i in range(0, n - 1):
            for j in range(1, m):
                a[i, j] = np.linspace(a[-1, j], a[i, 0], (n - i) + j + 1 - 1)[n - i - 1]
    elif up and down:
        for j in range(m):
            a[1:-1, j] = np.linspace(a[0, j], a[-1, j], n)[1:-1]
    elif left and right:
        for i in range(n):
            a[i, 1:-1] = np.linspace(a[i, 0], a[i, -1], m)[1:-1]
    else:
        print("no inerpolate =(")

    return a

def L2_norm(A, B):
    diff = A - B
    return np.sqrt(np.sum(diff ** 2))

def solve(x_range, y_range, h_x, h_y, method, theta, eps):
    if method == "exact":
        exact, iterss = solve_exact(x_range, y_range, h_x, h_y)
        return exact, iterss
    d = Diffur()
    obhod = "rd"

    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)
    n = len(x)
    m = len(y)
    res = np.zeros((n, m))

    for i in range(m):
        if d.alpha_1 == 0:
            res[0][i] = 1 / d.beta_1 * d.phi_1(y[i])
        if d.alpha_2 == 0:
            res[-1][i] = 1 / d.beta_2 * d.phi_2(y[i])
    for j in range(n):
        if d.alpha_3 == 0:
            res[j][0] = 1 / d.beta_3 * d.phi_3(x[j])
        if d.alpha_4 == 0:
            res[j][-1] = 1 / d.beta_4 * d.phi_4(x[j])
    
    n, m = res.shape
    interpolate_2d(n, m, res)

    iters = 1
    while True:
        res_prev = res
        res = np.zeros((n, m))
        for i in range(m):
            if d.alpha_1 == 0:
                res[0][i] = 1 / d.beta_1 * d.phi_1(y[i])
            if d.alpha_2 == 0:
                res[-1][i] = 1 / d.beta_2 * d.phi_2(y[i])
        for j in range(n):
            if d.alpha_3 == 0:
                res[j][0] = 1 / d.beta_3 * d.phi_3(x[j])
            if d.alpha_4 == 0:
                res[j][-1] = 1 / d.beta_4 * d.phi_4(x[j])

        obhods = ['rd', 'ld', 'lu', 'ru']
        obhods_info = [
            (1, n - 2, 1, m - 2),
            (1, n - 2, m - 2, 1),
            (n - 2, 1, m - 2, 1),
            (n - 2, 1, 1, m - 2),
        ]  # ОО - обратный обход
        obhod_idx = obhods.index(obhod)
        x_start, x_end, y_start, y_end = obhods_info[obhod_idx]
        x_dir = int((x_end - x_start) > 0) - int((x_end - x_start) < 0)
        y_dir = int((y_end - y_start) > 0) - int((y_end - y_start) < 0)

        uij_coeff = (d.c - 2 / h_x ** 2 - 2 / h_y ** 2)
        for i in range(x_start, x_end + x_dir, x_dir):
            for j in range(y_start, y_end + y_dir, y_dir):
                if method == "simple":
                    part_d2u_dx2 = 1 / h_x ** 2 * (res_prev[i + x_dir][j] + res_prev[i - x_dir][j])
                    part_d2u_dy2 = 1 / h_y ** 2 * (res_prev[i][j + y_dir] + res_prev[i][j - y_dir])
                    du_dx = d.a / (2 * h_x) * (res_prev[i + x_dir][j] - res_prev[i - x_dir][j])
                    du_dy = d.b / (2 * h_y) * (res_prev[i][j + y_dir] - res_prev[i][j - y_dir])
                    res[i][j] = 1 / uij_coeff * (d.f(x[i], y[j]) - (part_d2u_dx2 + part_d2u_dy2 + du_dx + du_dy))
                elif method == "zeidel" or method == "relaxation":
                    part_d2u_dx2 = 1 / h_x ** 2 * (res_prev[i + x_dir][j] + res[i - x_dir][j])
                    part_d2u_dy2 = 1 / h_y ** 2 * (res_prev[i][j + y_dir] + res[i][j - y_dir])
                    du_dx = d.a / (2 * h_x) * (res_prev[i + x_dir][j] - res[i - x_dir][j])
                    du_dy = d.b / (2 * h_y) * (res_prev[i][j + y_dir] - res[i][j - y_dir])
                    res[i][j] = (
                            theta * (1 / uij_coeff * (
                            d.f(x[i], y[j]) - (part_d2u_dx2 + part_d2u_dy2 + du_dx + du_dy))) +
                            (1 - theta) * res_prev[i][j]
                    )
        
        for i in range(1, m - 1):
            if d.alpha_1 != 0:
                u0j_coef = 2 * h_x * d.beta_1 - 3 * d.alpha_1
                res[0][i] = 1 / u0j_coef * (
                        2 * h_x * d.phi_1(y[i]) - d.alpha_1 * (4 * res[1][i] - res[2][i]))
            if d.alpha_2 != 0:
                unj_coeff = 2 * h_x * d.beta_2 + 3 * d.alpha_2
                res[-1][i] = 1 / unj_coeff * (
                        2 * h_x * d.phi_2(y[i]) + d.alpha_2 * (4 * res[-2][i] - res[-3][i]))
        
        for j in range(1, n - 1):
            if d.alpha_3 != 0:
                ui0_coef = 2 * h_y * d.beta_3 - 3 * d.alpha_3
                res[j][0] = 1 / ui0_coef * (
                        2 * h_y * d.phi_3(x[j]) - d.alpha_3 * (4 * res[j][1] - res[j][2]))
            if d.alpha_4 != 0:
                uim_coeff = 2 * h_y * d.beta_4 + 3 * d.alpha_4
                res[j][-1] = 1 / uim_coeff * (
                        2 * h_y * d.phi_4(x[j]) + d.alpha_4 * (4 * res[j][-2] - res[j][-3]))
        
        if L2_norm(res, res_prev) < eps:
            break

        iters += 1
    
    return res, iters

def solve_exact(x_range, y_range, h_x, h_y):
    d = Diffur()
    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)

    res = np.zeros((len(x), len(y)))
    for idx in range(len(x)):
        for idy in range(len(y)):
            res[idx][idy] = d.exact(x[idx], y[idy])

    return res, 1

def MAE(A, B):
    return abs(A - B).mean()

def maxAE(A, B):
    return abs(A - B).max()

def plot_results(x_end, y_end, x_steps, y_steps, theta, eps):
    x_range = (0, x_end)
    y_range = (0, y_end)
    h_x = x_end / x_steps
    h_y = y_end / y_steps

    exact, _ = solve(x_range, y_range, h_x, h_y, "exact", theta, eps)

    simple, simple_iters = solve(x_range, y_range, h_x, h_y, "simple", theta, eps)
    print("Метод простых итераций:")
    print(f'iters: {simple_iters}')
    print(f'Mean Abs Err: {MAE(simple, exact)}')
    print()

    zeidel, zeidel_iters = solve(x_range, y_range, h_x, h_y, "zeidel", 1, eps)
    print("Метод Зейделя:")
    print(f'iters: {zeidel_iters}')
    print(f'Mean Abs Err: {MAE(zeidel, exact)}')
    print()

    relaxation, relaxation_iters = solve(x_range, y_range, h_x, h_y, "relaxation", theta, eps)
    print("Метод верхней релаксации:")
    print(f'iters: {relaxation_iters}')
    print(f'Mean Abs Err: {MAE(relaxation, exact)}')
    print()

    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)

    fig, axs = plt.subplots(1, 2, figsize=(9, 3))

    lines_x = []
    lines_y = []
    solutions = {
        "exact": exact,
        "simple": simple,
        "zeidel": zeidel,
        "relaxation": relaxation,
    }
    for method_name, solution in solutions.items():
        line_x, = axs[0].plot(y, solution[1, :], label=method_name)
        lines_x.append(line_x)

        line_y, = axs[1].plot(x, solution[:, 1], label=method_name)
        lines_y.append(line_y)
    
    axs[0].set_title('u(x, y)')
    axs[0].set_xlabel('y')
    axs[0].set_ylabel('u(x, y)')
    axs[0].legend()

    axs[1].set_title('u(x, y)')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('u(x, y)')
    axs[1].legend()

    fig2, axs2 = plt.subplots(1, 2, figsize=(9, 3))
    h_x_s = np.arange(x_steps // 2, x_steps, 2)
    h_y_s = np.arange(y_steps // 2, y_steps, 2)

    MAES_x_simple = []
    MAES_x_zeidel = []
    MAES_x_relaxation = []
    for i in h_x_s:
        h_x = x_end / i
        h_y = y_end / y_steps
        MAES_x_simple.append(MAE(solve(x_range, y_range, h_x, h_y, "simple", theta, eps)[0], solve(x_range, y_range, h_x, h_y, "exact", theta, eps)[0]))
        MAES_x_zeidel.append(MAE(solve(x_range, y_range, h_x, h_y, "zeidel", 1, eps)[0], solve(x_range, y_range, h_x, h_y, "exact", theta, eps)[0]))
        MAES_x_relaxation.append(MAE(solve(x_range, y_range, h_x, h_y, "relaxation", theta, eps)[0], solve(x_range, y_range, h_x, h_y, "exact", theta, eps)[0]))
    MAES_y_simple = []
    MAES_y_zeidel = []
    MAES_y_relaxation = [] 
    for j in h_y_s:
        h_x = x_end / x_steps
        h_y = y_end / j
        MAES_y_simple.append(MAE(solve(x_range, y_range, h_x, h_y, "simple", theta, eps)[0], solve(x_range, y_range, h_x, h_y, "exact", theta, eps)[0]))
        MAES_y_zeidel.append(MAE(solve(x_range, y_range, h_x, h_y, "zeidel", 1, eps)[0], solve(x_range, y_range, h_x, h_y, "exact", theta, eps)[0]))
        MAES_y_relaxation.append(MAE(solve(x_range, y_range, h_x, h_y, "relaxation", theta, eps)[0], solve(x_range, y_range, h_x, h_y, "exact", theta, eps)[0]))

    print(f"интервалы h_x: {h_x_s}")
    print(f"интервалы h_y: {h_y_s}")

    solutions_MAES = {
        "simple": (MAES_x_simple, MAES_y_simple),
        "zeidel": (MAES_x_zeidel, MAES_y_zeidel),
        "relaxation": (MAES_x_relaxation , MAES_y_relaxation),
    }
    for method_name, solution in solutions_MAES.items():
        line_x, = axs2[0].plot(h_x_s, solution[0], label=method_name)
        lines_x.append(line_x)

        line_y, = axs2[1].plot(h_y_s, solution[1], label=method_name)
        lines_y.append(line_y)

    axs2[0].set_title('Ошибка по h_x')
    axs2[0].set_xlabel('h_x')
    axs2[0].set_ylabel('MAE')
    axs2[0].legend()

    axs2[1].set_title('Ошибка по h_y')
    axs2[1].set_xlabel('h_y')
    axs2[1].set_ylabel('MAE')
    axs2[1].legend()

    plt.show()


plot_results(pi / 2, pi / 2, 33, 33, 1.5, 1e-1)
