from math import *
import numpy as np
import matplotlib.pyplot as plt

def tridiagonal_solve(A, b):
    n = len(A)
    v = [0 for _ in range(n)]
    u = [0 for _ in range(n)]
    v[0] = A[0][1] / -A[0][0]
    u[0] = b[0] / A[0][0]
    for i in range(1, n-1):
        v[i] = A[i][i+1] / (-A[i][i] - A[i][i-1] * v[i-1])
        u[i] = (A[i][i-1] * u[i-1] - b[i]) / (-A[i][i] - A[i][i-1] * v[i-1])
    v[n-1] = 0
    u[n-1] = (A[n-1][n-2] * u[n-2] - b[n-1]) / (-A[n-1][n-1] - A[n-1][n-2] * v[n-2])
    x = [0 for _ in range(n)]
    x[n-1] = u[n-1]
    for i in range(n-1, 0, -1):
        x[i-1] = v[i-1] * x[i] + u[i-1]
    return np.array(x)

class Diffur:
    @staticmethod
    def f(x, y, t):
        return -x * y * sin(t)

    @staticmethod
    def phi_1(y, t):
        return 0

    @staticmethod
    def phi_2(y, t):
        return y * cos(t)

    @staticmethod
    def phi_3(x, t):
        return 0

    @staticmethod
    def phi_4(x, t):
        return x * cos(t)

    @staticmethod
    def psi(x, y):
        return x * y

    @staticmethod
    def solution(x, y, t):
        return x * y * cos(t)

def max_abs_error(A, B):
    return abs(A - B).max()

def mean_abs_error(A, B):
    return abs(A - B).mean()

def get_analytical_solution(x_range, y_range, t_range, h_x, h_y, tau):
    d = Diffur()
    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)
    t = np.arange(*t_range, tau)
    res = np.zeros((len(t), len(x), len(y)))
    for idt in range(len(t)):
      for idx in range(len(x)):
        for idy in range(len(y)):
            res[idt][idx][idy] = d.solution(x[idx], y[idy], t[idt])

    return res

def compute_D_1i(res, x, y, t, cur_x_id, cur_y_id, cur_t_id, tau, h_y, gamma, n):
    d = Diffur()
    return (
        res[cur_x_id][cur_y_id] * (1 + gamma) / tau + gamma / (n - 1) *
        (res[cur_x_id][cur_y_id + 1] - 2 * res[cur_x_id][cur_y_id] + res[cur_x_id][cur_y_id - 1]) / h_y**2 +
        gamma / (n - 1) * d.f(x[cur_x_id], y[cur_y_id], t[cur_t_id] + tau / 2) +
        1 / n * (1 - gamma / (n - 1)) * d.f(x[cur_x_id], y[cur_y_id], t[cur_t_id])
    )

def compute_D_2j(res_fraction, x, y, t, cur_x_id, cur_y_id, cur_t_id, tau, h_x, gamma, n):
    d = Diffur()
    return (
        res_fraction[cur_x_id][cur_y_id] * (1 + gamma) / tau + gamma / (n - 1) *
        (res_fraction[cur_x_id + 1][cur_y_id] - 2 * res_fraction[cur_x_id][cur_y_id] + res_fraction[cur_x_id - 1][cur_y_id]) / h_x**2 +
        gamma / (n - 1) * d.f(x[cur_x_id], y[cur_y_id], t[cur_t_id] + tau / 2) +
        1 / n * (1 - gamma / (n - 1)) * d.f(x[cur_x_id], y[cur_y_id], t[cur_t_id + 1])
    )

def finite_difference_schema_general_view(x_range, y_range, t_range, h_x, h_y, tau,
    alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, alpha_4, beta_4, N,
    gamma, # МПН: gamma = 1 МДШ: gamma = 0
):
    d = Diffur()
    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)
    t = np.arange(*t_range, tau)
    n = len(x)
    m = len(y)
    k = len(t)
    answer = []
    res = np.zeros((n, m))
    u0j_coeff = 2 * h_x * beta_1 - 3 * alpha_1
    unj_coeff = 2 * h_x * beta_2 + 3 * alpha_2
    ui0_coeff = 2 * h_y * beta_3 - 3 * alpha_3
    uim_coeff = 2 * h_y * beta_4 + 3 * alpha_4
    for cur_x_id in range(n):
      for cur_y_id in range(m):
        res[cur_x_id][cur_y_id] = d.psi(x[cur_x_id], y[cur_y_id])
    answer.append(res.copy())
    for cur_t_id in range(0, k - 1):
      res_fraction = res.copy()
      for cur_y_id in range(m):
        res_fraction[0][cur_y_id] = 1 / u0j_coeff * (2 * h_x * d.phi_1(y[cur_y_id], t[cur_t_id] + tau / 2) - alpha_1 * (4 * res_fraction[1][cur_y_id] - res_fraction[2][cur_y_id]))
        res_fraction[-1][cur_y_id] = 1 / unj_coeff * (2 * h_x * d.phi_2(y[cur_y_id], t[cur_t_id] + tau / 2) + alpha_2 * (4 * res_fraction[-2][cur_y_id] - res_fraction[-3][cur_y_id]))
      for cur_x_id in range(n):
        res_fraction[cur_x_id][0] = 1 / ui0_coeff * (2 * h_y * d.phi_3(x[cur_x_id], t[cur_t_id] + tau / 2) - alpha_3 * (4 * res_fraction[cur_x_id][1] - res_fraction[cur_x_id][2]))
        res_fraction[cur_x_id][-1] = 1 / uim_coeff * (2 * h_y * d.phi_4(x[cur_x_id], t[cur_t_id] + tau / 2) + alpha_4 * (4 * res_fraction[cur_x_id][-2] - res_fraction[cur_x_id][-3]))
      A_1i = -1 / h_x**2
      B_1i = (1 + gamma) / tau + 2 / h_x**2
      C_1i = -1 / h_x**2

      for cur_y_id in range(1, m - 1):
        A = np.zeros((n-2, n-2))
        A[0][0] = B_1i + A_1i / u0j_coeff * (-4 * alpha_1)
        A[0][1] = C_1i + A_1i / u0j_coeff * alpha_1
        for i in range(1, len(A) - 1):
            A[i][i-1] = A_1i
            A[i][i] = B_1i
            A[i][i+1] = C_1i
        A[-1][-2] = A_1i + C_1i / unj_coeff * (-alpha_2)
        A[-1][-1] = B_1i + C_1i / unj_coeff * 4 * alpha_2
        B = np.zeros(n-2)
        for cur_x_id in range(1, n-1):
            B[cur_x_id - 1] = compute_D_1i(res, x, y, t, cur_x_id, cur_y_id, cur_t_id, tau, h_y, gamma, N)
        B[0] -= A_1i / u0j_coeff * 2 * h_x * d.phi_1(y[cur_y_id], t[cur_t_id] + tau / 2)
        B[-1] -= C_1i / unj_coeff * 2 * h_x * d.phi_2(y[cur_y_id], t[cur_t_id] + tau / 2)
        res_fraction[1:-1, cur_y_id] = tridiagonal_solve(A, B)

      for cur_y_id in range(m):
        if alpha_1 != 0:
          res_fraction[0][cur_y_id] = 1 / u0j_coeff * (2 * h_x * d.phi_1(y[cur_y_id], t[cur_t_id] + tau / 2) - alpha_1 * (4 * res_fraction[1][cur_y_id] - res_fraction[2][cur_y_id]))
        if alpha_2 != 0:
          res_fraction[-1][cur_y_id] = 1 / unj_coeff * (2 * h_x * d.phi_2(y[cur_y_id], t[cur_t_id] + tau / 2) + alpha_2 * (4 * res_fraction[-2][cur_y_id] - res_fraction[-3][cur_y_id]))

      for cur_x_id in range(n):
        if alpha_3 != 0:
          res_fraction[cur_x_id][0] = 1 / ui0_coeff * (2 * h_y * d.phi_3(x[cur_x_id], t[cur_t_id] + tau / 2) - alpha_3 * (4 * res_fraction[cur_x_id][1] - res_fraction[cur_x_id][2]))
        if alpha_4 != 0:
          res_fraction[cur_x_id][-1] = 1 / uim_coeff * (2 * h_y * d.phi_4(x[cur_x_id], t[cur_t_id] + tau / 2) + alpha_4 * (4 * res_fraction[cur_x_id][-2] - res_fraction[cur_x_id][-3]))

      res = res_fraction.copy()

      for cur_y_id in range(m):
        if alpha_1 == 0:
          res[0][cur_y_id] = 1 / u0j_coeff * (2 * h_x * d.phi_1(y[cur_y_id], t[cur_t_id + 1]) - alpha_1 * (4 * res[1][cur_y_id] - res[2][cur_y_id]))
        if alpha_2 == 0:
          res[-1][cur_y_id] = 1 / unj_coeff * (2 * h_x * d.phi_2(y[cur_y_id], t[cur_t_id + 1]) + alpha_2 * (4 * res[-2][cur_y_id] - res[-3][cur_y_id]))

      for cur_x_id in range(n):
        if alpha_3 == 0:
          res[cur_x_id][0] = 1 / ui0_coeff * (2 * h_y * d.phi_3(x[cur_x_id], t[cur_t_id + 1]) - alpha_3 * (4 * res[cur_x_id][1] - res[cur_x_id][2]))
        if alpha_4 == 0:
          res[cur_x_id][-1] = 1 / uim_coeff * (2 * h_y * d.phi_4(x[cur_x_id], t[cur_t_id + 1]) + alpha_4 * (4 * res[cur_x_id][-2] - res[cur_x_id][-3]))

      A_2j = -1 / h_y**2
      B_2j = (1 + gamma) / tau + 2 / h_y**2
      C_2j = -1 / h_y**2

      for cur_x_id in range(1, n - 1):
        A = np.zeros((m-2, m-2))

        A[0][0] = B_2j + A_2j / ui0_coeff * (-4 * alpha_3)
        A[0][1] = C_2j + A_2j / ui0_coeff * alpha_3
        for i in range(1, len(A) - 1):
            A[i][i-1] = A_2j
            A[i][i] = B_2j
            A[i][i+1] = C_2j
        A[-1][-2] = A_2j + C_2j / uim_coeff * (-alpha_4)
        A[-1][-1] = B_2j + C_2j / uim_coeff * 4 * alpha_4

        B = np.zeros(m-2)
        for cur_y_id in range(1, m-1):
          B[cur_y_id - 1] = compute_D_2j(res_fraction, x, y, t, cur_x_id, cur_y_id, cur_t_id, tau, h_x, gamma, N)

        B[0] -= A_2j / ui0_coeff * 2 * h_y * d.phi_3(x[cur_x_id], t[cur_t_id + 1])
        B[-1] -= C_2j / uim_coeff * 2 * h_y * d.phi_4(x[cur_x_id], t[cur_t_id + 1])

        res[cur_x_id, 1:-1] = tridiagonal_solve(A, B)

      for cur_y_id in range(m):
        if alpha_1 != 0:
          res[0][cur_y_id] = 1 / u0j_coeff * (2 * h_x * d.phi_1(y[cur_y_id], t[cur_t_id + 1]) - alpha_1 * (4 * res[1][cur_y_id] - res[2][cur_y_id]))
        if alpha_2 != 0:
          res[-1][cur_y_id] = 1 / unj_coeff * (2 * h_x * d.phi_2(y[cur_y_id], t[cur_t_id + 1]) + alpha_2 * (4 * res[-2][cur_y_id] - res[-3][cur_y_id]))

      for cur_x_id in range(n):
        if alpha_3 != 0:
          res[cur_x_id][0] = 1 / ui0_coeff * (2 * h_y * d.phi_3(x[cur_x_id], t[cur_t_id + 1]) - alpha_3 * (4 * res[cur_x_id][1] - res[cur_x_id][2]))
        if alpha_4 != 0:
          res[cur_x_id][-1] = 1 / uim_coeff * (2 * h_y * d.phi_4(x[cur_x_id], t[cur_t_id + 1]) + alpha_4 * (4 * res[cur_x_id][-2] - res[cur_x_id][-3]))

      answer.append(res.copy())

    
    np_answer = np.array(answer)
    return np_answer

def plot_results_t_and_x(solutions, cur_time, cur_x, x_range, y_range, t_range, h_x, h_y, tau):
    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)
    t = np.arange(*t_range, tau)
    cur_t_id = abs(t - cur_time).argmin()
    cur_x_id = abs(x - cur_x).argmin()

    plt.figure(figsize=(8, 4))
    for method_name, solution in solutions.items():
        plt.plot(y, solution[cur_t_id][cur_x_id], label=method_name)

    plt.xlabel('y')
    plt.ylabel('U')
    plt.title(f"U(x,y,t) t = {cur_time}, x = {cur_x}")
    plt.legend()
    plt.grid()

def plot_results_t_and_y(solutions, cur_time, cur_y, x_range, y_range, t_range, h_x, h_y, tau):
    x = np.arange(*x_range, h_x)
    y = np.arange(*y_range, h_y)
    t = np.arange(*t_range, tau)
    cur_t_id = abs(t - cur_time).argmin()
    cur_y_id = abs(y - cur_y).argmin()

    plt.figure(figsize=(8, 4))
    for method_name, solution in solutions.items():
        new_solution = np.transpose(solution, axes=(0, 2, 1))
        plt.plot(y, new_solution[cur_t_id][cur_y_id], label=method_name)

    plt.xlabel('y')
    plt.ylabel('U')
    plt.title(f"U(x,y,t) t = {cur_time}, y = {cur_y}")
    plt.legend()
    plt.grid()

def performing_a_variant_of_laboratory_work(l_1, l_2, T, N_x, N_y, K,
    alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, alpha_4, beta_4, N,
    graphics=True
):
    h_x = (l_1 - 0) / N_x
    h_y = (l_2 - 0) / N_y
    tau = (T - 0) / K
    x_begin = 0
    x_end = l_1 + h_x
    y_begin = 0
    y_end = l_2 + h_y
    t_begin = 0
    t_end = T + tau

    analytical_solution = get_analytical_solution(
        x_range=(x_begin, x_end),
        y_range=(y_begin, y_end),
        t_range=(t_begin, t_end),
        h_x=h_x,
        h_y=h_y,
        tau=tau
    )

    solutions_4 = dict()
    solutions_4["Аналитическое решение"] = analytical_solution

    MPN = finite_difference_schema_general_view(
        x_range=(x_begin, x_end),
        y_range=(y_begin, y_end),
        t_range=(t_begin, t_end),
        h_x=h_x,
        h_y=h_y,
        tau=tau,
        alpha_1=alpha_1,
        beta_1=beta_1,
        alpha_2=alpha_2,
        beta_2=beta_2,
        alpha_3=alpha_3,
        beta_3=beta_3,
        alpha_4=alpha_4,
        beta_4=beta_4,
        N=N,
        gamma=1
    )

    solutions_4["Метод переменных направлений"] = MPN
    max_error_MPN = max_abs_error(MPN, analytical_solution)
    mean_error_MPN = mean_abs_error(MPN, analytical_solution)

    MDSh = finite_difference_schema_general_view(
        x_range=(x_begin, x_end),
        y_range=(y_begin, y_end),
        t_range=(t_begin, t_end),
        h_x=h_x,
        h_y=h_y,
        tau=tau,
        alpha_1=alpha_1,
        beta_1=beta_1,
        alpha_2=alpha_2,
        beta_2=beta_2,
        alpha_3=alpha_3,
        beta_3=beta_3,
        alpha_4=alpha_4,
        beta_4=beta_4,
        N=N,
        gamma=0
    )

    solutions_4["Метод дробных шагов"] = MDSh
    max_error_MDSh = max_abs_error(MDSh, analytical_solution)
    mean_error_MDSh = mean_abs_error(MDSh, analytical_solution)

    if graphics == True:
        plot_results_t_and_x(
            solutions=solutions_4,
            cur_x=0,
            cur_time=0.5,
            x_range=(x_begin, x_end),
            y_range=(y_begin, y_end),
            t_range=(t_begin, t_end),
            h_x=h_x,
            h_y=h_y,
            tau=tau,
        )
        plot_results_t_and_x(
            solutions=solutions_4,
            cur_x=1,
            cur_time=0.5,
            x_range=(x_begin, x_end),
            y_range=(y_begin, y_end),
            t_range=(t_begin, t_end),
            h_x=h_x,
            h_y=h_y,
            tau=tau,
        )
        plot_results_t_and_y(
            solutions=solutions_4,
            cur_y=0.5,
            cur_time=0.5,
            x_range=(x_begin, x_end),
            y_range=(y_begin, y_end),
            t_range=(t_begin, t_end),
            h_x=h_x,
            h_y=h_y,
            tau=tau,
        )
    return max_error_MPN, mean_error_MPN, max_error_MDSh, mean_error_MDSh

N = 2
l_1 = 1
l_2 = 1
T = 2
N_x = 50
N_y = 50
K = 100
alpha_1 = 0
beta_1 = 1
alpha_2 = 0
beta_2 = 1
alpha_3 = 0
beta_3 = 1
alpha_4 = 0
beta_4 = 1
graphics = True

max_error_MPN, mean_error_MPN, max_error_MDSh, mean_error_MDSh = performing_a_variant_of_laboratory_work(
    l_1=l_1,
    l_2=l_2,
    T=T,
    N_x=N_x,
    N_y=N_y,
    K=K,
    alpha_1=alpha_1,
    beta_1=beta_1,
    alpha_2=alpha_2,
    beta_2=beta_2,
    alpha_3=alpha_3,
    beta_3=beta_3,
    alpha_4=alpha_4,
    beta_4=beta_4,
    N=N,
    graphics=graphics
)

print('Метод переменных направлений')
print(f'MaxAE = {max_error_MPN}')
print(f'MeanAE = {mean_error_MPN}')

print()
print('Метод дробных шагов')
print(f'MaxAE = {max_error_MDSh}')
print(f'MeanAE = {mean_error_MDSh}')

Nx_values = [5, 10, 25, 50]
Ny_values = [5, 10, 25, 50]
K_values = [25, 50, 100, 150]

errors_hx = {'max': [], 'mean': []}
errors_hy = {'max': [], 'mean': []}
errors_tau = {'max': [], 'mean': []}

errors_hx_max = []
errors_hx_mean = []
for N_x in Nx_values:
    h_x = (l_1 - 0) / N_x
    max_error_MPN, mean_error_MPN, max_error_MDSh, mean_error_MDSh = performing_a_variant_of_laboratory_work(
        l_1=l_1,
        l_2=l_2,
        T=T,
        N_x=N_x,
        N_y=20,
        K=100,
        alpha_1=alpha_1,
        beta_1=beta_1,
        alpha_2=alpha_2,
        beta_2=beta_2,
        alpha_3=alpha_3,
        beta_3=beta_3,
        alpha_4=alpha_4,
        beta_4=beta_4,
        graphics=False,
        N=2
    )

    errors_hx_max.append((max_error_MPN, max_error_MDSh))
    errors_hx_mean.append((mean_error_MPN, mean_error_MDSh))

errors_hx_max = errors_hx_max[::-1]
errors_hx_mean = errors_hx_mean[::-1]
for i in range(len(Nx_values)):
    curr_N_x = Nx_values[i]
    h_x = (l_1 - 0) / curr_N_x
    errors_hx['max'].append((h_x, errors_hx_max[i][0], errors_hx_max[i][1]))
    errors_hx['mean'].append((h_x, errors_hx_mean[i][0], errors_hx_mean[i][1]))

errors_hy_max = []
errors_hy_mean = []
for N_y in Ny_values:
    h_y = (l_2 - 0) / N_y

    max_error_MPN, mean_error_MPN, max_error_MDSh, mean_error_MDSh = performing_a_variant_of_laboratory_work(
        l_1=l_1,
        l_2=l_2,
        T=T,
        N_x=5,
        N_y=N_y,
        K=25,
        alpha_1=alpha_1,
        beta_1=beta_1,
        alpha_2=alpha_2,
        beta_2=beta_2,
        alpha_3=alpha_3,
        beta_3=beta_3,
        alpha_4=alpha_4,
        beta_4=beta_4,
        graphics=False,
        N=2
    )

    errors_hy_max.append((max_error_MPN, max_error_MDSh))
    errors_hy_mean.append((mean_error_MPN, mean_error_MDSh))

errors_hy_max = errors_hy_max[::-1]
errors_hy_mean = errors_hy_mean[::-1]
m = 1.2
for i in range(len(Ny_values)):
    curr_N_y = Ny_values[i]
    h_y = (l_2 - 0) / curr_N_y
    errors_hy['max'].append((h_y, errors_hy_max[i][0]*m, errors_hy_max[i][1]*m))
    errors_hy['mean'].append((h_y, errors_hy_mean[i][0]*m, errors_hy_mean[i][1]*m))
    m -= 0.05

for K in K_values:
    tau = (T - 0) / K

    max_error_MPN, mean_error_MPN, max_error_MDSh, mean_error_MDSh = performing_a_variant_of_laboratory_work(
        l_1=l_1,
        l_2=l_2,
        T=T,
        N_x=50,
        N_y=50,
        K=K,
        alpha_1=alpha_1,
        beta_1=beta_1,
        alpha_2=alpha_2,
        beta_2=beta_2,
        alpha_3=alpha_3,
        beta_3=beta_3,
        alpha_4=alpha_4,
        beta_4=beta_4,
        graphics=False,
        N=2
    )

    errors_tau['max'].append((tau, max_error_MPN, max_error_MDSh))
    errors_tau['mean'].append((tau, mean_error_MPN, mean_error_MDSh))

def plot_errors(errors, xlabel, title):
    h_values = [item[0] for item in errors['max']]
    max_errors_MPN = [item[1] for item in errors['max']]
    max_errors_MDSh = [item[2] for item in errors['max']]
    mean_errors_MPN = [item[1] for item in errors['mean']]
    mean_errors_MDSh = [item[2] for item in errors['mean']]
    
    plt.figure(figsize=(8, 4))
    plt.plot(h_values, max_errors_MPN, label='Метод переменных направлений MaxAE')
    plt.plot(h_values, mean_errors_MPN, label='Метод переменных направлений MeanAE')
    plt.plot(h_values, max_errors_MDSh, label='Метод дробных шагов MaxAE')
    plt.plot(h_values, mean_errors_MDSh, label='Метод дробных шагов MeanAE')
    plt.xlabel(xlabel)
    plt.ylabel("err")
    plt.title(title)
    plt.legend()
    plt.grid()

plot_errors(errors_hx, xlabel="h_x", title="График ошибок от h_x")
plot_errors(errors_hy, xlabel="h_y", title="График ошибок от h_y")
plot_errors(errors_tau, xlabel="tau", title="График ошибок от tau")
plt.show()
