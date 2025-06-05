import numpy as np
from scipy.optimize import root
from scipy.integrate import quad
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def sphaleron_htm(rho1, rho2, rho3, rho4, rho5):
    """
    Computes the sphaleron energy in the Higgs Triplet Model.
    
    Parameters:
        rho1, rho2, rho3, rho4, rho5: float
            Model parameters.
            
    Returns:
        energy: float
            The sphaleron energy.
    """
    # Parameters
    beta = 1 + 2 * rho3
    n = 100  # number of grid points
    a = 30  # cutoff

    # Chebyshev grid points
    x = np.cos(np.pi * np.arange(n + 1) / n)
    xx = x[1:-1]  # inner grid points

    # Chebyshev differentiation matrix
    def chebyshev_diff_matrix(n):
        """
        Constructs the Chebyshev differentiation matrix.

        Parameters:
            n: int
            Number of Chebyshev grid points.

        Returns:
            D: np.ndarray
            The Chebyshev differentiation matrix of size (n+1, n+1).
        """
    # Initialize the differentiation matrix
        D = np.zeros((n + 1, n + 1))
    
    # Chebyshev grid points
        x = np.cos(np.pi * np.arange(n + 1) / n)
    
    # Fill the first and last rows and columns
        D[0, 0] = (2 * n**2 + 1) / 6
        D[0, n] = 0.5 * (-1)**n
        D[n, 0] = -0.5 * (-1)**n
        D[n, n] = -((2 * n**2 + 1) / 6)
    
        for j in range(1, n):  # For the first row
            D[0, j] = 2 * (-1)**j / (1 - x[j])
    
        for j in range(1, n):  # For the last row
            D[n, j] = -2 * (-1)**(n + j) / (1 + x[j])
    
        for i in range(1, n):  # For the first column
            D[i, 0] = -0.5 * (-1)**i / (1 - x[i])
    
        for i in range(1, n):  # For the last column
            D[i, n] = 0.5 * (-1)**(n + i) / (1 + x[i])
    
        for j in range(1, n):  # Diagonal elements
            D[j, j] = -x[j] / (2 * (1 - x[j]**2))
    
        for i in range(1, n):  # Off-diagonal elements
            for j in range(1, n):
                if i != j:
                    D[i, j] = (-1)**(i + j) / (x[i] - x[j])
    
        return D

    D = chebyshev_diff_matrix(n)
    D2 = np.dot(D, D)  # second derivative matrix
    D_core = D[1:-1, 1:-1]
    D2_core = D2[1:-1, 1:-1]

    # Initial guesses for the fields
    f_init = np.random.uniform(0.9, 1, n - 1)
    h_init = np.random.uniform(0.9, 1, n - 1)
    hd_init = np.random.uniform(0.9, 1, n - 1)
    var_init = np.concatenate([f_init, h_init, hd_init])

    # Equations of motion
    def equations(var):
        f = var[:n - 1]
        h = var[n - 1:2 * (n - 1)]
        hd = var[2 * (n - 1):]
        fp = np.dot(D_core, f)
        fpp = np.dot(D2_core, f)
        hp = np.dot(D_core, h)
        hpp = np.dot(D2_core, h)
        hdp = np.dot(D_core, hd)
        hdpp = np.dot(D2_core, hd)

        eq1 = 2 * (1 + xx)**2 * fpp - (2 * f + 1 + xx) * (2 * f - 1 + xx) * (2 * f + xx) \
              - (a**2 / (16 * beta)) * (1 + xx)**2 * (2 * f - 1 + xx) * (2 * h + 1 + xx)**2 \
              - ((a**2 * rho3) / (6 * beta)) * (1 + xx)**2 * (2 * f - 1 + xx) * (2 * hd + 1 + xx)**2

        eq2 = (1 + xx)**2 * hpp + (1 + xx) * (2 * hp + 1) - \
              (1 / 4) * (2 * f - 1 + xx)**2 * (2 * h + 1 + xx) + \
              (a**2 / (8 * beta)) * (1 + xx)**2 * (
                  (rho1 - rho2) * (2 * h + 1 + xx) * (4 - (2 * h + 1 + xx)**2) -
                  rho2 * (2 * h + 1 + xx) * ((2 * h + 1 + xx)**2 - 2 * (2 * hd + 1 + xx)) +
                  4 * (rho4 - rho1 + rho2) * (2 * h + 1 + xx) +
                  2 * (np.sqrt(2 * rho2 * rho3 * rho5) - rho2) * (2 * h + 1 + xx) * (2 * hd + 1 + xx) +
                  (rho1 - rho4 - np.sqrt(2 * rho2 * rho3 * rho5)) * (2 * h + 1 + xx) * (2 * hd + 1 + xx)**2)

        eq3 = rho3 * (1 + xx)**2 * hdpp + rho3 * (1 + xx) * (2 * hdp + 1) - \
              ((2 * rho3) / 3) * (2 * f - 1 + xx)**2 * (2 * hd + 1 + xx) + \
              ((a**2 * rho2) / (8 * beta)) * (1 + xx)**2 * ((2 * h + 1 + xx)**2 - 2 * (2 * hd + 1 + xx)) - \
              (a**2 / (8 * beta)) * (1 + xx)**2 * (
                  2 * (2 * rho3 * rho5 - rho2) * (2 * hd + 1 + xx) -
                  (np.sqrt(2 * rho2 * rho3 * rho5) - rho2) * (2 * h + 1 + xx)**2 -
                  (rho1 - rho4 - np.sqrt(2 * rho2 * rho3 * rho5)) * (2 * h + 1 + xx)**2 * (2 * hd + 1 + xx) +
                  (rho1 - rho4 - rho3 * rho5 - np.sqrt(rho2 * rho3 * rho5 / 2)) * (2 * hd + 1 + xx)**3)

        return np.concatenate([eq1, eq2, eq3])

    # Solve the equations
    sol = root(equations, var_init, method='hybr')
    if not sol.success:
        raise ValueError("Failed to solve equations of motion.")

    f_sol = np.concatenate([[0], sol.x[:n - 1], [0]])
    h_sol = np.concatenate([[0], sol.x[n - 1:2 * (n - 1)], [0]])
    hd_sol = np.concatenate([[0], sol.x[2 * (n - 1):], [0]])

    # 定义积分范围
    x_vals = np.linspace(-1, 1, 1000)

    # 将积分函数作用于数组
    def energy_integrand_array(x_vals):
        f0 = np.array([np.interp(xi, xx, f_sol[1:-1]) for xi in x_vals])
        h0 = np.array([np.interp(xi, xx, h_sol[1:-1]) for xi in x_vals])
        hd0 = np.array([np.interp(xi, xx, hd_sol[1:-1]) for xi in x_vals])
        
        f0p = np.gradient(f0, x_vals)
        h0p = np.gradient(h0, x_vals)
        hd0p = np.gradient(hd0, x_vals)

        return (4 / a**2) * f0p**2 + (8 / (a**2 * (1 + x_vals)**2)) * f0**2 * (1 - f0)**2 + \
               (1 / beta) * (1 - f0)**2 * h0**2 + (1 / (2 * beta)) * (1 + x_vals)**2 * h0p**2

    energy = np.trapz(energy_integrand_array(x_vals), x_vals)
    return (x+1)*a, f_sol+(1+x)/2, h_sol+(1+x)/2, hd_sol+(1+x)/2

# Example usage
plt.figure()
plt.subplots_adjust(left=0.15, bottom=0.15)
xx, fsol, hsol, hdsol = sphaleron_htm(0.6, 0.01, 1e-3, 0.59, 0.15)
plt.plot(xx,fsol,label=r"$f(\xi)$")
plt.plot(xx,hsol,label=r"$h(\xi)$")
plt.plot(xx,hdsol,label=r"$h_\Delta(\xi)$")
plt.xlabel(r"$\xi$",fontsize=18)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.xlim(0,20)
plt.legend(fontsize=15)
plt.show()