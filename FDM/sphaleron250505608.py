import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import openpyxl
import time
import pandas as pd

# Constants
MW = 1.0  # W-boson mass
MH = 1.556  # Higgs mass
MW2 = MW * MW
MH2 = MH * MH

# Radial domain
RMIN = 0.001
RMAX = 30.0
Nx = 3001
r_vals = np.linspace(RMIN, RMAX, Nx)
dr = (RMAX - RMIN) / (Nx - 1)

# Boundary condition functions
def fA_boundary(r):
    return 1.0 - 2 * (1 - 1.0 / np.cosh(1.150155 * r))


def KK_boundary(r):
    return np.tanh(1.0463536 * r)

def fB_boundary(r):
    return np.sqrt(0.01) * 10.2199 * np.exp(-2.74187 * r) * (r ** 1.6506)

def fC_boundary(r):
    return (np.sqrt(2) / r) * np.sqrt(0.01) * 7.48017 * np.exp(-2.86292 * r) * (r ** 1.67796)

def H_boundary(r):
    return (1 / np.sqrt(2)) * (1 / r) * np.sqrt(0.01) * 4.00337 * np.exp(-2.62178 * r) * (r ** 1.5913)

# Functions for derivatives
def first_derivative(f, i):
    if i == 0:
        return (f[i + 1] - f[i]) / dr
    elif i == Nx - 1:
        return (f[i] - f[i - 1]) / dr
    else:
        return (f[i + 1] - f[i - 1]) / (2.0 * dr)

def second_derivative(f, i):
    return (f[i + 1] - 2.0 * f[i] + f[i - 1]) / (dr * dr)

def second_deriv_r_times_field(f, i):
    return (r_vals[i + 1] * f[i + 1] - 2.0 * r_vals[i] * f[i] + r_vals[i - 1] * f[i - 1]) / (dr * dr)

# Set initial conditions
def set_initial_conditions(fA, vA, fB, vB, fC, vC, H, vH, KK, vK):
    for i in range(Nx):
        rr = r_vals[i]
        
        # fA, KK from boundary expressions across domain
        fA[i] = fA_boundary(rr)
        KK[i] = KK_boundary(rr)

        # fB, fC, H from exponent/r^power profiles
        fB[i] = np.sqrt(0.01) * 10.243 * np.exp(-2.7439 * rr) * (rr ** 1.652)
        fC[i] = (np.sqrt(2) / rr) * np.sqrt(0.01) * 7.52242 * np.exp(-2.8682 * rr) * (rr ** 1.68214)
        H[i] = (1 / np.sqrt(2)) * (1 / rr) * np.sqrt(0.01) * 4.0121 * np.exp(-2.62378 * rr) * (rr ** 1.59307)

        # Zero velocities at t=0
        vA[i] = 0.0
        vB[i] = 0.0
        vC[i] = 0.0
        vH[i] = 0.0
        vK[i] = 0.0

        return fA, vA, fB, vB, fC, vC, KK, vK, H, vH

# Compute the right-hand side of the PDE
def compute_pde_rhs(fA, vA, fB, vB, fC, vC, H, vH, KK, vK):


    # Dirichlet BC for fA, KK at r=RMIN, RMAX
    fA[0] = fA_boundary(r_vals[0])
    KK[0] = KK_boundary(r_vals[0])
    vA[0] = 0.0
    vK[0] = 0.0

    fA[Nx - 1] = fA_boundary(r_vals[Nx - 1])
    KK[Nx - 1] = KK_boundary(r_vals[Nx - 1])
    vA[Nx - 1] = 0.0
    vK[Nx - 1] = 0.0

    # For fB, fC, H => 0 at boundaries, also zero velocities.
    fB[0] = fB_boundary(r_vals[0])
    fC[0] = fC_boundary(r_vals[0])
    H[0] = H_boundary(r_vals[0])
    
    vB[0] = vC[0] = vH[0] = 0.0

    fB[Nx - 1] = fB_boundary(r_vals[Nx - 1])
    fC[Nx - 1] = fC_boundary(r_vals[Nx - 1])
    H[Nx - 1] = H_boundary(r_vals[Nx - 1])
    
    vB[Nx - 1] = vC[Nx - 1] = vH[Nx - 1] = 0.0

    # Initialize k arrays
    dvA = np.zeros(Nx)
    dvB = np.zeros(Nx)
    dvC = np.zeros(Nx)
    dvH = np.zeros(Nx)
    dvK = np.zeros(Nx)

    # PDE interior
    for i in range(1, Nx - 1):
        r = r_vals[i]

        # eq1 => dvA[i] = fA_tt
        fA_rr = second_derivative(fA, i)
        dfB_dr = first_derivative(fB, i)
        dfC_dr = first_derivative(fC, i)

        eq1RHS = (fA_rr
                   - ((fA[i] ** 2 + fB[i] ** 2 - 1.0) * fA[i]) / (r * r)
                   - (MW2 * ((H[i] ** 2 + KK[i] ** 2) * fA[i] + KK[i] ** 2 - H[i] ** 2))
                   - (fA[i] * fC[i] ** 2)
                   + 2.0 * dfB_dr * fC[i]
                   + fB[i] * dfC_dr)

        dvA[i] = eq1RHS

        # eq2 => dvB[i] = fB_tt
        fB_rr = second_derivative(fB, i)
        dfA_dr = first_derivative(fA, i)

        eq2RHS = (fB_rr
                   - ((fA[i] ** 2 + fB[i] ** 2 - 1.0) * fB[i]) / (r * r)
                   - (MW2 * ((H[i] ** 2 + KK[i] ** 2) * fB[i] - 2.0 * KK[i] * H[i]))
                   - (fB[i] * fC[i] ** 2)
                   - 2.0 * dfA_dr * fC[i]
                   - fA[i] * dfC_dr)

        dvB[i] = eq2RHS

        # eq3 => dvC[i] = fC_tt
        dH_dr = first_derivative(H, i)
        dK_dr = first_derivative(KK, i)

        eq3RHS = -((2.0 / (r * r)) * (fA[i] ** 2 + fB[i] ** 2) * fC[i]
                    + (MW2 * (H[i] ** 2 + KK[i] ** 2) * fC[i])
                    + (2.0 * MW2) * (dH_dr * KK[i] - H[i] * dK_dr)
                    + (2.0 / (r * r)) * (dfA_dr * fB[i] - fA[i] * dfB_dr))

        dvC[i] = eq3RHS

        # eq4 => dvH[i] = H_tt
        d2_rH = second_deriv_r_times_field(H, i)
        bracket_4 = ((1.0 / (2.0 * r * r)) * (fA[i] ** 2 + fB[i] ** 2 + 1.0) * H[i]
                      - (1.0 / (r * r)) * (H[i] * fA[i] + KK[i] * fB[i])
                      + 0.5 * MH2 * (H[i] ** 2 + KK[i] ** 2 - 1.0) * H[i]
                      - (1.0 / r) * KK[i] * fC[i]
                      + 0.25 * H[i] * fC[i] ** 2
                      - fC[i] * dK_dr
                      - 0.5 * KK[i] * dfC_dr)

        dvH[i] = (1.0 / r) * d2_rH - bracket_4

        # eq5 => dvK[i] = KK_tt
        d2_rK = second_deriv_r_times_field(KK, i)
        bracket_5 = ((1.0 / (2.0 * r * r)) * (fA[i] ** 2 + fB[i] ** 2 + 1.0) * KK[i]
                      + (1.0 / (r * r)) * (KK[i] * fA[i] - H[i] * fB[i])
                      + 0.5 * MH2 * (H[i] ** 2 + KK[i] ** 2 - 1.0) * KK[i]
                      + (1.0 / r) * H[i] * fC[i]
                      + 0.25 * KK[i] * fC[i] ** 2
                      + fC[i] * dH_dr
                      + 0.5 * H[i] * dfC_dr)

        dvK[i] = (1.0 / r) * d2_rK - bracket_5

    # Update derivatives
    dfA = vA.copy()
    dfB = vB.copy()
    dfC = vC.copy()
    dH = vH.copy()
    dKK = vK.copy()

    # Boundary conditions for derivatives
    for arr in [dfA, dfB, dfC, dH, dKK]:
        arr[0] = 0.0
        arr[Nx - 1] = 0.0

    return dfA, dvA, dfB, dvB, dfC, dvC, dH, dvH, dKK, dvK

# RK4 integration step
def rk4_step(fA, vA, fB, vB, fC, vC, H, vH, KK, vK, dt):

    # Compute k1
    k1fA, k1vA, k1fB, k1vB, k1fC, k1vC, k1H, k1vH, k1K, k1vK = compute_pde_rhs(fA, vA, fB, vB, fC, vC, H, vH, KK, vK)


    fA_temp = np.zeros(Nx)
    vA_temp = np.zeros(Nx)
    fB_temp = np.zeros(Nx)
    vB_temp = np.zeros(Nx)
    fC_temp = np.zeros(Nx)
    vC_temp = np.zeros(Nx)
    H_temp = np.zeros(Nx)
    vH_temp = np.zeros(Nx)
    KK_temp = np.zeros(Nx)
    vK_temp = np.zeros(Nx)
    # Compute k2
    for i in range(Nx):
        fA_temp[i] = fA[i] + 0.5 * dt * k1fA[i]
        vA_temp[i] = vA[i] + 0.5 * dt * k1vA[i]
        fB_temp[i] = fB[i] + 0.5 * dt * k1fB[i]
        vB_temp[i] = vB[i] + 0.5 * dt * k1vB[i]
        fC_temp[i] = fC[i] + 0.5 * dt * k1fC[i]
        vC_temp[i] = vC[i] + 0.5 * dt * k1vC[i]
        H_temp[i] = H[i] + 0.5 * dt * k1H[i]
        vH_temp[i] = vH[i] + 0.5 * dt * k1vH[i]
        KK_temp[i] = KK[i] + 0.5 * dt * k1K[i]
        vK_temp[i] = vK[i] + 0.5 * dt * k1vK[i]

    #print(f"A boundary value: {fA_boundary(r_vals[0])}, Type: {type(fA_boundary(r_vals[0]))}")
    #print(f"A boundary value: {fA_temp[0]}, Type: {type(fA_temp[0])}")

    # Store k2 results
    k2fA, k2vA, k2fB, k2vB, k2fC, k2vC, k2H, k2vH, k2K, k2vK = compute_pde_rhs(fA_temp, vA_temp, fB_temp, vB_temp, fC_temp, vC_temp, H_temp, vH_temp, KK_temp, vK_temp)

    fA_temp = fA + 0.5 * dt * k2fA
    vA_temp = vA + 0.5 * dt * k2vA
    fB_temp = fB + 0.5 * dt * k2fB
    vB_temp = vB + 0.5 * dt * k2vB
    fC_temp = fC + 0.5 * dt * k2fC
    vC_temp = vC + 0.5 * dt * k2vC
    H_temp = H + 0.5 * dt * k2H
    vH_temp = vH + 0.5 * dt * k2vH
    KK_temp = KK + 0.5 * dt * k2K
    vK_temp = vK + 0.5 * dt * k2vK

    k3fA, k3vA, k3fB, k3vB, k3fC, k3vC, k3H, k3vH, k3K, k3vK = compute_pde_rhs(fA_temp, vA_temp, fB_temp, vB_temp, fC_temp, vC_temp, H_temp, vH_temp, KK_temp, vK_temp)

    fA_temp = fA + dt * k3fA
    vA_temp = vA + dt * k3vA
    fB_temp = fB + dt * k3fB
    vB_temp = vB + dt * k3vB
    fC_temp = fC + dt * k3fC
    vC_temp = vC + dt * k3vC
    H_temp = H + dt * k3H
    vH_temp = vH + dt * k3vH
    KK_temp = KK + dt * k3K
    vK_temp = vK + dt * k3vK

    k4fA, k4vA, k4fB, k4vB, k4fC, k4vC, k4H, k4vH, k4K, k4vK = compute_pde_rhs(fA_temp, vA_temp, fB_temp, vB_temp, fC_temp, vC_temp, H_temp, vH_temp, KK_temp, vK_temp)

    fA[:] += (dt / 6.0) * (k1fA + 2.0 * k2fA + 2.0 * k3fA + k4fA)
    vA[:] += (dt / 6.0) * (k1vA + 2.0 * k2vA + 2.0 * k3vA + k4vA)
    fB[:] += (dt / 6.0) * (k1fB + 2.0 * k2fB + 2.0 * k3fB + k4fB)
    vB[:] += (dt / 6.0) * (k1vB + 2.0 * k2vB + 2.0 * k3vB + k4vB)
    fC[:] += (dt / 6.0) * (k1fC + 2.0 * k2fC + 2.0 * k3fC + k4fC)
    vC[:] += (dt / 6.0) * (k1vC + 2.0 * k2vC + 2.0 * k3vC + k4vC)
    H[:] += (dt / 6.0) * (k1H + 2.0 * k2H + 2.0 * k3H + k4H)
    vH[:] += (dt / 6.0) * (k1vH + 2.0 * k2vH + 2.0 * k3vH + k4vH)
    KK[:] += (dt / 6.0) * (k1K + 2.0 * k2K + 2.0 * k3K + k4K)
    vK[:] += (dt / 6.0) * (k1vK + 2.0 * k2vK + 2.0 * k3vK + k4vK)
    # Combine results
    #for i in range(Nx):
    #    fA[i] += (dt / 6.0) * (k1fA[i] + 2 * k2fA[i] + 2 * k3fA[i] + k4fA[i])
    #    vA[i] += (dt / 6.0) * (k1vA[i] + 2 * k2vA[i] + 2 * k3vA[i] + k4vA[i])
        # Repeat for fB, vB, fC, vC, H, vH, KK, vK

    return fA, vA, fB, vB, fC, vC, H, KK

# Main function to initialize and run the simulation
def main():
    # Initialize arrays
    fA = np.zeros(Nx)
    vA = np.zeros(Nx)
    fB = np.zeros(Nx)
    vB = np.zeros(Nx)
    fC = np.zeros(Nx)
    vC = np.zeros(Nx)
    H = np.zeros(Nx)
    vH = np.zeros(Nx)
    KK = np.zeros(Nx)
    vK = np.zeros(Nx)

    # Set initial conditions
    for i in range(Nx):
        rr = r_vals[i]
        fA[i] = fA_boundary(rr)
        KK[i] = KK_boundary(rr)
        fB[i] = fB_boundary(rr)
        fC[i] = fC_boundary(rr)
        H[i] = H_boundary(rr)
        vA[i] = 0.0
        vB[i] = 0.0
        vC[i] = 0.0
        vH[i] = 0.0
        vK[i] = 0.0

    # Initialize arrays
    fA2 = np.zeros(Nx)
    vA2 = np.zeros(Nx)
    fB2 = np.zeros(Nx)
    vB2 = np.zeros(Nx)
    fC2 = np.zeros(Nx)
    vC2 = np.zeros(Nx)
    H2 = np.zeros(Nx)
    vH2 = np.zeros(Nx)
    KK2 = np.zeros(Nx)
    vK2 = np.zeros(Nx)

    # Set initial conditions
    for i in range(Nx):
        rr = r_vals[i]
        fA2[i] = fA_boundary(rr)
        KK2[i] = KK_boundary(rr)
        fB2[i] = fB_boundary(rr)
        fC2[i] = fC_boundary(rr)
        H2[i] = H_boundary(rr)
        vA2[i] = 0.0
        vB2[i] = 0.0
        vC2[i] = 0.0
        vH2[i] = 0.0
        vK2[i] = 0.0

    # Initialize arrays
    fA3 = np.zeros(Nx)
    vA3 = np.zeros(Nx)
    fB3 = np.zeros(Nx)
    vB3 = np.zeros(Nx)
    fC3 = np.zeros(Nx)
    vC3 = np.zeros(Nx)
    H3 = np.zeros(Nx)
    vH3 = np.zeros(Nx)
    KK3 = np.zeros(Nx)
    vK3 = np.zeros(Nx)

    # Set initial conditions
    for i in range(Nx):
        rr = r_vals[i]
        fA3[i] = fA_boundary(rr)
        KK3[i] = KK_boundary(rr)
        fB3[i] = fB_boundary(rr)
        fC3[i] = fC_boundary(rr)
        H3[i] = H_boundary(rr)
        vA3[i] = 0.0
        vB3[i] = 0.0
        vC3[i] = 0.0
        vH3[i] = 0.0
        vK3[i] = 0.0

    # Time integration loop
    dt = 0.0005
    for step in range(4000):  # Adjust the number of steps as needed
        fA, vA, fB, vB, fC, vC, H, KK = rk4_step(fA, vA, fB, vB, fC, vC, H, vH, KK, vK, dt)
        print(step)

    for step in range(10000):  # Adjust the number of steps as needed
        fA2, vA2, fB2, vB2, fC2, vC2, H2, KK2 = rk4_step(fA2, vA2, fB2, vB2, fC2, vC2, H2, vH2, KK2, vK2, dt)
        print(step)

    for step in range(20000):  # Adjust the number of steps as needed
        fA3, vA3, fB3, vB3, fC3, vC3, H3, KK3 = rk4_step(fA3, vA3, fB3, vB3, fC3, vC3, H3, vH3, KK3, vK3, dt)
        print(step)
    
    plt.figure()
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.plot(fA,label=r"$t=2$")
    plt.plot(fA2,label=r"$t=5$")
    plt.plot(fA3,label=r"$t=10$")
    plt.legend(fontsize=15)
    plt.ylabel(r"$f_A(r,t)$",fontsize=18)
    plt.title(r"$Fig.3~\mathrm{of}~\mathrm{2505.05608}$",fontsize=18)
    plt.grid()
    plt.savefig('fA.pdf')
    plt.show()

    plt.figure()
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.plot(fB,label=r"$t=2$")
    plt.plot(fB2,label=r"$t=5$")
    plt.plot(fB3,label=r"$t=10$")
    plt.legend(fontsize=15)
    plt.ylabel(r"$f_B(r,t)$",fontsize=18)
    plt.grid()
    plt.savefig('fB.pdf')
    plt.show()

    plt.figure()
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.plot(fC,label=r"$t=2$")
    plt.plot(fC2,label=r"$t=5$")
    plt.plot(fC3,label=r"$t=10$")
    plt.legend(fontsize=15)
    plt.ylabel(r"$f_C(r,t)$",fontsize=18)
    plt.grid()
    plt.savefig('fC.pdf')
    plt.show()

    plt.figure()
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.plot(H,label=r"$t=2$")
    plt.plot(H2,label=r"$t=5$")
    plt.plot(H3,label=r"$t=10$")
    plt.legend(fontsize=15)
    plt.ylabel(r"$H(r,t)$",fontsize=18)
    plt.grid()
    plt.savefig('H.pdf')
    plt.show()

    plt.figure()
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.plot(KK,label=r"$t=2$")
    plt.plot(KK2,label=r"$t=5$")
    plt.plot(KK3,label=r"$t=10$")
    plt.legend(fontsize=15)
    plt.ylabel(r"$K(r,t)$",fontsize=18)
    plt.grid()
    plt.savefig('K.pdf')
    plt.show()

if __name__ == "__main__":
    main()