import numpy as np
import matplotlib.pyplot as plt


def dimer_equilibria(C2, KD, f, C_l):
    """
    Model function for normalized apparent diffusion time (tau_app) and
    species fractions for a monomer-dimer equilibrium.

    Parameters
    ----------
    C2 : float or array
        Total protein concentration expressed in terms of dimer (C0 = 2*C2).
    KD : float
        Equilibrium dissociation constant (nM) for monomer <-> dimer.
    f : float
        Labelling efficiency 
    C_l : float
        Concentration of labeled protein in terms of dimer (nM).
    Cl : float
        Concentration of labeled protein in terms of monomer (nM).

    Returns
    -------
    tau_app : array
        Apparent diffusion time, normalized to the dimer diffusion time
        (r2 = tauD1/tauD2 = 0.79).
    alpha1 : array
        Fraction of total protein present as monomer.
    alpha2 : array
        Fraction of total protein present as dimer.
    """
    C = 2 * C2
    r2 = 0.79
    alpha1 = (1 / (C * 4)) * (-KD + np.sqrt(8 * C * KD + KD**2))
    alpha2 = 1 - alpha1
    Cl = 2 * C_l
    lf    = (Cl * f) / C 
    c = 1 - 2 * (alpha1 / (1 + (1 - alpha1) * lf))
    tau_app = 0.5 * (c * (1 - r2) + np.sqrt(c**2 * (1 - r2)**2 + 4 * r2))
    return tau_app, alpha1, alpha2

def trimer_equilibria(C3, KDE, f, C_l):
    """
    Model function for normalized apparent diffusion time (tau_app) and
    species fractions for a monomer-trimer equilibrium.

    Parameters
    ----------
    C3 : float or array
        Total protein concentration expressed in terms of trimer (C0 = 3*C3).
    KDE : float
        Effective dissociation constant (nM); the trimer KD is computed
        internally as KD = (3/4) * KDE**2.
    f : float
        Labelling efficiency 
    Cl : float
        Concentration of labeled protein in terms of monomer (nM).
    C_l : float
        Concentration of labeled protein in terms of trimer (nM).

    Returns
    -------
    tau_app : array
        Apparent diffusion time, normalized to the trimer diffusion time
        (r3 = tauD1/tauD3 = 0.69).
    alpha1 : array
        Fraction of total protein present as monomer.
    alpha3 : array
        Fraction of total protein present as trimer.

    """
    r3 = 0.69
    C = 3 * C3
    KD = (3/4)*KDE**2
    a1 = 9 * C**4 * KD + np.sqrt(81 * C**8 * KD**2 + 4 * C**6 * KD**3)
    f1 = (2**(1/3) * KD) / (a1**(1/3))
    f2 = (a1**(1/3)) / (2**(1/3) * C**2)
    alpha1 = (1/3) * (-f1 + f2)
    alpha3 = 1 - alpha1
    Cl = 2 * C_l
    lf    = (Cl * f) / C 
    c = 1 - 2 * (alpha1 / (1 + (1 - alpha1) * (2 * lf)))
    tau_app = 0.5 * (c * (1 - r3) + np.sqrt(c**2 * (1 - r3)**2 + 4 * r3))
    return tau_app, alpha1, alpha3

def tetramer_equilibria(C4, KD1, KD2, f, C_l):
    """
   Model function for normalized apparent diffusion time (tau_app) and
   species fractions for a monomer-dimer-tetramer equilibrium system.

   Parameters
   ----------
   C4 : float or array
       Total protein concentration expressed in terms of tetramer 
       (C0 = 4*C4).
   KD1 : float
       Dissociation constant (nM) for dimer <-> tetramer.
   KD2 : float
       Dissociation constant (nM) for monomer <-> dimer.
   f : float
       Labelling efficiency 
   Cl : float
       Concentration of labeled protein in terms of monomer (nM).
   C_l : float
       Concentration of labeled protein in terms of tetramer (nM).

   Returns
   -------
   tau_app : array
       Apparent diffusion time, normalized to the tetramer diffusion time
       (r1 = tauD1/tauD4 = 0.79**2, 
        r2 = tauD2/tauD4 = 0.79).
   alpha1 : array
       Fraction of total protein present as monomer.
   alpha2 : array
       Fraction of total protein present as dimer.
   alpha4 : array
       Fraction of total protein present as tetramer.
   """
    C = np.atleast_1d(np.array(4 * C4, dtype=float))
    r1 = 0.79 * 0.79
    r2 = 0.79
    c1 = np.zeros(len(C))
    for i, Ctot in enumerate(C):
        pp = [4 / (KD2**2 * KD1), 0, 2 / KD2, 1, -Ctot]
        roots = np.roots(pp)
        valid = roots[(roots.real > 0) & (np.abs(roots.imag) < 1e-10)].real
        c1[i] = valid[0]
    c2 = c1**2 / KD2
    c4 = c2**2 / KD1
    alpha1 = c1      / C
    alpha2 = 2 * c2  / C
    alpha4 = 4 * c4  / C
    Cl = 4 * C_l
    lf    = (Cl* f) / C 
    denom = 1 + alpha2 * lf + alpha4 * 3 * lf
    a1 = alpha1                  / denom
    a2 = alpha2 * (1 +     lf)  / denom
    a3 = alpha4 * (1 + 3 * lf)  / denom
    b1 = r1 + r2 + 1
    b2 = a1 * r1 + a2 * r2 + a3
    b3 = r1 * r2 + r1 + r2
    b4 = a1 * r1 * (r2 + 1) + a2 * r2 * (r1 + 1) + a3 * (r2 + r1)
    b5 = r1 * r2
    tau_app = np.zeros(len(C))
    for i in range(len(C)):
        pp = [1, (b1 - 2 * b2[i]), (b3 - 2 * b4[i]), -b5]
        roots = np.roots(pp)
        valid = roots[(roots.real > 0) & (np.abs(roots.imag) < 1e-10)].real
        tau_app[i] = np.max(valid)
    return tau_app, alpha1, alpha2, alpha4


# Plot tau_app ratio vs concentration for each model

C_test = np.logspace(np.log10(1), np.log10(1e3), 100)

#dimer-monomer
plt.figure()
plt.semilogx(C_test, dimer_equilibria(C_test, 100, 0.5, 1)[0], label='dimer')
plt.xlabel('Conc (nM)')
plt.ylabel('tau_D')
plt.legend()
plt.show()

#trimer-monomer
plt.figure()
plt.semilogx(C_test, trimer_equilibria(C_test, 100, 0.5, 1)[0], label='trimer')
plt.xlabel('Conc (nM)')
plt.ylabel('tau_D')
plt.legend()
plt.show()

#tetramer-dimer-monomer
plt.figure()
plt.semilogx(C_test, tetramer_equilibria(C_test, 100, 10, 0.5, 1)[0], label='tetramer')
plt.xlabel('Conc (nM)')
plt.ylabel('tau_D')
plt.legend()
plt.show()