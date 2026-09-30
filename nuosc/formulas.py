"""
Three-flavour oscillation probabilities
=======================================

All functions return the full probability matrix

        P[alpha, beta, ...] = P(nu_alpha -> nu_beta)

with alpha, beta = 0 (e), 1 (mu), 2 (tau). The trailing dimensions have the
broadcast shape of the energy E [GeV] and baseline L [km] you pass in, so
you can use scalars, 1D arrays or 2D meshgrids (e.g. for oscillograms).

Three ways of computing P are provided:

1. prob_vacuum         -> textbook analytic formula in vacuum
2. prob_matter_nufast  -> NuFast algorithm for constant matter density
                          (Denton & Parke, arXiv:2405.02400) - fast & accurate
3. prob_matter_exact   -> brute-force numerical diagonalisation of the
                          Hamiltonian; slow, used only as a cross-check.
"""

import numpy as np

from .constants import K_PHASE, K_MATTER


# --------------------------------------------------------------------------- #
#  Building blocks                                                            #
# --------------------------------------------------------------------------- #
def pmns_matrix(s2_12, s2_13, s2_23, delta_cp):
    """
    PMNS matrix in the standard (PDG) parametrisation, U = U23 U13 U12.

    Parameters
    ----------
    s2_12, s2_13, s2_23 : sin^2 of the mixing angles
    delta_cp            : CP phase [degrees]

    Returns
    -------
    U : (3, 3) complex ndarray, rows = flavour (e, mu, tau), cols = mass (1,2,3)
    """
    s12, s13, s23 = np.sqrt(s2_12), np.sqrt(s2_13), np.sqrt(s2_23)
    c12, c13, c23 = np.sqrt(1 - s2_12), np.sqrt(1 - s2_13), np.sqrt(1 - s2_23)
    eid = np.exp(1j * np.deg2rad(delta_cp))

    U23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]], dtype=complex)
    U13 = np.array([[c13, 0, s13 / eid], [0, 1, 0], [-s13 * eid, 0, c13]],
                   dtype=complex)
    U12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]], dtype=complex)
    return U23 @ U13 @ U12


def matter_potential(E, rho, Ye=0.5, antineutrino=False):
    """
    Matter term A = 2 sqrt(2) G_F N_e E  [eV^2]   (E in GeV, rho in g/cm^3).
    The sign flips for antineutrinos.
    """
    A = K_MATTER * Ye * rho * np.asarray(E, dtype=float)
    return -A if antineutrino else A


def _broadcast(E, L):
    E, L = np.broadcast_arrays(np.asarray(E, dtype=float),
                               np.asarray(L, dtype=float))
    return E, L


# --------------------------------------------------------------------------- #
#  1. Vacuum                                                                  #
# --------------------------------------------------------------------------- #
def prob_vacuum(E, L, s2_12, s2_13, s2_23, delta_cp, dm2_21, dm2_31,
                antineutrino=False):
    r"""
    Vacuum oscillation probability for all 9 channels.

    P(nu_a -> nu_b) = delta_ab
        - 4 sum_{i>j} Re[U*_ai U_bi U_aj U*_bj] sin^2(Delta_ij)
        + 2 sum_{i>j} Im[U*_ai U_bi U_aj U*_bj] sin(2 Delta_ij)

    with Delta_ij = dm2_ij L / (4E) = 1.267 dm2_ij[eV^2] L[km] / E[GeV].
    For antineutrinos U -> U* (i.e. delta_cp -> -delta_cp).
    """
    E, L = _broadcast(E, L)
    U = pmns_matrix(s2_12, s2_13, s2_23, delta_cp)
    if antineutrino:
        U = U.conj()

    m2 = np.array([0.0, dm2_21, dm2_31])          # m_i^2 - m_1^2
    P = np.zeros((3, 3) + E.shape)

    for a in range(3):
        for b in range(3):
            P[a, b] = 1.0 if a == b else 0.0
            for i in range(3):
                for j in range(i):
                    # quartic product of mixing matrix elements
                    W = U[a, i].conj() * U[b, i] * U[a, j] * U[b, j].conj()
                    Delta = K_PHASE * (m2[i] - m2[j]) * L / E
                    P[a, b] += (-4.0 * W.real * np.sin(Delta) ** 2
                                + 2.0 * W.imag * np.sin(2.0 * Delta))
    return P


# --------------------------------------------------------------------------- #
#  2. Matter, constant density: NuFast                                        #
# --------------------------------------------------------------------------- #
def prob_matter_nufast(E, L, s2_12, s2_13, s2_23, delta_cp, dm2_21, dm2_31,
                       rho, Ye=0.5, antineutrino=False, N_Newton=0):
    """
    Oscillation probabilities in constant-density matter with NuFast
    (P. B. Denton and S. J. Parke, arXiv:2405.02400).

    Idea of the algorithm:
      * get the largest eigenvalue lambda3 of the Hamiltonian from the
        DMP approximation (optionally refined with N_Newton Newton steps);
      * the other two eigenvalues follow from the characteristic polynomial;
      * the |U_ai|^2 in matter follow from the "Eigenvector-Eigenvalue
        identity" (a.k.a. "Rosetta"), and the Jarlskog invariant from the
        Naumov-Harrison-Scott identity;
      * P_ee, P_mumu and P_mue are computed explicitly, the other six
        channels follow from unitarity.

    N_Newton = 0 already gives relative precision ~1e-4 or better at
    accelerator/atmospheric energies; use 1-2 for extreme precision.
    """
    E, L = _broadcast(E, L)

    # Antineutrinos: A -> -A and delta -> -delta
    delta = np.deg2rad(delta_cp)
    if antineutrino:
        delta = -delta
    Amatter = matter_potential(E, rho, Ye, antineutrino)   # [eV^2]
    Lover4E = K_PHASE * L / E                               # [1/eV^2]

    # ------------------------------------------------------------------ #
    # Useful functions of the vacuum parameters                          #
    # ------------------------------------------------------------------ #
    c13sq = 1.0 - s2_13

    Ue2sq = c13sq * s2_12
    Ue3sq = s2_13

    Um3sq = c13sq * s2_23
    # Um2sq and Ut2sq are temporary here, properly defined below
    Ut2sq = s2_13 * s2_12 * s2_23
    Um2sq = (1.0 - s2_12) * (1.0 - s2_23)

    Jrr = np.sqrt(Um2sq * Ut2sq)
    sind, cosd = np.sin(delta), np.cos(delta)

    Um2sq = Um2sq + Ut2sq - 2.0 * Jrr * cosd
    Jmatter = 8.0 * Jrr * c13sq * sind
    Dmsqee = dm2_31 - s2_12 * dm2_21

    # Coefficients of the characteristic polynomial
    #   lambda^3 - A lambda^2 + B lambda - C = 0
    A = dm2_21 + dm2_31                       # temporary
    See = A - dm2_21 * Ue2sq - dm2_31 * Ue3sq
    Tmm = dm2_21 * dm2_31                     # temporary
    Tee = Tmm * (1.0 - Ue3sq - Ue2sq)
    C = Amatter * Tee
    A = A + Amatter

    # ------------------------------------------------------------------ #
    # lambda3 from lambda+ of DMP                                        #
    # ------------------------------------------------------------------ #
    xmat = Amatter / Dmsqee
    tmp = 1.0 - xmat
    lambda3 = dm2_31 + 0.5 * Dmsqee * (xmat - 1.0 +
                                       np.sqrt(tmp * tmp + 4.0 * s2_13 * xmat))

    # Newton iterations (B only needed here)
    B = Tmm + Amatter * See
    for _ in range(N_Newton):
        lambda3 = ((lambda3 * lambda3 * (lambda3 + lambda3 - A) + C)
                   / (lambda3 * (2.0 * (lambda3 - A) + lambda3) + B))

    # ------------------------------------------------------------------ #
    # Eigenvalue differences                                             #
    # ------------------------------------------------------------------ #
    tmp = A - lambda3
    Dlambda21 = np.sqrt(tmp * tmp - 4.0 * C / lambda3)
    lambda2 = 0.5 * (A - lambda3 + Dlambda21)
    Dlambda32 = lambda3 - lambda2
    Dlambda31 = Dlambda32 + Dlambda21

    # ------------------------------------------------------------------ #
    # |U_ai|^2 in matter from the eigenvector-eigenvalue identity        #
    # ------------------------------------------------------------------ #
    PiDlambdaInv = 1.0 / (Dlambda31 * Dlambda32 * Dlambda21)
    Xp3 = PiDlambdaInv * Dlambda21
    Xp2 = -PiDlambdaInv * Dlambda31

    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2

    Smm = A - dm2_21 * Um2sq - dm2_31 * Um3sq
    Tmm = Tmm * (1.0 - Um3sq - Um2sq) + Amatter * (See + Smm - A)

    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2

    # Jarlskog in matter (Naumov-Harrison-Scott identity)
    Jmatter = Jmatter * dm2_21 * dm2_31 * (dm2_31 - dm2_21) * PiDlambdaInv

    # Remaining elements from unitarity
    Ue1sq = 1.0 - Ue3sq - Ue2sq
    Um1sq = 1.0 - Um3sq - Um2sq
    Ut3sq = 1.0 - Um3sq - Ue3sq
    Ut2sq = 1.0 - Um2sq - Ue2sq
    Ut1sq = 1.0 - Um1sq - Ue1sq

    # ------------------------------------------------------------------ #
    # Kinematic terms                                                    #
    # ------------------------------------------------------------------ #
    sinD21 = np.sin(Dlambda21 * Lover4E)
    sinD31 = np.sin(Dlambda31 * Lover4E)
    sinD32 = np.sin(Dlambda32 * Lover4E)

    triple_sin = sinD21 * sinD31 * sinD32
    sinsqD21_2 = 2.0 * sinD21 * sinD21
    sinsqD31_2 = 2.0 * sinD31 * sinD31
    sinsqD32_2 = 2.0 * sinD32 * sinD32

    # ------------------------------------------------------------------ #
    # The three independent probabilities (CP-conserving / violating)    #
    # ------------------------------------------------------------------ #
    Pme_CPC = ((Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
               + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
               + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2)
    Pme_CPV = -Jmatter * triple_sin

    Pmm = 1.0 - 2.0 * (Um2sq * Um1sq * sinsqD21_2
                       + Um3sq * Um1sq * sinsqD31_2
                       + Um3sq * Um2sq * sinsqD32_2)
    Pee = 1.0 - 2.0 * (Ue2sq * Ue1sq * sinsqD21_2
                       + Ue3sq * Ue1sq * sinsqD31_2
                       + Ue3sq * Ue2sq * sinsqD32_2)

    # ------------------------------------------------------------------ #
    # All nine channels (rows and columns of P sum to 1)                 #
    # ------------------------------------------------------------------ #
    P = np.empty((3, 3) + E.shape)
    P[0, 0] = Pee
    P[0, 1] = Pme_CPC - Pme_CPV            # e  -> mu
    P[0, 2] = 1.0 - P[0, 0] - P[0, 1]      # e  -> tau
    P[1, 0] = Pme_CPC + Pme_CPV            # mu -> e
    P[1, 1] = Pmm
    P[1, 2] = 1.0 - P[1, 0] - P[1, 1]      # mu -> tau
    P[2, 0] = 1.0 - P[0, 0] - P[1, 0]      # tau -> e
    P[2, 1] = 1.0 - P[0, 1] - P[1, 1]      # tau -> mu
    P[2, 2] = 1.0 - P[0, 2] - P[1, 2]      # tau -> tau
    return P


# --------------------------------------------------------------------------- #
#  3. Matter, constant density: exact numerical solution (cross-check)       #
# --------------------------------------------------------------------------- #
def prob_matter_exact(E, L, s2_12, s2_13, s2_23, delta_cp, dm2_21, dm2_31,
                      rho, Ye=0.5, antineutrino=False):
    """
    Solve i d/dx nu = H nu exactly for constant density by diagonalising

        2E H = U diag(0, dm2_21, dm2_31) U^dagger + diag(A, 0, 0)

    S = exp(-i H L) = V exp(-i lambda L / 2E) V^dagger,
    P(a -> b) = |S_ba|^2.
    """
    E, L = _broadcast(E, L)
    shape = E.shape
    Ef, Lf = E.ravel(), L.ravel()

    U = pmns_matrix(s2_12, s2_13, s2_23, delta_cp)
    if antineutrino:
        U = U.conj()
    M2 = U @ np.diag([0.0, dm2_21, dm2_31]) @ U.conj().T      # [eV^2]

    H = np.repeat(M2[None, :, :], Ef.size, axis=0)            # (N, 3, 3)
    H[:, 0, 0] += matter_potential(Ef, rho, Ye, antineutrino)

    lam, V = np.linalg.eigh(H)                                # (N,3), (N,3,3)
    # phase lambda L / (2E) = 2 * K_PHASE * lambda[eV^2] L[km] / E[GeV]
    phase = np.exp(-2j * K_PHASE * lam * (Lf / Ef)[:, None])
    S = np.einsum("nik,nk,njk->nij", V, phase, V.conj())      # S[n, b, a]

    P = np.abs(S) ** 2                                        # P[n, b, a]
    P = np.transpose(P, (2, 1, 0))                            # P[a, b, n]
    return P.reshape((3, 3) + shape)
