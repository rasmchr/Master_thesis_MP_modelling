import numpy as np

def B_funky(z_phys, t, u1, u2, alpha, Da, kappa, C, t_array):
    """Fully standalone B function. Computes M_vals internally.
        Parameters:
        z_phys  : physical z coordinate array (full domain 0 to 150)    
        t       : time
        u1, u2  : velocities for region 1 and 2
        alpha, Da, kappa : physical parameters
        t_array : time array (used for M interpolation)
    """
    # ── B1 ──────────────────────────────────────────────────────────────────────
    def B_func51(z, t_val, u1, alpha, lambda_, kappa, C):
        a  = kappa - lambda_ * u1
        b_ss  = (alpha / a) * C * np.exp(-lambda_ * z) * (1 - np.exp(-a * z / u1))
        b_tr  = (alpha / a) * np.exp(-lambda_ * z) * (1 - np.exp(-a * t_val))
        return np.where(z <= u1 * t_val, b_ss, b_tr)


    # ── B1 evaluated exactly at z = 50 ──────────────────────────────────────────
    def B1_at_50(t_array, u1, alpha, lambda_, kappa, C):
        """Simply evaluate B_func51 at z=50 for each t — no clamp."""
        return B_func51(50.0, t_array, u1, alpha, lambda_, kappa, C)


    # ── Interface ODE — convolution integral ─────────────────────────────────────
    def interface_M(t_array, u1, u2, alpha, lambda_, kappa, C):
        """
        M(t) = ∫₀ᵗ  u1 · B1(50, τ) · exp(−u2·(t−τ))  dτ

        Computed incrementally so we only do O(N) work instead of O(N²).
        Update rule (trapezoidal, step by step):
            M(t_{n+1}) = M(t_n)·exp(−u2·Δt)
                    + Δt·u1·½·[B1(50,t_n)·exp(−u2·Δt) + B1(50,t_{n+1})]
        This is the exact integral assuming B1 is linear between steps.
        """
        t_array = np.asarray(t_array)
        M       = np.zeros_like(t_array)
        B1_vals = B1_at_50(t_array, u1, alpha, lambda_, kappa, C)
        #print("These are the B1 vals: ", B1_vals, "This is the length: ", len(B1_vals))

        for i in range(1, len(t_array)):
            dt_i   = t_array[i] - t_array[i - 1]
            decay  = np.exp(-u2 * dt_i)
            # trapezoidal on the integrand f(τ) = u1·B1(50,τ) between t_{i-1} and t_i
            f_prev = u1 * B1_vals[i - 1]
            f_curr = u1 * B1_vals[i]
            M[i]   = M[i - 1] * decay + dt_i * 0.5 * (f_prev * decay + f_curr)

        return M


    # ── B2 ──────────────────────────────────────────────────────────────────────
    def B_func52(z2_phys, t_val, u2, alpha, lambda_, t_array, kappa, M_vals, C):
        a2 = kappa - lambda_ * u2
        z2 = z2_phys - 50.0        # shifted coordinate (0 at interface)

        # Characteristic travel time from interface to z
        tau_travel = z2 / u2

        # Retarded time (clipped to [0, t_val] and to array bounds)
        tau = np.clip(t_val - tau_travel, 0.0, t_array[-1])

        # M interpolated at the retarded time
        M_tau = np.interp(tau, t_array, M_vals)

        # Branch 1: characteristic has arrived (z2 ≤ u2·t_val)
        b1 = ( (alpha / a2) * C * np.exp(-lambda_ * z2_phys)
            + (M_tau - (alpha / a2) * C * np.exp(-lambda_ * 50.0))
                * np.exp(-kappa * tau_travel) )

        # Branch 2: characteristic has not yet arrived
        b2 = (alpha / a2) * C * np.exp(-lambda_ * z2_phys) * (1.0 - np.exp(-a2 * t_val))

        return np.where(z2 <= u2 * t_val, b1, b2)
    lambda_ = np.sqrt(alpha / Da)
    z1 = z_phys[z_phys <= 50]
    z2 = z_phys[z_phys > 50]
    B1 = B_func51(z1, t, u1, alpha, lambda_, kappa, C)
    M_vals = interface_M(t_array, u1, u2, alpha, lambda_, kappa, C)
    B2 = B_func52(z2, t, u2, alpha, lambda_, t_array, kappa, M_vals, C)
    B = np.concatenate((B1, B2))
    A = C * np.exp(-lambda_ * z_phys)
    full_solutiuon = B + A
    return full_solutiuon
    
