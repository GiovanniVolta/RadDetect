import warnings
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc

try:
    from iminuit import Minuit
except ImportError:
    Minuit = None

# 222Rn half-life in hours and decay constant in 1/h
RN222_HALF_LIFE_H = 3.823 * 24.0
RN222_LAMBDA_H = np.log(2) / RN222_HALF_LIFE_H

# Unit conversion constants
FROM_NM2_H_TO_CM2_SEC = 1e-14 / 3600.0
FROM_CM2_SEC_TO_NM2_H = 3600.0 / 1e-14


def calculate_A_volumetric(activity_bq, mu_x, sigma_x, lambda_decay, S_nm2):
    """
    Calculates the 3D volumetric amplitude (A) of the initial Gaussian distribution.

    Because mu_x > 0 (the peak is in the vacuum), the sample only captures the
    extreme left-hand tail of the distribution (x <= 0). This function ensures
    absolute mass conservation of that specific tail region.

    Parameters:
        activity_bq (float): Initial measured alpha activity in the sample [Bq = decays/s].
        mu_x (float): Implantation depth parameter [nm] (MUST be positive).
        sigma_x (float): Depth straggle standard deviation [nm].
        lambda_decay (float): Radioactive decay constant [1/h].
        S_nm2 (float): Exposed outgassing surface area of the sample [nm^2].

    Returns:
        float: Normalized volumetric peak amplitude (A) of the Gaussian profile [atoms / nm^3].
    """
    # Total initial atoms N0: [decays / h] / [1 / h] = [atoms]
    N0 = (activity_bq * 3600.0) / lambda_decay
    # Area under the 1D tail (from x = -inf up to the boundary at x = 0).
    # Since mu_x is positive, erfc(mu_x / ...) inherently calculates 
    # just the tail volume.
    integral_1D = (
        sigma_x * np.sqrt(np.pi / 2.0) * erfc(mu_x / (np.sqrt(2.0) * sigma_x))
    )
    # Amplitude A: [atoms] / ([nm^2] * [nm]) = [atoms / nm^3]
    A = N0 / (S_nm2 * integral_1D)
    
    return A


def calculate_activity_from_A(A, mu_x, sigma_x, lambda_decay, S_nm2):
    """
    Calculates the initial alpha activity from the volumetric amplitude (A).

    Inverts the normalization integral by integrating the Gaussian tail across
    the physical domain (x <= 0) and multiplying by the decay constant.

    Parameters:
        A (float): Volumetric amplitude [atoms / nm^3].
        mu_x (float): Implantation depth parameter [nm].
        sigma_x (float): Depth straggle [nm].
        lambda_decay (float): Radioactive decay constant [1/h].
        S_nm2 (float): Exposed surface area [nm^2].

    Returns:
        float: Total initial alpha activity of the sample [Bq = decays/s].
    """
    # 1. Re-calculate the integral of the Gaussian tail
    integral_1D = sigma_x * np.sqrt(np.pi / 2) * erfc(mu_x / (np.sqrt(2) * sigma_x))
    
    # 2. Calculate Total Atoms (N0) = A * Area * Integral
    N0 = A * S_nm2 * integral_1D
    
    # 3. Convert N0 atoms to Activity [Bq]
    activity_bq = (N0 * lambda_decay) / 3600.0
    
    return activity_bq


def C_x_t(x, t, D, mu_x, sigma_x, A, lambda_decay=0.0):
    """
    Calculates the concentration profile C(x,t) inside the substrate.

    Accounts for 1D Fickian diffusion, radioactive decay, and an absorbing
    boundary at x = 0 (vacuum chamber interface) solved via the Method of Images.

    Coordinate System:
        x <= 0 : Inside the metal sample (physical domain).
        x = 0  : The substrate surface boundary.
        x > 0  : The vacuum chamber (masked to exactly 0.0).

    Parameters:
        x (float or array-like): Spatial coordinate(s) [nm].
        t (float or array-like): Elapsed diffusion time [h].
        D (float): Diffusion coefficient [nm^2 / h].
        mu_x (float): Initial Gaussian center [nm].
        sigma_x (float): Initial Gaussian width [nm].
        A (float): Volumetric amplitude [atoms / nm^3].
        lambda_decay (float, optional): Decay constant [1/h]. Defaults to 0.0 (no decay).

    Returns:
        float or numpy.ndarray: Volumetric concentration C(x, t) at the specified
        position(s) and time [atoms / nm^3].
    """
    x_arr = np.asarray(x, dtype=float)
    t_val = float(t) if np.isscalar(t) else np.asarray(t, dtype=float)

    if np.isscalar(t_val) and t_val == 0:
        C_initial = (
            A
            * np.exp(-((x_arr - mu_x) ** 2) / (2.0 * sigma_x**2))
            * np.exp(-lambda_decay * t_val)
        )
        return np.where(x_arr >= 0, 0.0, C_initial)

    Sigma = np.sqrt(sigma_x**2 + 2.0 * D * t_val)
    sqrt_Dt = np.sqrt(D * t_val)

    # TERM 1: The REAL source
    term1_exp = np.exp(-((x_arr - mu_x) ** 2) / (2.0 * Sigma**2))
    term1_erfc = erfc(
        (x_arr * sigma_x**2 + 2.0 * mu_x * D * t_val)
        / (2.0 * sigma_x * Sigma * sqrt_Dt)
    )

    # TERM 2: The IMAGE sink
    term2_exp = np.exp(-((x_arr + mu_x) ** 2) / (2.0 * Sigma**2))
    term2_erfc = erfc(
        (-x_arr * sigma_x**2 + 2.0 * mu_x * D * t_val)
        / (2.0 * sigma_x * Sigma * sqrt_Dt)
    )

    prefactor = np.exp(-lambda_decay * t_val) * (A * sigma_x) / (2.0 * Sigma)

    C = prefactor * (term1_exp * term1_erfc - term2_exp * term2_erfc)
    return np.where(x_arr >= 0, 0.0, C)


def J_diff_t_reconstructed(t, D, mu_x, sigma_x, A):
    """
    Calculates the pure diffusion outgassing flux crossing the boundary at x = 0.

    Evaluates the spatial derivative J = D * (dC/dx)|_{x=0} analytically.
    Radioactive decay is excluded to allow for mass-balance integration.

    Parameters:
        t (float or array-like): Elapsed time point(s) [h].
        D (float): Diffusion coefficient [nm^2 / h].
        mu_x (float): Implantation depth parameter [nm].
        sigma_x (float): Depth straggle [nm].
        A (float): Volumetric amplitude [atoms / nm^3].

    Returns:
        float or numpy.ndarray: Instantaneous surface outgassing flux [atoms / (nm^2 * h)].
    """
    t_arr = np.asarray(t, dtype=float)
    is_scalar = t_arr.ndim == 0
    t_flat = t_arr.ravel()

    result = np.zeros_like(t_flat, dtype=float)
    mask = t_flat > 0

    if np.any(mask):
        t_v = t_flat[mask]
        Sigma2 = sigma_x**2 + 2.0 * D * t_v

        Vc = (sigma_x**2 * 2.0 * D * t_v) / Sigma2
        c = (2.0 * mu_x * D * t_v) / Sigma2
        
        prefactor = A / np.sqrt(4.0 * np.pi * D * t_v**3)
        global_exp = np.exp(-(mu_x**2) / (2.0 * Sigma2))
        
        term_A = Vc * np.exp(-(c**2) / (2.0 * Vc))
        term_B = c * np.sqrt(np.pi * Vc / 2.0) * erfc(c / np.sqrt(2.0 * Vc))
        
        result[mask] = prefactor * global_exp * (term_A - term_B)

    res_reshaped = result.reshape(t_arr.shape)
    return float(res_reshaped) if is_scalar else res_reshaped


def Q_t(t, D, mu_x, sigma_x, A):
    """
    Calculates the cumulative material outgassed per unit area by time t.

    Integrates the initial concentration profile convolved with the complementary
    error function kernel, assuming no radioactive decay.

    Parameters:
        t (float or array-like): Elapsed time point(s) [h].
        D (float): Diffusion coefficient [nm^2 / h].
        mu_x (float): Implantation depth parameter [nm].
        sigma_x (float): Depth straggle [nm].
        A (float): Volumetric amplitude [atoms / nm^3].

    Returns:
        float or numpy.ndarray: Cumulative outgassed areal density [atoms / nm^2].
    """
    t_arr = np.asarray(t, dtype=float)
    is_scalar = t_arr.ndim == 0
    t_flat = t_arr.ravel()
    results = np.zeros_like(t_flat, dtype=float)

    mask = t_flat > 0
    for i in np.where(mask)[0]:
        t_val = t_flat[i]

        def integrand(x):
            C0_x = C_x_t(x, 0, D, mu_x, sigma_x, A, lambda_decay=0.0)
            diff_term = erfc(-x / np.sqrt(4.0 * D * t_val))
            return C0_x * diff_term

        val, _ = quad(integrand, -np.inf, 0)
        results[i] = val

    res_reshaped = results.reshape(t_arr.shape)
    return float(res_reshaped) if is_scalar else res_reshaped


def N_chamber_t(t, D, mu_x, sigma_x, A, lambda_h, S_nm2):
    """Calculates the instantaneous inventory of Radon atoms in the vacuum chamber.

    Integrates the cumulative surface flux over the surface area S and attenuates
    by the radioactive decay factor exp(-lambda_h * t).

    Parameters:
        t (float or array-like): Elapsed time point(s) [h].
        D (float): Diffusion coefficient [nm^2 / h].
        mu_x (float): Implantation depth parameter [nm].
        sigma_x (float): Depth straggle [nm].
        A (float): Volumetric amplitude [atoms / nm^3].
        lambda_h (float): Radioactive decay constant [1/h].
        S_nm2 (float): Exposed surface area [nm^2].

    Returns:
        float or numpy.ndarray: Total number of active (undecayed) Radon-222 atoms
        present in the chamber [atoms] (dimensionless count).
    """
    t_arr = np.asarray(t, dtype=float)
    is_scalar = t_arr.ndim == 0
    t_flat = t_arr.ravel()
    N_values = np.zeros_like(t_flat, dtype=float)
    
    mask = t_flat > 0

    def inventory_integrand(tau):
        return J_diff_t_reconstructed(tau, D, mu_x, sigma_x, A) * S_nm2

    for i in np.where(mask)[0]:
        t_val = t_flat[i]
        integral, _ = quad(inventory_integrand, 0, t_val)
        N_values[i] = np.exp(-lambda_h * t_val) * integral

    res_reshaped = N_values.reshape(t_arr.shape)
    return float(res_reshaped) if is_scalar else res_reshaped


def model_rate(t, D, A, mu_z, sigma_z, lambda_h, S_nm2, eff):
    """
    Calculates the expected detector count rate at elapsed time t.

    Computes the instantaneous chamber activity A_ch = N(t) * lambda and
    scales by the detector efficiency: R(t) = eff * N(t) * (lambda_h / 3600).

    Parameters:
        t (float or array-like): Elapsed observation time [h].
        D (float): Diffusion coefficient [nm^2 / h].
        A (float): Volumetric amplitude [atoms / nm^3].
        mu_z (float): Implantation depth [nm].
        sigma_z (float): Depth straggle [nm].
        lambda_h (float): Radioactive decay constant [1/h].
        S_nm2 (float): Exposed surface area [nm^2].
        eff (float): Total detection efficiency [counts / decay] (dimensionless, 0 to 1).

    Returns:
        float or numpy.ndarray: Expected detector event rate [Hz = counts / s].
    """
    N = N_chamber_t(t, D, mu_z, sigma_z, A, lambda_h, S_nm2)
    return eff * N * (lambda_h / 3600.0)


def create_cost_function(times_h, raw_rates, raw_yerr, priors=None):
    """
    Creates a chi2 cost function with Gaussian pull terms for iminuit fitting.

    Parameters:
        times_h (array-like): Sample observation times [h].
        raw_rates (array-like): Measured detector event rates [Hz].
        raw_yerr (array-like): 1-sigma statistical uncertainties on event rates [Hz].
        priors (dict, optional): Dictionary mapping parameter names to (mean, std)
            tuples defining Gaussian penalty terms.

    Returns:
        callable: Cost function `cost(D, activity_bq, mu_z, sigma_z, lambda_h, S_nm2, eff)`
        returning the total penalized objective value [dimensionless] (chi2 + penalties),
        configured with `cost.errordef = 1.0` for least-squares optimization.
    """    
    _times = np.asarray(times_h, dtype=float)
    _rates = np.asarray(raw_rates, dtype=float)
    _yerr = np.asarray(raw_yerr, dtype=float)
    _priors = priors or {}

    def cost(D, activity_bq, mu_z, sigma_z, lambda_h, S_nm2, eff):
        _A = calculate_A_volumetric(activity_bq, mu_z, sigma_z, lambda_h, S_nm2)
        model = model_rate(_times, D, _A, mu_z, sigma_z, lambda_h, S_nm2, eff)
        chi2 = np.sum(((_rates - model) / _yerr) ** 2)

        penalties = 0.0
        param_values = {
            "D": D,
            "activity_bq": activity_bq,
            "mu_z": mu_z,
            "sigma_z": sigma_z,
            "lambda_h": lambda_h,
            "S_nm2": S_nm2,
            "eff": eff,
        }
        for pname, (p_mean, p_std) in _priors.items():
            if pname in param_values and p_std > 0:
                penalties += ((param_values[pname] - p_mean) / p_std) ** 2

        return chi2 + penalties

    # Explicitly tell iminuit this is a least-squares (chi2) cost function (errordef = 1.0)
    cost.errordef = 1.0
    return cost


class Radon222Diffusion1D:
    """1D Radon Diffusion Model class using Method of Images for Gaussian distribution."""

    RN222_HALF_LIFE_H = RN222_HALF_LIFE_H
    RN222_LAMBDA_H = RN222_LAMBDA_H
    FROM_NM2_H_TO_CM2_SEC = FROM_NM2_H_TO_CM2_SEC
    FROM_CM2_SEC_TO_NM2_H = FROM_CM2_SEC_TO_NM2_H

    # Expose mathematical functions as staticmethods
    calculate_A_volumetric = staticmethod(calculate_A_volumetric)
    calculate_activity_from_A = staticmethod(calculate_activity_from_A)
    concentration = staticmethod(C_x_t)
    C_x_t = staticmethod(C_x_t)
    flux = staticmethod(J_diff_t_reconstructed)
    J_diff_t_reconstructed = staticmethod(J_diff_t_reconstructed)
    cumulative_outgassed = staticmethod(Q_t)
    Q_t = staticmethod(Q_t)
    chamber_inventory = staticmethod(N_chamber_t)
    N_chamber_t = staticmethod(N_chamber_t)
    model_rate = staticmethod(model_rate)
    create_cost_function = staticmethod(create_cost_function)

    def __init__(
        self,
        mu_z,
        sigma_z,
        activity_bq,
        eff=0.30,
        material=None,
        lambda_h=RN222_LAMBDA_H,
        S_nm2=4e14,
        D_init=1.0,
    ):
        """Initializes Radon222Diffusion1D with depth, straggle, activity, and efficiency.

        Parameters:
            mu_z (float): Implantation depth [nm].
            sigma_z (float): Depth straggle [nm].
            activity_bq (float): Initial activity [Bq].
            eff (float, optional): Detection efficiency (0 to 1). Defaults to 0.30 (30%).
            material (str, optional): Material key ('Cu', 'Ti', 'SS') for reference.
            lambda_h (float, optional): Radioactive decay constant [1/h]. Defaults to 222Rn decay.
            S_nm2 (float, optional): Surface area [nm^2]. Defaults to 4e14 (4 cm^2).
            D_init (float, optional): Initial diffusion coefficient [nm^2 / h]. Defaults to 1.0.
        """
        self.mu_z = mu_z
        self.sigma_z = sigma_z
        self.activity_bq = activity_bq
        self.eff = eff
        self.material = material
        self.lambda_h = lambda_h
        self.S_nm2 = S_nm2
        self.D_init = D_init

    def calculate_A(self, activity_bq=None):
        """
        Calculates volumetric amplitude A from instance or specified activity.

        Parameters:
            activity_bq (float, optional): Activity [Bq]. Uses `self.activity_bq` if None.

        Returns:
            float: Volumetric amplitude [atoms / nm^3].
        """
        _act = activity_bq if activity_bq is not None else self.activity_bq
        return calculate_A_volumetric(
            _act, self.mu_z, self.sigma_z, self.lambda_h, self.S_nm2
        )

    def eval_rate(self, t, D=None, activity_bq=None, eff=None):
        """
        Evaluates model count rate at time points t.

        Parameters:
            t (float or array-like): Time point(s) [h].
            D (float, optional): Diffusion coefficient [nm^2 / h]. Defaults to `self.D_init`.
            activity_bq (float, optional): Initial activity [Bq]. Defaults to `self.activity_bq`.
            eff (float, optional): Detection efficiency. Defaults to `self.eff`.

        Returns:
            float or numpy.ndarray: Modeled detector rate [Hz].
        """
        _D = D if D is not None else self.D_init
        _eff = eff if eff is not None else self.eff
        _activity_bq = activity_bq if activity_bq is not None else self.activity_bq

        if _activity_bq is None:
            raise ValueError(
                "You need to pass the initial activity to compute the volumetric amplitude A, "
                "or set activity_bq in the instance."
            )

        _A = calculate_A_volumetric(
            _activity_bq, self.mu_z, self.sigma_z, self.lambda_h, self.S_nm2
        )
        return model_rate(
            t, _D, _A, self.mu_z, self.sigma_z, self.lambda_h, self.S_nm2, _eff
        )

    def get_minuit_fitting_object(
        self,
        times_h,
        rates_hz,
        rate_err,
        D_init=None,
        priors=None,
        fixed=None,
        limits=None,
    ):
        """
        Builds and configures the iminuit.Minuit optimization object.

        Parameters:
            times_h (array-like): Observation timestamps [h].
            rates_hz (array-like): Measured count rates [Hz].
            rate_err (array-like): Uncertainties on rates [Hz].
            D_init (float, optional): Initial guess for D [nm^2 / h].
            priors (dict, optional): Gaussian prior constraints `{param: std}` or `{param: (mean, std)}`.
            fixed (dict, optional): Parameter fixation overrides `{param: bool}`.
            limits (dict, optional): Parameter boundary overrides `{param: (min, max)}`.

        Returns:
            iminuit.Minuit: Configured Minuit instance initialized with starting values,
            limits, parameter errors, and cost function ready for `.migrad()`.
        """
        if Minuit is None:
            raise ImportError("iminuit package is required for fitting functionality.")

        _D = D_init if D_init is not None else self.D_init

        starting_values = {
            "D": _D,
            "activity_bq": self.activity_bq,
            "mu_z": self.mu_z,
            "sigma_z": self.sigma_z,
            "lambda_h": self.lambda_h,
            "S_nm2": self.S_nm2,
            "eff": self.eff,
        }

        # Normalize priors: allow either `std` (float) or `(mean, std)` (tuple)
        normalized_priors = {}
        if priors is not None:
            for param, prior_val in priors.items():
                if param not in starting_values:
                    raise ValueError(
                        f"Unknown parameter '{param}' in priors. "
                        f"Valid parameters are: {list(starting_values.keys())}"
                    )
                start_val = starting_values[param]

                if isinstance(prior_val, (int, float, np.number)):
                    prior_mean = start_val
                    prior_std = float(prior_val)
                else:
                    prior_mean, prior_std = prior_val
                    if not np.isclose(prior_mean, start_val, rtol=1e-5, atol=1e-12):
                        warnings.warn(
                            f"Mismatch in '{param}': The prior center ({prior_mean}) "
                            f"does not match the initialized starting value ({start_val}). "
                            "It is highly recommended to start the fit at the prior's center."
                        )

                if prior_std <= 0:
                    raise ValueError(
                        f"Prior standard deviation for '{param}' must be > 0, got {prior_std}."
                    )

                normalized_priors[param] = (prior_mean, prior_std)

        cost_func = create_cost_function(
            times_h, rates_hz, rate_err, priors=normalized_priors or None
        )

        m = Minuit(
            cost_func,
            D=_D,
            activity_bq=self.activity_bq,
            mu_z=self.mu_z,
            sigma_z=self.sigma_z,
            lambda_h=self.lambda_h,
            S_nm2=self.S_nm2,
            eff=self.eff,
        )

        # Fix all nuisance parameters by default; float only D
        default_fixed = {
            "D": False,
            "activity_bq": True,
            "mu_z": True,
            "sigma_z": True,
            "lambda_h": True,
            "S_nm2": True,
            "eff": True,
        }

        # Automatically release parameters that have a prior assigned
        for param in normalized_priors:
            if param in default_fixed:
                default_fixed[param] = False

        # Allow explicit user `fixed` dict to override
        if fixed:
            default_fixed.update(fixed)

        for param, is_fixed in default_fixed.items():
            if param in m.parameters:
                m.fixed[param] = is_fixed

        # Set physical boundaries on parameters so Minuit doesn't step into unphysical regions
        default_limits = {
            "D": (0, None),
            "activity_bq": (0, None),
            "mu_z": (0, None),
            "sigma_z": (1e-6, None),
            "lambda_h": (0, None),
            "S_nm2": (0, None),
            "eff": (0, 1),
        }
        if limits:
            default_limits.update(limits)

        for param, lim in default_limits.items():
            if param in m.parameters:
                m.limits[param] = lim

        # Set initial step sizes (m.errors) to the prior standard deviations (sigma_prior).
        # This seeds Minuit's initial Hessian estimate with the expected penalty curvature
        # (d^2chi2/dtheta^2 ~ 2/sigma^2) and ensures finite-difference gradient steps
        # match the physical parameter scale, preventing line-search overshoot in MIGRAD.
        for param, (_, prior_std) in normalized_priors.items():
            if param in m.parameters:
                m.errors[param] = prior_std

        # Bind references directly to the Minuit instance
        self.initial_priors = normalized_priors
        self.initial_fixed = default_fixed
        self.initial_limits = default_limits
        
        return m

    def fit(
        self,
        times_h,
        rates_hz,
        rate_err,
        D_init=None,
        priors=None,
        fixed=None,
        limits=None,
    ):
        """
        Convenience wrapper returning a configured iminuit.Minuit fitting object.

        Returns:
            iminuit.Minuit: Configured Minuit instance.
        """
        m = self.get_minuit_fitting_object(
            times_h,
            rates_hz,
            rate_err,
            D_init=D_init,
            priors=priors,
            fixed=fixed,
            limits=limits,
        )

        m.migrad()
        m.hesse()
        return m

    def plot_fit(
        self,
        times_sec,
        rates_hz,
        rate_err,
        minuit_obj,
        mask_fit=None,
        title=None,
        n_points=500,
    ):
        """
        Plots observed rates, model prediction, and pull residuals.

        Parameters:
            times_sec (array-like): Elapsed time series [s].
            rates_hz (array-like): Measured count rates [Hz].
            rate_err (array-like): Rate uncertainties [Hz].
            minuit_obj (iminuit.Minuit): Fitted Minuit object containing parameter values.
            mask_fit (array-like of bool, optional): Boolean mask of fitted data points.
            title (str, optional): Figure suptitle.
            n_points (int, optional): Number of points for continuous curve. Defaults to 500.

        Returns:
            tuple[matplotlib.figure.Figure, tuple[matplotlib.axes.Axes, matplotlib.axes.Axes]]:
            The generated figure and axis objects `(fig, (ax_data, ax_residuals))`.
        """
        times_sec = np.asarray(times_sec, dtype=float)
        rates_hz = np.asarray(rates_hz, dtype=float)
        rate_err = np.asarray(rate_err, dtype=float)

        if mask_fit is None:
            mask_fit = np.ones_like(times_sec, dtype=bool)
        else:
            mask_fit = np.asarray(mask_fit, dtype=bool)

        times_days = times_sec / 86400.0
        times_fine_sec = np.linspace(np.min(times_sec), np.max(times_sec), n_points)
        times_fine_days = times_fine_sec / 86400.0

        _A = calculate_A_volumetric(
            minuit_obj.values["activity_bq"],
            minuit_obj.values["mu_z"],
            minuit_obj.values["sigma_z"],
            minuit_obj.values["lambda_h"],
            minuit_obj.values["S_nm2"],
        )

        model_rates_fine = model_rate(
            times_fine_sec / 3600.0,
            minuit_obj.values["D"],
            _A,
            minuit_obj.values["mu_z"],
            minuit_obj.values["sigma_z"],
            minuit_obj.values["lambda_h"],
            minuit_obj.values["S_nm2"],
            minuit_obj.values["eff"],
        )

        d_cm2_s = minuit_obj.values["D"] * self.FROM_NM2_H_TO_CM2_SEC
        d_err_cm2_s = (
            minuit_obj.errors["D"] * self.FROM_NM2_H_TO_CM2_SEC
            if "D" in minuit_obj.errors
            else 0.0
        )

        label = f"Diffusion constant: {d_cm2_s:.3e} ± {d_err_cm2_s:.3e} cm$^2$ s$^{-1}$"

        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(8, 6),
            dpi=150,
            gridspec_kw={"height_ratios": [2, 1]},
            sharex=True,
        )

        if np.any(~mask_fit):
            unfit_times = times_days[~mask_fit]
            for ax in (ax1, ax2):
                ax.axvspan(
                    np.min(unfit_times) - 0.1,
                    np.max(unfit_times) + 0.1,
                    color="grey",
                    lw=0,
                    alpha=0.25,
                )

        ax1.errorbar(
            times_days,
            rates_hz,
            yerr=rate_err,
            lw=0,
            elinewidth=1,
            marker="o",
            ms=2,
            color="blue",
            alpha=0.6,
            label="Data",
        )
        ax1.plot(times_fine_days, model_rates_fine, color="red", lw=1.5, label=label)
        ax1.set_ylabel("Rate [Hz]")
        ax1.legend(loc="upper right")
        ax1.grid(True, linestyle="--", alpha=0.5)

        fitted_rates = model_rate(
            times_sec[mask_fit] / 3600.0,
            minuit_obj.values["D"],
            _A,
            minuit_obj.values["mu_z"],
            minuit_obj.values["sigma_z"],
            minuit_obj.values["lambda_h"],
            minuit_obj.values["S_nm2"],
            minuit_obj.values["eff"],
        )
        residuals = (rates_hz[mask_fit] - fitted_rates) / rate_err[mask_fit]

        # Draw 2-sigma band first so 1-sigma band renders cleanly on top
        ax2.fill_between(times_fine_days, -2, 2, color="yellow", alpha=0.5, lw=0)
        ax2.fill_between(times_fine_days, -1, 1, color="green", alpha=0.5, lw=0)

        ax2.errorbar(
            times_days[mask_fit],
            residuals,
            yerr=1,
            lw=0,
            elinewidth=1,
            marker="o",
            ms=2,
            color="blue",
            alpha=0.6,
        )

        ax2.set_xlabel("Elapsed time [days]")
        ax2.set_ylabel("Normalized\nResiduals")
        ax2.grid(True, linestyle="--", alpha=0.5)
        ax2.set_ylim(-3, 3)

        if title:
            fig.suptitle(title)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)

        # return fig, (ax1, ax2)
        
    
    def plot_pulls(
        self,
        minuit_obj,
    ):
        """
        Plots standardized pulls for nuisance parameters directly from the fit result.

        Parameters:
            minuit_obj (iminuit.Minuit): Fitted Minuit object containing .values and .errors.

        Returns:
            tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: (fig, ax)
        """
        priors = getattr(self, "initial_priors", None)
        if not priors:
            raise ValueError(
                "No `initial_priors` found on this Radon222Diffusion1D instance. "
                "Ensure you ran `fit()` with priors specified."
            )

        param_names = []
        pulls = []
        pull_errs = []

        for name, prior in priors.items():
            p_mean, p_std = prior if isinstance(prior, (tuple, list)) else (minuit_obj.values[name], float(prior))
            val = minuit_obj.values[name]
            err = minuit_obj.errors[name] if (name in minuit_obj.errors and not minuit_obj.fixed.get(name, False)) else 0.0

            param_names.append(name)
            pulls.append((val - p_mean) / p_std)
            pull_errs.append(err / p_std)

        
        fig, ax = plt.subplots(figsize=(7, 4), dpi=150)

        y_pos = np.arange(len(param_names))
        max_dev = max(3.0, np.max(np.abs(np.array(pulls) + np.array(pull_errs))) * 1.15)

        # +/- 1 sigma and +/- 2 sigma prior expectation bands
        ax.axvspan(-2.0, 2.0, color="#fff2a8", alpha=0.6, lw=0, label=r"$\pm 2\sigma$ prior band")
        ax.axvspan(-1.0, 1.0, color="#a1e8a1", alpha=0.6, lw=0, label=r"$\pm 1\sigma$ prior band")
        ax.axvline(0.0, color="black", linestyle="-", lw=1.0, alpha=0.7)
        for ref in (-2.0, -1.0, 1.0, 2.0):
            ax.axvline(ref, color="gray", linestyle="--", lw=0.8, alpha=0.5)

        # Post-fit pulls
        ax.errorbar(
            pulls,
            y_pos,
            xerr=pull_errs,
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.5,
            capsize=3,
            capthick=1.2,
            ms=5,
            zorder=5,
            label="Post-fit pull",
        )

        ax.set_yticks(y_pos)
        ax.set_yticklabels(param_names, fontsize=10)
        ax.invert_yaxis()
        ax.set_xlim(-max_dev, max_dev)
        ax.set_xlabel(r"Pull: $(\hat{\theta} - \theta_{\mathrm{prior}}) \,/\, \sigma_{\mathrm{prior}}$", fontsize=11)
        ax.set_title("Nuisance Parameter Pulls", fontsize=12, pad=10)
        ax.grid(True, axis="x", linestyle=":", alpha=0.5)
        ax.legend(loc="upper right", framealpha=0.85, fontsize=8)

        plt.tight_layout()
        # return fig, ax