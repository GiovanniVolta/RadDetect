import warnings

import numpy as np
from iminuit import Minuit
from scipy.stats import crystalball, norm, skewnorm


def _warn_if_multiple_modes(mode, peak_name):
    """Internal helper to check active modes and raise a warning if > 1.

    Parameters:
        mode (str or list): The mode(s) activated for the peak.
        peak_name (str): Identifier for the peak to display in the warning.
    """
    active_modes = [m for m in ["cb", "sg", "norm"] if m in mode]
    if len(active_modes) > 1:
        warnings.warn(
            f"Multiple modes {active_modes} are active for {peak_name}. "
            "Fitting multiple complex shapes simultaneously can lead to "
            "highly correlated parameters and failure to converge. "
            "Ensure this is intentional.",
            UserWarning,
            stacklevel=3,
        )


class SinglePeak:
    """Functor model for fitting a single alpha peak with an optional flat background.

    Because this is a callable class, you instantiate it with your configuration
    settings, and then pass the instance directly to iminuit.
    """

    def __init__(self, peak_mode="cb", use_bkg=True):
        """Configures the SinglePeak instance."""
        self.peak_mode = peak_mode
        self.use_bkg = use_bkg
        _warn_if_multiple_modes(self.peak_mode, "SinglePeak")

        # Route the evaluation functions ONCE during initialization
        self._calc_bkg = self._eval_bkg if self.use_bkg else self._eval_zero
        self._calc_cb = self._eval_cb if "cb" in self.peak_mode else self._eval_zero
        self._calc_sg = self._eval_sg if "sg" in self.peak_mode else self._eval_zero
        self._calc_n = self._eval_norm if "norm" in self.peak_mode else self._eval_zero

    def _eval_zero(self, x, *args):
        """Catch-all for disabled components; absorbs unused parameters via *args."""
        return np.zeros_like(x, dtype=float)

    def _eval_bkg(self, x, bkg):
        return np.full_like(x, bkg, dtype=float)

    def _eval_cb(self, x, beta, m, loc, scale, A):
        return A * crystalball.pdf(x, beta, m, loc=loc, scale=scale)

    def _eval_sg(self, x, a, loc, scale, A):
        return A * skewnorm.pdf(x, a, loc=loc, scale=scale)

    def _eval_norm(self, x, loc, scale, A):
        return A * norm.pdf(x, loc=loc, scale=scale)

    def total_model(
        self,
        x,
        bkg,
        beta_cb,
        m_cb,
        loc_cb,
        scale_cb,
        A_cb,
        a_sg,
        loc_sg,
        scale_sg,
        A_sg,
        loc_n,
        scale_n,
        A_n,
    ):
        """Evaluates the model array.

        This method is inspected and called by iminuit.
        """
        # Execute the routed functions
        self.bkg_term = self._calc_bkg(x, bkg)
        self.cb = self._calc_cb(x, beta_cb, m_cb, loc_cb, scale_cb, A_cb)
        self.sg = self._calc_sg(x, a_sg, loc_sg, scale_sg, A_sg)
        self.n = self._calc_n(x, loc_n, scale_n, A_n)

        self.peak = self.cb + self.sg + self.n

        # Sum them up
        self.total = self.bkg_term + self.peak

        return self.total

    # This is for iminuit
    __call__ = total_model

    def prepare_fit_args(self, p_cfg):
        """Parses a master parameter config dictionary and auto-fixes parameters for
        components that are disabled in this instance.

        Parameters:
            p_cfg (dict): Dictionary where values are (init, fixed, limits).

        Returns:
            tuple: (init_dict, fixed_dict, limits_dict)
        """
        init = {k: v[0] for k, v in p_cfg.items()}
        fixed = {k: v[1] for k, v in p_cfg.items()}
        limits = {k: v[2] for k, v in p_cfg.items()}

        if "cb" not in self.peak_mode:
            fixed.update({k: True for k in fixed if "_cb" in k})
        if "sg" not in self.peak_mode:
            fixed.update({k: True for k in fixed if "_sg" in k})
        if "norm" not in self.peak_mode:
            fixed.update({k: True for k in fixed if "_n" in k})

        if not self.use_bkg:
            fixed["bkg"] = True
            init["bkg"] = 0.0

        return init, fixed, limits

    def extract_params(self, m: Minuit):
        """Extracts the location and scale parameters for the active
        shape from a fitted
        Minuit object.

        Parameters:
            m (Minuit): The iminuit object containing the completed fit.

        Returns:
            tuple: (loc, loc_err, scale, scale_err)
        """
        if not m.valid:
            warnings.warn("Extracting parameters from an invalid fit.")

        if "cb" in self.peak_mode:
            return (
                m.values["loc_cb"],
                m.errors["loc_cb"],
                m.values["scale_cb"],
                m.errors["scale_cb"],
            )
        if "sg" in self.peak_mode:
            return (
                m.values["loc_sg"],
                m.errors["loc_sg"],
                m.values["scale_sg"],
                m.errors["scale_sg"],
            )
        if "norm" in self.peak_mode:
            return (
                m.values["loc_n"],
                m.errors["loc_n"],
                m.values["scale_n"],
                m.errors["scale_n"],
            )

        raise ValueError("No valid peak mode found to extract.")


class DoublePeak:
    """Functor model for fitting two adjacent alpha peaks with an optional flat
    background.

    Instantiate with specific configurations for Peak 1 and Peak 2,
    then pass the instance directly to iminuit.
    """

    def __init__(
        self,
        peak1_mode="cb",
        peak2_mode="cb",
        use_peak1=True,
        use_peak2=True,
        use_bkg=True,
    ):
        """Configures the DoublePeak instance."""
        self.peak1_mode = peak1_mode
        self.peak2_mode = peak2_mode
        self.use_peak1 = use_peak1
        self.use_peak2 = use_peak2
        self.use_bkg = use_bkg

        if self.use_peak1:
            _warn_if_multiple_modes(self.peak1_mode, "Peak 1")
        if self.use_peak2:
            _warn_if_multiple_modes(self.peak2_mode, "Peak 2")

        # Route evaluation functions ONCE during initialization
        self._calc_bkg = self._eval_bkg if self.use_bkg else self._eval_zero

        # Peak 1 routing
        self._calc_cb_1 = (
            self._eval_cb
            if ("cb" in self.peak1_mode and self.use_peak1)
            else self._eval_zero
        )
        self._calc_sg_1 = (
            self._eval_sg
            if ("sg" in self.peak1_mode and self.use_peak1)
            else self._eval_zero
        )
        self._calc_n_1 = (
            self._eval_norm
            if ("norm" in self.peak1_mode and self.use_peak1)
            else self._eval_zero
        )

        # Peak 2 routing
        self._calc_cb_2 = (
            self._eval_cb
            if ("cb" in self.peak2_mode and self.use_peak2)
            else self._eval_zero
        )
        self._calc_sg_2 = (
            self._eval_sg
            if ("sg" in self.peak2_mode and self.use_peak2)
            else self._eval_zero
        )
        self._calc_n_2 = (
            self._eval_norm
            if ("norm" in self.peak2_mode and self.use_peak2)
            else self._eval_zero
        )

    def _eval_zero(self, x, *args):
        """Catch-all for disabled components; absorbs extra parameters via *args."""
        return np.zeros_like(x, dtype=float)

    def _eval_bkg(self, x, bkg):
        return np.full_like(x, bkg, dtype=float)

    def _eval_cb(self, x, beta, m, loc, scale, A):
        return A * crystalball.pdf(x, beta, m, loc=loc, scale=scale)

    def _eval_sg(self, x, a, loc, scale, A):
        return A * skewnorm.pdf(x, a, loc=loc, scale=scale)

    def _eval_norm(self, x, loc, scale, A):
        return A * norm.pdf(x, loc=loc, scale=scale)

    def total_model(
        self,
        x,
        bkg,
        beta_cb_1,
        m_cb_1,
        loc_cb_1,
        scale_cb_1,
        A_cb_1,
        a_sg_1,
        loc_sg_1,
        scale_sg_1,
        A_sg_1,
        loc_n_1,
        scale_n_1,
        A_n_1,
        beta_cb_2,
        m_cb_2,
        loc_cb_2,
        scale_cb_2,
        A_cb_2,
        a_sg_2,
        loc_sg_2,
        scale_sg_2,
        A_sg_2,
        loc_n_2,
        scale_n_2,
        A_n_2,
    ):
        """Evaluates the double peak model array.

        Inspected and called by iminuit.
        """
        self.bkg_term = self._calc_bkg(x, bkg)

        # Evaluate Peak 1
        self.cb_1 = self._calc_cb_1(x, beta_cb_1, m_cb_1, loc_cb_1, scale_cb_1, A_cb_1)
        self.sg_1 = self._calc_sg_1(x, a_sg_1, loc_sg_1, scale_sg_1, A_sg_1)
        self.n_1 = self._calc_n_1(x, loc_n_1, scale_n_1, A_n_1)
        self.peak1 = self.cb_1 + self.sg_1 + self.n_1

        # Evaluate Peak 2
        self.cb_2 = self._calc_cb_2(x, beta_cb_2, m_cb_2, loc_cb_2, scale_cb_2, A_cb_2)
        self.sg_2 = self._calc_sg_2(x, a_sg_2, loc_sg_2, scale_sg_2, A_sg_2)
        self.n_2 = self._calc_n_2(x, loc_n_2, scale_n_2, A_n_2)
        self.peak2 = self.cb_2 + self.sg_2 + self.n_2

        # Sum total
        self.total = self.bkg_term + self.peak1 + self.peak2

        return self.total

    # This is for iminuit
    __call__ = total_model

    def prepare_fit_args(self, p_cfg):
        """Parses a master parameter config dictionary and auto-fixes parameters for
        components that are disabled in this instance.

        Parameters:
            p_cfg (dict): Dictionary where values are (init, fixed, limits).

        Returns:
            tuple: (init_dict, fixed_dict, limits_dict)
        """
        init = {k: v[0] for k, v in p_cfg.items()}
        fixed = {k: v[1] for k, v in p_cfg.items()}
        limits = {k: v[2] for k, v in p_cfg.items()}

        # Peak 1 components
        if "cb" not in self.peak1_mode or not self.use_peak1:
            fixed.update({k: True for k in fixed if "_cb_1" in k})
        if "sg" not in self.peak1_mode or not self.use_peak1:
            fixed.update({k: True for k in fixed if "_sg_1" in k})
        if "norm" not in self.peak1_mode or not self.use_peak1:
            fixed.update({k: True for k in fixed if "_n_1" in k})

        # Peak 2 components
        if "cb" not in self.peak2_mode or not self.use_peak2:
            fixed.update({k: True for k in fixed if "_cb_2" in k})
        if "sg" not in self.peak2_mode or not self.use_peak2:
            fixed.update({k: True for k in fixed if "_sg_2" in k})
        if "norm" not in self.peak2_mode or not self.use_peak2:
            fixed.update({k: True for k in fixed if "_n_2" in k})

        # Background
        if not self.use_bkg:
            fixed["bkg"] = True
            init["bkg"] = 0.0

        return init, fixed, limits

    def extract_peak1_params(self, m: Minuit):
        """Extracts loc and scale parameters for Peak 1 from a fitted Minuit object."""
        if not m.valid:
            warnings.warn("Extracting from invalid fit.")
        if not self.use_peak1:
            raise ValueError("Peak 1 is disabled.")

        if "cb" in self.peak1_mode:
            return (
                m.values["loc_cb_1"],
                m.errors["loc_cb_1"],
                m.values["scale_cb_1"],
                m.errors["scale_cb_1"],
            )
        if "sg" in self.peak1_mode:
            return (
                m.values["loc_sg_1"],
                m.errors["loc_sg_1"],
                m.values["scale_sg_1"],
                m.errors["scale_sg_1"],
            )
        if "norm" in self.peak1_mode:
            return (
                m.values["loc_n_1"],
                m.errors["loc_n_1"],
                m.values["scale_n_1"],
                m.errors["scale_n_1"],
            )

        raise ValueError("No valid peak1_mode found.")

    def extract_peak2_params(self, m: Minuit):
        """Extracts loc and scale parameters for Peak 2 from a fitted Minuit object."""
        if not m.valid:
            warnings.warn("Extracting from invalid fit.")
        if not self.use_peak2:
            raise ValueError("Peak 2 is disabled.")

        if "cb" in self.peak2_mode:
            return (
                m.values["loc_cb_2"],
                m.errors["loc_cb_2"],
                m.values["scale_cb_2"],
                m.errors["scale_cb_2"],
            )
        if "sg" in self.peak2_mode:
            return (
                m.values["loc_sg_2"],
                m.errors["loc_sg_2"],
                m.values["scale_sg_2"],
                m.errors["scale_sg_2"],
            )
        if "norm" in self.peak2_mode:
            return (
                m.values["loc_n_2"],
                m.errors["loc_n_2"],
                m.values["scale_n_2"],
                m.errors["scale_n_2"],
            )

        raise ValueError("No valid peak2_mode found.")


class TriplePeak:
    """Functor model for fitting three adjacent peaks with an optional flat background.

    Instantiate with specific configurations for Peak 1, Peak 2, and Peak 3, then pass
    the instance directly to iminuit.
    """

    def __init__(
        self,
        peak1_mode="cb",
        peak2_mode="cb",
        peak3_mode="cb",
        use_peak1=True,
        use_peak2=True,
        use_peak3=True,
        use_bkg=True,
    ):
        """Configures the TriplePeak instance."""
        self.peak1_mode = peak1_mode
        self.peak2_mode = peak2_mode
        self.peak3_mode = peak3_mode
        self.use_peak1 = use_peak1
        self.use_peak2 = use_peak2
        self.use_peak3 = use_peak3
        self.use_bkg = use_bkg

        # Note: Assuming _warn_if_multiple_modes is defined elsewhere in your module
        if self.use_peak1:
            _warn_if_multiple_modes(self.peak1_mode, "Peak 1")
        if self.use_peak2:
            _warn_if_multiple_modes(self.peak2_mode, "Peak 2")
        if self.use_peak3:
            _warn_if_multiple_modes(self.peak3_mode, "Peak 3")

        # Route evaluation functions ONCE during initialization
        self._calc_bkg = self._eval_bkg if self.use_bkg else self._eval_zero

        # Peak 1 routing
        self._calc_cb_1 = (
            self._eval_cb
            if ("cb" in self.peak1_mode and self.use_peak1)
            else self._eval_zero
        )
        self._calc_sg_1 = (
            self._eval_sg
            if ("sg" in self.peak1_mode and self.use_peak1)
            else self._eval_zero
        )
        self._calc_n_1 = (
            self._eval_norm
            if ("norm" in self.peak1_mode and self.use_peak1)
            else self._eval_zero
        )

        # Peak 2 routing
        self._calc_cb_2 = (
            self._eval_cb
            if ("cb" in self.peak2_mode and self.use_peak2)
            else self._eval_zero
        )
        self._calc_sg_2 = (
            self._eval_sg
            if ("sg" in self.peak2_mode and self.use_peak2)
            else self._eval_zero
        )
        self._calc_n_2 = (
            self._eval_norm
            if ("norm" in self.peak2_mode and self.use_peak2)
            else self._eval_zero
        )

        # Peak 3 routing
        self._calc_cb_3 = (
            self._eval_cb
            if ("cb" in self.peak3_mode and self.use_peak3)
            else self._eval_zero
        )
        self._calc_sg_3 = (
            self._eval_sg
            if ("sg" in self.peak3_mode and self.use_peak3)
            else self._eval_zero
        )
        self._calc_n_3 = (
            self._eval_norm
            if ("norm" in self.peak3_mode and self.use_peak3)
            else self._eval_zero
        )

    def _eval_zero(self, x, *args):
        """Catch-all for disabled components; absorbs extra parameters via *args."""
        return np.zeros_like(x, dtype=float)

    def _eval_bkg(self, x, bkg):
        return np.full_like(x, bkg, dtype=float)

    def _eval_cb(self, x, beta, m, loc, scale, A):
        return A * crystalball.pdf(x, beta, m, loc=loc, scale=scale)

    def _eval_sg(self, x, a, loc, scale, A):
        return A * skewnorm.pdf(x, a, loc=loc, scale=scale)

    def _eval_norm(self, x, loc, scale, A):
        return A * norm.pdf(x, loc=loc, scale=scale)

    def total_model(
        self,
        x,
        bkg,
        beta_cb_1,
        m_cb_1,
        loc_cb_1,
        scale_cb_1,
        A_cb_1,
        a_sg_1,
        loc_sg_1,
        scale_sg_1,
        A_sg_1,
        loc_n_1,
        scale_n_1,
        A_n_1,
        beta_cb_2,
        m_cb_2,
        loc_cb_2,
        scale_cb_2,
        A_cb_2,
        a_sg_2,
        loc_sg_2,
        scale_sg_2,
        A_sg_2,
        loc_n_2,
        scale_n_2,
        A_n_2,
        beta_cb_3,
        m_cb_3,
        loc_cb_3,
        scale_cb_3,
        A_cb_3,
        a_sg_3,
        loc_sg_3,
        scale_sg_3,
        A_sg_3,
        loc_n_3,
        scale_n_3,
        A_n_3,
    ):
        """Evaluates the triple peak model array.

        Inspected and called by iminuit.
        """
        self.bkg_term = self._calc_bkg(x, bkg)

        # Evaluate Peak 1
        self.cb_1 = self._calc_cb_1(x, beta_cb_1, m_cb_1, loc_cb_1, scale_cb_1, A_cb_1)
        self.sg_1 = self._calc_sg_1(x, a_sg_1, loc_sg_1, scale_sg_1, A_sg_1)
        self.n_1 = self._calc_n_1(x, loc_n_1, scale_n_1, A_n_1)
        self.peak1 = self.cb_1 + self.sg_1 + self.n_1

        # Evaluate Peak 2
        self.cb_2 = self._calc_cb_2(x, beta_cb_2, m_cb_2, loc_cb_2, scale_cb_2, A_cb_2)
        self.sg_2 = self._calc_sg_2(x, a_sg_2, loc_sg_2, scale_sg_2, A_sg_2)
        self.n_2 = self._calc_n_2(x, loc_n_2, scale_n_2, A_n_2)
        self.peak2 = self.cb_2 + self.sg_2 + self.n_2

        # Evaluate Peak 3
        self.cb_3 = self._calc_cb_3(x, beta_cb_3, m_cb_3, loc_cb_3, scale_cb_3, A_cb_3)
        self.sg_3 = self._calc_sg_3(x, a_sg_3, loc_sg_3, scale_sg_3, A_sg_3)
        self.n_3 = self._calc_n_3(x, loc_n_3, scale_n_3, A_n_3)
        self.peak3 = self.cb_3 + self.sg_3 + self.n_3

        # Sum total
        self.total = self.bkg_term + self.peak1 + self.peak2 + self.peak3

        return self.total

    # This is for iminuit
    __call__ = total_model

    def prepare_fit_args(self, p_cfg):
        """Parses a master parameter config dictionary and auto-fixes parameters for
        components that are disabled in this instance.

        Parameters:
            p_cfg (dict): Dictionary where values are (init, fixed, limits).

        Returns:
            tuple: (init_dict, fixed_dict, limits_dict)
        """
        init = {k: v[0] for k, v in p_cfg.items()}
        fixed = {k: v[1] for k, v in p_cfg.items()}
        limits = {k: v[2] for k, v in p_cfg.items()}

        # Peak 1 components
        if "cb" not in self.peak1_mode or not self.use_peak1:
            fixed.update({k: True for k in fixed if "_cb_1" in k})
        if "sg" not in self.peak1_mode or not self.use_peak1:
            fixed.update({k: True for k in fixed if "_sg_1" in k})
        if "norm" not in self.peak1_mode or not self.use_peak1:
            fixed.update({k: True for k in fixed if "_n_1" in k})

        # Peak 2 components
        if "cb" not in self.peak2_mode or not self.use_peak2:
            fixed.update({k: True for k in fixed if "_cb_2" in k})
        if "sg" not in self.peak2_mode or not self.use_peak2:
            fixed.update({k: True for k in fixed if "_sg_2" in k})
        if "norm" not in self.peak2_mode or not self.use_peak2:
            fixed.update({k: True for k in fixed if "_n_2" in k})

        # Peak 3 components
        if "cb" not in self.peak3_mode or not self.use_peak3:
            fixed.update({k: True for k in fixed if "_cb_3" in k})
        if "sg" not in self.peak3_mode or not self.use_peak3:
            fixed.update({k: True for k in fixed if "_sg_3" in k})
        if "norm" not in self.peak3_mode or not self.use_peak3:
            fixed.update({k: True for k in fixed if "_n_3" in k})

        # Background
        if not self.use_bkg:
            fixed["bkg"] = True
            init["bkg"] = 0.0

        return init, fixed, limits

    def extract_peak1_params(self, m: Minuit):
        """Extracts loc and scale parameters for Peak 1 from a fitted Minuit object."""
        if not m.valid:
            warnings.warn("Extracting from invalid fit.")
        if not self.use_peak1:
            raise ValueError("Peak 1 is disabled.")

        if "cb" in self.peak1_mode:
            return (
                m.values["loc_cb_1"],
                m.errors["loc_cb_1"],
                m.values["scale_cb_1"],
                m.errors["scale_cb_1"],
            )
        if "sg" in self.peak1_mode:
            return (
                m.values["loc_sg_1"],
                m.errors["loc_sg_1"],
                m.values["scale_sg_1"],
                m.errors["scale_sg_1"],
            )
        if "norm" in self.peak1_mode:
            return (
                m.values["loc_n_1"],
                m.errors["loc_n_1"],
                m.values["scale_n_1"],
                m.errors["scale_n_1"],
            )

        raise ValueError("No valid peak1_mode found.")

    def extract_peak2_params(self, m: Minuit):
        """Extracts loc and scale parameters for Peak 2 from a fitted Minuit object."""
        if not m.valid:
            warnings.warn("Extracting from invalid fit.")
        if not self.use_peak2:
            raise ValueError("Peak 2 is disabled.")

        if "cb" in self.peak2_mode:
            return (
                m.values["loc_cb_2"],
                m.errors["loc_cb_2"],
                m.values["scale_cb_2"],
                m.errors["scale_cb_2"],
            )
        if "sg" in self.peak2_mode:
            return (
                m.values["loc_sg_2"],
                m.errors["loc_sg_2"],
                m.values["scale_sg_2"],
                m.errors["scale_sg_2"],
            )
        if "norm" in self.peak2_mode:
            return (
                m.values["loc_n_2"],
                m.errors["loc_n_2"],
                m.values["scale_n_2"],
                m.errors["scale_n_2"],
            )

        raise ValueError("No valid peak2_mode found.")

    def extract_peak3_params(self, m: Minuit):
        """Extracts loc and scale parameters for Peak 3 from a fitted Minuit object."""
        if not m.valid:
            warnings.warn("Extracting from invalid fit.")
        if not self.use_peak3:
            raise ValueError("Peak 3 is disabled.")

        if "cb" in self.peak3_mode:
            return (
                m.values["loc_cb_3"],
                m.errors["loc_cb_3"],
                m.values["scale_cb_3"],
                m.errors["scale_cb_3"],
            )
        if "sg" in self.peak3_mode:
            return (
                m.values["loc_sg_3"],
                m.errors["loc_sg_3"],
                m.values["scale_sg_3"],
                m.errors["scale_sg_3"],
            )
        if "norm" in self.peak3_mode:
            return (
                m.values["loc_n_3"],
                m.errors["loc_n_3"],
                m.values["scale_n_3"],
                m.errors["scale_n_3"],
            )

        raise ValueError("No valid peak3_mode found.")


class AccumulationModel:
    def __init__(self):
        # Setup configuration here if needed in the future
        pass

    def total_model(self, x, A, tau, t0):
        """Evaluates the accumulation model.

        x: Time array
        A: Saturation amplitude
        tau: Mean lifetime
        t0: Start time of accumulation
        """
        self.exp = A * (1 - np.exp(-(x - t0) / tau))

        return self.exp

    __call__ = total_model


class ExpansionModel:
    def __init__(self):
        # Setup configuration here if needed in the future
        pass

    def total_model(self, x, A, tau, t0):
        """Evaluates the expansion model.

        x: Time array
        A: Saturation amplitude
        tau: Mean lifetime
        t0: Start time of accumulation
        """
        self.exp = A * np.exp(-(x - t0) / tau)

        return self.exp

    __call__ = total_model


def linear_model(x, m, q):
    """
    Linear model for energy calibration (Energy = m * channel + q).

    Parameters:
        x (array or float): MCA channel.
        m (float): Gain (Energy/channel).
        q (float): Offset (Energy).

    Returns:
        array or float: Calibrated energy.
    """
    return m * x + q


def res_model(E, a, b):
    """
    Inverse square root model for relative energy resolution (R = a / sqrt(E) + b).

    Parameters:
        E (array or float): Energy in MeV.
        a (float): Stochastic noise term.
        b (float): Constant noise term.

    Returns:
        array or float: Relative resolution (sigma/E or scale/loc).
    """
    return a / np.sqrt(E) + b


def calculate_fractions(
    model_instance,
    params,
    mca_range,
    x=np.arange(0, 1300 + 1, 0.1),
    flip_peaks=False,
    model_type="triple",
):
    """Calculates the primary peak fraction in an ROI and secondary peak contamination
    for a given set of parameters.

    Parameters:
    - model_type: "triple" (uses peak2 and peak3) or "double" (uses peak1 and peak2).
    - flip_peaks: If True, swaps the primary and secondary peaks in the calculation.
    """
    # Evaluate the model to populate the peak attributes
    model_instance(x, *params)

    # Apply mask
    mask = (x >= mca_range[0]) & (x <= mca_range[1])

    # Determine the peaks based on the model type
    if model_type.lower() == "double":
        peak_a = model_instance.peak1
        peak_b = model_instance.peak2
    elif model_type.lower() == "triple":
        peak_a = model_instance.peak2
        peak_b = model_instance.peak3
    else:
        raise ValueError("model_type must be either 'double' or 'triple'")

    # Swap the peaks if requested
    if flip_peaks:
        primary_peak = peak_b
        secondary_peak = peak_a
    else:
        primary_peak = peak_a
        secondary_peak = peak_b

    # Calculate sums and ratios
    sum_primary_tot = np.sum(primary_peak)
    fraction = (
        np.sum(primary_peak[mask]) / sum_primary_tot if sum_primary_tot > 0 else 0
    )

    sum_selected_tot = np.sum(primary_peak[mask] + secondary_peak[mask])
    ratio = (
        np.sum(secondary_peak[mask]) / sum_selected_tot if sum_selected_tot > 0 else 0
    )

    return fraction, ratio
