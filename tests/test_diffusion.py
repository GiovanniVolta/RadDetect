import numpy as np

import raddetect as rd
from raddetect.diffusion import Radon222Diffusion1D


def test_diffusion_imports():
    assert hasattr(rd, "diffusion")
    assert hasattr(rd, "Radon222Diffusion1D")


def test_calculate_A_and_activity_roundtrip():
    activity_in = 5.29  # Bq
    mu_x = 4.3  # nm
    sigma_x = 6.7  # nm
    lambda_h = rd.diffusion.RN222_LAMBDA_H
    S_nm2 = 1e14

    A = rd.diffusion.calculate_A_volumetric(
        activity_in, mu_x, sigma_x, lambda_h, S_nm2
    )
    activity_out = rd.diffusion.calculate_activity_from_A(
        A, mu_x, sigma_x, lambda_h, S_nm2
    )

    np.testing.assert_allclose(activity_out, activity_in, rtol=1e-10)


def test_concentration_boundary():
    A = 1e-4
    mu_x = 4.3
    sigma_x = 6.7
    D = 1.0

    # For x >= 0, concentration must be strictly 0
    x_pos = np.array([0.0, 1.0, 10.0])
    c_pos = rd.diffusion.C_x_t(x_pos, t=10.0, D=D, mu_x=mu_x, sigma_x=sigma_x, A=A)
    np.testing.assert_array_equal(c_pos, 0.0)

    # For x < 0 and t > 0, concentration should be positive
    x_neg = np.array([-1.0, -5.0, -10.0])
    c_neg = rd.diffusion.C_x_t(x_neg, t=10.0, D=D, mu_x=mu_x, sigma_x=sigma_x, A=A)
    assert np.all(c_neg > 0.0)


def test_flux_time_zero():
    A = 1e-4
    mu_x = 4.3
    sigma_x = 6.7
    D = 1.0

    # At t <= 0, outgassing flux should be 0
    assert rd.diffusion.J_diff_t_reconstructed(0.0, D, mu_x, sigma_x, A) == 0.0
    assert rd.diffusion.J_diff_t_reconstructed(-1.0, D, mu_x, sigma_x, A) == 0.0

    # At t > 0, flux should be positive
    flux_t10 = rd.diffusion.J_diff_t_reconstructed(10.0, D, mu_x, sigma_x, A)
    assert flux_t10 > 0.0


def test_radon_diffusion_1d_class():
    diff = Radon222Diffusion1D(material="Cu", activity_bq=5.29)
    assert diff.mu_z == 4.3
    assert diff.sigma_z == 6.7
    assert diff.A_init is not None

    times_h = np.linspace(1, 24, 10)
    rates = diff.eval_rate(times_h)
    assert len(rates) == 10
    assert np.all(rates > 0)


def test_minuit_fitting_object_creation():
    diff = Radon222Diffusion1D(material="SS", activity_bq=10.0)
    times_h = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    rates_hz = np.array([0.1, 0.2, 0.25, 0.3, 0.32])
    rate_err = np.array([0.01, 0.01, 0.01, 0.01, 0.01])

    m = diff.get_minuit_fitting_object(times_h, rates_hz, rate_err, D_init=1.0)
    assert m is not None
    assert "D" in m.parameters
    assert m.values["D"] == 1.0
    assert m.fixed["A"] is True
    assert m.fixed["D"] is False
