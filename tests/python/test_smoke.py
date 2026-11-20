"""Smoke tests of the Python module multiax (C++ core)."""
import multiax


def test_build_info():
    info = multiax.build_info()
    assert info['cplusplus'] >= 202002
    assert info['openmp_version'] > 0
    assert info['max_threads'] >= 1
    assert info['version'] == multiax.__version__


def test_woehler_curve_parameters():
    wc = multiax.WoehlerCurve(k=3.0, k_star=22.0, L_k=109.19, N_k=5e5)
    assert (wc.k, wc.k_star, wc.L_k, wc.N_k) == (3.0, 22.0, 109.19, 5e5)
    assert 'WoehlerCurve(' in repr(wc)
