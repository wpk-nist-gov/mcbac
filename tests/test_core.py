import numpy as np

import mcbac

# import pytest


def test_pymatgen(fixture):

    c = mcbac.CifHelper.from_pymatgen(fixture.path_cif)

    np.testing.assert_array_equal(c.nebr_stack().values, fixture.nebrs)

    out = c.nebr_stack_charge()

    np.testing.assert_allclose(out["val_total"], fixture.charge)
    np.testing.assert_allclose(
        out["val_zero"], fixture.charge_zero, rtol=1e-6, atol=1e-5
    )
    np.testing.assert_allclose(
        out["val_shift"], fixture.charge_shift, rtol=1e-6, atol=1e-5
    )


def test_gemmi(fixture):
    c = mcbac.CifHelper.from_gemmi(str(fixture.path_cif))

    np.testing.assert_array_equal(c.nebr_stack().values, fixture.nebrs)

    out = c.nebr_stack_charge()

    np.testing.assert_allclose(out["val_total"], fixture.charge)
    np.testing.assert_allclose(
        out["val_zero"], fixture.charge_zero, rtol=1e-6, atol=1e-5
    )
    np.testing.assert_allclose(
        out["val_shift"], fixture.charge_shift, rtol=1e-6, atol=1e-5
    )
