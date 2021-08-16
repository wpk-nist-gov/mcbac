from pathlib import Path

import pandas as pd
import pytest

DATA_DIR = Path(__file__).absolute().parent / "test_data"


class TestData:
    def __init__(
        self,
        path_cif,
        # path_nebrs,
        # path_charge,
        # path_charge_zero,
        # path_charge_shift,
    ):

        print("name", path_cif)
        path_nebrs = "{}.dat.total.txt".format(path_cif)
        path_charge = "Final_Charge_{}.txt".format(path_cif)
        path_charge_zero = "Final_Charge_Zeroed_{}.txt".format(path_cif)
        path_charge_shift = "Shift_Charge_{}.txt".format(path_cif)

        self.path_cif = DATA_DIR / path_cif

        self.nebrs = pd.read_csv(DATA_DIR / path_nebrs, header=None, sep=r"\s+").values
        self.charge = pd.read_csv(DATA_DIR / path_charge, header=None)[0].values
        self.charge_zero = pd.read_csv(DATA_DIR / path_charge_zero, header=None)[
            0
        ].values
        self.charge_shift = pd.read_csv(DATA_DIR / path_charge_shift, header=None)[
            0
        ].values


names = [x.name for x in DATA_DIR.glob("*.cif") if not x.name.startswith("FINAL_")]
params = [(x,) for x in names]


@pytest.fixture(
    scope="module",
    #             params=[
    # ("Example.gemmi.cif", "Example.cif.dat.total.txt", "Final_Charge_Example.cif.txt", "Final_Charge_Zeroed_Example.cif.txt", "Shift_Charge_Example.cif.txt")]
    # params=[('Example.cif',), ('jacs.6b08724_ja6b08724_si_016_clean.cif',)]
    params=params,
)
def fixture(request):
    return TestData(*request.param)
