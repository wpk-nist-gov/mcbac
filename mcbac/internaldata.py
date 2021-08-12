from functools import lru_cache
from pathlib import Path

import pandas as pd

DATA_PATH = Path(__file__).absolute().parent / "data"


@lru_cache(maxsize=10)
def get_radii(path=None):
    import json

    if path is None:
        path = DATA_PATH / "radii.json"
    with open(path, "r") as f:
        out = json.load(f)

    return out


def _get_database_frame(path, prefix="nebr_", val_name="val"):
    df = pd.read_csv(path, sep=r"\s+", header=None)
    if prefix is not None:
        df = df.rename(columns=lambda x: "{}{}".format(prefix, x))

    if val_name is not None:
        df.columns = list(df.columns[:-1]) + [val_name]
    return df


@lru_cache(maxsize=10)
def get_database_0(path=None, val_name="val_0"):
    if path is None:
        path = DATA_PATH / "Database_0th.txt"
    return _get_database_frame(path, val_name=val_name)


@lru_cache(maxsize=10)
def get_database_1(path=None, val_name="val_1"):
    if path is None:
        path = DATA_PATH / "Database_1st.txt"
    return _get_database_frame(path, val_name=val_name)


@lru_cache(maxsize=10)
def get_database_2(path=None, val_name="val_2"):
    if path is None:
        path = DATA_PATH / "Database_2nd.txt"
    return _get_database_frame(path, val_name=val_name)
