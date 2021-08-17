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
    """
    get charge database for zeroth neighbors

    Parameters
    ----------
    path : str or pathlib.Path, optional
        Defaults to internal path
    val_name : str, default="val_0"
        name of values column in output database

    Returns
    -------
    database_0 : pandas.DataFrame
        Dataframe with columns [nebr_0, ..., val_name]
    """
    if path is None:
        path = DATA_PATH / "Database_0th.txt"
    return _get_database_frame(path, val_name=val_name)


@lru_cache(maxsize=10)
def get_database_1(path=None, val_name="val_1"):
    """
    get charge database for first neighbors

    Parameters
    ----------
    path : str or pathlib.Path, optional
        Defaults to internal path
    val_name : str, default="val_1"
        name of values column in output database

    Returns
    -------
    database_1 : pandas.DataFrame
        Dataframe with columns [nebr_0, ..., val_name]


    """
    if path is None:
        path = DATA_PATH / "Database_1st.txt"
    return _get_database_frame(path, val_name=val_name)


@lru_cache(maxsize=10)
def get_database_2(path=None, val_name="val_2"):
    """
    get charge database for second neighbors

    Parameters
    ----------
    path : str or pathlib.Path, optional
        Defaults to internal path
    val_name : str, default="val_2"
        name of values column in output database

    Returns
    -------
    database_2 : pandas.DataFrame
        Dataframe with columns [nebr_0, ..., val_name]

    """
    if path is None:
        path = DATA_PATH / "Database_2nd.txt"
    return _get_database_frame(path, val_name=val_name)


def get_database(paths=None, val_names=None):

    if paths is None:
        paths = [None] * 3

    if val_names is None:
        val_names = ["val_0", "val_1", "val_2"]

    funcs = [get_database_0, get_database_1, get_database_2]

    return {
        i: func(path=path, val_name=val_name)
        for i, (path, val_name, func) in zip(paths, val_names, funcs)
    }
