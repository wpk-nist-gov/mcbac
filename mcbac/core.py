import itertools

import numpy as np
import pandas as pd
import xarray as xr

from .cached_decorators import gcached
from .internaldata import get_database_0, get_database_1, get_database_2, get_radii


# Decorators for class interaction
def _prop_wrapper(key):
    @property
    def new_prop(self):
        return self.data[key]

    return new_prop


def _prop_wrapper_float(key):
    @property
    def new_prop(self):
        return float(self.data[key])

    return new_prop


def _prop_wrapper_array(key, dtype=None):
    @property
    def new_prop(self):
        return np.array(self.data[key], dtype=dtype)

    return new_prop


def add_properties(float_mapping=None, array_mapping=None, mapping=None):
    if mapping is None:
        mapping = {}
    if float_mapping is None:
        float_mapping = {}
    if array_mapping is None:
        array_mapping = {}

    def wrapper(cls):
        for name, key in mapping.items():
            setattr(cls, name, _prop_wrapper(key))
        for name, key in float_mapping.items():
            setattr(cls, name, _prop_wrapper_float(key))
        for name, key_dtype in array_mapping.items():
            if isinstance(key_dtype, (tuple, list)):
                key, dtype = key_dtype
            else:
                key = key_dtype
                dtype = None
            setattr(cls, name, _prop_wrapper_array(key, dtype))

        return cls

    return wrapper


@add_properties(
    float_mapping={
        "a": "_cell_length_a",
        "b": "_cell_length_b",
        "c": "_cell_length_c",
        "alpha": "_cell_angle_alpha",
        "beta": "_cell_angle_beta",
        "gamma": "_cell_angle_gamma",
    },
    array_mapping={
        "label": "_atom_site_label",
        "frac_x": ("_atom_site_fract_x", float),
        "frac_y": ("_atom_site_fract_y", float),
        "frac_z": ("_atom_site_fract_z", float),
        "charge": ("_atom_site_charge", float),
    },
)
class CifHelper:
    """
    Helper class to work with cif file

    Paramters
    ---------
    data : dict
        dictionary from cif file
    radii_dict : dict, optional
        mapping from elements -> radii

    Attributes
    ----------
    a, b, c : float
        values from data['_cell_length_a'] etc
    alpha, beta, gamma : float
        value from data['_cell_angle_alpha], etc
    label : array of strings
        value from data['_atom_site_label]
    frac_x, frac_y, frac_z, charge : array of floats
        value from data['_atom_site_fract_x'], etc

    """

    def __init__(self, data, radii_dict=None):
        if radii_dict is None:
            radii_dict = get_radii()
        self._radii_dict = radii_dict
        self.data = data

    def __getitem__(self, index):
        return self.data[index]

    @property
    def alpha_rad(self):
        return np.deg2rad(self.alpha)

    @property
    def beta_rad(self):
        return np.deg2rad(self.beta)

    @property
    def gamma_rad(self):
        return np.deg2rad(self.gamma)

    @property
    def volume(self):
        cos_alpha, cos_beta, cos_gamma = [
            np.cos(x) for x in (self.alpha_rad, self.beta_rad, self.gamma_rad)
        ]
        return (
            self.a
            * self.b
            * self.c
            * (
                1
                - cos_alpha ** 2
                - cos_beta ** 2
                - cos_gamma ** 2
                + 2 * cos_alpha * cos_beta * cos_gamma
            )
            ** 0.5
        )

    @gcached()
    def atom_radii(self):
        return np.array([self._radii_dict[k] for k in self.label])

    @property
    def frac_all(self):
        return np.stack((self.frac_x, self.frac_y, self.frac_z), axis=-1)

    @property
    def transform_matrix(self):
        """transform to cartisian coords using frac_coords .dot. transform_matrix.T"""
        a1 = self.a
        b1 = self.b * np.cos(self.gamma_rad)
        c1 = self.c * np.cos(self.beta_rad)

        a2 = 0.0
        b2 = self.b * np.sin(self.gamma_rad)
        c2 = self.c * (
            (np.cos(self.alpha_rad) - np.cos(self.beta_rad) * np.cos(self.gamma_rad))
            / np.sin(self.gamma_rad)
        )

        a3 = 0.0
        b3 = 0.0
        c3 = self.volume / (self.a * self.b * np.sin(self.gamma_rad))

        return np.array([[a1, b1, c1], [a2, b2, c2], [a3, b3, c3]])

    @property
    def coords(self):
        return np.dot(self.frac_all, self.transform_matrix.T)

    def coords_replicated(self, stack=True):
        """
        create coordinates for replicate.

        Parameters
        ----------
        stack : bool, default = True
            If `False`, return all coordinates, and output shape = (27, n, 3)
            If `True`, stack all replicates into first dimension ans shape = (27 * n, 3)
        """

        offsets = np.array(list(itertools.product(*(range(-1, 2),) * 3)))
        out = np.dot(
            self.frac_all[None, :, :] + offsets[:, None, :], self.transform_matrix.T
        )
        if stack:
            out = out.reshape(-1, 3)
        return out

    def _coords_replicated_dumb(self):
        """for testing"""
        a1 = self.a
        b1 = self.b * np.cos(self.gamma_rad)
        c1 = self.c * np.cos(self.beta_rad)

        a2 = 0.0
        b2 = self.b * np.sin(self.gamma_rad)
        c2 = self.c * (
            (np.cos(self.alpha_rad) - np.cos(self.beta_rad) * np.cos(self.gamma_rad))
            / np.sin(self.gamma_rad)
        )

        a3 = 0.0
        b3 = 0.0
        c3 = self.volume / (self.a * self.b * np.sin(self.gamma_rad))

        X = []
        Y = []
        Z = []

        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):

                    for xx, yy, zz in zip(self.frac_x, self.frac_y, self.frac_z):

                        X.append((xx + i) * a1 + (yy + j) * b1 + (zz + k) * c1)
                        Y.append((xx + i) * a2 + (yy + j) * b2 + (zz + k) * c2)
                        Z.append((xx + i) * a3 + (yy + j) * b3 + (zz + k) * c3)

        return np.stack([X, Y, Z], axis=-1)

    def _coords_dumb(self):
        """for testing"""
        a1 = self.a
        b1 = self.b * np.cos(self.gamma_rad)
        c1 = self.c * np.cos(self.beta_rad)

        a2 = 0.0
        b2 = self.b * np.sin(self.gamma_rad)
        c2 = self.c * (
            (np.cos(self.alpha_rad) - np.cos(self.beta_rad) * np.cos(self.gamma_rad))
            / np.sin(self.gamma_rad)
        )

        a3 = 0.0
        b3 = 0.0
        c3 = self.volume / (self.a * self.b * np.sin(self.gamma_rad))

        X = []
        Y = []
        Z = []

        for xx, yy, zz in zip(self.frac_x, self.frac_y, self.frac_z):

            X.append((xx) * a1 + (yy) * b1 + (zz) * c1)
            Y.append((xx) * a2 + (yy) * b2 + (zz) * c2)
            Z.append((xx) * a3 + (yy) * b3 + (zz) * c3)

        return np.stack([X, Y, Z], axis=-1)

    def __len__(self):
        return len(self.frac_x)

    def distance_matrix(self, stack=True):
        """
        compute distance matrix

        Paramters
        ----------
        stack : bool, default = True
            Whether to stack or not stack the replicate dimension

        Returns
        -------
        distance_matrix : array
            if not stack, then shape = (n, 27, n)
            and distance_matrix[c, r, n] = dist(coords[c,:], coords_repl[r, n, :]).
            if stack, then shape = (n, 27 * n) and
            distance_matrix[c, i] = dist(coords[c, :], coords_repl[i, :]).
        """

        from scipy.spatial.distance import cdist

        coords_cntr = self.coords
        coords_repl = self.coords_replicated(stack=True)

        distances = cdist(coords_cntr, coords_repl)

        if not stack:
            distances = distances.reshape(len(self), 27, len(self))
        return distances

    @property
    def distance_dataarray(self) -> xr.DataArray:

        radii = self.atom_radii
        label = self.label

        da = xr.DataArray(
            self.distance_matrix(stack=False),
            dims=["cntr", "repl", "nebr"],
            coords={
                "cntr_radii": ("cntr", radii),
                "nebr_radii": ("nebr", radii),
                "cntr_name": ("cntr", label),
                "nebr_name": ("nebr", label),
            },
            name="distance",
        )
        return da

    @gcached()
    def distance_frame(self) -> pd.DataFrame:
        return self.distance_dataarray.to_dataframe()

    @gcached(prop=False)
    def nebr_frame(self, rcut_fac=1.25):
        """neighbor frame"""
        df = self.distance_frame

        if rcut_fac is not None:
            df = df.assign(
                rcut=lambda x: (x["cntr_radii"] + x["nebr_radii"]) * rcut_fac
            ).query("distance < rcut")

        df = (
            df.reset_index()
            .sort_values(["cntr", "distance"])
            # assign ranks
            .groupby("cntr", as_index=False)
            .apply(lambda x: x.assign(rank=lambda x: np.arange(len(x))))
            .reset_index(drop=True)
        )
        return df

    @gcached(prop=False)
    def nebr_nebr_frame(self, rcut_fac=1.25):
        """neighbors of neighbors"""
        nebr_frame = self.nebr_frame(rcut_fac=rcut_fac)
        df = nebr_frame.query("rank > 0")
        out = pd.merge(
            df, df, left_on="nebr", right_on="cntr", how="left", suffixes=["", "_2nd"]
        )
        return out

    @gcached(prop=False)
    def nebr_stack(self, rcut_fac=1.25):
        """list of neighbors to use as keys in database"""
        nebr_frame = self.nebr_frame(rcut_fac=rcut_fac)
        nebr_nebr_frame = self.nebr_nebr_frame(rcut_fac=rcut_fac)

        def _sum_names(df, col):
            return (
                df[["cntr", col]].sort_values(["cntr", col]).groupby("cntr")[col].sum()
            )

        # 0th
        df_0 = pd.DataFrame(
            {"cntr": np.arange(len(self)), "nebr_0": self.label}
        ).set_index("cntr")["nebr_0"]
        df_1 = nebr_frame.query("cntr != nebr").pipe(_sum_names, col="nebr_name")

        # TODO : consider changing this?
        # this includes self and double counts
        # i.e., if nebr(0) -> [1, 2, 3], nebr(1) -> [0, 2], nebr(3) -> (0,1)
        # then nebr_of_nebr[0] -> nebr(1) + nebr(2) + nebr(3)  -> [0,2, 0, 1]
        df_2 = nebr_nebr_frame.query(
            "cntr != nebr and cntr != cntr_2nd"
        ).pipe(  # and cntr != nebr_2nd')
            _sum_names, col="nebr_name_2nd"
        )
        return pd.concat({"nebr_0": df_0, "nebr_1": df_1, "nebr_2": df_2}, axis=1)

    def nebr_stack_charge(
        self,
        rcut_fac=1.25,
        db_0=None,
        db_1=None,
        db_2=None,
        val_total="val_total",
        val_shift="val_shift",
        val_zero="val_zero",
        nebr_stack=None,
        fill_missing="offset",
    ):
        """
        merge charge values into `self.nebr_stack()`

        This finds the highest order value

        Paramters
        ---------
        db_0, db_1, db_2 : pandas.DataFrame
            db_i is the  i-th neighbor database.  It should have colums `[nebr_0, ..., nebr_i, val_i]`
        fill_missing :
        """

        if nebr_stack is None:
            nebr_stack = self.nebr_stack(rcut_fac=rcut_fac)

        if db_0 is None:
            db_0 = get_database_0()
        if db_1 is None:
            db_1 = get_database_1()
        if db_2 is None:
            db_2 = get_database_2()

        out = (
            nebr_stack.merge(db_0, how="left")
            .merge(db_1, how="left")
            .merge(db_2, how="left")
        )

        # pick out the highest order non nan value
        out.loc[:, val_total] = out[["val_0", "val_1", "val_2"]].fillna(
            method="ffill", axis=1
        )["val_2"]

        if fill_missing == "offset":
            # fill with charge to get zero total charge
            # sum of current charges / number of missing charges
            n_missing = out[val_total].isnull().sum()
            if n_missing > 0:
                value = -out[val_total].sum() / out[val_total].isnull().sum()
                out.loc[:, val_total] = out[val_total].fillna(value=value)

        if val_shift is not None or val_zero is not None:
            sum_abs = np.abs(out["val_total"]).sum()
            sum_real = out["val_total"].sum()
            fac = sum_real / sum_abs

            if val_shift is not None:
                out.loc[:, val_shift] = np.abs(out[val_total]) * fac
            if val_zero is not None:
                out.loc[:, val_zero] = out["val_total"] - np.abs(out[val_total]) * fac

        return out

    @classmethod
    def from_dict(cls, d, **kws):
        """
        for dict of multiple cif things {lable : {cif stuff}, ....}
        create a dict of CifHelpers
        """
        return {k: cls(data, **kws) for k, data in d.items()}

    @classmethod
    def from_pymatgen(cls, path, parser_kws=None, **kws):
        from pymatgen.io.cif import CifParser

        if parser_kws is None:
            parser_kws = {}
        p = CifParser(path, **parser_kws)

        return cls.from_dict(p.as_dict(), **kws)

    @classmethod
    def from_gemmi(cls, path, json_kws=None, **kws):
        import json

        from gemmi import cif

        p = cif.read_file(path)
        if json_kws is None:
            json_kws = {}

        d = json.loads(p.as_json())

        return cls.from_dict(d, **kws)
