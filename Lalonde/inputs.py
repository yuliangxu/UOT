"""Public LaLonde inputs and the manuscript's sample, cost, and tuning grid."""
from hashlib import sha256
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np
from scipy.spatial.distance import cdist

DATA = Path(__file__).resolve().parent / "data"


SOURCE = "https://users.nber.org/~rdehejia/data/"
SHA256 = {
    "nsw_control.txt": "8fd745ed2c3426bb77e34b395fb84d456d346ba545af58a65b4963160b0699fd",
    "nsw_treated.txt": "ab65cd58de17a78b692e66e4d714219259ac180428f24c42ddbb928cfb1820fe",
    "psid_controls.txt": "73ffcdd8feb7f965f8e1de9bfdb80a095542c7a95ce7f56b2e2b2a9716e858cc",
    "psid2_controls.txt": "b5754be1d42e59d8e8e7067a71b68dab6e8efeb24bb654b286f18d02a809c0a1",
    "psid3_controls.txt": "1150a17122674eb3b8b54be9b122c8c1d4245acdf85cb0e8b53345496ea8600e",
}


def verify(path, expected):
    actual = sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"SHA-256 mismatch for {path}: expected {expected}, found {actual}")


def prepare_inputs():
    """Download missing files and verify the bytes of every input before fitting."""
    DATA.mkdir(parents=True, exist_ok=True)
    for name, expected in SHA256.items():
        target = DATA / name
        if not target.exists():
            with TemporaryDirectory(dir=DATA) as temp:
                pending = Path(temp) / name
                with urlopen(SOURCE + name, timeout=60) as response:
                    pending.write_bytes(response.read())
                verify(pending, expected)
                # Publish only complete, validated bytes; concurrent runs may share this cache.
                try:
                    os.link(pending, target)
                except FileExistsError:
                    pass
        verify(target, expected)


def load_inputs():
    treated = np.loadtxt(DATA / "nsw_treated.txt")
    randomized_controls = np.loadtxt(DATA / "nsw_control.txt")
    psid = np.vstack([
        np.loadtxt(DATA / name) for name in
        ("psid_controls.txt", "psid2_controls.txt", "psid3_controls.txt")
    ])
    controls = np.delete(psid, 7, axis=1)  # NSW files omit re74
    n, m = len(treated), len(controls)
    predictors = np.vstack([treated[:, 1:-1], controls[:, 1:-1]])
    predictors = (predictors - predictors.mean(axis=0)) / predictors.std(axis=0, ddof=1)
    cost = cdist(predictors[:n], predictors[n:], metric="euclidean")
    cost /= cost.max()
    benchmark = treated[:, -1].mean() - randomized_controls[:, -1].mean()
    return cost, np.full(n, 1/n), np.full(m, 1/m), treated[:, -1], controls[:, -1], benchmark


def grids():
    def grid(low, high, size):
        return [float(f"{v:.3g}") for v in np.logspace(low, high, size)]
    return grid(-4, -2, 15), grid(-4, 1, 20)
