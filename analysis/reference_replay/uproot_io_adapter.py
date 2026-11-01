"""Explicit minimal uproot I/O adapter for the bundled, unchanged legacy scripts.

This is NOT ROOT. It implements only the exact reader calls used by the frozen
scripts. All selection, normalization, flux, binning, uncertainty and plotting
algorithms remain in those scripts. Unknown attributes/calls raise normally.
"""
from functools import lru_cache
import re
import types
import numpy as np
import uproot


@lru_cache(maxsize=8)
def arrays(path, tree):
    with uproot.open(path) as f:
        return f[tree].arrays(library="np")


class _Tree:
    def __init__(self, tree):
        self.tree = tree

    def GetBranch(self, name):
        return True if name in self.tree.keys() else None


class _File:
    def __init__(self, path):
        self.file = uproot.open(path)

    def IsZombie(self):
        return False

    def Get(self, name):
        return _Tree(self.file[name]) if name in self.file else None

    def Close(self):
        self.file.close()


class _Value:
    def __init__(self, value):
        self.value = value

    def GetValue(self):
        return self.value


class RDataFrame:
    def __init__(self, tree, path, mask=None):
        self.tree, self.path = tree, str(path)
        self.data = arrays(self.path, tree)
        self.mask = np.ones(len(next(iter(self.data.values()))), dtype=bool) if mask is None else mask

    def Filter(self, expression):
        # Only conjunctions of scalar comparisons and abs(column-constant) are
        # accepted. The expressions originate in the byte-preserved scripts.
        mask = self.mask.copy()
        for clause in expression.split("&&"):
            if not re.fullmatch(r"[A-Za-z0-9_\s.+\-<>=()]+", clause):
                raise ValueError(f"Unsupported legacy filter: {expression}")
            mask &= eval(clause.strip(), {"__builtins__": {}, "abs": np.abs}, self.data)
        return RDataFrame(self.tree, self.path, mask)

    def Sum(self, column):
        return _Value(float(np.sum(self.data[column][self.mask], dtype=np.float64)))

    def Count(self):
        return _Value(int(np.count_nonzero(self.mask)))

    def AsNumpy(self, columns):
        return {name: np.array(self.data[name][self.mask], copy=True) for name in columns}


def module():
    result = types.ModuleType("ROOT")
    result.__doc__ = __doc__
    result.RDataFrame = RDataFrame
    result.TFile = types.SimpleNamespace(Open=_File)
    result.__version__ = "EXPLICIT_UPROOT_IO_ADAPTER_NOT_ROOT"
    return result
