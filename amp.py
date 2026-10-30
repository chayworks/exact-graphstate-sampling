"""Exact amplitudes <b| W (x) |S> of a circuit = Clifford body + per-qubit single-qubit tail."""
import re, numpy as np
from tableau import clifford_body, Tableau
from graphform import graph_form
from plan import plan

def tail_matrix(ops):
    M = np.eye(2, dtype=complex)
    for name, _ in ops:
        if name == "h": G = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        elif name == "tdg": G = np.diag([1, np.exp(-1j * np.pi / 4)])
        elif name == "t": G = np.diag([1, np.exp(1j * np.pi / 4)])
        elif name.startswith("rz("):
            th = float(re.match(r"rz\((.*)\)", name).group(1)); G = np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])
        else: raise ValueError(f"unsupported tail gate {name}")
        M = G @ M
    return M

class Amplitudes:
    def __init__(self, n, T, tail, order=None):
        self.n = n; self.A, U = graph_form(T)
        per = [[] for _ in range(n)]
        for name, qs in tail:
            assert len(qs) == 1, "tail must be single-qubit"; per[qs[0]].append((name, qs))
        self.W = [tail_matrix(per[q]) for q in range(n)]
        self.M = [self.W[q] @ U[q].conj().T for q in range(n)]
        self.order = order or list(range(n)); self.steps = plan(self.A, self.order)
    def weights(self, b):
        return [(self.M[q][b[q], 0], self.M[q][b[q], 1]) for q in range(self.n)]
