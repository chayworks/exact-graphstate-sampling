"""Exact probabilities and exact samples for IQP circuits with one- and two-qubit generators on a sparse graph.

State: |psi> = exp(i sum_j theta_j X_{g_j}) |0^n>, measured in the Z basis (the IQPopt convention; since X is real and
|0> is real, theta -> -theta conjugates the amplitudes and leaves the distribution unchanged). With H X H = Z,
    psi(y) = 2^-n sum_x (-1)^{x.y} c(x),   c(x) = prod_j exp(i theta_j (-1)^{g_j . x}).
For generators of weight 1 and 2 on a graph, c is a product of vertex and edge factors, so psi is a matrix product
over a vertex order whose bond at step k is the frontier F_k (placed vertices with an unplaced neighbour); the bond
dimension is 2^|F_k|. Summing |psi|^2 over an output bit y_v forces x_v = x'_v, so marginals need only a right
environment E_k[r, r'] of size 4^|F_k| that is updated by a gather, and sampling is the chain rule left to right.
Cost: n 4^w for the environments and n 4^w per sample, w = max |F_k|; memory 4^w complex numbers.

Usage as a module: model = IQPModel(n, gates, thetas, order); model.prob(y); model.sample(rng, count).
Run as a script for the self-test against a dense statevector (python iqp_exact.py --selftest)."""
import sys, numpy as np

class IQPModel:
    def __init__(self, n, gates, thetas, order=None):
        """gates: list of tuples of qubit indices (weight 1 or 2); thetas: matching angles. Repeated generators add."""
        self.n = n; order = list(range(n)) if order is None else list(order); assert sorted(order) == list(range(n))
        h = np.zeros(n); J = {}
        for g, t in zip(gates, thetas):
            g = tuple(sorted(int(q) for q in g))
            if len(g) == 1: h[g[0]] += t
            elif len(g) == 2: J[g] = J.get(g, 0.0) + t
            else: raise ValueError("only weight-1 and weight-2 generators")
        self.h, self.J, self.order = h, J, order
        nb = [set() for _ in range(n)]
        for (a, b) in J: nb[a].add(b); nb[b].add(a)
        pos = {v: k for k, v in enumerate(order)}
        # frontier after placing step k: placed vertices with an unplaced neighbour, kept in a fixed list order
        self.F = []; cur = []
        for k, v in enumerate(order):
            cur = [u for u in cur + [v] if any(pos[w] > k for w in nb[u])]
            self.F.append(cur)
        self.width = max(len(f) for f in self.F)
        # site tables: for site k, left bond F_{k-1}, right bond F_k; for each left state l and x in {0,1}:
        # weight W[l, x] and the right state R[l, x]
        self.W, self.R = [], []
        prev = []
        for k, v in enumerate(order):
            Fl, Fr = prev, self.F[k]; nl = 1 << len(Fl); L = np.arange(nl)
            bits = {u: (L >> i) & 1 for i, u in enumerate(Fl)}
            W = np.empty((nl, 2), complex); R = np.zeros((nl, 2), np.int64)
            for x in (0, 1):
                ph = h[v] * (1 - 2 * x) * np.ones(nl)
                for u in Fl:
                    e = (min(u, v), max(u, v))
                    if e in J: ph = ph + J[e] * (1 - 2 * (bits[u] ^ x))
                W[:, x] = np.exp(1j * ph)
                r = np.zeros(nl, np.int64)
                for i, u in enumerate(Fr):
                    b = (x * np.ones(nl, np.int64)) if u == v else bits[u]
                    r |= b.astype(np.int64) << i
                R[:, x] = r
            self.W.append(W); self.R.append(R); prev = Fr
        self._env = None

    def amplitude(self, y_by_qubit):
        """psi(y) for y indexed by qubit."""
        Lv = np.ones(1, complex)
        for k, v in enumerate(self.order):
            W, R = self.W[k], self.R[k]; nr = 1 << len(self.F[k]); out = np.zeros(nr, complex); yv = int(y_by_qubit[v])
            for x in (0, 1):
                np.add.at(out, R[:, x], Lv * W[:, x] * (0.5 * (1 - 2 * (x & yv))))
            Lv = out
        return Lv[0]

    def prob(self, y_by_qubit):
        return abs(self.amplitude(y_by_qubit)) ** 2

    def envs(self):
        """E_k for k = 0..n: E_n = [[1]]; E_{k-1}[l, l'] = 1/2 sum_x W[l,x] conj(W[l',x]) E_k[R[l,x], R[l',x]]."""
        if self._env is not None: return self._env
        E = [None] * (self.n + 1); E[self.n] = np.ones((1, 1), complex)
        for k in range(self.n - 1, -1, -1):
            W, R, Ek = self.W[k], self.R[k], E[k + 1]
            acc = 0
            for x in (0, 1):
                acc = acc + 0.5 * np.outer(W[:, x], W[:, x].conj()) * Ek[np.ix_(R[:, x], R[:, x])]
            E[k] = acc
        self._env = E; return E

    def sample(self, rng, count):
        """Exact samples (chain rule over the order), returned as strings indexed by qubit, with their probabilities."""
        E = self.envs(); out = []
        for _ in range(count):
            Lv = np.ones(1, complex); y = np.zeros(self.n, np.int8)
            for k, v in enumerate(self.order):
                W, R = self.W[k], self.R[k]; nr = 1 << len(self.F[k]); vs = []
                for yv in (0, 1):
                    o = np.zeros(nr, complex)
                    for x in (0, 1):
                        np.add.at(o, R[:, x], Lv * W[:, x] * (0.5 * (1 - 2 * (x & yv))))
                    vs.append(o)
                p = np.array([np.real(vv @ E[k + 1] @ vv.conj()) for vv in vs]); p = np.maximum(p, 0)
                yv = int(rng.random() * p.sum() >= p[0]); y[v] = yv; Lv = vs[yv]
            out.append((y.copy(), abs(Lv[0]) ** 2))
        return out

def dense_probs(n, gates, thetas):
    """Brute force: build exp(i sum theta X_g)|0> by diagonalising in the X basis (H^n), return Born probs by qubit."""
    N = 1 << n; x = np.arange(N); ph = np.zeros(N)
    for g, t in zip(gates, thetas):
        par = np.zeros(N, np.int64)
        for q in g: par ^= (x >> q) & 1
        ph += t * (1 - 2 * par)
    c = np.exp(1j * ph)                      # coefficients in the x basis of H^n|psi>, up to normalisation
    # psi(y) = 2^-n sum_x (-1)^{x.y} c(x): fast Walsh-Hadamard transform
    a = c.astype(complex).copy(); hlen = 1
    while hlen < N:
        a = a.reshape(-1, 2, hlen); a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], 1).reshape(-1); hlen *= 2
    return np.abs(a / N) ** 2                # index y has bit q = qubit q

def statevector_probs(n, gates, thetas):
    """Independent check of dense_probs: apply each exp(i t X_g) as a matrix exponential on the statevector."""
    N = 1 << n; psi = np.zeros(N, complex); psi[0] = 1; idx = np.arange(N)
    for g, t in zip(gates, thetas):
        mask = 0
        for q in g: mask |= 1 << q
        psi = np.cos(t) * psi + 1j * np.sin(t) * psi[idx ^ mask]          # exp(i t P) = cos t + i sin t P, P = X_g
    return np.abs(psi) ** 2

def _selftest():
    rng = np.random.default_rng(7); worst = 0.0
    for trial in range(40):
        n = int(rng.integers(3, 12)); edges = set()
        for _ in range(int(rng.integers(n - 1, 2 * n))):
            a, b = rng.choice(n, 2, replace=False); edges.add((int(min(a, b)), int(max(a, b))))
        gates = [(i,) for i in range(n)] + sorted(edges); thetas = rng.uniform(-np.pi, np.pi, len(gates))
        if trial % 3 == 0: gates = gates + gates[:3]; thetas = np.concatenate([thetas, rng.uniform(-1, 1, 3)])  # repeats add
        order = list(rng.permutation(n))
        P1 = dense_probs(n, gates, thetas); P2 = statevector_probs(n, gates, thetas)
        assert abs(P1.sum() - 1) < 1e-12 and np.max(np.abs(P1 - P2)) < 1e-12, "dense references disagree"
        m = IQPModel(n, gates, thetas, order)
        Pm = np.array([m.prob([(y >> q) & 1 for q in range(n)]) for y in range(1 << n)])
        worst = max(worst, np.max(np.abs(Pm - P1)))
        E0 = m.envs()[0]; assert abs(E0[0, 0] - 1) < 1e-10, E0          # total probability through the environments
    print(f"selftest probabilities: 40 random graphs/orders up to 11 qubits, max |dP| = {worst:.2e}")
    assert worst < 1e-12
    # control: mis-wiring one edge must be detected (the check is not vacuous)
    n = 6; gates = [(i,) for i in range(n)] + [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)]
    thetas = rng.uniform(-np.pi, np.pi, len(gates)); P1 = dense_probs(n, gates, thetas)
    bad = IQPModel(n, gates[:-1] + [(0, 4)], thetas)
    db = np.max(np.abs(np.array([bad.prob([(y >> q) & 1 for q in range(n)]) for y in range(1 << n)]) - P1))
    print(f"control (one edge moved): max |dP| = {db:.2e}"); assert db > 1e-3
    # sampler: TV distance of 40,000 samples against the exact law, vs the TV of 40,000 ideal i.i.d. draws
    n = 8; gates = [(i,) for i in range(n)] + [(i, i + 1) for i in range(n - 1)] + [(0, 4), (2, 7)]
    thetas = rng.uniform(-np.pi, np.pi, len(gates)); P = dense_probs(n, gates, thetas)
    m = IQPModel(n, gates, thetas, [3, 1, 0, 2, 4, 5, 7, 6]); S = m.sample(rng, 40000)
    idx = np.array([sum(int(b) << q for q, b in enumerate(y)) for y, _ in S])
    emp = np.bincount(idx, minlength=1 << n) / len(idx); tv = 0.5 * np.abs(emp - P).sum()
    ideal = np.bincount(rng.choice(1 << n, 40000, p=P / P.sum()), minlength=1 << n) / 40000; tvi = 0.5 * np.abs(ideal - P).sum()
    pr = max(abs(p - P[i]) for (y, p), i in zip(S[:200], idx[:200]))
    print(f"sampler: TV {tv:.4f} vs ideal i.i.d. draws {tvi:.4f}; reported probabilities match to {pr:.1e}")
    assert tv < 1.5 * tvi + 0.005 and pr < 1e-12
    print("selftest PASS")

if __name__ == "__main__":
    if "--selftest" in sys.argv: _selftest()
