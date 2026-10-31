"""In-place cut-rank DP: memory is one array of 2^R entries (R = max cut rank).

Each step maps old coordinates (c, x_v) to new ones by a GF(2)-linear map Lam with kernel K. Let Q be the
span of the old-coordinate parts of K, and pick coset representatives c with the Q-pivot bits zero. Every
old coset c+Q (2^q entries) feeds m = 2^{q+1-d} new entries (d = dim K), and no other coset feeds them.
So each coset is a butterfly that may overwrite its own slots.

The physical address of a logical index is a linear map T (R x r). New logical coordinates are chosen as
(non-Q-pivot old coordinates, then tau), where tau selects which of the m outputs. T_new keeps the old
columns for the kept coordinates and maps tau to g: a vector T(e) with e in Q when q >= 1, or a free
physical bit when q = 0 (growth).
"""

def _reduce(vecs):
    """Row-echelon basis of span(vecs) (python ints); returns list of (vec, pivot_bit)."""
    basis = []
    for v in vecs:
        for b, p in basis:
            if (v >> p) & 1: v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(b ^ v if (b >> p) & 1 else b, pb) for b, pb in basis]
            basis.append((v, p))
    return basis

def _in_span(v, basis):
    for b, p in basis:
        if (v >> p) & 1: v ^= b
    return v == 0

def plan_inplace(A, order, R):
    n = len(A); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]
    B = []; T = []; done = 0; steps = []
    for v in order:
        r = len(B); vb = 1 << v; done |= vb
        rho = sum(1 << j for j, b in enumerate(B) if b & vb)
        gens = [b & ~vb for b in B] + [nb[v] & ~done]            # index r is x_v
        # kernel of Lam: combinations of gens that vanish
        red = []; K = []
        for i, g in enumerate(gens):
            comb = 1 << i; t = g
            for (bv, bc) in red:
                if (t >> (bv.bit_length() - 1)) & 1: t ^= bv; comb ^= bc
            if t: red.append((t, comb)); red.sort(key=lambda p: -p[0])
            else: K.append(comb)
        low = (1 << r) - 1
        Qb = _reduce([k & low for k in K if k & low])            # basis of Q with pivot bits
        q = len(Qb); qpiv = [p for _, p in Qb]
        keep = [i for i in range(r) if i not in qpiv]             # coset-rep coordinates (non-Q-pivot)
        # enumerate Q x {0,1}: element (qv, x) -> u = qv | x<<r ; group by image under Lam (= same coset of K)
        Qall = [0]
        for b, _ in Qb: Qall += [a ^ b for a in Qall]
        Kspan = [0]
        for k in K: Kspan += [a ^ k for a in Kspan]
        Kset = set(Kspan)
        reps = []; terms = []                                     # terms: (qv, x, tau)
        for qv in Qall:
            for x in (0, 1):
                u = qv | (x << r); tau = None
                for ti, rep in enumerate(reps):
                    if (u ^ rep) in Kset: tau = ti; break
                if tau is None: reps.append(u); tau = len(reps) - 1
                terms.append((qv, x, tau))
        m = len(reps); assert m in (1, 2) and reps[0] == 0
        # new basis: field vectors of kept coordinates, plus that of reps[1] if m == 2
        def field(u):
            f = 0
            for i in range(r + 1):
                if (u >> i) & 1: f ^= gens[i]
            return f
        newB = [gens[i] for i in keep] + ([field(reps[1])] if m == 2 else [])
        # physical map: q >= 1 -> second output goes to another slot of the same coset (T of a Q basis vector);
        # q == 0 (growth) -> lowest physical bit outside span(T), which exists because rn <= R.
        if m == 2:
            if q >= 1:
                e = Qb[0][0]; g = 0
                for i in range(r):
                    if (e >> i) & 1: g ^= T[i]
            else:
                tb = _reduce(list(T)); g = next(1 << j for j in range(R) if not _in_span(1 << j, tb))
            newT = [T[i] for i in keep] + [g]
        else:
            newT = [T[i] for i in keep]
        assert all(t < (1 << R) for t in newT) and len(_reduce(list(newT))) == len(newT), "layout not injective"
        # sanity: new basis independent and spanning the true field space
        nbasis = _reduce(newB); assert len(nbasis) == len(newB)
        assert all(_in_span(field(u), nbasis) for u in [1 << i for i in range(r + 1)])
        # physical address of slot (rep c in keep-coords, qv): T(c) ^ T(qv)
        Tq = {}
        for qv in Qall:
            a = 0
            for i in range(r):
                if (qv >> i) & 1: a ^= T[i]
            Tq[qv] = a
        steps.append(dict(v=v, r=r, rn=len(newB), keep=keep, Tkeep=[T[i] for i in keep], rho=rho, terms=terms,
                          Tq=Tq, m=m, g=(newT[-1] if m == 2 else None)))
        B, T = newB, newT
    assert not B
    return steps

def sample_chain_probs(steps, R, M, bits):
    """Exact conditional probabilities along a forced bit path, with the in-place layout (reference)."""
    import numpy as np
    mem = np.zeros(1 << R, complex); mem[0] = 1.0; out = []
    for st, bv in zip(steps, bits):
        r, keep, Tk, rho, terms, Tq, m, g = st["r"], st["keep"], st["Tkeep"], st["rho"], st["terms"], st["Tq"], st["m"], st["g"]
        Mv = M[st["v"]]; nk = len(keep)
        def outputs(c, b):
            base = 0
            for j in range(nk):
                if (c >> j) & 1: base ^= Tk[j]
            cl = 0
            for j in range(nk):
                if (c >> j) & 1: cl |= 1 << keep[j]
            vals = [0j] * m
            for qv, x, tau in terms:
                full = cl ^ qv; sgn = -1 if (x and bin(full & rho).count("1") & 1) else 1
                vals[tau] += Mv[b, x] * sgn * mem[base ^ Tq[qv]]
            return base, vals
        nrm = [0.0, 0.0]
        for c in range(1 << nk):
            for b in (0, 1):
                _, vals = outputs(c, b); nrm[b] += sum(abs(z) ** 2 for z in vals)
        p1 = nrm[1] / (nrm[0] + nrm[1]); out.append(p1 if bv else 1 - p1)
        new_writes = []; written = set(); oldslots = set()
        for c in range(1 << nk):
            base, vals = outputs(c, bv); new_writes.append((base, vals))
            own = {base ^ a for a in Tq.values()}; oldslots |= own
            tgt = [base] + ([base ^ g] if m == 2 else [])
            for a in tgt:
                assert a not in written, "two cosets write the same address"
                written.add(a)
        for base, vals in new_writes:
            own = {base ^ a for a in Tq.values()}
            tgt = [base] + ([base ^ g] if m == 2 else [])
            assert all(a in own or a not in oldslots for a in tgt), "write would clobber another coset's input"
        for base, vals in new_writes:                              # (reference: compute all, then write)
            mem[base] = vals[0] / np.sqrt(nrm[bv])
            if m == 2: mem[base ^ g] = vals[1] / np.sqrt(nrm[bv])
    return out

def write_inplace_plan(path, steps, R, M, n):
    def pext(x, mask):
        out = 0; j = 0
        for i in range(64):
            if (mask >> i) & 1:
                if (x >> i) & 1: out |= 1 << j
                j += 1
        return out
    with open(path, "w") as f:
        f.write(f"{n} {len(steps)} {R}\n")
        for s in steps:
            keepmask = sum(1 << i for i in s["keep"]); rhoc = pext(s["rho"], keepmask)
            Mv = M[s["v"]]
            f.write(f"{s['v']} {len(s['keep'])} {s['m']} {s['g'] or 0} {rhoc} {len(s['terms'])} ")
            f.write(" ".join(str(t) for t in s["Tkeep"]) + " ")
            f.write(" ".join(f"{s['Tq'][qv]} {x} {tau} {bin(qv & s['rho']).count('1') & 1}" for qv, x, tau in s["terms"]) + " ")
            f.write(" ".join("%.17g %.17g" % (complex(Mv[b, x]).real, complex(Mv[b, x]).imag) for b in (0, 1) for x in (0, 1)) + "\n")
