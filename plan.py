"""Plan the cut-rank dynamic program for  Z(w) = sum_x prod_j w_j(x_j) (-1)^{sum_{i<j} A_ij x_i x_j}.
This is the linear rank-width method (Van den Nest, Duer, Vidal, Briegel, PRA 75, 012337 (2007)): the largest cut rank
r_k along `order` is the width of the linear rank decomposition the order defines, an upper bound on the linear
rank-width of A, and the table at step k has 2^{r_k} entries.

Processing vertices in `order`, the processed vertices act on the future ones only through the field
s = A[F, P] x_P, which lives in an r_k-dimensional GF(2) space (r_k = cut rank). States are coordinates
c in a basis of that space. Each step is the surjective linear map Lam:(c, x_v) -> c' and we evaluate it
in gather form: new[c'] = sum_{u in R c' + ker Lam} old[u_c] w_v(u_x) (-1)^{u_x <rho, u_c>}.
"""

def plan(A, order):
    n = len(A); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]
    B = []; done = 0; steps = []
    for v in order:
        r = len(B); vb = 1 << v
        rho = sum(1 << j for j, b in enumerate(B) if b & vb)
        done |= vb
        gens = [b & ~vb for b in B] + [nb[v] & ~done]          # index r is x_v
        basis, piv, red = [], [], []                            # red: (reduced vector, combination mask over gens)
        lam = []
        for i, g in enumerate(gens):
            comb = 1 << i; t = g
            for (bv, bc) in red:
                if t & (1 << (bv.bit_length() - 1)): t ^= bv; comb ^= bc
            if t:
                red.append((t, comb)); red.sort(key=lambda p: -p[0]); piv.append(i)
            else:
                lam.append(comb)                                   # comb is a kernel vector of Lam
        # new basis = gens at pivot positions (coordinate j <-> gens[piv[j]])
        steps.append(dict(v=v, r=r, rn=len(piv), rho=rho, piv=piv, ker=lam))
        B = [gens[i] for i in piv]
        B_red_ok = True
    assert not B, "state space must be empty at the end"
    return steps
