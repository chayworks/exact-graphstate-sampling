"""Exact minimum stabilizer weight of a graph state by integer programming (scipy HiGHS).
Variables: v_j (chosen vertices), y_j = (Av)_j mod 2 via y_j = sum_i A_ji v_i - 2 t_j, s_j >= v_j, s_j >= y_j.
Minimise sum s_j subject to sum v_j >= 1. A proven optimum (status 0) is d_min; a time-limited run gives bounds.
Usage: python -P minweight_ilp.py <A.txt> [time_limit_s]"""
import sys, time, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
A = np.array([[int(c) for c in l.strip()] for l in open(sys.argv[1]) if l.strip()]); n = len(A)
tl = float(sys.argv[2]) if len(sys.argv) > 2 else 600
# variable layout: v[0:n], y[n:2n], s[2n:3n], t[3n:4n]
N = 4 * n; c = np.zeros(N); c[2 * n:3 * n] = 1
rows, lo, hi = [], [], []
for j in range(n):                       # sum_i A_ji v_i - 2 t_j - y_j = 0
    r = np.zeros(N); r[:n] = A[j]; r[3 * n + j] = -2; r[n + j] = -1; rows.append(r); lo.append(0); hi.append(0)
    r = np.zeros(N); r[2 * n + j] = 1; r[j] = -1; rows.append(r); lo.append(0); hi.append(np.inf)          # s_j >= v_j
    r = np.zeros(N); r[2 * n + j] = 1; r[n + j] = -1; rows.append(r); lo.append(0); hi.append(np.inf)      # s_j >= y_j
r = np.zeros(N); r[:n] = 1; rows.append(r); lo.append(1); hi.append(np.inf)                                # v != 0
ub = np.ones(N); ub[3 * n:] = n
t0 = time.time()
res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=np.ones(N), bounds=Bounds(np.zeros(N), ub),
           options={"time_limit": tl, "disp": False})
print(f"status {res.status} ({res.message}); {time.time() - t0:.1f} s")
if res.x is not None:
    v = np.nonzero(np.round(res.x[:n]))[0]; w = int(round(res.fun))
    chk = int(np.count_nonzero(np.r_[np.zeros(n, int)].astype(int) | 0) ) if False else None
    vv = np.zeros(n, int); vv[v] = 1; wd = int(np.count_nonzero(vv | (A @ vv % 2)))
    print(f"weight {w} (recomputed from v: {wd}) at v = {v.tolist()}; MIP gap {getattr(res, 'mip_gap', float('nan'))}, dual bound {getattr(res, 'mip_dual_bound', float('nan'))}")
    print("PROVEN d_min = %d" % w if res.status == 0 else "NOT PROVEN (time limit): d_min in [%s, %d]" % (getattr(res, 'mip_dual_bound', '?'), w))
