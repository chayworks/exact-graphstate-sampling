"""Statistical certificate for the hybrid sampler from its own proposals, without a Porter-Thomas assumption.

The sampler draws a prefix y uniformly, computes z_y = D * p(y) (the exact prefix marginal in units of the uniform
distribution, D = 2^(number of prefix bits)), accepts with probability min(1, z_y / c) and draws the suffix exactly.
The accepted prefix law is q(y) = min(z_y, c) / (D m) with m = E_u[min(z, c)], so with eps = E_u[(z - c)^+] = 1 - m:
    TV(q, p) <= 1/2 (E_u|w/m - w| + E_u|w - z|) = 1/2 ((1 - m) + eps) = eps,      w = min(z, c),
and the suffix is exact given the prefix, so the same bound holds for the full strings. For the classical fidelity
F = (sum sqrt(q p))^2 (the Bhattacharyya coefficient squared; for a prefix law with an exact suffix this is the
quantity #263 bounds), sqrt(w z) >= w gives sum sqrt(q p) = E_u sqrt(w z) / sqrt(m) >= sqrt(m), i.e. F >= m. Every proposal (accepted or not) gives an i.i.d. draw
of w in [0, c] under the uniform law, so a one-sided lower confidence bound on m bounds eps. We use the empirical
Bernstein bound of Maurer & Pontil (2009, Thm 4) for X = w / c in [0, 1]:
    P( E X >= mean - sqrt(2 V ln(2/delta) / n) - 7 ln(2/delta) / (3 (n - 1)) ) >= 1 - delta,  V the sample variance.
Floating-point error in z (complex64 tables; the gate agrees with #249 to 4e-6 in log2) is not included.

Usage: python hybrid_certificate.py <samples.txt> <cap> [delta]     (reads z= from every ACC/rej line)
       python hybrid_certificate.py --selftest"""
import sys, math, random

def bound(ws, c, delta):
    n = len(ws); x = [w / c for w in ws]; mu = sum(x) / n; v = sum((t - mu) ** 2 for t in x) / (n - 1)
    L = math.log(2 / delta)
    lo = mu - math.sqrt(2 * v * L / n) - 7 * L / (3 * (n - 1))
    hoeff = mu - math.sqrt(math.log(1 / delta) / (2 * n))
    return c * lo, c * hoeff, mu * c

def report(zs, c, delta):
    ws = [min(z, c) for z in zs]; mlo, mh, mhat = bound(ws, c, delta)
    tv = max(0.0, 1 - mlo); F = min(1.0, max(0.0, mlo))
    print(f"n = {len(zs)} proposals, cap {c}, max z {max(zs):.6f}, above cap {sum(z > c for z in zs)}")
    print(f"mean min(z,c) = {mhat:.6f}; {100 * (1 - delta):.0f}% lower bound (empirical Bernstein) {mlo:.4f} "
          f"(Hoeffding {mh:.4f})")
    print(f"=> TV(sampled, ideal) <= {tv:.4f} and classical fidelity >= {F:.4f} at {100 * (1 - delta):.0f}% confidence")
    return tv, F

def selftest():
    # 1. the inequality TV <= eps on exact small distributions, including heavy excluded mass
    rnd = random.Random(1)
    for c in (1.03, 1.5, 3.0):
        D = 4096; z = [rnd.expovariate(1.0) for _ in range(D)]; s = sum(z); z = [t * D / s for t in z]   # E_u z = 1
        w = [min(t, c) for t in z]; m = sum(w) / D; eps = 1 - m
        tv = 0.5 * sum(abs(wi / m - zi) for wi, zi in zip(w, z)) / D
        F = (sum(math.sqrt(wi / m * zi) for wi, zi in zip(w, z)) / D) ** 2
        assert tv <= eps + 1e-12, (c, tv, eps); assert F >= m - 1e-12, (c, F, m)
        print(f"selftest c={c}: TV {tv:.4f} <= eps {eps:.4f}; F {F:.4f} >= m {m:.4f}")
    # 2. the TV bound is not vacuous in the other direction: with no mass above c the true TV is 0
    z = [1.0] * 100; tv, _ = report(z, 1.03, 0.05); assert tv < 0.2
    # 3. coverage: the 95% bound on m must fail in at most ~5% of repetitions (boundary of the confidence claim)
    D = 1 << 14; z = [rnd.expovariate(1.0) for _ in range(D)]; s = sum(z); z = [t * D / s for t in z]; c = 1.5
    m = sum(min(t, c) for t in z) / D; fails = 0; R = 400
    for _ in range(R):
        ws = [min(z[rnd.randrange(D)], c) for _ in range(102)]
        if bound(ws, c, 0.05)[0] > m: fails += 1
    print(f"selftest coverage: bound above the true mean in {fails}/{R} repetitions (allowed about 5%)")
    assert fails <= 0.05 * R * 1.5
    print("selftest PASS")

if __name__ == "__main__":
    if sys.argv[1] == "--selftest": selftest(); raise SystemExit
    path, c = sys.argv[1], float(sys.argv[2]); delta = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    zs = [float(l.split("z=")[1].split()[0]) for l in open(path) if l.startswith(("ACC", "rej"))]
    report(zs, c, delta)
