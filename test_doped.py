"""End-to-end check of the T-gadget pipeline: a random small circuit with mid-circuit T gates is (a) simulated densely,
(b) converted by doped_to_qasm.convert and evaluated with the cut-rank amplitude code (amp.Amplitudes) on the gadget
circuit with all ancilla outcomes 0. The two must agree as |amp_gadget|^2 * 2^k = |amp_T|^2 for every data string.
Control: dropping the ancilla T tail (measuring the ancilla after H alone) must break the agreement."""
import numpy as np, itertools, os, tempfile
from doped_gadget import run_dense, ONE
from doped_to_qasm import convert
from tableau import clifford_body
from amp import Amplitudes
from ref import run
rng = np.random.default_rng(3); worst = 0.0; ctrl_bad = 0; trials = 25
for tr in range(trials):
    n = int(rng.integers(2, 5)); ops = [("h", [q]) for q in range(n)]
    for _ in range(int(rng.integers(8, 16))):
        g = rng.choice(["s", "sx", "sxdg", "h", "cz", "rz(pi/4)"], p=[.15, .2, .1, .15, .3, .1])
        if g == "cz": a, b = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(b)]))
        else: ops.append((str(g), [int(rng.integers(n))]))
    for _ in range(int(rng.integers(1, 4))): ops.insert(int(rng.integers(1, len(ops))), ("rz(pi/4)", [int(rng.integers(n))]))
    if not any(g == "rz(pi/4)" for g, _ in ops): ops.append(("rz(pi/4)", [0]))
    dense = run_dense([("t", qs) if g == "rz(pi/4)" else (g, qs) for g, qs in ops], n)
    N, body, tail, anc = convert(n, ops); k = len(anc)
    with tempfile.NamedTemporaryFile("w", suffix=".qasm", delete=False) as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
        path = f.name
    nn, T, tl = clifford_body(path); os.unlink(path); assert nn == N and len(tl) == 2 * k
    amp = Amplitudes(N, T, tl)
    for bits in itertools.product([0, 1], repeat=n):
        b = list(bits) + [0] * k; a = run(amp.steps, amp.weights(b)) / 2 ** (N / 2)
        i = int("".join(map(str, bits)), 2)                      # dense index: qubit 0 is the most significant axis
        worst = max(worst, abs(abs(a) ** 2 * 2 ** k - abs(dense[i]) ** 2))
    # control: wrong tail (H only)
    with tempfile.NamedTemporaryFile("w", suffix=".qasm", delete=False) as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + [op for op in tail if op[0] != "t"]: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
        path = f.name
    nn, T2, tl2 = clifford_body(path); os.unlink(path); amp2 = Amplitudes(N, T2, tl2)
    bad = max(abs(abs(run(amp2.steps, amp2.weights(list(bits) + [0] * k)) / 2 ** (N / 2)) ** 2 * 2 ** k - abs(dense[int("".join(map(str, bits)), 2)]) ** 2) for bits in itertools.product([0, 1], repeat=n))
    ctrl_bad += bad > 1e-3
print(f"{trials} random T-circuits (2-4 qubits, 1-4 T gates): max |2^k |amp_gadget|^2 - |amp_T|^2| = {worst:.2e}; control (no T in the ancilla tail) differs in {ctrl_bad}/{trials}")
# the control can agree when a T gate's phase is invisible in the probabilities (e.g. a T just before a Z measurement)
assert worst < 1e-10 and ctrl_bad >= trials // 4
