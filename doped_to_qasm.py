"""Write the T-gadget form of a doped circuit as a plain Clifford-plus-tail QASM so that the whole toolchain
(clifford_body -> graph_form -> plan/inplace -> gsamp2/dp.cpp) applies unchanged.

Each mid-circuit rz(pi/4) on q becomes  h a; cz q,a; h a;  (CNOT q->a onto a fresh ancilla a) in the Clifford body,
and the ancilla's tail is  t a; h a;  (measure after T then H; outcome 0 is the projection <phi| that implements T,
see doped_gadget.py). Data qubits keep whatever tail the original circuit has after its last Clifford gate (none for
the doped instance: Z measurement). Also writes the spatial-sweep vertex order and the ancilla table.
Usage: python doped_to_qasm.py <doped.qasm> <out_prefix>   -> <out_prefix>.qasm, <out_prefix>_sweep_order.json, <out_prefix>_anc.json"""
import sys, re, json
CLIFF = {"h", "s", "sdg", "sx", "sxdg", "x", "z", "cz"}
def parse(path):
    n = None; ops = []
    for line in open(path):
        line = line.strip().rstrip(";")
        if not line or line.startswith(("OPENQASM", "include")): continue
        m = re.match(r"qreg q\[(\d+)\]", line)
        if m: n = int(m.group(1)); continue
        name, args = line.split(" ", 1); qs = [int(a) for a in re.findall(r"q\[(\d+)\]", args)]
        ops.append((name, qs))
    return n, ops

def convert(n, ops):
    """Returns (N, body ops, tail ops, ancillas [(a, q, gate_index)]). Everything after the last mid-circuit T or
    Clifford gate that is not Clifford is treated as the data tail."""
    last_cliff = max(i for i, (g, qs) in enumerate(ops) if g in CLIFF or g == "rz(pi/4)")
    body_src, tail_src = ops[:last_cliff + 1], ops[last_cliff + 1:]
    body = []; anc = []; k = 0
    for i, (g, qs) in enumerate(body_src):
        if g == "rz(pi/4)":
            a = n + k; k += 1; anc.append((a, qs[0], i)); body += [("h", [a]), ("cz", [qs[0], a]), ("h", [a])]
        elif g in CLIFF: body.append((g, qs))
        else: raise ValueError(f"unexpected mid-circuit gate {g}")
    tail = list(tail_src) + [op for a, _, _ in anc for op in (("t", [a]), ("h", [a]))]
    return n + k, body, tail, anc

if __name__ == "__main__":
    n, ops = parse(sys.argv[1]); out = sys.argv[2]
    N, body, tail, anc = convert(n, ops)
    with open(out + ".qasm", "w") as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
    by_q = {q: [] for q in range(n)}
    for a, q, _ in anc: by_q[q].append(a)
    sweep = []
    for q in range(n): sweep += by_q[q] + [q]
    json.dump(sweep, open(out + "_sweep_order.json", "w"))
    json.dump([{"vertex": a, "qubit": q, "gate_index": i} for a, q, i in anc], open(out + "_anc.json", "w"))
    print(f"{N} qubits ({n} data + {len(anc)} ancillas), {len(body)} Clifford gates, {len(tail)} tail gates -> {out}.qasm")
