"""The ibm_aachen coupling map as written in IQP-NISQ's src/utils/hardware.py (MIT licence), and the gate list their
efficient_connectivity_gates builds: one X_i per qubit 0..N-1, then X_iX_j for every coupler with both ends < N,
sorted by (i, j). A vertex order of small frontier for the exact sampler: sweep the chip column by column."""
import numpy as np

def aachen_edges():
    E = []
    for r in range(8):
        b = 20 * r; E += [(b + i, b + i + 1) for i in range(15)]
    E += [(3, 16), (7, 17), (11, 18), (15, 19), (16, 23), (17, 27), (18, 31), (19, 35),
          (21, 36), (25, 37), (29, 38), (33, 39), (36, 41), (37, 45), (38, 49), (39, 53),
          (43, 56), (47, 57), (51, 58), (55, 59), (56, 63), (57, 67), (58, 71), (59, 75),
          (61, 76), (65, 77), (69, 78), (73, 79), (76, 81), (77, 85), (78, 89), (79, 93),
          (83, 96), (87, 97), (91, 98), (95, 99), (96, 103), (97, 107), (98, 111), (99, 115),
          (101, 116), (105, 117), (109, 118), (113, 119), (116, 121), (117, 125), (118, 129), (119, 133),
          (123, 136), (127, 137), (131, 138), (135, 139), (136, 143), (137, 147), (138, 151), (139, 155)]
    return E

def gate_list(N, layers=1):
    one = [(i,) for i in range(N)]
    two = sorted({tuple(sorted(e)) for e in aachen_edges() if e[0] < N and e[1] < N})
    return (one + two) * layers

def coords():
    """(row, column) of every qubit: rows of 16 at r = 2*row, bridge qubits between rows at r = 2*row + 1."""
    c = {}
    for r in range(8):
        for i in range(16): c[20 * r + i] = (2 * r, i)
    for (a, b) in aachen_edges():                       # bridge qubits: index 16..19 mod 20, between rows
        for q, other in ((a, b), (b, a)):
            if q % 20 >= 16 and q not in c: c[q] = (2 * (q // 20) + 1, other % 20)
    assert len(c) == 156, len(c)
    return c

def column_order(N):
    c = coords(); qs = [q for q in range(N)]
    return sorted(qs, key=lambda q: (c[q][1], c[q][0]))
