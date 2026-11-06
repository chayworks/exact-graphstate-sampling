#!/bin/bash
# Cloud run 14: IBM's current doped instance (64 qubits on a ring, depth 73, 314 T gates; tracker #245, the circuit of
# #249, MD5 054441020510c19c628284de146f5f34). T gates as projected ancillas -> 378-vertex graph state; Renaud's
# published order (width 36) with the last 20 data vertices moved to the end (width unchanged), sliced at ancilla
# vertices 317 and 211 -> four in-place plans of width 34 (2^34 complex64 = 137 GB; no twist needed). One 256 GiB node.
# Gate: IBM's shots 0, 1, 2 of Experiment 1 re-scored and compared with Renaud's published probabilities (Zenodo
# 22882048) to |d log2 Pr| < 1e-3. Then hybrid samples (44 data bits proposed, 20 drawn exactly, cap 1.03)
# until the budget ends. Auto-terminates.
set -u
BOOT=$(date +%s); BUDGET_MIN=__BUDGET_MIN__
shutdown -h +$(( BUDGET_MIN + 20 ))
B=__BUCKET__; P=ring64-run14; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude gsamp2 --exclude '*.bin' --exclude 'r64h_s*.txt' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
left() { echo $(( BUDGET_MIN * 60 - ($(date +%s) - BOOT) )); }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp2 gsamp2.cpp 2> build.log || finish "build gsamp2 failed"
export OMP_NUM_THREADS=$(nproc)
SL=r64h_s00.txt,r64h_s01.txt,r64h_s10.txt,r64h_s11.txt; K=358; FLAGS="--nd 64 --na 314 --slices $SL --K $K --bin ./gsamp2"
S=$(date +%s); /usr/bin/time -v python3 doped_sampler.py score r64h_s00.txt order_m20.json gate_shots.txt gate.txt --qiskit $FLAGS 2> time_gate.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_result.txt
import json
got = [float(l.split()[1]) for l in open("/opt/gs/gate.txt") if l.strip()]; exp = json.load(open("/opt/gs/gate_expected.json"))
d = [g - e for g, e in zip(got, exp)]
print("PASS" if len(got) == 3 and all(abs(x) < 1e-3 for x in d) else "FAIL", "ours", got, "Renaud", exp, "diff", d)
PY
log "gate (3 IBM shots vs Renaud's #249 probabilities): $(cat gate_result.txt); ${T}s; peak $(grep Maximum time_gate.txt)"
grep -q PASS gate_result.txt || finish "gate failed"
up
( while sleep 300; do up; log "progress samples=$(grep -c '^ACC' samples.txt 2>/dev/null) proposals=$(grep -c '^[Ar]' samples.txt 2>/dev/null) above_cap=$(grep -c Z_ABOVE_CAP samples.txt 2>/dev/null)"; done ) &
S=$(date +%s); /usr/bin/time -v timeout $(( $(left) - 10 * 60 ))s python3 doped_sampler.py sample r64h_s00.txt order_m20.json 100000 20261104 1.03 samples.txt $FLAGS 2> time_sample.txt
log "sampling ended: $(grep -c '^ACC' samples.txt) accepted of $(grep -c '^[Ar]' samples.txt) proposals in $(( $(date +%s) - S ))s; max z $(grep -o 'z=[0-9.]*' samples.txt | cut -d= -f2 | sort -g | tail -1); peak $(grep Maximum time_sample.txt)"
finish "complete"
