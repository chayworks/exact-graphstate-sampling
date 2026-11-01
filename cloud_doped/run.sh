#!/bin/bash
# Cloud run 7: IBM's doped 70-qubit instance (468 T gates) as a 538-vertex graph state with projected ancillas.
# In-place sampler at width 34 (one table of 2^34 complex64 = 128 GiB): needs a node with >= 200 GiB RAM
# (x2iedn.4xlarge 512 GiB / r7a.8xlarge 256 GiB). Steps: build; gate = two IBM shots scored with two different
# plans (m20 and m10 suffix layouts, different step structure) must agree to 1e-5; then score IBM's published shots
# (exact log2 Pr, full depth) for ~55% of the budget, then hybrid rejection sampling for the rest.
# Powers off (= terminate) after BUDGET_MIN + 20 minutes whatever happens.
set -u
BOOT=$(date +%s); BUDGET_MIN=__BUDGET_MIN__
shutdown -h +$(( BUDGET_MIN + 20 ))
B=__BUCKET__; P=doped-run7; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude gsamp2 --exclude '*_plan.txt' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
left() { echo $(( BUDGET_MIN * 60 - ($(date +%s) - BOOT) )); }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp2 gsamp2.cpp 2> build.log || finish "build gsamp2 failed"
export OMP_NUM_THREADS=$(nproc)
head -1 ibm_doped_shots_qiskit.txt > gate_shots.txt
# Gate: the same two hardware shots through two plans with different step structure (m20: last 20 data vertices
# deferred; m10: last 10). Also times one pass.
S=$(date +%s); /usr/bin/time -v python3 doped_sampler.py score doped_m20_inplace_plan.txt doped_m20_step_vertices.json gate_shots.txt gate_m20.txt --qiskit 2> time_gate.txt; T=$(( $(date +%s) - S ))
python3 doped_sampler.py score doped_m10_inplace_plan.txt doped_m10_step_vertices.json gate_shots.txt gate_m10.txt --qiskit
python3 - <<'PY' > gate_result.txt
a = [float(l.split()[1]) for l in open("/opt/gs/gate_m20.txt")]; b = [float(l.split()[1]) for l in open("/opt/gs/gate_m10.txt")]
print("PASS" if len(a) == 1 and len(b) == 1 and all(abs(x - y) < 1e-5 for x, y in zip(a, b)) else "FAIL", a, b)
PY
log "gate (cross-order, 1 shot): $(cat gate_result.txt); ${T}s for the m20 pass; peak $(grep Maximum time_gate.txt)"
grep -q PASS gate_result.txt || finish "gate failed"
up
( while sleep 300; do up; log "progress shots=$(wc -l < shots_scored.txt 2>/dev/null) samples=$(grep -c '^ACC' samples.txt 2>/dev/null) proposals=$(grep -c '^[Ar]' samples.txt 2>/dev/null)"; done ) &
# 1. Score IBM's shots (exact full-depth probabilities), interleaved with 16 uniform control strings first.
python3 - <<'PY'
import random; r = random.Random(20261101); shots = [l.strip() for l in open("/opt/gs/ibm_doped_shots_qiskit.txt") if l.strip()]
uni = ["".join(r.choice("01") for _ in range(70)) for _ in range(16)]
open("/opt/gs/score_list.txt", "w").write("\n".join(uni + shots) + "\n")
PY
S=$(date +%s); timeout $(( $(left) * 55 / 100 ))s python3 doped_sampler.py score doped_m20_inplace_plan.txt doped_m20_step_vertices.json score_list.txt shots_scored.txt --qiskit
log "scoring ended: $(wc -l < shots_scored.txt) strings (first 16 are uniform controls) in $(( $(date +%s) - S ))s"
up
# 2. Hybrid rejection sampling: 50 bits proposed uniformly, accepted on the exact marginal, last 20 bits drawn exactly.
S=$(date +%s); timeout $(( $(left) - 8 * 60 ))s python3 doped_sampler.py sample doped_m20_inplace_plan.txt doped_m20_step_vertices.json 100000 20261101 1.03 samples.txt
log "sampling ended: $(grep -c '^ACC' samples.txt) accepted of $(grep -c '^[Ar]' samples.txt) proposals in $(( $(date +%s) - S ))s; max z $(grep -o 'z=[0-9.]*' samples.txt | cut -d= -f2 | sort -g | tail -1)"
finish "complete"
