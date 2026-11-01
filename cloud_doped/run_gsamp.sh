#!/bin/bash
# Cloud run 11: the doped instance on the two-table sampler gsamp.cpp (streaming writes; measured 2.7x the in-place
# engine's rate per update on the fstate instance), width-34 order sliced at ancilla vertices 135 and 104 -> four
# plans of width 33, two tables of 2^33 complex64 = 137 GB: needs a 256 GiB node (spot r7a/r7i/r7iz.8xlarge).
# Per proposal: 'prefix' on each slice plan to step K = 518, then 'resume' on the summed dumps. Gate: IBM doped shot #1
# must reproduce the run-7/9/10 value -68.29533 (+/- 1e-5). Then score shots, then sample. Auto-terminates.
set -u
BOOT=$(date +%s); BUDGET_MIN=__BUDGET_MIN__
shutdown -h +$(( BUDGET_MIN + 20 ))
B=__BUCKET__; P=doped-run11; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude gsamp --exclude '*.bin' --exclude '*_s[01][01].txt' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
left() { echo $(( BUDGET_MIN * 60 - ($(date +%s) - BOOT) )); }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp gsamp.cpp 2> build.log || finish "build gsamp failed"
export OMP_NUM_THREADS=$(nproc)
SL=doped_m20g_s00.txt,doped_m20g_s01.txt,doped_m20g_s10.txt,doped_m20g_s11.txt; K=518; BIN=./gsamp
head -1 ibm_doped_shots_qiskit.txt > gate_shots.txt
S=$(date +%s); /usr/bin/time -v python3 doped_sampler.py score doped_m20g_s00.txt doped_m20_step_vertices.json gate_shots.txt gate.txt --qiskit --slices $SL --K $K --bin $BIN 2> time_gate.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_result.txt
a = [float(l.split()[1]) for l in open("/opt/gs/gate.txt")]
print("PASS" if len(a) == 1 and abs(a[0] - (-68.295331)) < 1e-5 else "FAIL", a, "reference -68.295331 (runs 7, 9, 10)")
PY
log "gate (shot #1 vs runs 7/9/10): $(cat gate_result.txt); ${T}s for one sliced pass (4 slice passes); peak $(grep Maximum time_gate.txt)"
grep -q PASS gate_result.txt || finish "gate failed"
up
( while sleep 300; do up; log "progress shots=$(wc -l < shots_scored.txt 2>/dev/null) samples=$(grep -c '^ACC' samples.txt 2>/dev/null) proposals=$(grep -c '^[Ar]' samples.txt 2>/dev/null)"; done ) &
python3 - <<'PY'
import random; r = random.Random(20261101); shots = [l.strip() for l in open("/opt/gs/ibm_doped_shots_qiskit.txt") if l.strip()]
uni = ["".join(r.choice("01") for _ in range(70)) for _ in range(8)]
open("/opt/gs/score_list.txt", "w").write("\n".join(uni + shots[1:]) + "\n")
PY
S=$(date +%s); timeout $(( $(left) * 55 / 100 ))s python3 doped_sampler.py score doped_m20g_s00.txt doped_m20_step_vertices.json score_list.txt shots_scored.txt --qiskit --slices $SL --K $K --bin $BIN
log "scoring ended: $(wc -l < shots_scored.txt) strings (first 8 are uniform controls) in $(( $(date +%s) - S ))s"
up
S=$(date +%s); timeout $(( $(left) - 8 * 60 ))s python3 doped_sampler.py sample doped_m20g_s00.txt doped_m20_step_vertices.json 100000 20261101 1.03 samples.txt --slices $SL --K $K --bin $BIN
log "sampling ended: $(grep -c '^ACC' samples.txt) accepted of $(grep -c '^[Ar]' samples.txt) proposals in $(( $(date +%s) - S ))s; max z $(grep -o 'z=[0-9.]*' samples.txt | cut -d= -f2 | sort -g | tail -1)"
finish "complete"
