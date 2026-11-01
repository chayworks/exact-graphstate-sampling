#!/bin/bash
# Kaggle (TPU VM, 330 GB RAM, 96 vCPU) runner for the doped instance. Expects the inputs next to this script
# (a Kaggle dataset mounted at /kaggle/input/<name>/); writes everything to /kaggle/working (saved with the notebook).
# Usage: bash run_kaggle.sh <input_dir> <budget_minutes>
set -u
IN=$1; BUDGET_MIN=$2; W=/kaggle/working/doped; mkdir -p $W; cd $W; BOOT=$(date +%s)
cp $IN/* $W/ 2>/dev/null; sha256sum $W/*.txt $W/*.json $W/*.cpp > inputs.sha256
left() { echo $(( BUDGET_MIN * 60 - ($(date +%s) - BOOT) )); }
log() { echo "$(date -u +%FT%TZ) $*" | tee -a run.log; }
log "nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G"
g++ -O3 -march=native -fopenmp -o gsamp2 gsamp2.cpp || { log "build failed"; exit 1; }
export OMP_NUM_THREADS=$(nproc)
head -2 ibm_doped_shots_qiskit.txt > gate_shots.txt
S=$(date +%s); python3 doped_sampler.py score doped_m20_inplace_plan.txt doped_m20_step_vertices.json gate_shots.txt gate_m20.txt --qiskit; T=$(( $(date +%s) - S ))
python3 doped_sampler.py score doped_m10_inplace_plan.txt doped_m10_step_vertices.json gate_shots.txt gate_m10.txt --qiskit
python3 - <<'PY' | tee gate_result.txt
a = [float(l.split()[1]) for l in open("gate_m20.txt")]; b = [float(l.split()[1]) for l in open("gate_m10.txt")]
print("PASS" if len(a) == 2 and len(b) == 2 and all(abs(x - y) < 1e-5 for x, y in zip(a, b)) else "FAIL", a, b)
PY
log "gate: $(cat gate_result.txt); ${T}s for two m20 passes"
grep -q PASS gate_result.txt || exit 1
python3 - <<'PY'
import random; r = random.Random(20261101); shots = [l.strip() for l in open("ibm_doped_shots_qiskit.txt") if l.strip()]
uni = ["".join(r.choice("01") for _ in range(70)) for _ in range(16)]
open("score_list.txt", "w").write("\n".join(uni + shots) + "\n")
PY
S=$(date +%s); timeout $(( $(left) * 55 / 100 ))s python3 doped_sampler.py score doped_m20_inplace_plan.txt doped_m20_step_vertices.json score_list.txt shots_scored.txt --qiskit
log "scoring ended: $(wc -l < shots_scored.txt) strings in $(( $(date +%s) - S ))s"
S=$(date +%s); timeout $(( $(left) - 5 * 60 ))s python3 doped_sampler.py sample doped_m20_inplace_plan.txt doped_m20_step_vertices.json 100000 20261101 1.03 samples.txt
log "sampling ended: $(grep -c '^ACC' samples.txt) accepted of $(grep -c '^[Ar]' samples.txt) in $(( $(date +%s) - S ))s"
log "done"
