#!/bin/bash
# Exact XEB of IBM's hardware samples (#151) with the in-place sampler. Powers off (= terminate) after 330 min.
set -u
BOOT=$(date +%s)
shutdown -h +330
B=graphsamp-096563427953-use1; P=ibmxeb-run1; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude gsamp2 --exclude '*.txt.in' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp2 gsamp2.cpp 2> build.log || finish "build failed"
export OMP_NUM_THREADS=$(nproc)
# Gate: the two cloud samples' exact log2 Pr, recomputed with this plan, must match the sampler run of 2026-10-29.
S=$(date +%s); ./gsamp2 nq70_inplace_plan.txt prob gate_strings.txt > gate_out.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_result.txt
exp = [-69.508530214929, -68.876165103181]
got = [float(l) for l in open("/opt/gs/gate_out.txt") if l.strip()]
ok = len(got) == 2 and all(abs(a - b) < 1e-5 for a, b in zip(exp, got))
print("PASS" if ok else "FAIL", got)
PY
log "gate: $(cat gate_result.txt) (${T}s for 2 strings)"
grep -q PASS gate_result.txt || finish "gate failed"
( while sleep 300; do up; log "progress uni=$(wc -l < uniform_out.txt 2>/dev/null) rev=$(wc -l < reversed_out.txt 2>/dev/null) ibm=$(wc -l < ibm_out.txt 2>/dev/null)"; done ) &
./gsamp2 nq70_inplace_plan.txt prob uniform.txt > uniform_out.txt; log "uniform done"
./gsamp2 nq70_inplace_plan.txt prob ibm_reversed100.txt > reversed_out.txt; log "reversed done"
/usr/bin/time -v timeout $(( (330 - 15) * 60 - ($(date +%s) - BOOT) ))s ./gsamp2 nq70_inplace_plan.txt prob ibm_qiskit.txt > ibm_out.txt 2> time_ibm.txt
log "ibm scoring ended: $(wc -l < ibm_out.txt) of 1726"
finish "complete"
