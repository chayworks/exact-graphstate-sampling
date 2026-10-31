#!/bin/bash
# Cloud run 6: independent recomputation (dp.cpp built with -DDP_FLOAT: complex64 tables, double accumulation) of the
# first two 76-qubit samples drawn in run 5, with one sliced vertex (dp width 32: two tables of 2^32 complex64 = 69 GB).
# Gate: the same binary must reproduce the two known 70-qubit probabilities (width-30 plan, unsliced).
# Powers off (= terminate) after 150 min whatever happens.
set -u
BOOT=$(date +%s)
shutdown -h +150
B=__BUCKET__; P=nq76-run6; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude dp --exclude '*_plan.txt' --exclude '*.pyc' --exclude '*.qasm' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time python3-numpy >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
aws s3 cp s3://$B/nq76-run5/out/samples_nq76.txt $W/samples_nq76.txt --only-show-errors || finish "samples download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -DDP_FLOAT -o dp dp.cpp 2> build.log || finish "build dp failed"
export OMP_NUM_THREADS=$(nproc)
python3 xcheck_make.py nq70.qasm nq70_order.json nq70_slices.json gate_samples.txt 2 g70 > gate_make.txt 2>&1 || finish "gate make failed"
S=$(date +%s); ./dp g70_plan.txt g70_jobs.txt > gate_dp.txt 2> gate_err.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_result.txt
import math
exp = [-69.508530214929, -68.876165103181]
got = [math.log2(float(l.split()[0]) ** 2 + float(l.split()[1]) ** 2) - 70 for l in open("/opt/gs/gate_dp.txt") if l.strip()]
print("PASS" if len(got) == 2 and all(abs(a - b) < 1e-5 for a, b in zip(exp, got)) else "FAIL", got)
PY
log "gate70 (dp float, width 30): $(cat gate_result.txt) (${T}s for 2 strings)"
grep -q PASS gate_result.txt || finish "gate failed"
python3 xcheck_make.py fi76.qasm fi76_order.json fi76_slices.json samples_nq76.txt 2 x76 > xcheck_make.txt 2>&1 || finish "xcheck_make failed"
log "plan: $(head -1 x76_plan.txt) ; $(cut -c1-160 xcheck_make.txt)"
( while sleep 300; do up; log "progress xcheck=$(wc -l < xcheck_out.txt 2>/dev/null) rss_kB=$(ps -o rss= -C dp 2>/dev/null | head -1)"; done ) &
S=$(date +%s); /usr/bin/time -v timeout $(( (150 - 10) * 60 - ($(date +%s) - BOOT) ))s ./dp x76_plan.txt x76_jobs.txt > xcheck_out.txt 2> time_dp.txt; RC=$?
log "xcheck ended: $(wc -l < xcheck_out.txt) lines in $(( $(date +%s) - S ))s, exit $RC"
finish "complete"
