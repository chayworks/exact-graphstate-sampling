#!/bin/bash
# Cloud run 5: exact Born samples of an 76-qubit, CZ-depth-76 instance of the #151 construction (seed 1000+76), with the
# in-place sampler at width 33 (69 GB table) on an r7a.4xlarge (16 vCPU, 128 GiB), then an independent
# recomputation of the first two samples with the sliced amplitude engine dp.cpp.
# Powers off (= terminate) after 180 min whatever happens.
set -u
BOOT=$(date +%s)
shutdown -h +180
B=__BUCKET__; P=nq76-run5; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude gsamp2 --exclude dp --exclude '*_plan.txt' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time python3-numpy >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp2 gsamp2.cpp 2> build.log || finish "build gsamp2 failed"
g++ -O3 -march=native -fopenmp -o dp dp.cpp 2>> build.log || finish "build dp failed"
export OMP_NUM_THREADS=$(nproc)
# Gate 1: the 70-qubit width-30 plan must reproduce the two reference samples (same binary, same machine).
S=$(date +%s); ./gsamp2 nq70_inplace_plan.txt prob gate_strings.txt > gate_out.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_result.txt
exp = [-69.508530214929, -68.876165103181]
got = [float(l) for l in open("/opt/gs/gate_out.txt") if l.strip()]
print("PASS" if len(got) == 2 and all(abs(a - b) < 1e-5 for a, b in zip(exp, got)) else "FAIL", got)
PY
log "gate70: $(cat gate_result.txt) (${T}s for 2 strings)"
grep -q PASS gate_result.txt || finish "gate failed"
( while sleep 300; do up; log "progress samples=$(wc -l < samples_nq76.txt 2>/dev/null) xcheck=$(wc -l < xcheck_out.txt 2>/dev/null)"; done ) &
# Samples: as many as fit in 60 min; the rest of the window goes to the independent check of sample 1.
S=$(date +%s)
/usr/bin/time -v timeout 60m ./gsamp2 fi76_inplace_plan.txt sample 1000 8000000 > samples_nq76.txt 2> time_sample.txt
log "sampling ended: $(wc -l < samples_nq76.txt) samples in $(( $(date +%s) - S ))s"
up
# Independent check of the first two samples with dp.cpp (different algorithm, sliced plan; 2 tables of 2^32 complex64).
python3 xcheck_make.py fi76.qasm fi76_order.json fi76_slices.json samples_nq76.txt 2 x76 > xcheck_make.txt 2>&1 || finish "xcheck_make failed"
S=$(date +%s); timeout $(( (180 - 10) * 60 - ($(date +%s) - BOOT) ))s ./dp x76_plan.txt x76_jobs.txt > xcheck_out.txt 2> xcheck_err.txt
log "xcheck ended: $(wc -l < xcheck_out.txt) lines in $(( $(date +%s) - S ))s"
finish "complete"
