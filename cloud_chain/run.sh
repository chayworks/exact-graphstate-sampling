#!/bin/bash
# Cloud run 3: exact per-step conditionals (gsamp chain mode) of IBM's 1,726 hardware shots and of our 1,012 exact
# samples, in the width-30 order and in its reverse (same width; covers the qubits that sit early in the forward
# order). Powers off (= terminate) after 420 min whatever happens.
set -u
BOOT=$(date +%s)
shutdown -h +420
B=__BUCKET__; P=chain-run3; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude gsamp --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp gsamp.cpp 2> build.log || finish "build failed"
export OMP_NUM_THREADS=$(nproc)
# Gate: log2 Pr of the two reference samples, forward and reversed plan, must match the earlier runs.
S=$(date +%s); ./gsamp plan_fwd.txt chain gate_fwd.txt | cut -d' ' -f1 > gate_fwd_out.txt; ./gsamp plan_rev.txt chain gate_rev.txt | cut -d' ' -f1 > gate_rev_out.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_result.txt
exp = [-69.508530214929, -68.876165103181]
ok = True
for f in ("/opt/gs/gate_fwd_out.txt", "/opt/gs/gate_rev_out.txt"):
    got = [float(l) for l in open(f) if l.strip()]
    ok = ok and len(got) == 2 and all(abs(a - b) < 1e-5 for a, b in zip(exp, got)); print(f, got)
print("PASS" if ok else "FAIL")
PY
log "gate: $(tail -1 gate_result.txt) (${T}s for 4 strings)"
grep -q PASS gate_result.txt || finish "gate failed"
( while sleep 300; do up; log "progress $(for f in ibm_fwd ibm_rev ideal_fwd ideal_rev; do echo -n "$f=$(wc -l < ${f}_chain.txt 2>/dev/null) "; done)"; done ) &
./gsamp plan_fwd.txt chain ibm_fwd.txt > ibm_fwd_chain.txt; log "ibm_fwd done"
./gsamp plan_rev.txt chain ibm_rev.txt > ibm_rev_chain.txt; log "ibm_rev done"
./gsamp plan_fwd.txt chain ideal_fwd.txt > ideal_fwd_chain.txt; log "ideal_fwd done"
/usr/bin/time -v timeout $(( (420 - 12) * 60 - ($(date +%s) - BOOT) ))s ./gsamp plan_rev.txt chain ideal_rev.txt > ideal_rev_chain.txt 2> time_last.txt
log "ideal_rev ended: $(wc -l < ideal_rev_chain.txt) of 1012"
finish "complete"
