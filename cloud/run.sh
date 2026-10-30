#!/bin/bash
# EC2 user-data for the nq70 exact-sampling run. Placeholders: __BUCKET__, __PREFIX__, __MAXMIN__, __NSAMP__.
# Safety: the instance powers off after __MAXMIN__ minutes no matter what; launched with
# InstanceInitiatedShutdownBehavior=terminate, so power-off = terminate (billing stops).
set -u
shutdown -h +__MAXMIN__
B=__BUCKET__; P=__PREFIX__; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
finish() { log "finish: $1"; aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude gsamp --only-show-errors; shutdown -h now; }
log "boot $(curl -s http://169.254.169.254/latest/meta-data/instance-type) nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G"
log "cpu $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp gsamp.cpp 2> build.log || finish "build failed"
log "built: $(g++ --version | head -1)"

# 1. Reproduce the three 55-qubit samples generated locally (same seeds). Bits must match exactly and
#    log2 Pr within 1e-6; otherwise stop before spending on the big run.
OMP_NUM_THREADS=$(nproc) ./gsamp check55_plan.txt sample 3 1 > check55_cloud.txt
python3 - <<'PY' > check55_result.txt
a = [l.split() for l in open("/opt/gs/check55_expected.txt") if l.strip()]
b = [l.split() for l in open("/opt/gs/check55_cloud.txt") if l.strip()]
ok = len(a) == len(b) == 3 and all(x[0] == y[0] and abs(float(x[1]) - float(y[1])) < 1e-6 and x[2] == y[2] for x, y in zip(a, b))
print("PASS" if ok else "FAIL"); [print(" ".join(x), "|", " ".join(y)) for x, y in zip(a, b)]
PY
log "check55: $(head -1 check55_result.txt)"
grep -q PASS check55_result.txt || finish "check55 failed"

# 2. Full 70-qubit exact samples, streamed; uploaded every 5 minutes.
( while sleep 300; do aws s3 cp $W/samples_nq70.txt s3://$B/$P/out/samples_nq70.txt --only-show-errors; log "progress $(wc -l < $W/samples_nq70.txt) samples"; done ) &
START=$(date +%s)
/usr/bin/time -v timeout $(( (__MAXMIN__ - 20) * 60 ))s env OMP_NUM_THREADS=$(nproc) ./gsamp nq70_plan.txt sample __NSAMP__ 1000000 > samples_nq70.txt 2> time_nq70.txt
log "sampling done in $(( $(date +%s) - START )) s, $(wc -l < samples_nq70.txt) samples"
finish "complete"
