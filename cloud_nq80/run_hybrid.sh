#!/bin/bash
# Cloud run 12: exact hybrid samples of the 80-qubit, depth-80 instance of the #151 construction (seed 1080, width 35)
# on one 256 GiB node. The width-35 order is sliced at vertex 34 (width 34, one in-place table of 2^34 complex64 =
# 137 GB); vertex 34 has neighbours in the sampled suffix (last 20 steps), handled by the translation 'twist' of
# sliced_plans.py / gsamp2 resume --twist. Per proposal: 'prefix' on both slice plans to K = 60, twisted sum,
# exact draw of the last 20 bits; accept with z/1.03. Gates: (1) the binary reproduces two known 70-qubit values;
# (2) at 70 qubits, hybrid unsliced == twisted two-slice (bits and log2 Pr) on 3 proposals. Optional cross-check:
# the first 2 accepted samples re-scored with a different slicing (pair plans) if those files are present.
set -u
BOOT=$(date +%s); BUDGET_MIN=__BUDGET_MIN__
shutdown -h +$(( BUDGET_MIN + 20 ))
B=__BUCKET__; P=nq80-run12; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude gsamp2 --exclude '*.bin' --exclude '*_s[01]*.txt' --exclude '*_plan.txt' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
left() { echo $(( BUDGET_MIN * 60 - ($(date +%s) - BOOT) )); }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -o gsamp2 gsamp2.cpp 2> build.log || finish "build gsamp2 failed"
export OMP_NUM_THREADS=$(nproc)
# Gate 1: two known 70-qubit values (width-30 plan)
S=$(date +%s); ./gsamp2 nq70_inplace_plan.txt prob gate_strings.txt > gate70_out.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate70_result.txt
exp = [-69.508530214929, -68.876165103181]
got = [float(l) for l in open("/opt/gs/gate70_out.txt") if l.strip()]
print("PASS" if len(got) == 2 and all(abs(a - b) < 1e-5 for a, b in zip(exp, got)) else "FAIL", got)
PY
log "gate70: $(cat gate70_result.txt) (${T}s for 2 strings)"
grep -q PASS gate70_result.txt || finish "gate70 failed"
# Gate 2: the twist at 70 qubits: unsliced hybrid vs twisted two-slice resume on 3 proposals (seed 17), K = 50
S=$(date +%s)
./gsamp2 nq70_inplace_plan.txt hybrid tw70_props.txt 17 > gate_tw_hybrid.txt
./gsamp2 tw70_s0.txt prefix tw70_props.txt 50 tw70_d0 > /dev/null; ./gsamp2 tw70_s1.txt prefix tw70_props.txt 50 tw70_d1 > /dev/null
./gsamp2 tw70_s0.txt resume tw70_props.txt 50 17 tw70_d0 tw70_d1 --twist tw70_twist.txt > gate_tw_resume.txt; T=$(( $(date +%s) - S ))
python3 - <<'PY' > gate_tw_result.txt
hy = [l.split() for l in open("/opt/gs/gate_tw_hybrid.txt") if l.strip()]; re = [l.split() for l in open("/opt/gs/gate_tw_resume.txt") if len(l.split()[0]) == 70]
ok = len(hy) == 3 == len(re) and all(h[0] == r[0] for h, r in zip(hy, re)) and max(abs(float(h[1]) - float(r[1])) for h, r in zip(hy, re)) < 1e-5
print("PASS" if ok else "FAIL", [h[1] for h in hy], [r[1] for r in re], "bits equal:", [h[0] == r[0] for h, r in zip(hy, re)])
PY
log "gate twist70: $(cat gate_tw_result.txt) (${T}s)"
grep -q PASS gate_tw_result.txt || finish "gate twist70 failed"
up
( while sleep 300; do up; log "progress samples=$(grep -c '^ACC' samples_nq80.txt 2>/dev/null) proposals=$(grep -c '^[Ar]' samples_nq80.txt 2>/dev/null) xcheck=$(wc -l < xcheck_nq80.txt 2>/dev/null)"; done ) &
# Sampling: 80 qubits, slice vertex 34, K = 60, 20 bits drawn exactly; one proposal = 2 passes of 2^34 + the suffix
XC=0; [ -f fi80p_s00.txt ] && XC=1
S=$(date +%s); /usr/bin/time -v timeout $(( $(left) - (XC ? 150 : 15) * 60 ))s python3 doped_sampler.py sample fi80h_s0.txt fi80_order.json 100000 20261102 1.03 samples_nq80.txt --slices fi80h_s0.txt,fi80h_s1.txt --K 60 --twist fi80h_twist.txt --nd 80 --na 0 --m 20 2> time_sample.txt
log "sampling ended: $(grep -c '^ACC' samples_nq80.txt) accepted of $(grep -c '^[Ar]' samples_nq80.txt) proposals in $(( $(date +%s) - S ))s; max z $(grep -o 'z=[0-9.]*' samples_nq80.txt | cut -d= -f2 | sort -g | tail -1); peak $(grep Maximum time_sample.txt)"
up
# Cross-check with a different slicing (pair plans fi80p_s00..s11, their own K and twist), if present
if [ "$XC" = 1 ]; then
  grep '^ACC' samples_nq80.txt | head -2 | cut -d' ' -f2 > xcheck_strings.txt
  KP=$(python3 -c "import json; print(json.load(open('/opt/gs/fi80p_slices.json'))['K_suffix_start'])")
  TW=""; [ -f fi80p_twist.txt ] && TW="--twist fi80p_twist.txt"
  S=$(date +%s); timeout $(( $(left) - 10 * 60 ))s python3 doped_sampler.py score fi80p_s00.txt fi80_order.json xcheck_strings.txt xcheck_nq80.txt --slices fi80p_s00.txt,fi80p_s01.txt,fi80p_s10.txt,fi80p_s11.txt --K $KP $TW --nd 80 --na 0
  log "xcheck ended: $(wc -l < xcheck_nq80.txt) strings re-scored with the pair slicing in $(( $(date +%s) - S ))s"
fi
finish "complete"
