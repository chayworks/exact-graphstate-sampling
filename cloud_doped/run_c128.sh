#!/bin/bash
# Cloud run 15: settle the common factor 1 - 1.8e-4 between our complex64 doped amplitudes and #242's.
# gsamp2 built with -DGS_DOUBLE (complex128 tables and dumps, validated locally: test_gsamp2_doped 1.2e-13 against
# dense, sliced driver test_slice_adjacent_c128). Run 10's two-slice plans (ancillas 135, 104; four plans of width 33),
# so one 2^33 complex128 table = 128 GiB on an r7a.8xlarge (256 GiB). Scores IBM doped shot 1 (complex64 value
# -68.295330880 from run 7, seven ways to 1.4e-6; #242: -68.295074) and then the last shot (run 11, complex64
# -74.73 region) if time allows. No gate is needed: the point is the value itself, compared afterwards.
# Powers off (= terminate) after BUDGET_MIN + 20 minutes whatever happens.
set -u
BOOT=$(date +%s); BUDGET_MIN=__BUDGET_MIN__
shutdown -h +$(( BUDGET_MIN + 20 ))
B=__BUCKET__; P=doped-run15-c128; W=/opt/gs; mkdir -p $W; cd $W
log() { echo "$(date -u +%FT%TZ) $*" >> $W/run.log; aws s3 cp $W/run.log s3://$B/$P/run.log --only-show-errors; }
up() { aws s3 cp $W/ s3://$B/$P/out/ --recursive --exclude '*.cpp' --exclude '*.py' --exclude gsamp2 --exclude '*_plan*.txt' --exclude '*.bin' --exclude '*_s[01][01].txt' --only-show-errors; }
finish() { log "finish: $1"; up; shutdown -h now; }
left() { echo $(( BUDGET_MIN * 60 - ($(date +%s) - BOOT) )); }
log "boot nproc=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
dnf install -y gcc-c++ libgomp time >> $W/dnf.log 2>&1 || finish "dnf failed"
aws s3 cp s3://$B/$P/in/ $W/ --recursive --only-show-errors || finish "download failed"
sha256sum $W/* > $W/inputs.sha256 2>/dev/null
g++ -O3 -march=native -fopenmp -DGS_DOUBLE -o gsamp2 gsamp2.cpp 2> build.log || finish "build gsamp2 failed"
export OMP_NUM_THREADS=$(nproc)
SL=doped_m20sl_s00.txt,doped_m20sl_s01.txt,doped_m20sl_s10.txt,doped_m20sl_s11.txt; K=518
head -1 ibm_doped_shots_qiskit.txt > shot1.txt; grep . ibm_doped_shots_qiskit.txt | tail -1 > shotlast.txt
( while sleep 300; do up; log "progress"; done ) &
S=$(date +%s); /usr/bin/time -v python3 doped_sampler.py score doped_m20sl_s00.txt doped_m20_step_vertices.json shot1.txt c128_shot1.txt --qiskit --slices $SL --K $K 2> time_shot1.txt; T=$(( $(date +%s) - S ))
log "shot 1 complex128: $(cat c128_shot1.txt) (complex64 -68.295330880; #242 -68.295074); ${T}s; peak $(grep Maximum time_shot1.txt)"
up
if [ $(left) -gt $(( T + 15 * 60 )) ]; then
  S=$(date +%s); /usr/bin/time -v python3 doped_sampler.py score doped_m20sl_s00.txt doped_m20_step_vertices.json shotlast.txt c128_shotlast.txt --qiskit --slices $SL --K $K 2> time_shotlast.txt
  log "last shot complex128: $(cat c128_shotlast.txt); $(( $(date +%s) - S ))s"
else
  log "last shot skipped: $(left)s left"
fi
finish "complete"
