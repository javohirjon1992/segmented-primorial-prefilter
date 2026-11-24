#!/usr/bin/env bash
# Runs the C++ experiments of the paper and collects the raw per-run CSV files in results/.
#   ./run_all.sh         full run with the settings of the paper (about 40-70 min on a laptop)
#   ./run_all.sh quick   short smoke test (2-3 min)
#   ./run_all.sh new     steps 11-15 only
#   ./run_all.sh ossl    step 13 only: OpenSSL algorithm with and without the word prefilter (5-10 min)
#   ./run_all.sh confirm step 14 only: confirmatory experiment, docs/confirmatory_protocol_EN.txt (15-25 min)
#   ./run_all.sh replicate  replicated seed-dependent counts, 22,860 provable primes, complete generators (about 1-2 h)
#   ./run_all.sh audit      ctgrind analysis of the compiled constant-time variants (ctgrind/run_ctgrind.sh; Valgrind)
#   ./run_all.sh plan       confirmatory V0-V3 experiment at 2048/3072 bits (docs/protocol_V0-V3_2048_3072_EN.txt; hours)
set -e
cd "$(dirname "$0")"
Q=${1:-full}
R=results; mkdir -p $R
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a $R/log.txt; }
if [ "$Q" = "quick" ]; then SBN=512; SBB=1000000; RK2=1; RK3=1; OK1=3; OK2=2; OKR=200; CK=4; else SBN=4096; SBB=30000000; RK2=12; RK3=4; OK1=200; OK2=100; OKR=5000; CK=1000; fi
if [ "$Q" = "quick" ]; then N=200; S=150; K1=2; K2=1; G1=3; G2=1; GA=40; RS=200; CT1=2000; CT2=500; SW=500
else N=2000; S=6000; K1=10; K2=3; G1=60; G2=12; GA=2000; RS=10000; CT1=30000; CT2=4000; SW=600; fi
log "Mode: $Q"
./system_info.sh > $R/system.txt
make -s
if [ "$Q" = "audit" ]; then ./ctgrind/run_ctgrind.sh; python3 ../../analysis/ctgrind_summary.py ctgrind/logs; exit 0; fi
if [ "$Q" = "replicate" ]; then
  log "replicate: provable primes at 18 odd bit lengths, replicated counts, complete generators"
  ./review_exp provodd 1270 33 35 39 47 63 65 97 127 129 193 255 385 511 769 1023 1025 1279 1535 > $R/provodd.csv
  ./review_exp genrep 1280 200 100000 40 1 > $R/genrep_1280.csv; ./review_exp genrep 2560 60 300000 40 1 > $R/genrep_2560.csv
  ./review_exp saferep 768 40 100000 > $R/saferep_768.csv;     ./review_exp saferep 1280 40 100000 > $R/saferep_1280.csv
  for b in 1024 1536 2048; do K=60; [ $b = 1536 ] && K=40; [ $b = 2048 ] && K=20
    ./fullgen prov $b $K 100000 > $R/prov_${b}_geom.csv; ./fullgen prov $b $K 0 > $R/prov_${b}_td.csv
    ./fullgen prob $b $K 0 > $R/prob_${b}_tdlp.csv;      ./fullgen prob $b $K 64 > $R/prob_${b}_t64.csv; done
  ./fullgen rsa 2048 30 > $R/rsa_2048.csv; ./fullgen rsa 3072 10 > $R/rsa_3072.csv
  log "DONE. Analyze with: python3 ../../analysis/replication_analysis.py $R"; exit 0
fi
if [ "$Q" = "plan" ]; then
  log "confirmatory V0-V3 experiment (docs/protocol_V0-V3_2048_3072_EN.txt); platform name: ${PLATFORM:=platform}"
  ./ossl_seg 2048 2000 30000 confirm4 20 > $R/conf4_${PLATFORM}_2048.csv 2> $R/openssl_version.txt
  ./ossl_seg 3072 2000 60000 confirm4 20 > $R/conf4_${PLATFORM}_3072.csv 2>/dev/null
  log "DONE. Analyze with: python3 ../../analysis/confirm4_analysis.py conf4_<platform>_<bits>.csv ..."; exit 0
fi
if [ "$Q" != "new" ] && [ "$Q" != "ossl" ] && [ "$Q" != "confirm" ]; then
log "1/15 platform constants";       ./diss_bench consts 768 1280 1536 2560 3584 > $R/consts.csv
log "2/15 sweep over B";              ./diss_bench sweep 768 $N > $R/sweep_768.csv; ./diss_bench sweep 1280 $N > $R/sweep_1280.csv
                                       ./diss_bench sweep 2560 $((N/2)) > $R/sweep_2560.csv; ./diss_bench sweep 3584 $SW > $R/sweep_3584.csv
log "3/15 strategies and s0";         ./diss_bench strat 768 > $R/strat_768.csv; ./diss_bench strat 1536 > $R/strat_1536.csv; ./diss_bench s0 1280 > $R/s0_1280.csv
log "4/15 safe primes: sweep";        ./diss_bench safesweep 768 $S > $R/safesweep_768.csv; ./diss_bench safesweep 1280 $S > $R/safesweep_1280.csv
log "5/15 safe-prime generation (longest step)"
                                       ./diss_bench safe 768 $K1 100000 > $R/safe_768.csv; ./diss_bench safe 1280 $K2 100000 > $R/safe_1280.csv
log "6/15 prime generation";          ./diss_bench gen 1280 $G1 100000 40 > $R/gen_1280.csv; ./diss_bench gen 2560 $G2 300000 40 > $R/gen_2560.csv
                                       ./diss_bench gen 1280 $((G1/2)) 100000 64 > $R/gen64_1280.csv; ./diss_bench gen 2560 $((G2*2/3+1)) 300000 128 > $R/gen128_2560.csv
log "7/15 comparison with OpenSSL";    { ./ossl_bench gen 1280 $((G1/2)) ; ./ossl_bench gen 2560 $((G2*2/3+1)) | tail -n +2; } > $R/ossl_gen.csv 2> $R/openssl_version.txt
log "8/15 prime gaps (exploratory)";   ./diss_bench gaps 768 $GA > $R/gaps_768.csv; ./diss_bench gaps 1536 $((GA*3/10)) > $R/gaps_1536.csv
log "9/15 residue distribution";      ./diss_bench resid 448 $RS > $R/resid_448.csv; ./ossl_bench resid 448 $RS > $R/ossl_resid_448.csv 2>/dev/null
log "10/15 timing leakage";           ./ct_experiment 768 1000 $CT1 $R/ct_768_1000.csv > $R/ct_768_1000.txt
                                       ./ct_experiment 1280 3000 $CT2 $R/ct_1280_3000.csv > $R/ct_1280_3000.txt
fi
if [ "$Q" != "ossl" ] && [ "$Q" != "confirm" ]; then
log "11/15 safe primes: batch versus segmented screening"
                                       ./safe_batch 1280 $SBN $SBB > $R/safe_batch_1280.csv
log "12/15 complete RSA key generation (spp and OpenSSL)"
                                       ./rsa_bench 2048 $RK2 > $R/rsa_2048.csv; ./rsa_bench 3072 $RK3 > $R/rsa_3072.csv
fi
if [ "$Q" != "confirm" ]; then
log "13/15 OpenSSL algorithm: stock, replica, with word prefilter (exploratory)"
                                       ./ossl_seg 1024 $OK1 3000,10000,30000 > $R/ossl_seg_1024.csv 2> $R/openssl_version.txt
                                       ./ossl_seg 1536 $OK2 3000,10000,30000 > $R/ossl_seg_1536.csv 2>/dev/null
                                       ./ossl_seg 448 $OKR 3000 resid > $R/ossl_seg_resid448.csv 2>/dev/null
fi
log "14/15 confirmatory experiment (pre-specified protocol): V0 and V3, B = 10000"
                                       ./ossl_seg 1024 $CK 10000 confirm > $R/confirm_1024.csv 2> $R/openssl_version.txt
                                       ./ossl_seg 1536 $CK 10000 confirm > $R/confirm_1536.csv 2>/dev/null
if [ "$Q" = "full" ] || [ "$Q" = "new" ] || [ "$Q" = "complete" ]; then
log "15/15 complete generators (GenPrime, GenProvablePrime, GenRSAPrimes)"
                                       { ./fullgen prob 1024 60 6; ./fullgen prob 1024 60 64 | tail -n +2; ./fullgen prob 1536 40 4 | tail -n +2; ./fullgen prob 1536 40 64 | tail -n +2;
                                         ./fullgen prob 2048 20 3 | tail -n +2; ./fullgen prob 2048 20 64 | tail -n +2; } > $R/complete_prob.csv
                                       { ./fullgen prov 1024 60 100000; ./fullgen prov 1024 60 0 | tail -n +2; ./fullgen prov 1536 40 100000 | tail -n +2; ./fullgen prov 1536 40 0 | tail -n +2;
                                         ./fullgen prov 2048 20 100000 | tail -n +2; ./fullgen prov 2048 20 0 | tail -n +2; } > $R/complete_prov.csv
                                       { ./fullgen rsa 2048 30; ./fullgen rsa 3072 10 | tail -n +2; } > $R/complete_rsa.csv
fi
tar -czf results.tar.gz $R
log "DONE. Raw results: $R/ (archive: results.tar.gz)"
