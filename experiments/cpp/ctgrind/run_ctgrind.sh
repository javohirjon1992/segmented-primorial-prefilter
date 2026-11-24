#!/bin/bash
# ctgrind matrix: GCC/Clang x -O2/-O3 x (without / with inline-assembly barrier) x 5 variants x 2 sizes.
# Usage (requires Valgrind, g++, clang++, GMP):  ./run_ctgrind.sh   -> ctgrind_results.csv and logs/
# Code locations as reported in the paper: python3 ../../../analysis/ctgrind_summary.py logs
cd "$(dirname "$0")"
SPP=${SPP:-../../../spp}
for cc in g++ clang++; do for O in O2 O3; do for bar in 0 1; do
  D=""; [ $bar = 1 ] && D="-DCT_ASM_BARRIER"; out=ctg_${cc%%+*}_${O}_b$bar
  $cc -std=c++17 -$O -g $D -I $SPP/include ctgrind_check.cpp $SPP/src/spp.cpp -lgmp -o $out
done; done; done
mkdir -p logs; echo "compiler,opt,asm_barrier,variant,bits,B,N,errors,contexts,cond_contexts,addr_contexts" > ctgrind_results.csv
for cfg in "768 1000" "1280 3000"; do set -- $cfg; bits=$1; B=$2
 for bin in ctg_g_O2_b0 ctg_g_O3_b0 ctg_clang_O2_b0 ctg_clang_O3_b0 ctg_g_O2_b1 ctg_g_O3_b1 ctg_clang_O2_b1 ctg_clang_O3_b1; do
  for v in TD_early SEG_early SEG_full CT_TD CT_TREE; do
   log=logs/${bin}_${v}_${bits}.log
   valgrind --error-limit=no --num-callers=4 ./$bin $v $bits $B 30 > $log 2>&1
   e=$(grep "ERROR SUMMARY" $log | sed -E 's/.*SUMMARY: ([0-9]+) errors from ([0-9]+) contexts.*/\1,\2/')
   cc=$(grep -c "Conditional jump or move" $log); ac=$(grep -c "Use of uninitialised value of size" $log)
   IFS=_ read -r _ comp opt b <<< "$bin"
   echo "$comp,$opt,${b#b},$v,$bits,$B,30,$e,$cc,$ac" >> ctgrind_results.csv
  done; done; done
