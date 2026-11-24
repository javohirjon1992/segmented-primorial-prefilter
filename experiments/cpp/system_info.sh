#!/usr/bin/env bash
# Records the hardware and software parameters of a platform (Table "platforms" of the paper).
#   ./system_info.sh > results/system.txt
# Works on Linux and on Windows with MSYS2; missing tools are skipped.
echo "date:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "os:          $(uname -srmo 2>/dev/null || uname -a)"
if command -v lscpu >/dev/null; then
  lscpu | grep -E "^(Model name|Architecture|CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|CPU max MHz|CPU MHz|L1d cache|L2 cache|L3 cache)" | sed 's/^/cpu:         /'
elif [ -r /proc/cpuinfo ]; then
  echo "cpu:         $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2-)"; echo "cpu:         logical CPUs $(grep -c ^processor /proc/cpuinfo)"
fi
command -v powershell.exe >/dev/null && powershell.exe -NoProfile -Command \
  "Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors,MaxClockSpeed | Format-List" 2>/dev/null | sed '/^\s*$/d;s/^/cpu:         /'
if [ -r /proc/meminfo ]; then echo "memory:      $(grep MemTotal /proc/meminfo | awk '{printf "%.1f GB", $2/1048576}')"; fi
command -v g++ >/dev/null && echo "g++:         $(g++ --version | head -1)"
command -v clang++ >/dev/null && echo "clang++:     $(clang++ --version | head -1)"
printf '#include <gmp.h>\n#include <stdio.h>\nint main(void){printf("%%s\\n", gmp_version);return 0;}\n' > /tmp/gmpv_$$.c
if gcc /tmp/gmpv_$$.c -lgmp -o /tmp/gmpv_$$ 2>/dev/null; then echo "gmp:         $(/tmp/gmpv_$$)"; fi
rm -f /tmp/gmpv_$$ /tmp/gmpv_$$.c
command -v openssl >/dev/null && echo "openssl:     $(openssl version)"
command -v python3 >/dev/null && echo "python:      $(python3 --version 2>&1)" && python3 - <<'PY' 2>/dev/null
import importlib
for m in ("numpy", "scipy", "pandas", "matplotlib", "gmpy2"):
    try:
        print(f"python pkg:  {m} {importlib.import_module(m).__version__}")
    except Exception:
        pass
PY
command -v cpupower >/dev/null && cpupower frequency-info -p 2>/dev/null | sed 's/^/governor:    /'
exit 0
