"""Summarizes the Valgrind logs of experiments/cpp/ctgrind/run_ctgrind.sh as in the paper: for every build and variant, the
number of distinct code locations (instruction addresses of the innermost frame) with a secret-dependent conditional
branch, conditional move, or memory address.

    python3 analysis/ctgrind_summary.py [logs_dir] [output.csv]
Default logs: data/raw/ctgrind/logs (unpack data/raw/ctgrind/valgrind_logs.tar.gz first);
default output: data/tables/ctgrind_code_locations.csv
"""
import csv, os, re, sys, tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'data', 'raw', 'ctgrind', 'logs')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'data', 'tables', 'ctgrind_code_locations.csv')
if not os.path.isdir(L):
    tgz = os.path.join(ROOT, 'data', 'raw', 'ctgrind', 'valgrind_logs.tar.gz')
    with tarfile.open(tgz) as t:
        try:
            t.extractall(os.path.dirname(L), filter='data')
        except TypeError:          # Python < 3.12
            t.extractall(os.path.dirname(L))
NAME = re.compile(r'ctg_(g|clang)_(O[23])_b([01])_(TD_early|SEG_early|SEG_full|CT_TD|CT_TREE)_(\d+)\.log$')
rows = []
for f in sorted(os.listdir(L)):
    m = NAME.match(f)
    if not m:
        continue
    lines = open(os.path.join(L, f)).read().split('\n')
    branch, addr = set(), set()
    for i, line in enumerate(lines[:-1]):
        kind = branch if 'Conditional jump or move' in line else addr if 'Use of uninitialised value of size' in line else None
        if kind is None:
            continue
        a = re.search(r'(?:at|by) 0x([0-9A-F]+):', lines[i + 1])
        if a:
            kind.add(a.group(1))
    comp = 'GCC 13.3' if m.group(1) == 'g' else 'Clang 18.1'
    rows.append(dict(compiler=comp, opt='-' + m.group(2), asm_barrier=int(m.group(3)), variant=m.group(4), bits=int(m.group(5)),
                     branch_locations=len(branch), address_locations=len(addr), locations=len(branch | addr)))
order = ['TD_early', 'SEG_early', 'SEG_full', 'CT_TD', 'CT_TREE']
rows.sort(key=lambda r: (order.index(r['variant']), r['asm_barrier'], r['compiler'], r['opt'], r['bits']))
with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"{'variant':10s} " + ' '.join(f'{c:>16s}' for c in ('GCC -O2', 'GCC -O3', 'Clang -O2', 'Clang -O3', 'barrier (all 4)')))
for v in order:
    cells = []
    for comp, opt in (('GCC 13.3', '-O2'), ('GCC 13.3', '-O3'), ('Clang 18.1', '-O2'), ('Clang 18.1', '-O3')):
        x = {r['bits']: r['locations'] for r in rows if r['variant'] == v and r['compiler'] == comp and r['opt'] == opt and r['asm_barrier'] == 0}
        cells.append(f"{x.get(768, '-')} / {x.get(1280, '-')}")
    b = {(r['bits']): set() for r in rows}
    bar = [r for r in rows if r['variant'] == v and r['asm_barrier'] == 1]
    vals = {bits: sorted({r['locations'] for r in bar if r['bits'] == bits}) for bits in (768, 1280)}
    cells.append(' / '.join(','.join(map(str, vals[k])) for k in (768, 1280)))
    print(f'{v:10s} ' + ' '.join(f'{c:>16s}' for c in cells))
print(f'(768 bits, B = 1000 / 1280 bits, B = 3000; 30 candidates each) -> {OUT}')
