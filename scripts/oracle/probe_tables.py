"""Compile every generated constant and table cell through the Cyrius compiler and check
each value against tests/golden. This proves src/constants.cyr and src/particle_data.cyr
reproduce the oracle's bits.

Usage, from the repo root:
    python3 -I scripts/oracle/probe_tables.py .

It writes build/oracle-probe/ (gitignored), builds there with the project's cyrius.cyml
(which supplies the stdlib), and exits 1 on any mismatch or any compiler warning.
Extend `checks` as data modules arrive.
"""
import re, subprocess, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
work = root / "build/oracle-probe"
work.mkdir(parents=True, exist_ok=True)

def snake(c): return re.sub(r"(?<!^)([A-Z])", r"_\1", c).upper()

def rows(name):
    for line in (root / f"tests/golden/{name}.tsv").read_text().splitlines():
        if not line.startswith("#"):
            fn, cls, ins, outs, _ = line.split("\t")
            yield fn, cls, ins.split(" "), outs.split(" ")

checks = []  # (cyrius expression, kind, expected token, label)
for fn, cls, ins, outs in rows("constants"):
    name = ins[0][2:]
    checks.append((f"TM_{name}", "f", outs[0], f"const {name}"))

table_of = {
    "mass_mev": "MASS_MEV", "charge": "CHARGE", "generation": "GENERATION",
    "lifetime_seconds": "LIFETIME_SECONDS", "is_stable": "IS_STABLE", "spin": "SPIN",
    "decay_width_mev": "DECAY_WIDTH_MEV", "relative_strength": "RELATIVE_STRENGTH",
    "range_meters": "RANGE_METERS",
}
computed = {"Lepton::decay_width_mev", "Boson::lifetime_seconds"}
for fn, cls, ins, outs in rows("particle"):
    if fn in computed:
        continue
    enum, method = fn.split("::")
    variant = ins[0].split(".", 1)[1]
    E = snake(enum)
    idx = f"TM_{E}_{snake(variant)}"
    if method == "mediator":
        n = int(outs[0][2:])
        checks.append((f"TM_{E}_MEDIATOR_COUNT[{idx}]", "i", f"i:{n}", f"{fn} {variant} count"))
        for k, tok in enumerate(outs[1:]):
            want = "TM_BOSON_" + snake(tok.split(".", 1)[1])
            expr = f"TM_{E}_MEDIATORS[TM_{E}_MEDIATOR_START[{idx}] + {k}]"
            checks.append((expr, "e", want, f"{fn} {variant} [{k}]"))
        continue
    kind = outs[0][0]
    checks.append((f"TM_{E}_{table_of[method]}[{idx}]", kind, outs[0], f"{fn} {variant}"))

src = ['include "src/constants.cyr"', 'include "src/particle_data.cyr"', "",
       "fn main() {", "    var x: f64 = 0.0;", "    var n = 0;"]
for k, (expr, kind, want, label) in enumerate(checks):
    if kind == "f":
        src.append(f"    x = {expr}; fmt_hex(x); println(\" {k}\");")
    else:
        src.append(f"    n = {expr}; fmt_int(n); println(\" {k}\");")
src += ["    return 0;", "}", "", "var r = main();", "sys_exit_group(r);", ""]
probe = work / "tables_probe.cyr"
probe.write_text("\n".join(src))

# Enum member values for `e` checks, from the generated source.
enum_vals = {}
for m in re.finditer(r"^\s+(TM_[A-Z0-9_]+) = (\d+);", (root / "src/particle_data.cyr").read_text(), re.M):
    enum_vals[m.group(1)] = int(m.group(2))

b = subprocess.run(["cyrius", "build", str(probe), str(work / "tables_probe")],
                   capture_output=True, text=True, cwd=root)
diag = [l for l in (b.stdout + b.stderr).splitlines() if re.search(r"warn|error|refus", l, re.I)]
if b.returncode != 0 or diag:
    print("BUILD:", b.returncode, *diag, sep="\n")
    sys.exit(2)
r = subprocess.run([str(work / "tables_probe")], capture_output=True, text=True)
got = {}
for line in r.stdout.splitlines():
    val, k = line.rsplit(" ", 1)
    got[int(k)] = val

fails = 0
for k, (expr, kind, want, label) in enumerate(checks):
    g = got.get(k)
    if kind == "f":
        ok = g is not None and int(g, 16) == int(want[2:], 16)
    elif kind == "i":
        ok = g is not None and int(g) == int(want[2:])
    else:  # enum member
        ok = g is not None and int(g) == enum_vals[want]
    if not ok:
        fails += 1
        print(f"FAIL {label}: {expr} -> {g}, want {want}")
print(f"checked {len(checks)} values through the Cyrius compiler; failures: {fails}")
sys.exit(1 if fails else 0)
