"""Independent check of tests/golden/*.tsv (format 1). Usage, from the repo root: python3 -I scripts/oracle/validate.py .

1. Every row parses: 5 tab fields, a known class, well-formed tokens.
2. constants.tsv matches the decimal literals in rust-old/src/constants.rs.
3. Class-B relativity rows recomputed in Python with Rust's operation order are bit-identical.
4. libm.tsv: sqrt vs math.sqrt (correctly rounded); cbrt of perfect cubes is exact.
5. atomic.tsv: the class-B Rydberg lines and fine-structure levels recomputed in Python with
   Rust's operation order are bit-identical; every configuration holds Z electrons, never
   overfills a subshell, and is the Madelung filling except at the 20 NIST exceptions.
"""
import math, re, struct, sys
from pathlib import Path

root = Path(sys.argv[1])
gold = root / "tests/golden"

def bits_to_f(h): return struct.unpack("<d", struct.pack("<Q", int(h, 16)))[0]
def f_to_bits(x): return struct.unpack("<Q", struct.pack("<d", x))[0]

def parse_tok(t):
    kind, _, val = t.partition(":")
    if kind == "f":
        assert re.fullmatch(r"[0-9a-f]{16}", val), t
        return ("f", bits_to_f(val), int(val, 16))
    if kind == "i":
        return ("i", int(val), None)
    if kind in ("e", "s", "E", "v", "n"):
        return (kind, val, None)
    raise AssertionError(f"bad token {t!r}")

def rows(name):
    out = []
    for n, line in enumerate((gold / f"{name}.tsv").read_text().splitlines(), 1):
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        assert len(cols) == 5, (name, n, line)
        fn, cls, ins, outs, _note = cols
        assert cls in "ABCD", (name, n)
        pi = [] if ins == "-" else [parse_tok(t) for t in ins.split(" ")]
        po = [] if outs == "-" else [parse_tok(t) for t in outs.split(" ")]
        out.append((fn, cls, pi, po))
    return out

def same(x, want_bits):
    if math.isnan(x):
        return math.isnan(bits_to_f(f"{want_bits:016x}"))
    return f_to_bits(x) == want_bits

fails = 0
def check(ok, msg):
    global fails
    if not ok:
        fails += 1
        if fails <= 20:
            print("FAIL", msg)

counts = {}
for p in sorted(gold.glob("*.tsv")):
    counts[p.stem] = len(rows(p.stem))
print("parsed:", counts)

# 2. constants vs the Rust source literals
src = (root / "rust-old/src/constants.rs").read_text()
lits = {}
for m in re.finditer(r"pub const ([A-Z0-9_]+): f64 = ([^;]+);", src):
    expr = m.group(2).replace("_", "")
    lits[m.group(1)] = eval(expr, {"__builtins__": {}})  # literals and 1.0 / x only
for fn, cls, ins, outs in rows("constants"):
    name = ins[0][1]
    check(name in lits, f"constant {name} not in source")
    check(same(lits[name], outs[0][2]), f"constant {name}: {lits[name]!r} vs {outs[0][1]!r}")
print("constants checked:", len(lits))

# 3. relativity class-B rows
HBAR_C = lits["HBAR_C_MEV_FM"]; C = lits["C"]
def sqrt(x): return math.nan if x < 0 else math.sqrt(x)
def lorentz_gamma(b):
    b = abs(b)
    return math.inf if b >= 1.0 else 1.0 / sqrt((1.0 - b) * (1.0 + b))
def gamma_to_beta(g):
    if g <= 1.0: return 0.0
    if math.isinf(g): return 1.0
    return sqrt((g - 1.0) * (g + 1.0)) / g
def velocity_addition(b1, b2):
    d = 1.0 + b1 * b2
    if abs(b1) > 1.0 or abs(b2) > 1.0 or d == 0.0: return math.nan
    return (b1 + b2) / d
def rel_e(m, p): return sqrt(m * m + p * p)
def rel_p(m, t): return sqrt(t * t + 2.0 * t * m)
def debroglie(p): return math.inf if p <= 0.0 else 2.0 * math.pi * HBAR_C / p
impls = {
    "lorentz_gamma": lambda a: [lorentz_gamma(*a)],
    "gamma_to_beta": lambda a: [gamma_to_beta(*a)],
    "velocity_addition": lambda a: [velocity_addition(*a)],
    "relativistic_energy": lambda a: [rel_e(*a)],
    "relativistic_momentum": lambda a: [rel_p(*a)],
    "de_broglie_wavelength_fm": lambda a: [debroglie(*a)],
    "velocity_to_beta": lambda a: [a[0] / C],
    "beta_to_velocity": lambda a: [a[0] * C],
    "FourMomentum::from_mass_and_momentum": lambda a: [rel_e(*a), 0.0, 0.0, a[1]],
    "FourMomentum::at_rest": lambda a: [a[0], 0.0, 0.0, 0.0],
    "FourMomentum::new": lambda a: list(a),
    "FourMomentum::add": lambda a: [a[0] + a[4], a[1] + a[5], a[2] + a[6], a[3] + a[7]],
}
n = 0
for fn, cls, ins, outs in rows("relativity"):
    if fn in impls:
        got = impls[fn]([t[1] for t in ins])
        ok = len(got) == len(outs) and all(same(g, o[2]) for g, o in zip(got, outs))
        check(ok, f"{fn}{tuple(t[1] for t in ins)} -> {got} vs {[o[1] for o in outs]}")
        n += 1
print("relativity rows recomputed:", n)

# 4. libm spot checks
n = 0
for fn, cls, ins, outs in rows("libm"):
    x = ins[0][1]
    if fn == "libm::sqrt":
        check(same(sqrt(x) if not math.isnan(x) else math.nan, outs[0][2]), f"sqrt({x!r})")
        n += 1
    if fn == "libm::cbrt" and math.isfinite(x) and x == int(x) and abs(x) < 2**53:
        r = round(abs(x) ** (1 / 3))
        for k in (r - 1, r, r + 1):
            if k ** 3 == abs(x):
                check(outs[0][1] == math.copysign(k, x), f"cbrt({x!r}) = {outs[0][1]!r}, want {k}")
                n += 1
print("libm rows checked:", n)

# 5. atomic
RY = lits["RYDBERG"]; RY_EV = lits["RYDBERG_EV"]; FS = lits["FINE_STRUCTURE"]
def line_nm(z, n1, n2):
    if z == 0: return "E:InvalidAtomicNumber"
    if n1 == 0 or n2 == 0 or n1 == n2: return "E:InvalidQuantumNumbers"
    lo, up = min(n1, n2), max(n1, n2)
    zf, lf, uf = float(z), float(lo), float(up)
    return 1e9 / (((RY * zf) * zf) * ((1.0 / (lf * lf)) - (1.0 / (uf * uf))))
def level_ev(z, n, tj):
    if n == 0 or tj % 2 == 0 or tj > 2 * n - 1: return "E:InvalidQuantumNumbers"
    zf, nf, jph = float(z), float(n), (float(tj) + 1.0) / 2.0
    e0 = ((-RY_EV * zf) * zf) / (nf * nf)
    az = FS * zf
    return e0 * (1.0 + (((az * az) / (nf * nf)) * ((nf / jph) - 0.75)))
n = 0
for fn, cls, ins, outs in rows("atomic"):
    a = [t[1] for t in ins]
    impl = {"spectral_line_nm": line_nm, "hydrogen_level_energy_ev": level_ev}.get(fn)
    if impl is None:
        continue
    got = impl(*a)
    if isinstance(got, str):
        check(outs[0][0] == "E" and got == "E:" + outs[0][1], f"{fn}{tuple(a)} -> {got}")
    else:
        check(outs[0][0] == "f" and same(got, outs[0][2]), f"{fn}{tuple(a)} -> {got!r} vs {outs[0][1]!r}")
    n += 1
print("atomic line/level rows recomputed:", n)

MAX = {"S": 2, "P": 6, "D": 10, "F": 14}
L = {"S": 0, "P": 1, "D": 2, "F": 3}
# Madelung: by n + l, then n.
order = sorted(((n_, o) for n_ in range(1, 9) for o in "SPDF" if L[o] < n_), key=lambda x: (x[0] + L[x[1]], x[0]))
def madelung(z):
    out, left = [], z
    for n_, o in order:
        if left == 0: break
        k = min(left, MAX[o]); out.append((n_, o, k)); left -= k
    return out
NIST = {24, 29, 41, 42, 44, 45, 46, 47, 57, 58, 64, 78, 79, 89, 90, 91, 92, 93, 96, 103}
n = 0
for fn, cls, ins, outs in rows("atomic"):
    if fn != "electron_configuration":
        continue
    z = ins[0][1]
    if not 1 <= z <= 118:
        check(outs[0][0] == "E", f"electron_configuration({z}) should be an Err")
        continue
    k = outs[0][1]
    cfg = [(outs[1 + 3 * j][1], outs[2 + 3 * j][1].split(".")[1], outs[3 + 3 * j][1]) for j in range(int(k))]
    check(sum(e for _, _, e in cfg) == z, f"Z={z}: {sum(e for _, _, e in cfg)} electrons")
    check(all(0 < e <= MAX[o] for _, o, e in cfg), f"Z={z}: a subshell is empty or overfilled")
    check((cfg == madelung(z)) == (z not in NIST), f"Z={z}: Madelung filling iff not a NIST exception")
    n += 1
print("atomic configurations checked:", n)
print("FAILURES:", fails)
sys.exit(1 if fails else 0)
