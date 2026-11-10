"""Generate tests/nist_levels_data.cyr from the NIST ASD snapshots in tests/reference/nist_asd/.

Usage, from the repo root:
    python3 -I scripts/reference/nist_levels.py .           # write the file
    python3 -I scripts/reference/nist_levels.py . --check   # exit 1 if it is out of date

The snapshots are NIST ASD v5.12 energy levels of He II (4He+) and Li III (7Li2+); see
tests/reference/nist_asd/README.md for the source and the query. Every row is checked before
anything is written:
- the file's SHA-256 matches the one recorded below (and in the README);
- the configuration "<n><l letter>" agrees with the term "2<L>[*]": the same L, l < n, and the
  odd-parity star exactly when l is odd;
- 2J is odd and J = l +- 1/2;
- each (n, l, J) appears once, and every (l, J) of n is present for n <= 11;
- the first row is the ground level 1s 2S1/2 at 0, and the last is the ionization limit.

Each value is written as an exact positional decimal of the double Python parses from NIST's
text (generated .cyr files never use exponent notation), together with that double's IEEE-754
bit pattern as a decimal integer. tests/hydrogenic_nist.tcyr compares the two, so every run
proves the compiler read each literal to the same bits.
"""
import hashlib, re, struct, sys
from decimal import Decimal
from pathlib import Path

IONS = [
    # (spectrum, file, Z, A, SHA-256 of the snapshot)
    ("He II", "he_ii_levels.tsv", 2, 4, "b5af0e773fc79e6af11c76046e85c3419ad0b149dd3f85c8b6d933610a1e870a"),
    ("Li III", "li_iii_levels.tsv", 3, 7, "cd1eec6c1d617b114fcbe50e97ac9b42ac54242b4134b5ec0d9191a05b0554a2"),
]
L_LETTERS = "SPDFGHIKLMNOQRTUVWXY"  # spectroscopic notation: J is skipped, and P and S are not reused
FULL_N_MAX = 11
OUT = "tests/nist_levels_data.cyr"
HEADER = ["Configuration", "Term", "J", "Prefix", "Level (cm-1)", "Suffix", "Uncertainty (cm-1)", "Reference"]


def fail(msg):
    print(f"nist_levels: {msg}", file=sys.stderr)
    sys.exit(1)


def positional(x):
    """The exact positional decimal of repr(x), the shortest string that round-trips x."""
    s = format(Decimal(repr(x)), "f")
    if "." not in s:
        s += ".0"
    return s


def bits(x):
    return struct.unpack("<q", struct.pack("<d", x))[0]


def parse(root, spectrum, name, z, a, sha):
    path = root / "tests/reference/nist_asd" / name
    raw = path.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha:
        fail(f"{name}: SHA-256 {got}, expected {sha}. A new snapshot needs its hash updated here and in the README")
    lines = raw.decode("utf-8").splitlines()
    if lines[0].split("\t") != HEADER:
        fail(f"{name}: unexpected header {lines[0]!r}")
    levels, limit, seen = [], None, set()
    for k, line in enumerate(lines[1:], start=2):
        cols = [c.strip('"') for c in line.split("\t")]
        if len(cols) != len(HEADER):
            fail(f"{name}:{k}: {len(cols)} columns")
        conf, term, jtxt, prefix, level, suffix, unc, ref = cols
        if term == "Limit":
            if limit is not None or k != len(lines):
                fail(f"{name}:{k}: the limit must be the last row, once")
            limit = (float(level), float(unc), ref)
            continue
        m = re.fullmatch(r"(\d+)([a-z])", conf)
        t = re.fullmatch(r"2([A-Z])(\*?)", term)
        j = re.fullmatch(r"(\d+)/2", jtxt)
        if not (m and t and j):
            fail(f"{name}:{k}: can't parse {conf!r} {term!r} {jtxt!r}")
        n, l = int(m.group(1)), L_LETTERS.index(m.group(2).upper())
        if L_LETTERS.index(t.group(1)) != l:
            fail(f"{name}:{k}: configuration {conf} and term {term} disagree on L")
        if (t.group(2) == "*") != (l % 2 == 1):
            fail(f"{name}:{k}: parity star wrong for l = {l}")
        two_j = int(j.group(1))
        if l >= n or abs(two_j - 2 * l) != 1:
            fail(f"{name}:{k}: (n, l, 2J) = ({n}, {l}, {two_j}) is not a one-electron level")
        if (n, l, two_j) in seen:
            fail(f"{name}:{k}: ({n}, {l}, {two_j}) repeated")
        seen.add((n, l, two_j))
        if (prefix, suffix) not in [("", ""), ("(", ")")]:
            fail(f"{name}:{k}: unexpected brackets {prefix!r} {suffix!r}")
        levels.append((n, l, two_j, float(level), float(unc), conf, term, jtxt, ref))
    if limit is None:
        fail(f"{name}: no ionization limit")
    if levels[0][:4] != (1, 0, 1, 0.0):
        fail(f"{name}: the first row is not the ground level 1s 2S1/2 at 0")
    for n in range(1, FULL_N_MAX + 1):
        want = {(n, l, 2 * l + s) for l in range(n) for s in (-1, 1) if 2 * l + s > 0}
        if not want <= seen:
            fail(f"{name}: n = {n} is missing {sorted(want - seen)}")
    return spectrum, z, a, levels, limit


def table(lines, kind, name, values, comments):
    lines.append(f"var {name}: {kind}[{len(values)}] = {{")
    for v, c in zip(values, comments):
        lines.append(f"    {v},  # {c}")
    lines.append("};")
    lines.append("")


def render(ions):
    out = [
        "# nist_levels_data.cyr — NIST ASD energy levels of ⁴He⁺ (He II) and ⁷Li²⁺ (Li III), the reference",
        "# data of tests/hydrogenic_nist.tcyr.",
        "#",
        "# GENERATED by scripts/reference/nist_levels.py from tests/reference/nist_asd/*.tsv. Do not edit.",
        "# Regenerate: python3 -I scripts/reference/nist_levels.py .",
        "#",
        "# Source: Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024). NIST Atomic Spectra",
        "# Database (ver. 5.12), https://physics.nist.gov/asd, retrieved 2026-11-09. n = 1, 2 and the limit:",
        "# Yerokhin & Shabaev, J. Phys. Chem. Ref. Data 44, 033103 (2015). n >= 3: Erickson, J. Phys. Chem.",
        "# Ref. Data 6, 831 (1977), scaled to CODATA 2018 by NIST. tests/reference/nist_asd/README.md has",
        "# the details.",
        "#",
        "# Levels are cm^-1 above the ground level 1s 2S1/2. Each f64 column has a _BITS twin, the IEEE-754",
        "# bit pattern of the same double as a decimal integer, so the test can prove the literal was read",
        "# exactly.",
        "",
        f"const _REF_ION_COUNT = {len(ions)};",
        f"const _REF_FULL_N_MAX = {FULL_N_MAX};",
        "",
    ]
    names = [spec for spec, *_ in ions]
    starts, k = [], 0
    for _, _, _, levels, _ in ions:
        starts.append(k)
        k += len(levels)
    out.append(f"const _REF_LEVEL_COUNT = {k};")
    out.append("")
    out.append("# Per ion, in the order He II, Li III.")
    table(out, "u8", "_REF_ION_Z", [z for _, z, *_ in ions], names)
    table(out, "u8", "_REF_ION_A", [a for _, _, a, *_ in ions], names)
    table(out, "i64", "_REF_ION_START", starts, names)
    table(out, "i64", "_REF_ION_LEVELS", [len(lv) for *_, lv, _ in ions], names)
    lims = [lim[0] for *_, lim in ions]
    table(out, "f64", "_REF_ION_LIMIT_CM", [positional(x) for x in lims],
          [f"{s} ionization limit ({lim[2]})" for s, *_, lim in ions])
    table(out, "i64", "_REF_ION_LIMIT_CM_BITS", [bits(x) for x in lims], names)

    rows = [(spec, lv) for spec, _, _, levels, _ in ions for lv in levels]
    label = [f"{spec} {lv[5]} {lv[6]} {lv[7]}" + (f" ({lv[8]})" if lv[8] else "") for spec, lv in rows]
    out.append("# Per level, ion by ion, in NIST's order.")
    table(out, "u8", "_REF_LEVEL_N", [lv[0] for _, lv in rows], label)
    table(out, "u8", "_REF_LEVEL_L", [lv[1] for _, lv in rows], label)
    table(out, "u8", "_REF_LEVEL_TWO_J", [lv[2] for _, lv in rows], label)
    table(out, "f64", "_REF_LEVEL_CM", [positional(lv[3]) for _, lv in rows], label)
    table(out, "i64", "_REF_LEVEL_CM_BITS", [bits(lv[3]) for _, lv in rows], label)
    table(out, "f64", "_REF_LEVEL_UNC_CM", [positional(lv[4]) for _, lv in rows], label)
    table(out, "i64", "_REF_LEVEL_UNC_CM_BITS", [bits(lv[4]) for _, lv in rows], label)
    while out[-1] == "":
        out.pop()
    return "\n".join(out) + "\n"


def main():
    if len(sys.argv) < 2 or sys.argv[2:] not in ([], ["--check"]):
        fail("usage: nist_levels.py <repo root> [--check]")
    root = Path(sys.argv[1]).resolve()
    ions = [parse(root, *ion) for ion in IONS]
    text = render(ions)
    path = root / OUT
    if sys.argv[2:] == ["--check"]:
        if not path.exists() or path.read_text() != text:
            fail(f"{OUT} is out of date: run python3 -I scripts/reference/nist_levels.py .")
        print(f"nist_levels: {OUT} is current ({sum(len(i[3]) for i in ions)} levels)")
        return
    path.write_text(text)
    print(f"nist_levels: wrote {OUT} ({sum(len(i[3]) for i in ions)} levels)")


main()
