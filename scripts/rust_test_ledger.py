"""The Rust-test ledger: every one of the oracle's tests is ported, or recorded with a reason.

Port plan M7 (the mela pattern). Usage, from the repo root:
    python3 -I scripts/rust_test_ledger.py .           # write docs/development/rust-test-ledger.md
    python3 -I scripts/rust_test_ledger.py . --check   # exit 1 if the ledger is stale or a test is unaccounted

It lists each `#[test]` in rust-old/src/**.rs and rust-old/tests/*.rs, and each doc test, and
looks for its port: a `fn _t_<name>(` in the Cyrius test file for that Rust file. A Rust test with
no port must have an entry in REASONS, and a REASONS entry for a test that is ported, or that no
longer exists, is an error too.
"""
import re, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
check = "--check" in sys.argv[2:]
out_path = root / "docs/development/rust-test-ledger.md"

# Rust file -> the Cyrius test file its tests are ported into (None: nothing is ported).
FILES = {
    "src/particle.rs": "tests/particle.tcyr",
    "src/relativity.rs": "tests/relativity.tcyr",
    "src/error.rs": None,
    "src/nucleus.rs": "tests/nucleus.tcyr",
    "src/decay.rs": "tests/decay.tcyr",
    "src/timekeeping.rs": "tests/timekeeping.tcyr",
    "src/atomic.rs": "tests/atomic.tcyr",
    "src/reaction.rs": "tests/reaction.tcyr",
    "src/scattering.rs": "tests/scattering.tcyr",
    "src/optics.rs": "tests/optics.tcyr",
    "src/bridge.rs": None,
    "src/integration/soorat.rs": None,
    "tests/integration.rs": "tests/integration.tcyr",
}

# The two doc tests: (Rust file, what it is) -> (Cyrius file, test name).
DOC_TESTS = {
    ("src/lib.rs", "the Quick Start example"): ("tests/integration.tcyr", "doc_lib_quick_start"),
    ("src/timekeeping.rs", "AtomicInstant's example"): ("tests/integration.tcyr", "doc_atomic_instant"),
}

SERDE = "n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7)"
EXCHANGE = ("deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius "
            "consumer exists yet (port plan M7)")
BRIDGE = ("deferred: bridge.rs waits on its consumers (bijli, kimiya) being in Cyrius "
          "(port plan §1, *Deferred until the consumers exist*)")
SOORAT = "deferred: the soorat integration waits on soorat (port plan §1)"

# (Rust file, test) -> why it has no port. Bridge and soorat tests are covered by file below.
REASONS = {
    ("src/particle.rs", "serde_roundtrip_quark"): SERDE,
    ("src/particle.rs", "serde_roundtrip_lepton"): SERDE,
    ("src/particle.rs", "serde_roundtrip_boson"): SERDE,
    ("src/particle.rs", "serde_roundtrip_force"): SERDE,
    ("src/relativity.rs", "serde_roundtrip_four_momentum"): EXCHANGE,
    ("src/error.rs", "serde_roundtrip_error"): SERDE,
    ("src/error.rs", "serde_roundtrip_error_with_string"): SERDE,
    ("src/nucleus.rs", "serde_roundtrip"): EXCHANGE,
    ("src/nucleus.rs", "serde_roundtrip_shell_level"): SERDE,
    ("src/nucleus.rs", "serde_roundtrip_nuclear_moments"): SERDE,
    ("src/nucleus.rs", "serde_roundtrip_superallowed_decay"): SERDE,
    ("src/nucleus.rs", "corrected_ft_value_returns_average"): (
        "partial: its `corrected_ft_value` assert is n/a (deprecated in 1.3.0, not carried over, port plan §1); "
        "its `superallowed_average_ft` assert is ported as "
        "`_t_corrected_ft_value_returns_average_superallowed_part`"),
    ("src/decay.rs", "serde_roundtrip_decay_mode"): SERDE,
    ("src/decay.rs", "serde_roundtrip_isotope"): EXCHANGE,
    ("src/decay.rs", "serde_roundtrip_isomeric_transition"): SERDE,
    ("src/timekeeping.rs", "quality_factors_reasonable"): (
        "n/a: every assert calls `FrequencyStandard::quality_factor`, deprecated in 1.3.0 and not carried "
        "over (port plan §1); `quality_factor_for_linewidth` is ported and golden-checked"),
    ("src/timekeeping.rs", "serde_roundtrip_frequency_standard"): SERDE,
    ("src/timekeeping.rs", "serde_roundtrip_time_scale"): SERDE,
    ("src/timekeeping.rs", "serde_roundtrip_atomic_instant"): EXCHANGE,
    ("src/timekeeping.rs", "atomic_instant_serde_normalizes"): EXCHANGE,
    ("src/atomic.rs", "serde_roundtrip_orbital_type"): SERDE,
    ("src/atomic.rs", "serde_roundtrip_quantum_numbers"): SERDE,
    ("src/atomic.rs", "serde_roundtrip_orbital_filling"): SERDE,
    ("src/atomic.rs", "serde_roundtrip_transition_type"): SERDE,
    ("src/reaction.rs", "serde_roundtrip_reaction"): EXCHANGE,
    ("src/reaction.rs", "serde_roundtrip_thermal_cross_section"): SERDE,
    ("src/reaction.rs", "serde_roundtrip_nucleosynthesis"): SERDE,
    ("src/reaction.rs", "serde_roundtrip_fission_yield"): SERDE,
    ("src/optics.rs", "serde_roundtrip_spectral_line"): SERDE,
    ("tests/integration.rs", "serde_roundtrip_nucleus"): EXCHANGE,
    ("tests/integration.rs", "serde_roundtrip_quark"): SERDE,
    ("tests/integration.rs", "serde_roundtrip_decay_mode"): SERDE,
    ("tests/integration.rs", "serde_roundtrip_quantum_numbers"): SERDE,
}
WHOLE_FILE = {"src/bridge.rs": BRIDGE, "src/integration/soorat.rs": SOORAT}


def rust_tests(rel):
    lines = (root / "rust-old" / rel).read_text().split("\n")
    names = []
    for i, line in enumerate(lines):
        if line.strip().startswith("#[test]"):
            for nxt in lines[i + 1 : i + 6]:
                m = re.match(r"\s*(?:pub )?fn (\w+)", nxt)
                if m:
                    names.append(m.group(1))
                    break
    return names


def cyrius_tests(rel):
    if rel is None:
        return set()
    return set(re.findall(r"^fn _t_(\w+)\(", (root / rel).read_text(), re.M))


errors = []
rows = []        # (rust file, cyrius file, tests, ported, n/a, deferred, partial)
notes = []       # (rust file, test, reason)
for rel, cy in FILES.items():
    tests = rust_tests(rel)
    ported_names = cyrius_tests(cy)
    counts = {"ported": 0, "n/a": 0, "deferred": 0, "partial": 0}
    for t in tests:
        if t in ported_names:
            counts["ported"] += 1
            if (rel, t) in REASONS:
                errors.append(f"{rel}: {t} is ported but has a REASONS entry")
            continue
        reason = REASONS.get((rel, t)) or WHOLE_FILE.get(rel)
        if reason is None:
            errors.append(f"{rel}: {t} is neither ported nor recorded")
            continue
        kind = reason.split(":")[0]
        counts[kind] += 1
        if rel not in WHOLE_FILE:
            notes.append((rel, t, reason))
    for (r, t) in REASONS:
        if r == rel and t not in tests:
            errors.append(f"{rel}: REASONS names {t}, which is not a test")
    rows.append((rel, cy, len(tests), counts))
doc_rows = []
for (rel, what), (cy, name) in DOC_TESTS.items():
    if name not in cyrius_tests(cy):
        errors.append(f"doc test {rel} ({what}): {cy} has no _t_{name}")
    doc_rows.append((rel, what, cy, name))

total = sum(r[2] for r in rows) + len(DOC_TESTS)
sums = {k: sum(r[3][k] for r in rows) for k in ("ported", "n/a", "deferred", "partial")}
sums["ported"] += len(DOC_TESTS)

md = []
md.append("# Rust-test ledger")
md.append("")
md.append("> **Generated** by `python3 -I scripts/rust_test_ledger.py .`; `--check` gates it (port plan M7, the mela")
md.append("> pattern). Every test of the frozen Rust oracle (`rust-old/`, tanmatra 1.3.0) is ported to a Cyrius test,")
md.append("> or recorded here with its reason. The golden corpus (`tests/golden/`) checks every ported function on its")
md.append("> grid besides; this ledger is about the oracle's own assertions.")
md.append("")
md.append(f"**{total} Rust tests:** {sums['ported']} ported, {sums['partial']} partly ported, {sums['n/a']} not "
          f"applicable, {sums['deferred']} deferred.")
md.append("")
md.append("| Rust file | Cyrius tests | Rust tests | Ported | Partial | N/A | Deferred |")
md.append("|---|---|---:|---:|---:|---:|---:|")
for rel, cy, n, c in rows:
    md.append(f"| `{rel}` | {'`' + cy + '`' if cy else '—'} | {n} | {c['ported']} | {c['partial']} | {c['n/a']} | "
              f"{c['deferred']} |")
md.append(f"| doc tests (`src/lib.rs`, `src/timekeeping.rs`) | `tests/integration.tcyr` | {len(DOC_TESTS)} | "
          f"{len(DOC_TESTS)} | 0 | 0 | 0 |")
md.append(f"| **Total** | | **{total}** | **{sums['ported']}** | **{sums['partial']}** | **{sums['n/a']}** | "
          f"**{sums['deferred']}** |")
md.append("")
md.append("## Doc tests")
md.append("")
for rel, what, cy, name in doc_rows:
    md.append(f"- `rust-old/{rel}`, {what}: `{cy}` `_t_{name}`.")
md.append("")
md.append("## Tests without a full port")
md.append("")
md.append("| Rust file | Test | Status |")
md.append("|---|---|---|")
for rel, t, reason in notes:
    md.append(f"| `{rel}` | `{t}` | {reason} |")
for rel, reason in WHOLE_FILE.items():
    md.append(f"| `{rel}` | all {len(rust_tests(rel))} | {reason} |")
md.append("")
text = "\n".join(md)

for e in errors:
    print("LEDGER", e)
if check:
    old = out_path.read_text() if out_path.exists() else ""
    if old != text:
        print(f"STALE    {out_path.relative_to(root)} (regenerate: python3 -I scripts/rust_test_ledger.py .)")
        errors.append("stale")
else:
    out_path.write_text(text)
print(f"rust-test ledger: {total} tests — {sums['ported']} ported, {sums['partial']} partial, {sums['n/a']} n/a, "
      f"{sums['deferred']} deferred; {len(errors)} problem(s)")
sys.exit(1 if errors else 0)
