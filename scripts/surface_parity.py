"""Surface parity: every in-scope public item of the Rust oracle has a Cyrius counterpart.

Port plan §3 (a per-milestone gate) and M7 (the release gate). Usage, from the repo root:
    python3 -I scripts/surface_parity.py .

For each ported module it lists rust-old/src/<module>.rs's public items before
`#[cfg(test)]` and checks that src/ defines the Cyrius name the porting rules give:
- free `pub fn f`                    -> `fn f(`
- `impl Value { pub fn m }`          -> `fn snake(Value)_m(`       (quark_mass_mev, nucleus_binding_energy)
  where Value is an enum or a type the port packs into a word (Nucleus)
- `impl Struct { pub fn m }`         -> `fn m(` inside `impl Tm<Struct> {`
- `impl Add for Struct`              -> `fn add(` inside `impl Tm<Struct> {` (the `+` operator)
- `pub const X` (constants.rs)       -> `const TM_X =`
- `TanmatraError::Variant`           -> `TM_ERR_<SNAKE>` in enum TmError
Anything else needs an entry in EXCEPTIONS with its reason. Every file under rust-old/src must be in
MODULES, DEFERRED_MODULES or NOT_APPLICABLE_MODULES (the full diff, port plan M7). Exits 1 on any item
missing or any file unaccounted.
"""
import re, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()

# The ported modules, with the Rust types whose methods become free functions ("values":
# enums and word-packed types) or impl methods ("structs").
MODULES = {
    "error": {"values": set(), "structs": set()},
    "constants": {"values": set(), "structs": set()},
    "particle": {"values": {"Quark", "Lepton", "Boson", "FundamentalForce"}, "structs": set()},
    "relativity": {"values": set(), "structs": {"FourMomentum"}},
    "nucleus": {"values": {"Nucleus"}, "structs": {"ShellLevel", "SuperallowedDecay"}},
    "decay": {"values": set(), "structs": set()},
    "timekeeping": {"values": {"FrequencyStandard"}, "structs": {"AtomicInstant"}},
    "atomic": {"values": {"OrbitalType"}, "structs": {"QuantumNumbers"}},
    "reaction": {"values": set(), "structs": {"NuclearReaction"}},
    "scattering": {"values": set(), "structs": set()},
    "optics": {"values": set(), "structs": {"SpectralLine"}},
}

# Whole Rust files the port does not carry, with why; their public items are counted, not checked.
DEFERRED_MODULES = {
    "bridge": "deferred: waits on its consumers (bijli, kimiya) being in Cyrius (port plan §1)",
    "integration/soorat": "deferred: waits on soorat (port plan §1)",
}
NOT_APPLICABLE_MODULES = {
    "lib": "n/a: module declarations and the prelude's re-exports; Cyrius has one flat namespace",
    "integration/mod": "n/a: the `soorat` module declaration",
}

# Rust item -> why it has no direct Cyrius name.
EXCEPTIONS = {
    "FundamentalForce::mediator": "a slice: fundamental_force_mediator_count + fundamental_force_mediator(f, k)",
    "shell_model_levels": "a slice: shell_model_level_count + shell_model_level(k)",
    "superallowed_ft_values": "a Vec: superallowed_ft_value_count + superallowed_ft_value(k)",
    "corrected_ft_value": "deprecated in 1.3.0, not carried over (port plan §1)",
    "known_isotopes": "a Vec: known_isotope_count + known_isotope(k), read from the baked table",
    "all_frequency_standards": "a slice: all_frequency_standard_count + all_frequency_standard(k)",
    "FrequencyStandard::quality_factor": "deprecated in 1.3.0, not carried over (port plan §1)",
    "const LEAP_SECOND_TABLE_VALID_UNTIL": "a (year, month, day) tuple: TM_LEAP_SECOND_TABLE_VALID_UNTIL_YEAR / _MONTH / _DAY",
    "thermal_neutron_cross_sections": "a Vec: thermal_neutron_cross_section_count + thermal_neutron_cross_section(k)",
    "u235_fission_yields": "a Vec: u235_fission_yield_count + u235_fission_yield(k)",
    "pu239_fission_yields": "a Vec: pu239_fission_yield_count + pu239_fission_yield(k)",
}


def snake(camel):
    return re.sub(r"(?<!^)([A-Z])", r"_\1", camel).lower()


def rust_items(module):
    src = (root / f"rust-old/src/{module}.rs").read_text()
    cut = src.find("#[cfg(test)]")
    src = src if cut < 0 else src[:cut]
    items = []
    impl = None
    depth = 0
    for line in src.splitlines():
        s = line.strip()
        m = re.match(r"impl(?:<[^>]*>)?\s+(?:(\w+(?:::\w+)*)\s+for\s+)?(\w+)\s*\{", s)
        if m and depth == 0:
            trait, ty = m.group(1), m.group(2)
            impl = (ty, trait)
            if trait and trait.split("::")[-1] == "Add":
                items.append((f"{ty}::add", ty, "add"))
        fm = re.match(r"pub (?:const )?fn (\w+)", s)
        if fm:
            if impl and depth >= 1:
                items.append((f"{impl[0]}::{fm.group(1)}", impl[0], fm.group(1)))
            elif depth == 0:
                items.append((fm.group(1), None, fm.group(1)))
        cm = re.match(r"pub const ([A-Z][A-Z0-9_]*)\s*:", s)
        if cm and depth == 0:
            items.append((f"const {cm.group(1)}", None, cm.group(1)))
        depth += line.count("{") - line.count("}")
        if depth == 0:
            impl = None
    if module == "error":
        body = src[src.index("pub enum TanmatraError"):]
        body = body[: body.index("\n}")]
        for v in re.findall(r"^\s{4}([A-Z]\w+)", body, re.M):
            items.append((f"TanmatraError::{v}", "TanmatraError", v))
    return items


cyr = "\n".join(p.read_text() for p in sorted((root / "src").glob("*.cyr")))


def impl_methods(struct):
    m = re.search(rf"impl Tm{struct} \{{(.*?)\n\}}", cyr, re.S)
    return set(re.findall(r"^\s+fn (\w+)\(", m.group(1), re.M)) if m else set()


missing = 0
for module, kinds in MODULES.items():
    found = 0
    items = rust_items(module)
    for label, ty, name in items:
        if label in EXCEPTIONS:
            found += 1
            continue
        if label.startswith("const "):
            ok = re.search(rf"^const TM_{name} =", cyr, re.M) is not None
        elif ty == "TanmatraError":
            ok = re.search(rf"\bTM_ERR_{snake(name).upper()} =", cyr) is not None
        elif ty in kinds["structs"]:
            ok = name in impl_methods(ty)
        elif ty in kinds["values"]:
            ok = re.search(rf"^fn {snake(ty)}_{name}\(", cyr, re.M) is not None
        else:
            ok = re.search(rf"^fn {name}\(", cyr, re.M) is not None
        if ok:
            found += 1
        else:
            missing += 1
            print(f"MISSING  {module}: {label}")
    print(f"{module:<12} {found}/{len(items)} Rust public items have a Cyrius counterpart")
# Every file under rust-old/src must be ported, deferred or not applicable.
known = set(MODULES) | set(DEFERRED_MODULES) | set(NOT_APPLICABLE_MODULES)
for path in sorted((root / "rust-old/src").rglob("*.rs")):
    mod = str(path.relative_to(root / "rust-old/src"))[:-3]
    if mod not in known:
        missing += 1
        print(f"UNKNOWN  rust-old/src/{mod}.rs is not ported, deferred or marked not applicable")
for mod, why in {**DEFERRED_MODULES, **NOT_APPLICABLE_MODULES}.items():
    print(f"{mod:<12} {len(rust_items(mod))} Rust public items, not ported: {why}")
print(f"surface parity: {missing} missing")
sys.exit(1 if missing else 0)
