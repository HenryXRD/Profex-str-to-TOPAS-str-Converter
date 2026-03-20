#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
ProfexStr_to_TOPAS_v2_stable
===============================================================================

Deterministic Profex/BGMN .str -> TOPAS .str converter 
    Author: Henry Spratt + copilot

If profex is installed, the str files can be converted to TOPAS format str files 

Design goals:
- Faithful preservation of crystallographic intent in profex.str files
- Strict enforcement of SPACEGRP.DAT symmetry definitions
    - Profex SPACEGRP.DAT used to insert missing fixed coordinates / check the structure 
- No invention of free positional or lattice parameters
- Deterministic, reproducible output suitable for batch conversion
- Explicit handling of ambiguity via rejection 

Fractional coordinates are treated modulo 1

Fixed-by-symmetry constants may be inserted when uniquely determined; free
parameters are never defaulted or guessed

The converter favours correctness and transparency over convenience. If a
structure is ambiguous or underspecified, it will fail clearly.

-------------------------------------------------------------------------------
EXPECTED SKIPS
-------------------------------------------------------------------------------
Some Profex models are out of scope and are expected to produce skipped sites, 
for reasons including but not limited to:
- Complex disorder models (e.g. turbostratic disorder / layer stacking)
- Wyckoff/coordinate inconsistencies (e.g. given coord values don't match letter, 
  or profex .str says setting 1 origin 1 but coords given in setting 2)

BGMN Folder
    Dataset size: 415 Profex .str files
    Expected files with skipped sites: 4
    Expected offenders:
        nontronite15a.str
        smectitedi2wfix1.str
        Arsenolite.str
        CRYOLITE.STR

Cement Folder
    Dataset size: 12 Profex .str files
    Expected files with skipped sites: 1
    Expected offender:
        CSH-0625.str

Ceramics Folder
    Dataset size: 101 Profex .str files
    Expected files with skipped sites: 3
    Expected offenders:
        CoFe2O4.str
        Cu2O.str
        ZnAl2O4.str

===============================================================================
===============================================================================
"""

import re
import argparse
import csv
#schorl fix 
import itertools
from pathlib import Path
from fractions import Fraction

DEFAULT_SPACEGRP_DAT = r"C:\Program Files\Profex5\BGMNwin\SPACEGRP.DAT"

FLOAT = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?"
FRAC  = r"[+-]?\d+/\d+"
# FLOAT must not match the leading part of a fraction (e.g. "1" in "1/4"),
# otherwise fixed Wyckoff coordinates can be misparsed.
FLOAT_NOFRAC = rf"(?:{FLOAT})(?!/)"
NUM_OR_FRAC = rf"(?:{FRAC}|{FLOAT_NOFRAC})"
NM_TO_A = 10.0
NM2_TO_A2 = 100.0

# v2_stable default bounds for unit cell when Profex omits min/max
DEFAULT_CELL_RELATIVE_TOL = 0.01  # ±1% for a, b, c (lengths)
DEFAULT_ANGLE_RELATIVE_TOL = 0.01  # ±1% for al, be, ga (angles)

print("RUNNING SCRIPT:", __file__)
#print("NUM_OR_FRAC:", NUM_OR_FRAC)

#from fractions import Fraction

def float_to_topas_rational(v: float, max_den: int = 96, tol: float = 1e-6) -> str:
    """
    Convert a float to a short rational string (e.g. 1/3) when it is very close
    to a simple fraction; otherwise return a trimmed decimal.
    
    max_den: largest denominator allowed (96 is a good default for crystallography:
             handles 1/2, 1/3, 1/4, 1/6, 1/8, 1/12, 1/16, 1/24, 1/32, etc.)
    tol: absolute tolerance for deciding if the float is 'equal' to the fraction
    """
    if v is None:
        return "0"

    # wrap into [0,1) like your internal matching does
    v = v % 1.0
    if v < 0:
        v += 1.0

    f = Fraction(v).limit_denominator(max_den)
    approx = float(f)

    if abs(approx - v) <= tol:
        # Prefer integers like "0" or "1" if they occur
        if f.denominator == 1:
            return str(f.numerator)
        return f"{f.numerator}/{f.denominator}"

    # fallback: print enough precision to avoid ugly repeating 6-decimal rounding
    s = f"{v:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"

def coord_value(expr):
    """
    Normalize a coordinate expression for TOPAS equation form:
        x = <value>;
    Ensures:
      - no leading '='
      - no trailing ';'
      - trims whitespace
    Accepts: "0.25", "0.25;", "=0.25", "= 0.25;", "1/4", "1/4;", "x+1/2;", etc.
    Returns: "0.25", "1/4", "x+1/2"
    """
    if expr is None:
        return "0"  # defensive fallback; valid sites should not hit this
    s = str(expr).strip()

    # Strip trailing semicolons (one or many)
    while s.endswith(";"):
        s = s[:-1].rstrip()

    # Strip one leading '=' if present
    if s.startswith("="):
        s = s[1:].lstrip()

    return s

def topas_xyz_clause(xo, yo, zo):
    """
    Returns: 'x = <...>; y = <...>; z = <...>;'
    with no double-semicolons regardless of xo/yo/zo content.
    """
    return f"x = {coord_value(xo)}; y = {coord_value(yo)}; z = {coord_value(zo)};"

def norm_el(sym: str) -> str:
    """
    Normalize element symbol to canonical chemical form:
    first letter uppercase, second (if any) lowercase.
    """
    sym = sym.strip()
    if len(sym) == 1:
        return sym.upper()
    return sym[0].upper() + sym[1].lower()
    
def use_numeric_spacegroup(lattice: str) -> bool:
    """
    Decide whether TOPAS space_group should use numeric SpacegroupNo
    instead of Hermann–Mauguin symbol, based on crystal system.
    """
    if not lattice:
        return False
    return lattice.lower() in {
        "cubic",
        # you may extend later:
        "tetragonal",
        "orthorhombic",
    }

# ---------------- Phase-name normalization (for TOPAS-safe identifiers) ----------------
_GREEK_MAP = {
    'α':'alpha','β':'beta','γ':'gamma','δ':'delta','ε':'epsilon','ζ':'zeta','η':'eta','θ':'theta',
    'ι':'iota','κ':'kappa','λ':'lambda','μ':'mu','ν':'nu','ξ':'xi','ο':'omicron','π':'pi','ρ':'rho',
    'σ':'sigma','ς':'sigma','τ':'tau','υ':'upsilon','φ':'phi','χ':'chi','ψ':'psi','ω':'omega',
    'Α':'alpha','Β':'beta','Γ':'gamma','Δ':'delta','Ε':'epsilon','Ζ':'zeta','Η':'eta','Θ':'theta',
    'Ι':'iota','Κ':'kappa','Λ':'lambda','Μ':'mu','Ν':'nu','Ξ':'xi','Ο':'omicron','Π':'pi','Ρ':'rho',
    'Σ':'sigma','Τ':'tau','Υ':'upsilon','Φ':'phi','Χ':'chi','Ψ':'psi','Ω':'omega',
}
_PRIME_MAP = {"′":"_prime","’":"_prime","'":"_prime"}

def normalize_phase_name(name: str) -> str:
    s = ''.join(_GREEK_MAP.get(ch, ch) for ch in name)
    for k, v in _PRIME_MAP.items():
        s = s.replace(k, v)
    s = s.replace('-', '_')
    s = re.sub(r'[^A-Za-z0-9_]', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    return s or "Phase"

# ---------------- SPACEGRP.DAT parsing ----------------
# If SPACEGRP.DAT includes origin choice in the header line, capture it.
RE_DAT_ORIGIN = re.compile(r"\bOriginChoice\s*=\s*(\d+)", re.I)

RE_BLOCK = re.compile(
    r"SpacegroupNo=(?P<sg>\d+)\s+HermannMauguin=(?P<hm>\S+)\s+Setting=(?P<set>\d+).*?"
    r"Lattice=(?P<lat>\S+)(?:\s+UniqueAxis=(?P<ua>[abc]))?",
    re.I
)
RE_WYCK = re.compile(r"Wyckoff=(?P<w>[a-z])\s+N=", re.I)
RE_COORD_LINE = re.compile(r"^\s*([xyz0-9/\-+.*\s]+)\s*$", re.I)

'''  #old version, was working. 
def parse_spacegrp_dat(p: Path):
    lines = p.read_text(errors="ignore").splitlines()
    hm_info = {}
    sg_default = {}

    # Additional reverse lookup maps for robust resolution
    hm_by_sg_setting = {}                # (sgno, setting) -> hm
    hm_by_sg_setting_origin = {}         # (sgno, setting, origin) -> hm

    i = 0
    while i < len(lines):
        m = RE_BLOCK.search(lines[i])
        if not m:
            i += 1
            continue
        sgno = int(m.group("sg"))
        hm = m.group("hm").strip()
        setting = int(m.group("set"))
        lattice = (m.group("lat") or "").strip()
        raw_header = lines[i].strip()
        mo = RE_DAT_ORIGIN.search(raw_header)
        origin_choice = int(mo.group(1)) if mo else None
        ua = (m.group("ua") or None)
        ua = ua.lower() if ua else None

        hm_info.setdefault(hm, {
            "sgno": sgno,
            "setting": setting,
            "lattice": lattice,
            "unique_axis": ua,
            "wyckoff_proto": {},
            "wyckoff_ops": {},   # all equivalent-position triplets for each Wyckoff letter
            "origin_choice": origin_choice,
            "raw_header": raw_header
        })
        if setting == 1 and sgno not in sg_default:
            sg_default[sgno] = hm

        # Reverse indices for ALL settings (first seen wins)
        hm_by_sg_setting.setdefault((sgno, setting), hm)
        if origin_choice is not None:
            hm_by_sg_setting_origin.setdefault((sgno, setting, origin_choice), hm)

        i += 1
        while i < len(lines) and not RE_BLOCK.search(lines[i]):
            mw = RE_WYCK.search(lines[i])
            if not mw:
                i += 1
                continue
            w = mw.group("w").lower()
            i += 1

            ops = []
            proto = None

            while i < len(lines):
                s = lines[i].strip()
                if not s:
                    i += 1
                    continue
                if RE_WYCK.search(s) or RE_BLOCK.search(s):
                    break
                if RE_COORD_LINE.match(s):
                    toks = s.split()
                    if len(toks) >= 3:
                        trip = [t.lower() for t in toks[:3]]
                        ops.append(trip)
                        if proto is None:
                            proto = trip
                i += 1

            if proto:
                hm_info[hm]["wyckoff_proto"][w] = proto
            if ops:
                hm_info[hm]["wyckoff_ops"][w] = ops

    return hm_info, sg_default, hm_by_sg_setting, hm_by_sg_setting_origin
'''

def parse_spacegrp_dat(p: Path):
    lines = p.read_text(errors="ignore").splitlines()

    hm_info = {}        # hm -> [entry, entry, ...]
    sg_default = {}     # sgno -> hm (default = first Setting=1 seen)

    hm_by_sg_setting = {}               # (sgno, setting) -> hm
    hm_by_sg_setting_origin = {}        # (sgno, setting, origin) -> hm

    info_by_sg_setting = {}             # (sgno, setting) -> entry
    info_by_sg_setting_origin = {}      # (sgno, setting, origin) -> entry

    i = 0
    while i < len(lines):
        m = RE_BLOCK.search(lines[i])
        if not m:
            i += 1
            continue

        sgno = int(m.group("sg"))
        hm = m.group("hm").strip()
        setting = int(m.group("set"))
        lattice = (m.group("lat") or "").strip()
        raw_header = lines[i].strip()

        mo = RE_DAT_ORIGIN.search(raw_header)
        origin_choice = int(mo.group(1)) if mo else None

        ua = (m.group("ua") or None)
        ua = ua.lower() if ua else None

        entry = {
            "sgno": sgno,
            "hm": hm,
            "setting": setting,
            "origin_choice": origin_choice,
            "lattice": lattice,
            "unique_axis": ua,
            "raw_header": raw_header,
            "wyckoff_proto": {},  # w -> [expr,expr,expr] (first triplet)
            "wyckoff_ops": {},    # w -> [[...],[...],...] (all triplets)
        }

        # establish default HM for sgno from the first Setting=1 encountered
        if setting == 1 and sgno not in sg_default:
            sg_default[sgno] = hm

        # reverse indices (first seen wins)
        hm_by_sg_setting.setdefault((sgno, setting), hm)
        if origin_choice is not None:
            hm_by_sg_setting_origin.setdefault((sgno, setting, origin_choice), hm)

        i += 1
        while i < len(lines) and not RE_BLOCK.search(lines[i]):
            mw = RE_WYCK.search(lines[i])
            if not mw:
                i += 1
                continue

            w = mw.group("w").lower()
            i += 1

            ops = []
            proto = None

            while i < len(lines):
                s = lines[i].strip()
                if not s:
                    i += 1
                    continue
                if RE_WYCK.search(s) or RE_BLOCK.search(s):
                    break
                if RE_COORD_LINE.match(s):
                    toks = s.split()
                    if len(toks) >= 3:
                        trip = [t.lower() for t in toks[:3]]
                        ops.append(trip)
                        if proto is None:
                            proto = trip
                i += 1

            if proto:
                entry["wyckoff_proto"][w] = proto
            if ops:
                entry["wyckoff_ops"][w] = ops

        # store entry keyed by sg/setting(/origin)
        info_by_sg_setting[(sgno, setting)] = entry
        if origin_choice is not None:
            info_by_sg_setting_origin[(sgno, setting, origin_choice)] = entry

        hm_info.setdefault(hm, []).append(entry)

    return (
        hm_info,
        sg_default,
        hm_by_sg_setting,
        hm_by_sg_setting_origin,
        info_by_sg_setting,
        info_by_sg_setting_origin,
    )

def hm_from_sg_setting(hm_info, sgno, setting_no, unique_axis_profex=None):
    for hm_key, info in hm_info.items():
        if (
            info.get("sgno") == sgno
            and info.get("setting") == setting_no
            and ua_matches(info.get("unique_axis"), unique_axis_profex)
        ):
            return hm_key
    return None
    
def ua_matches(dat_ua, profex_ua):
    if profex_ua is None:
        return True
    if dat_ua is None:
        return True
    return str(dat_ua).lower() == str(profex_ua).lower()

# ---------------- Profex parsing ----------------
RE_PHASE   = re.compile(r"^\s*PHASE\s*=\s*([^\s/]+)", re.I)
RE_SGNO_STRICT = re.compile(r"\bSpacegroupNo\s*=\s*(\d+)", re.I)
RE_SGNO_FUZZY  = re.compile(r"\bSpaceg\w{3,8}No\s*=\s*(\d+)", re.I)
RE_SETTING = re.compile(r"\bSetting\s*=\s*(\d+)", re.I)
RE_UA_PROFEX = re.compile(r"\bUniqueAxis\s*=\s*([abc])", re.I)
RE_ORIGIN  = re.compile(r"\bOriginChoice\s*=\s*(\d+)", re.I)

# IMPORTANT: allow / and - etc in HM by stopping only at whitespace
RE_HM      = re.compile(r"\bHermannMauguin\s*=\s*([^\s]+)", re.I)
RE_FORMULA = re.compile(r"\bFormula\s*=\s*(.+)$", re.I)

RE_PARAM = re.compile(r"\bPARAM\s*=\s*([A-Z]+)\s*=\s*([^\s/]+)", re.I)

# PARAM parsing inside atom lines (supports lowercase names like p, fe, xMg, etc.)
RE_PARAM_IN_LINE = re.compile(
    r"\bPARAM\s*=\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^\s/]+)",
    re.I
)

# Bare (non-PARAM=) cell parameters, accepted ONLY in the header window
RE_CELL_BARE = re.compile(
    r"\b(A|B|C|ALPHA|BETA|GAMMA)\s*=\s*([^\s/]+)",
    re.I
)

RE_ATOM_SINGLE = re.compile(
    rf"\bE=(?P<el>[A-Za-z][A-Za-z]?)(?P<charge>[+-]\d+)?"
    rf"(?:\((?P<occ>{FLOAT})\))?"
    rf"\s+.*?\bWyckoff=(?P<w>[a-z])"
    rf"(?:\s+x=(?P<x>{NUM_OR_FRAC}))?"
    rf"(?:\s+y=(?P<y>{NUM_OR_FRAC}))?"
    rf"(?:\s+z=(?P<z>{NUM_OR_FRAC}))?"
    rf"(?:\s+TDS=(?P<tds>{NUM_OR_FRAC}))?",
    re.I
)

RE_ATOM_MIXED = re.compile(
    rf"\bE=\((?P<inner>.+?)\)\s+.*?\bWyckoff=(?P<w>[a-z])"
    rf"(?:\s+x=(?P<x>{NUM_OR_FRAC}))?"
    rf"(?:\s+y=(?P<y>{NUM_OR_FRAC}))?"
    rf"(?:\s+z=(?P<z>{NUM_OR_FRAC}))?"
    rf"(?:\s+TDS=(?P<tds>{NUM_OR_FRAC}))?",
    re.I
)

RE_MIX_COMP = re.compile(
    rf"(?P<el>[A-Za-z][A-Za-z]?)(?P<charge>[+-]\d+)?"
    rf"(?:\((?P<occ>[^)]+)\))?"
    rf"(?:\s*\*TDS=(?P<tds>{NUM_OR_FRAC}))?",
    re.I
)

def extract_xyz_tds_from_line(raw: str):
    """
    Extract x,y,z and TDS from a Profex atom line regardless of token order.
    Returns:
      given: {"x": (float|None, text|None), "y": (...), "z": (...)}
      site_tds: float|None
    """
    given = {"x": (None, None), "y": (None, None), "z": (None, None)}
    site_tds = None

    # Find x=..., y=..., z=..., TDS=... in ANY order
    for key, val in re.findall(
        r"\b(x|y|z|TDS)\s*=\s*(" + NUM_OR_FRAC + r")",
        raw,
        flags=re.I
    ):
        k = key.lower()

        if k == "tds":
            if "/" in val:
                n, d = val.split("/", 1)
                site_tds = float(Fraction(int(n), int(d)))
            else:
                site_tds = float(val)
            continue

        # x / y / z
        if "/" in val:
            n, d = val.split("/", 1)
            fv = float(Fraction(int(n), int(d)))
            txt = val
        else:
            fv = float(val)
            txt = float_to_topas_rational(fv)

        given[k] = (fv, txt)

    return given, site_tds

def parse_bgmn_mid_min_max(tok):
    if any(c.isalpha() for c in tok):
        raise ValueError(f"Non-numeric PARAM value: {tok}")
    if "_" in tok and "^" in tok:
        mid, rest = tok.split("_", 1)
        lo, hi = rest.split("^", 1)
        return float(mid), float(lo), float(hi)
    return float(tok), None, None

def cell_key_map(k):
    return {"A":"a","B":"b","C":"c","ALPHA":"al","BETA":"be","GAMMA":"ga"}[k]

def sanitize_param_name(name: str) -> str:
    # TOPAS-safe parameter name token
    s = re.sub(r'[^A-Za-z0-9_]', '_', name.strip())
    s = re.sub(r'_+', '_', s).strip('_')
    return s or "p"

def parse_occ_token(tok: str):
    """
    Returns (occ_val, occ_expr)
    - occ_val: float when numeric
    - occ_expr: string when symbolic (e.g. 'p', 'fe')
    """
    if tok is None:
        return (None, None)
    s = tok.strip()
    # numeric?
    try:
        if "/" in s:
            n, d = s.split("/", 1)
            return (float(Fraction(int(n), int(d))), None)
        return (float(s), None)
    except ValueError:
        return (None, s)

def rewrite_expr(expr: str, name_map: dict) -> str:
    """
    Replace whole-word param identifiers in expr using name_map.
    """
    out = expr
    for k, v in name_map.items():
        out = re.sub(rf"\b{re.escape(k)}\b", v, out)
    return out

def infer_rhombohedral_axes_from_cell(cell_ranges: dict) -> bool:
    """
    Heuristic: Profex files written in rhombohedral (primitive) axes often specify:
      - a and alpha
      - NOT c and NOT gamma
    i.e. (a, α) only, implying a=b=c and α=β=γ.
    """
    return ("a" in cell_ranges and "al" in cell_ranges and "c" not in cell_ranges and "ga" not in cell_ranges)

# ---------------- Wyckoff coordinate handling ----------------
'''
def token_to_float(tok):
    tok = tok.lower().strip()
    if any(c.isalpha() for c in tok):
        raise ValueError(f"Internal error: token_to_float called on non-numeric token '{tok}'")
    if "/" in tok:
        n, d = tok.split("/", 1)
        return float(Fraction(int(n), int(d)))
    return float(tok)
'''    
_VAR_RE = re.compile(
    r"^(?P<sign>-)?(?:(?P<coef>\d+)\*)?(?P<var>[xyz])(?P<off>([+-].+))?$",
    re.I
)

_CONST_DEC_RE  = re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$')
_CONST_FRAC_RE = re.compile(r'^[+-]?\d+/\d+$')

def _wrap01(v: float) -> float:
    v = v % 1.0
    return v + 1.0 if v < 0 else v

#schorl fix
def _solve_mod_linear(gv: float, coeff: float, off: float, max_k: int = 12):
    """
    Solve coeff*x + off = gv (mod 1) for x when coeff is a small integer.
    Returns candidate x values in [0,1).
    """
    if coeff == 0:
        return []
    k = int(round(coeff))
    # only support small integer coefficients deterministically
    if abs(k) < 1 or abs(k) > max_k or abs(coeff - k) > 1e-9:
        return []
    base = _wrap01((gv - off) / k)
    n = abs(k)
    return [_wrap01(base + m / n) for m in range(n)]

def centring_vectors_from_hm(hm: str):
    """
    Return lattice centring translation vectors for the space group symbol.
    Always includes (0,0,0). Values are fractional coordinates modulo 1.
    """
    if not hm:
        return [(0.0, 0.0, 0.0)]
    c = hm.strip()[0].upper()

    # Primitive
    if c == "P":
        return [(0.0, 0.0, 0.0)]

    # Base-centred
    if c == "A":
        return [(0.0, 0.0, 0.0), (0.0, 0.5, 0.5)]
    if c == "B":
        return [(0.0, 0.0, 0.0), (0.5, 0.0, 0.5)]
    if c == "C":
        return [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]

    # Body-centred
    if c == "I":
        return [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]

    # Face-centred
    if c == "F":
        return [
            (0.0, 0.0, 0.0),
            (0.0, 0.5, 0.5),
            (0.5, 0.0, 0.5),
            (0.5, 0.5, 0.0),
        ]

    # Rhombohedral (hexagonal axes setting; common in HM symbols starting with R)
    # The conventional R-centring vectors in hex axes are (2/3,1/3,1/3) and (1/3,2/3,2/3) plus (0,0,0).
    if c == "R":
        return [
            (0.0, 0.0, 0.0),
            (2.0/3.0, 1.0/3.0, 1.0/3.0),
            (1.0/3.0, 2.0/3.0, 2.0/3.0),
        ]

    # Unknown / default to primitive
    return [(0.0, 0.0, 0.0)]


def shift_given_xyz(given_xyz, t):
    """Shift provided numeric coordinates by -t (mod 1). None stays None."""
    tx, ty, tz = t
    out = {}
    for ax, tv in zip(("x", "y", "z"), (tx, ty, tz)):
        gv = given_xyz[ax]
        out[ax] = None if gv is None else _wrap01(gv - tv)
    return out


def unshift_out_xyz(out_xyz, t):
    """Shift inferred output coordinates by +t (mod 1). None stays None."""
    tx, ty, tz = t
    out = {}
    for ax, tv in zip(("x", "y", "z"), (tx, ty, tz)):
        v = out_xyz[ax]
        out[ax] = None if v is None else _wrap01(v + tv)
    return out

def _parse_const(s: str) -> float:
    s = s.strip()
    if any(c.isalpha() for c in s):
        raise ValueError(f"Non-numeric constant in Wyckoff expression: '{s}'")
    if "/" in s:
        n, d = s.split("/", 1)
        return float(Fraction(int(n), int(d)))
    return float(s)

def _parse_offset(off: str) -> float:
    if not off:
        return 0.0
    total = 0.0
    for part in re.findall(r"[+-][^+-]+", off.replace(" ", "")):
        sign = 1.0 if part[0] == "+" else -1.0
        total += sign * _parse_const(part[1:])
    return total

def _parse_expr(expr: str):
    expr = expr.strip().lower()

    if expr in ("x", "y", "z"):
        return (expr, 1.0, 0.0)

    # constant (only if the ENTIRE token is numeric or a simple fraction)
    if _CONST_DEC_RE.match(expr) or _CONST_FRAC_RE.match(expr):
        return (None, 0.0, _wrap01(_parse_const(expr)))

    m = _VAR_RE.match(expr)
    if not m:
        raise ValueError(f"Unsupported Wyckoff expression: {expr}")

    var = m.group("var").lower()
    sign = -1.0 if m.group("sign") else 1.0
    coef = int(m.group("coef")) if m.group("coef") else 1
    coeff = sign * coef

    off = _wrap01(_parse_offset(m.group("off") or ""))
    return (var, coeff, off)

def _eval_expr(parsed, vals):
    var, coeff, off = parsed
    if var is None:
        return off
    return _wrap01(coeff * vals[var] + off)

def _diff_mod1(a, b):
    d = abs(_wrap01(a - b))
    return min(d, 1.0 - d)

#old working version before schorl fix 
'''
def try_match_triplet(triplet, given_xyz, tol=5e-4):
    """
    triplet: [e1,e2,e3] strings from SPACEGRP.DAT
    given_xyz: {'x': float|None, 'y': float|None, 'z': float|None}
    Returns (matched: bool, out_xyz: dict|None, inferred_vars: dict)
      out_xyz contains computed x/y/z where determinable, else None
      inferred_vars are solved free parameters (x/y/z) used internally
    """
    try:
        parsed = [_parse_expr(e) for e in triplet]
    except ValueError:
        # Unsupported Wyckoff form in SPACEGRP.DAT (e.g. x-y, -x+y, etc.)
        # Treat as "no match" rather than aborting the whole conversion.
        return (False, None, {})
    vals = {}

    # Solve parameters from provided coords
    for axis, p in zip(("x","y","z"), parsed):
        gv = given_xyz[axis]
        var, coeff, off = p
        if gv is None:
            continue

        if var is None:
            if _diff_mod1(gv, off) > tol:
                return (False, None, {})
            continue

        # If coefficient is not ±1 and the underlying variable is already provided,
        # validate directly rather than trying to invert modulo 1.
        if abs(coeff) != 1 and given_xyz.get(var) is not None:
            test = _wrap01(coeff * given_xyz[var] + off)
            if _diff_mod1(gv, test) > tol:
                return (False, None, {})
            candidate = _wrap01(given_xyz[var])
        else:
            candidate = _wrap01((gv - off) / coeff)

        if var in vals and _diff_mod1(vals[var], candidate) > tol:
            return (False, None, {})
        vals[var] = candidate

    # Evaluate coordinates
    out = {}
    for axis, p in zip(("x","y","z"), parsed):
        var, sign, off = p
        out[axis] = off if var is None else (_eval_expr(p, vals) if var in vals else None)

    # Verify provided coords match evaluated coords
    for axis in ("x","y","z"):
        gv = given_xyz[axis]
        if gv is None:
            continue
        if out[axis] is None:
            return (False, None, {})
        if _diff_mod1(gv, out[axis]) > tol:
            return (False, None, {})

    return (True, out, vals)
'''

#schorl fix version: -2x -x coord forms etc
def try_match_triplet(triplet, given_xyz, tol=5e-4):
    """
    triplet: [e1,e2,e3] strings from SPACEGRP.DAT
    given_xyz: {'x': float|None, 'y': float|None, 'z': float|None}
    Returns (matched: bool, out_xyz: dict|None, inferred_vars: dict)
      out_xyz contains computed x/y/z where determinable, else None
      inferred_vars are solved free parameters (x/y/z) used internally
    """
    try:
        parsed = [_parse_expr(e) for e in triplet]  # each: (var|None, coeff, off)
    except ValueError:
        return (False, None, {})

    # Build candidate sets for internal variables (x/y/z) implied by provided coords
    cand_map = {}  # var -> list of candidates in [0,1)

    for axis, p in zip(("x", "y", "z"), parsed):
        gv = given_xyz.get(axis)
        var, coeff, off = p

        if gv is None:
            continue

        if var is None:
            # axis is constant; must match
            if _diff_mod1(gv, off) > tol:
                return (False, None, {})
            continue

        # Solve coeff*var + off = gv (mod 1)
        if abs(coeff) == 1:
            candidates = [_wrap01((gv - off) / coeff)]
        else:
            candidates = _solve_mod_linear(gv, coeff, off)
            if not candidates:
                return (False, None, {})  # unsupported coefficient form

        # Intersect with any existing candidate set for this variable
        if var not in cand_map:
            cand_map[var] = candidates
        else:
            new = []
            for a in cand_map[var]:
                if any(_diff_mod1(a, b) <= tol for b in candidates):
                    new.append(a)
            # de-dup within tol
            dedup = []
            for v in new:
                if not any(_diff_mod1(v, u) <= tol for u in dedup):
                    dedup.append(v)
            cand_map[var] = dedup
            if not cand_map[var]:
                return (False, None, {})

    # Evaluate helper: returns value if evaluable, else None
    def eval_axis(p, vals):
        var, coeff, off = p
        if var is None:
            return off
        if var not in vals:
            return None
        return _wrap01(coeff * vals[var] + off)

    # If no variables were constrained, we can still succeed only if provided coords
    # were constants (handled above). Otherwise we cannot verify.
    vars_in_play = sorted(cand_map.keys())
    if not vars_in_play:
        out = {}
        for axis, p in zip(("x", "y", "z"), parsed):
            out[axis] = eval_axis(p, {})
        # verify provided coords match evaluable axes
        for axis in ("x", "y", "z"):
            gv = given_xyz.get(axis)
            if gv is None:
                continue
            if out[axis] is None or _diff_mod1(gv, out[axis]) > tol:
                return (False, None, {})
        return (True, out, {})

    # Try all combinations (tiny in practice: coeff=2 => 2 candidates)
    candidate_lists = [cand_map[v] for v in vars_in_play]
    for combo in itertools.product(*candidate_lists):
        vals = dict(zip(vars_in_play, combo))

        # Verify all provided coords match
        ok = True
        for axis, p in zip(("x", "y", "z"), parsed):
            gv = given_xyz.get(axis)
            if gv is None:
                continue
            ev = eval_axis(p, vals)
            if ev is None or _diff_mod1(gv, ev) > tol:
                ok = False
                break
        if not ok:
            continue

        # Build out_xyz (constants and solved axes evaluated; unsolved axes -> None)
        out = {}
        for axis, p in zip(("x", "y", "z"), parsed):
            out[axis] = eval_axis(p, vals)
        return (True, out, vals)

    return (False, None, {})

'''
def format_coord_from_proto(proto_token, given_val, given_txt, axis_label, tol=5e-4):
    """
    Returns (out_text, mismatch_flag, mismatch_detail)
    - If proto_token is x/y/z: free variable. If missing in Profex, assume 0 and flag mismatch.
    - Else proto_token is constant/fraction: enforce it; if Profex differs, flag mismatch.
    """
    p = proto_token.lower().strip()

    # treat variable and affine-variable tokens as "free-form"
    if p in ("x", "y", "z") or p.startswith(("x+", "x-", "y+", "y-", "z+", "z-", "-x", "-y", "-z")):
        if given_txt is None:
            # follow your current policy here:
            # either return None (strict) OR keep old behavior
            return None, True, f"missing free {axis_label}"
        return given_txt, False, None

    # constant/fraction
    out_txt = p
    if given_val is None or given_txt is None:
        return out_txt, False, None

    if abs(given_val - token_to_float(p)) > tol:
        return out_txt, True, f"{axis_label} {given_txt}->{out_txt}"
    return out_txt, False, None

def resolve_coords(proto3, given_xyz):
    (xv, xt) = given_xyz["x"]
    (yv, yt) = given_xyz["y"]
    (zv, zt) = given_xyz["z"]

    xo, mx, dx = format_coord_from_proto(proto3[0], xv, xt, "x")
    yo, my, dy = format_coord_from_proto(proto3[1], yv, yt, "y")
    zo, mz, dz = format_coord_from_proto(proto3[2], zv, zt, "z")

    mismatch = mx or my or mz
    details = [d for d in (dx, dy, dz) if d]
    return xo, yo, zo, mismatch, details
'''

def pick_setting_from_profex_hm(info_by_sg_setting, sgno: int, hm_profex: str):
    """
    If Profex provides a Hermann–Mauguin symbol but not Setting,
    try to select the SPACEGRP.DAT entry for this SpacegroupNo whose HM matches.
    Returns: setting (int) or None if no unique match.
    """
    if not hm_profex:
        return None
    hm_profex = hm_profex.strip()

    matches = []
    for (g, s), entry in info_by_sg_setting.items():
        if g != sgno:
            continue
        if str(entry.get("hm", "")).strip() == hm_profex:
            matches.append(s)

    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    return None

# ---------------- Cell constraints (lattice-based) ----------------
def apply_lattice_constraints(out, phase_topas, cell_ranges, lattice, unique_axis, rhombo_axes=False):
    lat = (lattice or "").strip().lower()
    ua = (unique_axis or "b").strip().lower()

    def pname(axis):  # parameter names follow same scheme
        return f"{axis}_{phase_topas}"

    def have(axis):
        return axis in cell_ranges

    # Cubic: a=b=c, angles 90
    if lat == "cubic":
        if have("a"):
            out.append(f"  b = {pname('a')};")
            out.append(f"  c = {pname('a')};")
        out.append("  al 90")
        out.append("  be 90")
        out.append("  ga 90")
        return

    # Tetragonal: a=b != c, angles 90
    if lat == "tetragonal":
        if have("a"):
            out.append(f"  b = {pname('a')};")
        out.append("  al 90")
        out.append("  be 90")
        out.append("  ga 90")
        return

    # Hexagonal / Trigonal:
    # - default: hexagonal axes constraints (a=b, α=β=90, γ=120)
    # - rhombo_axes=True: primitive rhombohedral axes constraints (a=b=c, α=β=γ)
    if lat in ("hexagonal", "trigonal"):
        if rhombo_axes:
            # a=b=c
            if have("a"):
                out.append(f" b = {pname('a')};")
                out.append(f" c = {pname('a')};")
            # α=β=γ (link β and γ to α if α is present)
            if have("al"):
                out.append(f" be = {pname('al')};")
                out.append(f" ga = {pname('al')};")
            return
        else:
            # hex axes
            if have("a"):
                out.append(f" b = {pname('a')};")
                out.append(" al 90")
                out.append(" be 90")
                out.append(" ga 120")
            return

    # Orthorhombic: all angles 90
    if lat == "orthorhombic":
        out.append("  al 90")
        out.append("  be 90")
        out.append("  ga 90")
        return

    # Monoclinic: two angles 90 depending on unique axis
    if lat == "monoclinic":
        if ua == "b":
            out.append("  al 90")
            out.append("  ga 90")
        elif ua == "a":
            out.append("  be 90")
            out.append("  ga 90")
        elif ua == "c":
            out.append("  al 90")
            out.append("  be 90")
        return

    # Triclinic: no automatic fixing
    return

# ---------------- Main conversion ----------------
def convert_one_str(path: Path, hm_info, sg_default, hm_by_sg_setting, hm_by_sg_setting_origin, info_by_sg_setting, info_by_sg_setting_origin, phase_counts):
    lines = path.read_text(errors="ignore").splitlines()

    phase_raw = None
    sgno = None
    setting_profex = None
    hm_profex = None
    formula_line = None
    unique_axis_profex = None
    origin_profex = None
    profex_sym_line = None  # full line containing SpacegroupNo=... (for provenance)
    sg_fuzzy_used = False


    cell_ranges = {}  # only what appears in Profex: ax -> (mid, lo, hi)
    atoms = []        # expanded atoms
    site_key_counter = 0

    in_cell_header_window = False
    cell_window_closed = False

    for ln in lines:
        if phase_raw is None:
            m = RE_PHASE.search(ln)
            if m:
                phase_raw = m.group(1).strip()

        # Open window when we encounter the spacegroup line
        if (not in_cell_header_window) and RE_SGNO_STRICT.search(ln):
            in_cell_header_window = True

        # Close window when we encounter the first RP=... line
        if in_cell_header_window and (not cell_window_closed) and re.search(r"\bRP\s*=", ln, flags=re.I):
            cell_window_closed = True
            in_cell_header_window = False
        # --- END STEP 3 ---

        m_sg = RE_SGNO_STRICT.search(ln)
        if m_sg:
            if sgno is None:
                sgno = int(m_sg.group(1))
            if profex_sym_line is None:
                profex_sym_line = ln.strip()
        else:
            # Fuzzy fallback ONLY if strict match failed
            m_fz = RE_SGNO_FUZZY.search(ln)
            if m_fz and sgno is None:
                sgno = int(m_fz.group(1))
                if profex_sym_line is None:
                    profex_sym_line = ln.strip()
                sg_fuzzy_used = True

        if setting_profex is None:
            m = RE_SETTING.search(ln)
            if m:
                setting_profex = int(m.group(1))
                
        if unique_axis_profex is None:
            m = RE_UA_PROFEX.search(ln)
            if m:
                unique_axis_profex = m.group(1).lower()
                
        if origin_profex is None:
            m = RE_ORIGIN.search(ln)
            if m:
                origin_profex = int(m.group(1))

        if hm_profex is None:
            m = RE_HM.search(ln)
            if m:
                hm_profex = m.group(1).strip()

        if formula_line is None:
            m = RE_FORMULA.search(ln)
            if m:
                formula_line = ln.strip()

        for mp in RE_PARAM.finditer(ln):
            k = mp.group(1).upper()
            v = mp.group(2)
            if k in ("A", "B", "C", "ALPHA", "BETA", "GAMMA"):
                try:
                    cell_ranges[cell_key_map(k)] = parse_bgmn_mid_min_max(v)
                except ValueError:
                    pass

        # Also accept bare A/B/C/ALPHA/BETA/GAMMA in the header window ONLY
        if in_cell_header_window and not cell_window_closed:
            for mp in RE_CELL_BARE.finditer(ln):
                k = mp.group(1).upper()
                v = mp.group(2)

                # Skip if PARAM= already handled this key (spacing-safe)
                if re.search(rf"\bPARAM\s*=\s*{k}\s*=", ln, flags=re.I):
                    continue

                cell_key = cell_key_map(k)

                # Do not overwrite a value already parsed via PARAM=
                if cell_key in cell_ranges:
                    continue

                if k in ("A", "B", "C", "ALPHA", "BETA", "GAMMA"):
                    try:
                        cell_ranges[cell_key] = parse_bgmn_mid_min_max(v)
                    except ValueError:
                        pass

        # Mixed occupancy first
        mm = RE_ATOM_MIXED.search(ln)
        if mm:
            raw = ln.strip()
            w = mm.group("w").lower()
            given, site_tds = extract_xyz_tds_from_line(raw)

            # site-local PARAM= definitions on this same atom line
            local_params = {}
            for pm in RE_PARAM_IN_LINE.finditer(raw):
                pname = pm.group(1)
                try:
                    local_params[pname] = parse_bgmn_mid_min_max(pm.group(2))
                except ValueError:
                    pass

            site_key_counter += 1
            site_key = site_key_counter

            inner = mm.group("inner")
            parts = [p.strip() for p in inner.split(",")]

            comps = []
            for part in parts:
                mc = RE_MIX_COMP.search(part)
                if not mc:
                    continue
                el = norm_el(mc.group("el"))
                if len(el) not in (1, 2):
                    raise ValueError(f"Invalid element symbol after normalization: {el}")

                occ_val, occ_expr = parse_occ_token(mc.group("occ"))
                comps.append({
                    "el": el,
                    "occ": occ_val,
                    "occ_expr": occ_expr,
                    "tds": float(mc.group("tds")) if mc.group("tds") else None
                })

            # infer remainder occupancy for 2-component sites:
            # if one comp has a symbolic occ (e.g. p) and the other has none -> set other to (1-p)
            exprs = [c for c in comps if c["occ_expr"] is not None]
            missing = [c for c in comps if c["occ_expr"] is None and c["occ"] is None]
            if len(comps) == 2 and len(exprs) == 1 and len(missing) == 1:
                p = exprs[0]["occ_expr"]
                missing[0]["occ_expr"] = f"(1-{p})"

            # emit expanded atoms
            for c in comps:
                occ = 1.0 if (c["occ"] is None and c["occ_expr"] is None) else c["occ"]
                occ_expr = c["occ_expr"]

                tds = c["tds"] if c["tds"] is not None else site_tds
                beq = (tds * NM2_TO_A2) if tds is not None else None

                atoms.append({
                    "el": c["el"],
                    "occ": occ if occ_expr is None else None,
                    "occ_expr": occ_expr,
                    "wyck": w,
                    "given": given,
                    "beq": beq,
                    "raw": raw,
                    "mixed": True,
                    "site_key": site_key,
                    "local_params": local_params,
                })
            continue

        ms = RE_ATOM_SINGLE.search(ln)
        if ms:
            raw = ln.strip()
            w = ms.group("w").lower()
            given, site_tds = extract_xyz_tds_from_line(raw)

            # site-local PARAM= definitions on this same atom line
            local_params = {}
            for pm in RE_PARAM_IN_LINE.finditer(raw):
                pname = pm.group(1)
                try:
                    local_params[pname] = parse_bgmn_mid_min_max(pm.group(2))
                except ValueError:
                    pass

            site_key_counter += 1
            site_key = site_key_counter

            el = norm_el(ms.group("el"))
            if len(el) not in (1, 2):
                raise ValueError(f"Invalid element symbol after normalization: {el}")

            occ = float(ms.group("occ")) if ms.group("occ") else 1.0
            tds = site_tds
            beq = (tds * NM2_TO_A2) if tds is not None else None

            atoms.append({
                "el": el,
                "occ": occ,
                "occ_expr": None,
                "wyck": w,
                "given": given,
                "beq": beq,
                "raw": raw,
                "mixed": False,
                "site_key": site_key,
                "local_params": local_params,
            })

    if phase_raw is None:
        phase_raw = path.stem
    if phase_counts.get(phase_raw, 0) > 1:
        # Duplicate PHASE= → use filename
        phase_topas = normalize_phase_name(path.stem)
    else:
        # Unique PHASE= → use it directly
        phase_topas = normalize_phase_name(phase_raw)

    # ---------------- Space-group resolution (FIX 1 / Option 1: select entry first) ----------------

    if sgno is None:
        raise ValueError(f"{path.name}: Profex file does not specify SpacegroupNo")

    # Determine which Setting to use:
    # 1) If Profex provides Setting=, use it.
    # 2) Else prefer Setting=1 if SPACEGRP.DAT has it.
    # 3) Else use the first available setting for this SpacegroupNo.
    if setting_profex is not None:
        setting_use = setting_profex
        sg_setting_note = None
    else:
        # 1) If Profex HM uniquely identifies a SPACEGRP.DAT entry for this sgno, use it.
        inferred_from_hm = pick_setting_from_profex_hm(info_by_sg_setting, sgno, hm_profex)

        if inferred_from_hm is not None:
            setting_use = inferred_from_hm
            sg_setting_note = None

        else:
            # 2) Otherwise, fall back to your existing rhombohedral-axes inference (only relevant when trigonal/R)
            want_rhombo_axes = infer_rhombohedral_axes_from_cell(cell_ranges)

            if want_rhombo_axes and (sgno, 2) in info_by_sg_setting:
                setting_use = 2
                sg_setting_note = None

            elif (sgno, 1) in info_by_sg_setting:
                setting_use = 1
                sg_setting_note = " /* NOTE: Profex did not specify Setting; using SPACEGRP.DAT Setting=1 default. */"

            else:
                avail = sorted({s for (g, s) in info_by_sg_setting.keys() if g == sgno})
                if not avail:
                    raise ValueError(f"{path.name}: No SPACEGRP.DAT entries found for SpacegroupNo={sgno}")
                setting_use = avail[0]
                sg_setting_note = None

    # Determine OriginChoice to use (only if Profex provides it; otherwise leave None)
    origin_use = origin_profex

    # Select the specific SPACEGRP.DAT entry (OriginChoice-aware if possible)
    info = (
        info_by_sg_setting_origin.get((sgno, setting_use, origin_use))
        if origin_use is not None
        else None
    ) or info_by_sg_setting.get((sgno, setting_use))

    if info is None:
        raise ValueError(
            f"{path.name}: No SPACEGRP.DAT entry found for SpacegroupNo={sgno}, Setting={setting_use}, OriginChoice={origin_use}"
        )

    # Optional: if Profex provides UniqueAxis, warn if it doesn't match the chosen entry
    if unique_axis_profex is not None and not ua_matches(info.get("unique_axis"), unique_axis_profex):
        sg_setting_note += (
            f"\n  /* WARNING: Profex UniqueAxis={unique_axis_profex} does not match SPACEGRP.DAT UniqueAxis={info.get('unique_axis')}. */"
        )

    # HM is provenance-only (keep for comments), but output HM comes from the selected SPACEGRP.DAT entry
    hm = info["hm"]
    centring_vecs = centring_vectors_from_hm(hm)
    hm_mismatch = (hm_profex is not None and hm_profex.strip() != hm)
    sg_hm_note = None
    if hm_profex is not None:
        sg_hm_note = None

    # Pull all symmetry/Wyckoff data from the selected entry (this is the point of fix 1)
    wy = info["wyckoff_proto"]
    wy_ops = info.get("wyckoff_ops", {})
    lattice = info.get("lattice", "")
    unique_axis = info.get("unique_axis", None)

    # Decide how to emit the space_group for TOPAS
    # Numeric output is ONLY safe when the resolved setting is canonical (Setting=1)
    # and no HM-based setting inference was required.

    if (
        use_numeric_spacegroup(lattice)
        and setting_use == 1
        and (
            setting_profex is not None
            or inferred_from_hm is None
        )
    ):
        # Safe numeric output (unique representation)
        space_group_out = str(sgno)
        space_group_mode = "numeric"
    else:
        # Force HM output to preserve setting / axis / origin conventions
        space_group_out = hm
        space_group_mode = "hermann_mauguin"

    # ---------------- Build TOPAS output ----------------
    any_mixed_expansion = any(a["mixed"] for a in atoms)

    out = [
        "/* Generated using ProfexStr_to_TOPAS_v2.py by Henry Spratt */",
        "",
        "str",
        f'  phase_name "{phase_topas}"',
        f'  space_group "{space_group_out}"'
    ]
    warnings = []
    
    # Safety/provenance note: ask user to confirm TOPAS accepts the symbol and intent matches Profex
    out.append("  /*")
    out.append("     Space-group entry for TOPAS was generated automatically.")
    out.append(f"     Output mode: {space_group_mode}.")
    out.append("     Please confirm that TOPAS accepts this symbol and that it")
    out.append("     matches the intended symmetry from the original Profex .str.")
    out.append("  */")

    # --- Profex str provenance (raw first, then interpreted) ---

    # 1) EXACT symmetry line as it appeared in the Profex .str
    if profex_sym_line:
        out.append(f"  /* Profex symmetry line (verbatim): {profex_sym_line} */")

    # 2) Parsed symmetry fields
    out.append(f"  /* Parsed SpacegroupNo={sgno} */")

    if setting_profex is not None:
        out.append(f"  /* Parsed Setting={setting_profex} */")

    if hm_profex is not None:
        out.append(f"  /* Parsed Profex HermannMauguin={hm_profex} */")

    if origin_profex is not None:
        out.append(f"  /* Parsed OriginChoice={origin_profex} */")

    # 3) Warnings and resolution notes
    if sg_fuzzy_used:
        out.append(
            "  /* WARNING: SpacegroupNo appears to be misspelled in Profex; "
            "parsed using fuzzy matching. Please correct the profex .str file. */"
        )

    if hm_mismatch:
        out.append(
            f"  /* WARNING: Profex HermannMauguin='{hm_profex}' does not match "
            f"SPACEGRP.DAT HermannMauguin='{hm}'. SPACEGRP.DAT definition was used. */"
        )

    # 4) Other metadata
    if formula_line:
        out.append(f"  /* {formula_line} */")

    out.append("")

    # Track whether we applied default bounds (for a single phase-level note)
    applied_default_cell_bounds = False

    for ax in ("a", "b", "c"):
        if ax in cell_ranges:
            mid, lo, hi = cell_ranges[ax]

            if lo is not None:
                # Profex explicitly provided bounds → respect them
                out.append(
                    f" {ax} {ax}_{phase_topas} {mid * NM_TO_A:.6f} "
                    f"min={lo * NM_TO_A:.6f}; max={hi * NM_TO_A:.6f};"
                )
            else:
                # v2_stable: apply default ±1% bounds
                lo_def = mid * (1.0 - DEFAULT_CELL_RELATIVE_TOL)
                hi_def = mid * (1.0 + DEFAULT_CELL_RELATIVE_TOL)
                applied_default_cell_bounds = True

                out.append(
                    f" {ax} {ax}_{phase_topas} {mid * NM_TO_A:.6f} "
                    f"min={lo_def * NM_TO_A:.6f}; max={hi_def * NM_TO_A:.6f};"
                )

    for ax in ("al", "be", "ga"):
        if ax in cell_ranges:
            mid, lo, hi = cell_ranges[ax]

            if lo is not None:
                out.append(
                    f" {ax} {ax}_{phase_topas} {mid:.6f} "
                    f"min={lo:.6f}; max={hi:.6f};"
                )
            else:
                lo_def = mid * (1.0 - DEFAULT_ANGLE_RELATIVE_TOL)
                hi_def = mid * (1.0 + DEFAULT_ANGLE_RELATIVE_TOL)
                applied_default_cell_bounds = True

                out.append(
                    f" {ax} {ax}_{phase_topas} {mid:.6f} "
                    f"min={lo_def:.6f}; max={hi_def:.6f};"
                )

    out.append("")
    out.append("  /* Unit cell constraints from SPACEGRP.DAT lattice system */")
    rhombo_axes = infer_rhombohedral_axes_from_cell(cell_ranges) and (setting_use == 2)
    apply_lattice_constraints(out, phase_topas, cell_ranges, lattice, unique_axis, rhombo_axes=rhombo_axes)

    out.append("")
    out.append("  /* Phase scale, MVW and crystallite size. Suggested CS_L Limits for FP */")
    out.append(f"  scale scale_{phase_topas} 0.001")
    out.append(f"  MVW(ZM_{phase_topas} 0, V_{phase_topas} 0, wp_{phase_topas} 0)")
    out.append(f"  CS_L(csl_{phase_topas}, 200 min 32 max 5000)")

    # -------- site-local occupancy parameter definitions (avoid collisions) --------
    site_param_name_map = {}   # (site_key, base) -> unique_name
    prm_defs = []              # list of (unique_name, mid, lo, hi)

    for a in atoms:
        sk = a.get("site_key")
        lp = a.get("local_params") or {}
        for base, (mid, lo, hi) in lp.items():
            key = (sk, base)
            if key in site_param_name_map:
                continue
            base_s = sanitize_param_name(base)
            # unique per Profex atom line (site_key) AND wyckoff letter for readability
            uniq = f"{base_s}_{a['wyck']}{sk}_{phase_topas}"
            site_param_name_map[key] = uniq
            prm_defs.append((uniq, mid, lo, hi))

    # attach per-atom map for expression rewriting
    for a in atoms:
        sk = a.get("site_key")
        a["_param_map"] = {}
        for base in (a.get("local_params") or {}).keys():
            a["_param_map"][base] = site_param_name_map[(sk, base)]

    # emit prm defs
    if prm_defs:
        out.append("")
        out.append(" /* Occupancy parameters from Profex (site-local) */")
        for name, mid, lo, hi in prm_defs:
            if lo is not None:
                out.append(f" prm {name} {mid:g} min {lo:g} max {hi:g}")
            else:
                out.append(f" prm {name} {mid:g}")
        out.append("")

    out.append("")
    out.append("  /* Atom sites */")
    out.append("  /* fractional coords output as n/d only when inserted from SPACEGRP.DAT or profex supplied high precision decimals */")

    # Create TOPAS site IDs per element
    counts = {}
    skipped_sites = 0
    any_sites_skipped = False
    any_coords_transformed = False
    any_constants_filled = False
    skip_reasons = {}

    for a in atoms:
        proto = wy.get(a["wyck"])
        if proto is None:
            available = "".join(sorted(wy.keys()))
            raise ValueError(
                f"{path.name}: Wyckoff letter '{a['wyck']}' not found for resolved space_group='{hm}'. "
                f"Available letters: {available}. Profex atom line: {a['raw']}"
            )

        ops = wy_ops.get(a["wyck"], None)

        given_num = {
            "x": a["given"]["x"][0],
            "y": a["given"]["y"][0],
            "z": a["given"]["z"][0],
        }
        given_txt = {
            "x": a["given"]["x"][1],
            "y": a["given"]["y"][1],
            "z": a["given"]["z"][1],
        }

        #if a["wyck"] == "c" and "z=1/4" in a["raw"]:
          #  out.append(f" /* DEBUG GIVEN (raw): {a['raw']} */")
           # out.append(f" /* DEBUG GIVEN_NUM: {given_num} */")
           # out.append(f" /* DEBUG GIVEN_TXT: {given_txt} */")

        is_general = (proto == ["x","y","z"])

        # Count skipped sites for a header warning
        # (Add `skipped_sites = 0` before the loop.)
        # skipped_sites = 0

        def _emit_skipped(reason: str):
            nonlocal skipped_sites, any_sites_skipped, skip_reasons
            skipped_sites += 1
            any_sites_skipped = True
            skip_reasons[reason] = skip_reasons.get(reason, 0) + 1
            out.append(f"  /* ERROR: {reason}. No TOPAS site was written. Original Profex atom line: */")
            out.append(f"  /* {a['raw']} */")

        def _mark_transformed():
            nonlocal any_coords_transformed
            any_coords_transformed = True

        def _mark_constant_filled():
            nonlocal any_constants_filled
            any_constants_filled = True

        # --- General position rule: must be complete ---
        if is_general:
            if (given_num["x"] is None or given_num["y"] is None or given_num["z"] is None):
                _emit_skipped("General Wyckoff position requires x,y,z but Profex is missing one or more coordinates")
                continue
            # General position: keep Profex coords as-is
            xo, yo, zo = given_txt["x"], given_txt["y"], given_txt["z"]
            mm = False
            details = []
        else:
            # --- Special position: validate against any triplet + infer dependent missing coords ---
            if ops is None:
                # Cannot validate or infer special-position coordinates without equivalent-position triplets.
                # Do NOT fall back to resolve_coords(), because proto tokens may include variables (e.g. y, -y+1/2)
                # and would cause float('y') errors or unintended enforcement.
                _emit_skipped(
                    "No SPACEGRP.DAT equivalent-position triplets available for this Wyckoff letter; "
                    "special-position validation/inference cannot be performed"
                )
                continue

            else:
                matches = []
                for trip in ops:
                    # Try each centring coset by shifting the given coords into the base set
                    for t in centring_vecs:
                        shifted_given = shift_given_xyz(given_num, t)
                        ok, out_xyz, vars_used = try_match_triplet(trip, shifted_given)

                        #if a["wyck"] == "c" and "z=1/4" in a["raw"]:
                          #  out.append(f" /* DEBUG trying triplet: {trip} shifted_given: {shifted_given} */")
                           # out.append(f" /* DEBUG try_match_triplet ok={ok} out_xyz={out_xyz} vars_used={vars_used} */")

                        if ok:
                            # Convert inferred coords back into original frame
                            out_xyz2 = unshift_out_xyz(out_xyz, t)
                            matches.append(out_xyz2)
                            break  # stop after first t that works for this triplet

                if not matches:
                    out.append(f"  /* DEBUG: SPACEGRP.DAT triplets for Wyckoff={a['wyck']} (first 5): {ops[:5]} */")
                    _emit_skipped("Special-position coordinates do not match any equivalent-position triplet in SPACEGRP.DAT")
                    continue

                # Accept the FIRST matching equivalent-position triplet as canonical.
                # Other symmetry mates (e.g. y vs -y) are algebraically equivalent and
                # must not cause rejection or conflict for special positions.
                filled = matches[0]

                # Missing coordinate policy: we only fill if it is determinable (not introducing a new free var)
                invalid = False
                for ax in ("x","y","z"):
                    if given_num[ax] is None:
                        idx = ("x","y","z").index(ax)

                        # If proto token is free and cannot be determined, error
                        if proto[idx] in ("x","y","z") and filled[ax] is None:
                            _emit_skipped(f"Missing free coordinate '{ax}' in special position cannot be inferred uniquely")
                            invalid = True
                            break

                if invalid:
                    continue  # <-- CRITICAL: do NOT emit a site line

                # Build output coords
                def fmt(ax):
                    idx = ("x", "y", "z").index(ax)

                    # --- FIXED-BY-SYMMETRY OVERRIDE ---
                    # If this coordinate is fixed by Wyckoff symmetry (proto token not x/y/z),
                    # then output the symmetry-consistent value from SPACEGRP.DAT matching (filled[ax]),
                    # not Profex's supplied numeric approximation.
                    if proto[idx] not in ("x", "y", "z"):
                        _mark_transformed()
                        _mark_constant_filled()
                        return float_to_topas_rational(filled[ax])  # uses your tol/max_den settings

                    # --- FREE COORDINATE (x/y/z) ---
                    # If Profex provided it, keep it (your parse-time rationalization already
                    # converts 0.333333 -> 1/3 etc., so this remains clean).
                    if given_txt[ax] is not None:
                        return given_txt[ax]

                    # Profex did not provide it. For free coordinates, we only allow filling
                    # if it is determinable; otherwise the caller should have already skipped.
                    _mark_transformed()
                    return float_to_topas_rational(filled[ax])

                xo, yo, zo = fmt("x"), fmt("y"), fmt("z")
                mm = False
                details = []

        # If this site got here, it is valid and we emit a TOPAS site line.
        counts[a["el"]] = counts.get(a["el"], 0) + 1
        site_id = f"{a['el']}{counts[a['el']]}"
        beq_txt = "1" if a["beq"] is None else f"{a['beq']:.4f}"

        # ---- occupancy output supports expressions (e.g. p, (1-p)) ----
        if a.get("occ_expr"):
            # a["_param_map"] is created earlier when you emit site-local prm definitions
            expr = rewrite_expr(a["occ_expr"], a.get("_param_map", {}))
            occ_clause = f"occ {a['el']} = {expr};"
        else:
            occ_clause = f"occ {a['el']} {a['occ']:g}"

        line = (
            f" site {site_id} "
            f"{topas_xyz_clause(xo, yo, zo)} "
            f"{occ_clause} "
            f"beq {beq_txt}"
        )

        if a["beq"] is None:
            line += " /* beq 1 inserted because this is the TOPAS GUI default */"

        # Always preserve exact Profex atom line
        line += f" /* {a['raw']} */"
        out.append(line)

    # Only warn globally if sites were skipped
    if skipped_sites > 0:
        warnings.append(
            f" /* WARNING: {skipped_sites} atom site(s) were skipped because the given Profex coordinates did not match "
            f"any equivalent-position triplet in SPACEGRP.DAT for the resolved space group entry. "
            f"The structure is ambiguous; verify with the Bilbao Crystallographic Server or other information. */"
        )

    # DO NOT emit global warnings for coordinate transforms or fixed constants
    # These are documented locally at each affected site

    # Header notes/warnings insertion (single block, appears near top)
    hdr = []

    if sg_setting_note:
        hdr += [sg_setting_note, ""]
    if sg_hm_note:
        hdr += [sg_hm_note, ""]

    if any_mixed_expansion:
        hdr += [
            "  /* NOTE: One or more mixed-occupancy atomic sites from Profex were expanded",
            "     into multiple TOPAS site entries with shared coordinates. */",
            ""
        ]

    if warnings:
        hdr = warnings + [""] + hdr

    if hdr:
        out = out[:3] + hdr + out[3:]

    metrics = {
        "input_file": str(path),
        "phase_name": phase_topas,
        "sgno": sgno,
        "setting_use": setting_use,
        "origin_use": origin_profex,
        "hm_used": hm,
        "lattice": lattice,
        "unique_axis": unique_axis,
        "skipped_sites": skipped_sites,
        "coords_transformed": any_coords_transformed,
        "constants_filled": any_constants_filled,
        "mixed_expansion": any_mixed_expansion,
        "skip_reasons": "; ".join(f"{k}={v}" for k, v in skip_reasons.items())
    }
    return phase_topas, "\n".join(out), metrics

# ---------------- CLI ----------------

def main():
    ap = argparse.ArgumentParser(description="Convert Profex/BGMN .str files to TOPAS .str files")
    ap.add_argument("--dat", default=DEFAULT_SPACEGRP_DAT, help="Path to SPACEGRP.DAT")
    ap.add_argument("--in", dest="in_dir", required=True, help="Input folder with Profex .str files")
    ap.add_argument("--out", dest="out_dir", required=True, help="Output folder for TOPAS .str files")
    args = ap.parse_args()

    dat = Path(args.dat)
    if not dat.exists():
        raise FileNotFoundError(f"SPACEGRP.DAT not found: {dat}")

    # hm_info, sg_default, hm_by_sg_setting, hm_by_sg_setting_origin = parse_spacegrp_dat(dat)  # old version of parse_spacegrp_dat

    (
        hm_info,
        sg_default,
        hm_by_sg_setting,
        hm_by_sg_setting_origin,
        info_by_sg_setting,
        info_by_sg_setting_origin,
    ) = parse_spacegrp_dat(dat)

    in_dir = Path(args.in_dir)

    # First pass: count PHASE= occurrences
    phase_counts = {}

    for f in sorted(in_dir.glob("*.str")):
        lines = f.read_text(errors="ignore").splitlines()
        for ln in lines:
            m = RE_PHASE.search(ln)
            if m:
                phase = m.group(1).strip()
                phase_counts[phase] = phase_counts.get(phase, 0) + 1
                break

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_metrics = []

    for f in sorted(in_dir.glob("*.str")):
        try:
            phase_topas, text, metrics = convert_one_str(
                f,
                hm_info,
                sg_default,
                hm_by_sg_setting,
                hm_by_sg_setting_origin,
                info_by_sg_setting,
                info_by_sg_setting_origin,
                phase_counts,
            )
            all_metrics.append(metrics)
            
            out_file = out_dir / f.name

            if out_file.exists():
                stem = out_file.stem
                suffix = out_file.suffix  # ".str"
                i = 1
                while True:
                    candidate = out_dir / f"{stem}_{i}{suffix}"
                    if not candidate.exists():
                        out_file = candidate
                        break
                    i += 1

            out_file.write_text(text, encoding="utf-8")
        except Exception as e:
            # Write per-file error report and continue
            (out_dir / f"{f.stem}_ERROR.txt").write_text(str(e), encoding="utf-8")

        # Write canary/regression report
        from datetime import datetime

        report_path = out_dir / f"_canary_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        with report_path.open("w", newline="", encoding="utf-8") as fp:
            fieldnames = [
                "input_file","phase_name","sgno","setting_use","origin_use","hm_used",
                "lattice","unique_axis","skipped_sites","coords_transformed",
                "constants_filled","mixed_expansion","skip_reasons"
            ]
            writer = csv.DictWriter(fp, fieldnames=fieldnames)
            writer.writeheader()
            for row in all_metrics:
                # Ensure missing keys don't crash CSV writing
                writer.writerow({k: row.get(k, "") for k in fieldnames})

        print(f"Wrote canary report: {report_path}")

        # Console summary
        n_total = len(all_metrics)
        n_skips = sum(1 for r in all_metrics if int(r.get("skipped_sites", 0) or 0) > 0)
        n_trans = sum(1 for r in all_metrics if r.get("coords_transformed"))
        n_const = sum(1 for r in all_metrics if r.get("constants_filled"))

        print(f"Canary summary: total={n_total}, files_with_skips={n_skips}, files_with_transforms={n_trans}, files_with_constants={n_const}")

        worst = sorted(all_metrics, key=lambda r: int(r.get("skipped_sites", 0) or 0), reverse=True)[:10]
        print("Top offenders (by skipped_sites):")
        for r in worst:
            if int(r.get("skipped_sites", 0) or 0) <= 0:
                break
            print(f"  {r.get('input_file')}  skipped_sites={r.get('skipped_sites')}  sgno={r.get('sgno')}  setting={r.get('setting_use')}  hm={r.get('hm_used')}")


if __name__ == "__main__":
    main()