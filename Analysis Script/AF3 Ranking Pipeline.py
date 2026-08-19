#!/usr/bin/env python3
"""
AF3 candidate QC + ranking metric pipeline for Rx/Sr35 vLRR screens.

VERSION 10 REVIEWED:
  - Retains the v9 median-based steric-clash QC.
  - Adds an AF3 summary-matrix order safeguard for chain_pair_* metrics.
  - Adds model-support fraction to the custom interface-consistency tiebreaker,
    so two agreeing interface-bearing seeds cannot look equivalent to five.
  - Author/native metrics remain delegated to official AFM-LIS / IPSAE or read
    directly from native AF3 JSON; project-derived metrics remain labelled.

VERSION 8 VERIFIED-OFFICIAL:
  - LIS/cLIS/iLIS/LIA/cLIA/iLIA/iLISA/LIR/cLIR/pLDDT come from the
    authors' official AFM-LIS lis.py, not a local reimplementation.
  - ipSAE/pDockQ/pDockQ2 come from DunbrackLab's official ipsae.py.
  - Native AF3 confidence fields are read directly from AF3 JSON.
  - Derived/custom metrics (BSA, 4.5A geometry, contact-probability
    summaries, vLRR fraction, Jaccard consistency, final ranking) are
    explicitly labelled as derived analyses; they are not claimed to be
    official AlphaFold/AFM-LIS/ipSAE metrics.
  - Windows UTF-8 and AFM-LIS output-directory fixes are retained.

=============================== VERSION 2 ===================================
Patched from the original pipeline. Every change is tagged with  # [v2]  so
you can diff against your friend's original.

WHAT CHANGED
------------
CRITICAL (affects your conclusions)
  1. contact_probs channel added. AF3's distogram-derived contact probability
     map lives in the SAME full-data JSON the script already loads. It was
     never read. New columns: contact_prob_max, n_contact_pairs_gt0.3,
     n_contact_pairs_gt0.5.
  2. Ranking now sorts on MEDIAN not MEAN (medians were computed but unused),
     and the second sort key is now the contact channel, which is independent
     of PAE, instead of ipSAE which is ~collinear with iLIS.
  3. vLRR localisation can now be COMPUTED from cLIR indices instead of read
     from a manual Yes/Partially/No annotation. Supply --vlrr-ranges. The
     continuous fraction is binned back into Yes/Partially/No so it behaves
     like the manual key, and manual-vs-computed agreement is reported.

SILENT-FAILURE FIXES
  4. mean_pairwise_jaccard returns None (not 0.0) when there is not enough
     data. 0.0 was indistinguishable from genuine seed disagreement.
  5. choose_lis_row now warns when it cannot match a model.
  6. compute_pae_qc now warns on shape mismatch instead of silently
     returning {}. PAE columns are now part of the critical-missing check.
  7. ipSAE output is parsed BY HEADER NAME, with the positional layout only
     as a fallback. Tool URLs can be pinned to a git ref and verified against
     an expected SHA256.
  8. SASA: if the three subset calculations do not all use the same backend,
     BSA is set to None instead of subtracting two different SASA
     definitions from each other.

DESIGN / CONSISTENCY
  9. Chain assignment is still positional (chain1=NLR, chain2=effector,
     chain3=ATP) but token counts per chain are now reported and sanity
     checked, with a warning if the NLR chain is not the largest.
 10. rank_within_effector added alongside rank_within_backbone. Comparing
     absolute iLIS across different effectors is confounded; the
     within-effector rank is the defensible one.
 11. annotation_warning (lookup effector mismatch) now surfaces in QC_notes.
 12. BSA convention is stated explicitly in tool_provenance.txt.
     BSA_total_A2        = SASA_A + SASA_B - SASA_AB   (total buried)
     BSA_interface_area_A2 = BSA_total / 2             (per-side convention)
     Make sure whichever one you quote matches your Rx benchmark table.

IMPORTANT ASSUMPTION (unchanged, matching this project):
  chain 1 = NLR
  chain 2 = effector
  chain 3 = ATP
All inter-protein ranking metrics are restricted to chain1-chain2.

Recalculation sources:
  * iLIS/LIS/cLIS/LIR/cLIR/per-chain pLDDT: official AFM-LIS lis.py
    https://github.com/flyark/AFM-LIS
  * ipSAE: official DunbrackLab IPSAE ipsae.py
    https://github.com/DunbrackLab/IPSAE
  * pTM/ipTM/chain-pair ipTM/has_clash/PAE/contact_probs: native AF3 JSON
    https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md
  * BSA: FreeSASA when available; Bio.PDB ShrakeRupley fallback otherwise.

Outputs:
  raw_data.csv                one row per AF3 model/sample
  summary_mean.csv            one row per candidate/backbone/effector
  qc.csv                      QC-focused table
  AF3_candidate_metrics.xlsx  same three tables as sheets
  tool_provenance.txt         tool refs, SHA256s, all cutoffs, BSA convention

Example (Windows):
  python af3_candidate_ranking_pipeline_v2.py ^
      --roots "D:\\AF3\\Rx" "D:\\AF3\\Sr35" ^
      --interface-lookup "D:\\AF3\\vlrr_interface_lookup.csv" ^
      --out "D:\\AF3\\ranking_analysis" ^
      --vlrr-ranges "Rx=489-593" "Sr35=..." ^
      --lis-workers 2

Install dependencies:
  python -m pip install numpy scipy biopython xlsxwriter
  python -m pip install freesasa   # optional but preferred for BSA
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
from scipy.spatial import cKDTree

try:
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Polypeptide import is_aa
except ImportError as exc:
    raise SystemExit("Biopython is required: python -m pip install biopython") from exc


# ------------------------------- configuration ------------------------------

# [v2] Tool sources are now pinnable. "main" is a REPRODUCIBILITY HAZARD:
# re-running in six months can silently give you different numbers because
# upstream changed. Replace these with a commit SHA once you are happy, e.g.
#   AFM_LIS_REF = "3f9c1ab..."
# and/or lock the file hash with --expect-lis-sha256 / --expect-ipsae-sha256.
AFM_LIS_REF_DEFAULT = "main"
IPSAE_REF_DEFAULT = "main"
AFM_LIS_URL_TMPL = "https://raw.githubusercontent.com/flyark/AFM-LIS/{ref}/lis.py"
IPSAE_URL_TMPL = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/{ref}/ipsae.py"

LIS_PAE_CUTOFF = 12.0   # native iLIS definition; also what Rx_benchmarking_table.py used
LIS_CB_CUTOFF = 8.0     # native cLIS/cLIR definition
# [v3] WAS 15.0. Rx_benchmarking_table.py used IPSAE_CUTOFF = 10.0, so 15 would
# have made the two tables incomparable. Overridable with --ipsae-pae-cutoff.
IPSAE_PAE_CUTOFF = 10.0
IPSAE_DIST_CUTOFF = 10.0
HEAVY_ATOM_CLASH_CUTOFF = 1.5
EXPECTED_MODELS_DEFAULT = 5

# [v3] SASA settings copied verbatim from Rx_benchmarking_table.py so the BSA
# column of the candidate table is numerically comparable with the benchmark
# table. Shrake-Rupley, 1.4 A probe, 256-point Fibonacci sphere, Bondi radii,
# heavy atoms only.
BONDI_VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
             "F": 1.47, "CL": 1.75, "MG": 1.73}
SASA_PROBE = 1.4
SASA_SPHERE_POINTS = 256
SASA_DEFAULT_RADIUS = 1.70

# [v2] contact-probability thresholds. 0.3 reproduces the column used in the
# Rx benchmark table ("AF3 contact pairs > 0.3").
CONTACT_PROB_THRESHOLDS = (0.3, 0.5)

# [v2] vLRR fraction -> Yes / Partially / No bin edges (overridable on CLI).
VLRR_BIN_LOW_DEFAULT = 0.20
VLRR_BIN_HIGH_DEFAULT = 0.60

# [v2] Fallback column layout for ipsae.py output, used ONLY if no header
# line is found. Header-name parsing is preferred.
IPSAE_FALLBACK_COLUMNS = [
    "Chn1", "Chn2", "PAE", "Dist", "Type",
    "ipSAE", "ipSAE_d0chn", "ipSAE_d0dom",
    "ipTM_af", "ipTM_d0chn",
    "pDockQ", "pDockQ2", "LIS",
    "n0res", "n0chn", "n0dom",
    "d0res", "d0chn", "d0dom",
    "nres1", "nres2", "dist1", "dist2", "Model",
]

# Numeric columns for model -> candidate mean aggregation.
SUMMARY_NUMERIC_COLUMNS = [
    "af3_global_iptm", "af3_global_ptm", "af3_ranking_score",
    "pairwise_iptm_nlr_effector", "pairwise_pae_min_nlr_to_effector",
    "pairwise_pae_min_effector_to_nlr",
    "iLIS", "LIS", "cLIS", "LIA", "cLIA", "iLIA", "iLISA",
    "afm_lis_ipSAE_crosscheck",
    "ipSAE_official_max", "ipSAE_official_d0dom_max", "ipSAE_official_d0chn_max",
    "ipSAE_official_nlr_to_effector", "ipSAE_official_effector_to_nlr",
    "pDockQ_official", "pDockQ2_official",                      # [v2]
    "contact_prob_max", "contact_prob_mean",                     # [v2]
    "n_contact_pairs_gt0.3", "n_contact_pairs_gt0.5",            # [v2]
    "vLRR_fraction", "vLRR_cLIR_n_in_range",                     # [v2]
    "nlr_pLDDT", "effector_pLDDT",
    "LIR_nlr", "LIR_effector", "cLIR_nlr", "cLIR_effector",
    "LIpLDDT_nlr", "LIpLDDT_effector", "cLIpLDDT_nlr", "cLIpLDDT_effector",
    "nlr_intra_PAE_mean", "nlr_intra_PAE_median", "nlr_intra_PAE_p90",
    "effector_intra_PAE_mean", "effector_intra_PAE_median", "effector_intra_PAE_p90",
    "PAE_nlr_to_effector_mean", "PAE_nlr_to_effector_min",
    "PAE_effector_to_nlr_mean", "PAE_effector_to_nlr_min",
    "AB_heavy_atom_clash_pairs_lt1p5A", "AB_min_heavy_atom_distance_A",
    "SASA_nlr_A2", "SASA_effector_A2", "SASA_AB_A2",
    "BSA_total_A2", "BSA_interface_area_A2", "BSA_fraction_combined",
    "n_tokens_nlr", "n_tokens_effector", "n_tokens_atp",          # [v2]
    "iface_res_nlr_n", "iface_res_effector_n",                    # [v4]
    "geometric_contact_pairs_lt4p5A",                             # [v4]
]

# [v2] Columns that also get median + SD, because ranking uses medians.
MEDIAN_SD_COLUMNS = [
    "iLIS",
    "ipSAE_official_max",
    "n_contact_pairs_gt0.3",
    "contact_prob_max",
    "vLRR_fraction",
    "AB_heavy_atom_clash_pairs_lt1p5A",   # [v9] clash QC is now median-based
]


# -------------------------------- utilities ---------------------------------

def eprint(*args: Any, **kwargs: Any) -> None:
    print(*args, file=sys.stderr, **kwargs)


def norm_candidate(value: Any) -> str:
    s = str(value).strip()
    if not s:
        return ""
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except ValueError:
        pass
    return s.upper() if s.lower() == "wt" else s


def norm_effector(value: Any) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(value)).upper()


def norm_backbone(value: Any) -> str:
    s = str(value).strip().lower()
    if s == "rx":
        return "Rx"
    if s == "sr35":
        return "Sr35"
    return str(value).strip()


def as_float(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        f = float(value)
        if math.isnan(f):
            return None
        return f
    except (TypeError, ValueError):
        return None


def safe_mean(values: Iterable[Any]) -> Optional[float]:
    xs = [as_float(x) for x in values]
    xs = [x for x in xs if x is not None]
    return float(np.mean(xs)) if xs else None


def safe_median(values: Iterable[Any]) -> Optional[float]:
    xs = [as_float(x) for x in values]
    xs = [x for x in xs if x is not None]
    return float(np.median(xs)) if xs else None


def safe_sd(values: Iterable[Any]) -> Optional[float]:
    xs = [as_float(x) for x in values]
    xs = [x for x in xs if x is not None]
    return float(np.std(xs, ddof=1)) if len(xs) >= 2 else (0.0 if len(xs) == 1 else None)


def first_seen(seq: Iterable[str]) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def child_env() -> Dict[str, str]:
    """
    [v5] Environment for the lis.py / ipsae.py subprocesses.

    On a Chinese-locale Windows the console codepage is GBK (cp936). A child
    Python writing to a pipe picks its stdout encoding from that locale, so
    lis.py's progress bar (which uses U+2591 block characters) raises
    UnicodeEncodeError and the whole tool crashes -- taking iLIS, cLIR and the
    per-chain pLDDT columns with it. Forcing UTF-8 in the child fixes it
    without touching the user's system settings.
    """
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    return env


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_tool(url: str, dest: Path, refresh: bool = False,
                expect_sha256: Optional[str] = None) -> Path:
    """[v2] Download if needed, then optionally verify against a pinned hash."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if refresh or not dest.exists():
        eprint(f"Downloading official tool: {url}")
        urllib.request.urlretrieve(url, dest)
    if expect_sha256:
        actual = sha256_file(dest)
        if actual.lower() != expect_sha256.strip().lower():
            raise SystemExit(
                f"TOOL HASH MISMATCH for {dest.name}\n"
                f"  expected: {expect_sha256}\n"
                f"  actual:   {actual}\n"
                f"Upstream changed, or you pinned the wrong hash. Refusing to run: "
                f"a silently-changed tool would give you numbers that are not "
                f"comparable with your previous analysis."
            )
    return dest


def parse_job_name(folder_name: str) -> Optional[Tuple[str, str, str]]:
    """Parse examples: rx52_and_sre19_atp, sr35_52_and_sre19_atp."""
    # [v3] accept _adp too, so the Rx benchmark models can be run through the
    # SAME pipeline as the candidates instead of a separate hand-written script.
    m = re.match(r"^(rx|sr35)_?(.*?)_and_(.*?)_(?:atp|adp)$", folder_name, flags=re.I)
    if not m:
        return None
    backbone = norm_backbone(m.group(1))
    candidate = norm_candidate(m.group(2))
    effector = m.group(3).strip()
    if not candidate or not effector:
        return None
    return backbone, candidate, effector


def parse_range_spec(spec: str) -> Set[int]:
    """[v2] '489-593,700,810-820' -> set of ints."""
    out: Set[int] = set()
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            a, _, b = part.partition("-")
            try:
                aa, bb = int(a), int(b)
            except ValueError:
                raise SystemExit(f"Bad residue range: {part!r}")
            out.update(range(min(aa, bb), max(aa, bb) + 1))
        else:
            try:
                out.add(int(part))
            except ValueError:
                raise SystemExit(f"Bad residue number: {part!r}")
    return out


def parse_vlrr_ranges(specs: Optional[Sequence[str]]) -> Dict[str, Set[int]]:
    """[v2] --vlrr-ranges "Rx=489-593" "Sr35=500-620,700-740" """
    out: Dict[str, Set[int]] = {}
    for spec in specs or []:
        if "=" not in spec:
            raise SystemExit(
                f"--vlrr-ranges entry must look like Backbone=start-end[,start-end]: {spec!r}"
            )
        bb, _, rng = spec.partition("=")
        bb = norm_backbone(bb)
        residues = parse_range_spec(rng)
        if not residues:
            raise SystemExit(f"--vlrr-ranges entry produced no residues: {spec!r}")
        out[bb] = residues
    return out


def load_interface_lookup(path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    lookup: Dict[Tuple[str, str], Dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            backbone = norm_backbone(row.get("backbone", ""))
            candidate = norm_candidate(row.get("candidate_id", ""))
            effector = row.get("effector", "").strip()
            cls = row.get("vLRR_Interface", "").strip()
            if backbone and candidate:
                lookup[(backbone, candidate)] = {
                    "effector": effector,
                    "vLRR_Interface": cls,
                }
    return lookup


def discover_jobs(roots: Sequence[Path]) -> List[Dict[str, Any]]:
    jobs: List[Dict[str, Any]] = []
    seen: Set[Path] = set()
    for root in roots:
        root = root.resolve()
        if not root.exists():
            eprint(f"WARNING: root not found: {root}")
            continue
        candidates = [root] if parse_job_name(root.name) else []
        candidates.extend(p for p in root.rglob("*") if p.is_dir() and parse_job_name(p.name))
        for p in candidates:
            rp = p.resolve()
            if rp in seen:
                continue
            seen.add(rp)
            parsed = parse_job_name(p.name)
            if parsed:
                backbone, candidate, effector = parsed
                jobs.append({
                    "job_dir": rp,
                    "folder_name": p.name,
                    "backbone": backbone,
                    "candidate_id": candidate,
                    "effector": effector,
                })
    jobs.sort(key=lambda j: (j["backbone"],
                             (int(j["candidate_id"]) if j["candidate_id"].isdigit() else 10**9),
                             j["candidate_id"]))
    return jobs


def find_model_files(job_dir: Path) -> List[Dict[str, Path]]:
    """Support AF3 Server flat output and official local seed/sample output."""
    cifs = list(job_dir.rglob("*.cif"))
    cifs = [p for p in cifs if "_tools" not in p.parts and "ranking_analysis" not in p.parts]

    nested = [p for p in cifs if re.search(r"seed-\d+_sample-\d+", p.as_posix())]
    flat_numbered = [p for p in cifs if re.search(r"_model_\d+\.cif$", p.name)]
    selected = flat_numbered if flat_numbered else (nested if nested else cifs)

    models: List[Dict[str, Path]] = []
    for cif in selected:
        full: Optional[Path] = None
        summary: Optional[Path] = None
        model_id = ""

        m = re.match(r"^(.*)_model_(\d+)\.cif$", cif.name)
        if m:
            prefix, idx = m.group(1), m.group(2)
            model_id = f"model_{idx}"
            full_cands = [
                cif.with_name(f"{prefix}_full_data_{idx}.json"),
                cif.with_name(f"{prefix}_confidences_{idx}.json"),
            ]
            summary_cands = [cif.with_name(f"{prefix}_summary_confidences_{idx}.json")]
            full = next((x for x in full_cands if x.exists()), None)
            summary = next((x for x in summary_cands if x.exists()), None)
        else:
            m2 = re.match(r"^(.*)_model\.cif$", cif.name)
            if m2:
                prefix = m2.group(1)
                full = cif.with_name(f"{prefix}_confidences.json")
                summary = cif.with_name(f"{prefix}_summary_confidences.json")
                if not full.exists():
                    full = None
                if not summary.exists():
                    summary = None
            elif cif.name == "model.cif":
                full = cif.with_name("confidences.json")
                summary = cif.with_name("summary_confidences.json")
                if not full.exists():
                    full = None
                if not summary.exists():
                    summary = None
            sm = re.search(r"seed-(\d+)_sample-(\d+)", cif.as_posix())
            model_id = f"seed-{sm.group(1)}_sample-{sm.group(2)}" if sm else cif.stem

        if full and summary:
            models.append({"cif": cif, "full": full, "summary": summary, "model_id": Path(model_id)})
        else:
            eprint(f"WARNING: incomplete AF3 model triplet, skipping: {cif}")

    def model_sort_key(m: Dict[str, Path]) -> Tuple[int, str]:
        s = str(m["model_id"])
        q = re.search(r"(?:model_|sample-)(\d+)", s)
        return (int(q.group(1)) if q else 999999, s)

    models.sort(key=model_sort_key)
    return models


def read_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def chain_order_from_json(summary: Dict[str, Any], full: Dict[str, Any]) -> List[str]:
    """
    [v5] Distinct chain IDs in input/matrix order.

    MY BUG, and a bad one: the previous version returned summary["chain_ids"]
    verbatim. In this AF3 output that key is a PER-TOKEN list, so the result
    was ["A","A","A",...] and chain_order[0] == chain_order[1] == "A" -- i.e.
    the NLR and the "effector" were the same chain. Every interface number in
    such a run is meaningless, not merely missing.

    token_chain_ids is now the primary source: de-duplicated it gives exactly
    the chain order that summary's chain_pair_* matrices are indexed by. Every
    branch is de-duplicated, so a per-token list can never leak through again.
    """
    token_ids = full.get("token_chain_ids", [])
    if token_ids:
        return first_seen(str(x) for x in token_ids)
    ids = summary.get("chain_ids")
    if isinstance(ids, list) and ids:
        return first_seen(str(x) for x in ids)
    atom_ids = full.get("atom_chain_ids", [])
    if atom_ids:
        return first_seen(str(x) for x in atom_ids)
    return []


def matrix_chain_order_from_summary(summary: Dict[str, Any],
                                    full: Dict[str, Any]) -> List[str]:
    """
    [v10] Chain order used ONLY to index AF3 summary chain_pair_* matrices.

    Current AF3 documentation defines summary['chain_ids'] as the chain order
    for chain-level and chain-pair matrices. Some AF3 Server exports encountered
    in this project exposed repeated IDs instead, so we accept summary chain_ids
    only when it is already one ID per distinct chain; otherwise we fall back to
    de-duplicated token_chain_ids.
    """
    ids = summary.get("chain_ids")
    if isinstance(ids, list) and ids:
        vals = [str(x) for x in ids]
        uniq = first_seen(vals)
        if len(vals) == len(uniq):
            return vals
    token_ids = full.get("token_chain_ids", [])
    if token_ids:
        return first_seen(str(x) for x in token_ids)
    atom_ids = full.get("atom_chain_ids", [])
    if atom_ids:
        return first_seen(str(x) for x in atom_ids)
    return []


def chain_token_counts(full: Dict[str, Any]) -> Dict[str, int]:
    """[v2] Tokens per chain, for the chain-assignment sanity check."""
    tci = full.get("token_chain_ids", [])
    return dict(Counter(str(x) for x in tci))


def pae_block_stats(pae: np.ndarray, row_mask: np.ndarray, col_mask: np.ndarray,
                    offdiag: bool = False) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    block = pae[np.ix_(row_mask, col_mask)]
    if block.size == 0:
        return None, None, None
    if offdiag and block.shape[0] == block.shape[1]:
        mask = ~np.eye(block.shape[0], dtype=bool)
        vals = block[mask]
    else:
        vals = block.reshape(-1)
    if vals.size == 0:
        return None, None, None
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return None, None, None
    return float(np.mean(vals)), float(np.median(vals)), float(np.percentile(vals, 90))


def compute_pae_qc(full: Dict[str, Any], nlr_chain: str, eff_chain: str,
                   label: str = "") -> Dict[str, Optional[float]]:
    pae = np.asarray(full.get("pae", []), dtype=float)
    token_chain_ids = np.asarray(full.get("token_chain_ids", []), dtype=object)
    # [v2] used to return {} silently on any mismatch
    if pae.ndim != 2:
        eprint(f"WARNING: PAE matrix missing or not 2-D ({label}); PAE columns will be blank")
        return {}
    if token_chain_ids.size != pae.shape[0]:
        eprint(f"WARNING: token_chain_ids ({token_chain_ids.size}) does not match PAE "
               f"dimension ({pae.shape[0]}) ({label}); PAE columns will be blank")
        return {}
    m_n = token_chain_ids == nlr_chain
    m_e = token_chain_ids == eff_chain
    if not m_n.any() or not m_e.any():
        eprint(f"WARNING: chain {nlr_chain!r} or {eff_chain!r} has no tokens ({label})")
        return {}
    n_mean, n_med, n_p90 = pae_block_stats(pae, m_n, m_n, offdiag=True)
    e_mean, e_med, e_p90 = pae_block_stats(pae, m_e, m_e, offdiag=True)
    ne_mean, ne_med, _ = pae_block_stats(pae, m_n, m_e)
    en_mean, en_med, _ = pae_block_stats(pae, m_e, m_n)
    ne = pae[np.ix_(m_n, m_e)]
    en = pae[np.ix_(m_e, m_n)]
    return {
        "nlr_intra_PAE_mean": n_mean,
        "nlr_intra_PAE_median": n_med,
        "nlr_intra_PAE_p90": n_p90,
        "effector_intra_PAE_mean": e_mean,
        "effector_intra_PAE_median": e_med,
        "effector_intra_PAE_p90": e_p90,
        "PAE_nlr_to_effector_mean": ne_mean,
        "PAE_nlr_to_effector_median": ne_med,
        "PAE_nlr_to_effector_min": float(np.min(ne)) if ne.size else None,
        "PAE_effector_to_nlr_mean": en_mean,
        "PAE_effector_to_nlr_median": en_med,
        "PAE_effector_to_nlr_min": float(np.min(en)) if en.size else None,
    }


def compute_contact_metrics(full: Dict[str, Any], nlr_chain: str, eff_chain: str,
                            thresholds: Sequence[float] = CONTACT_PROB_THRESHOLDS,
                            label: str = "") -> Dict[str, Any]:
    """
    [v2] AF3 distogram-derived contact probability, restricted to the
    NLR-effector block. This is the channel that separated Rx+PVX-CP
    (max 0.76, 36 pairs > 0.3) from Rx+SRE23 (max 0.14, 0 pairs > 0.3).
    It is INDEPENDENT of PAE, so it is not covered by iLIS/ipSAE.
    """
    cp = np.asarray(full.get("contact_probs", []), dtype=float)
    tci = np.asarray(full.get("token_chain_ids", []), dtype=object)
    if cp.ndim != 2:
        eprint(f"WARNING: contact_probs missing from full-data JSON ({label}); "
               f"contact columns will be blank")
        return {}
    if tci.size != cp.shape[0]:
        eprint(f"WARNING: token_chain_ids ({tci.size}) does not match contact_probs "
               f"dimension ({cp.shape[0]}) ({label})")
        return {}
    m_n = tci == nlr_chain
    m_e = tci == eff_chain
    if not m_n.any() or not m_e.any():
        return {}
    block = cp[np.ix_(m_n, m_e)]
    if block.size == 0:
        return {}
    block = block[np.isfinite(block)] if block.ndim == 1 else np.where(np.isfinite(block), block, 0.0)
    out: Dict[str, Any] = {
        "contact_prob_max": float(np.max(block)),
        "contact_prob_mean": float(np.mean(block)),
    }
    for t in thresholds:
        out[f"n_contact_pairs_gt{t}"] = int(np.count_nonzero(block > t))
    return out


def matrix_value(summary: Dict[str, Any], key: str, i: int, j: int) -> Optional[float]:
    try:
        return float(summary[key][i][j])
    except (KeyError, IndexError, TypeError, ValueError):
        return None


def run_lis_for_job(lis_script: Path, job_dir: Path, out_csv: Path, workers: int,
                    extra_args: Optional[Sequence[str]] = None) -> List[Dict[str, str]]:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(lis_script), str(job_dir),
        "-o", out_csv.name,
        "-d", str(out_csv.parent),
        "--pae-cutoff", str(LIS_PAE_CUTOFF),
        "--cb-cutoff", str(LIS_CB_CUTOFF),
        "--no-skip-existing",
        "-w", str(max(1, workers)),
    ]
    cmd.extend(extra_args or [])   # [v5] e.g. --lis-extra-args --platform af3
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          encoding="utf-8", errors="replace", env=child_env())  # [v5]
    if proc.returncode != 0:
        tail = (proc.stderr or "")[-2000:]
        eprint(f"WARNING: AFM-LIS failed for {job_dir.name}:\n{tail}")
        if "UnicodeEncodeError" in tail or "codec can't encode" in tail:
            eprint("HINT: that is a console-encoding crash inside lis.py, not a data problem. "
                   "This build already forces PYTHONIOENCODING=utf-8 for the child process; "
                   "if you still see it, run `chcp 65001` in the terminal first.")
        return []
    if not out_csv.exists():
        eprint(f"WARNING: AFM-LIS produced no CSV for {job_dir.name}")
        return []
    with out_csv.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    # [v2] fail loudly if the expected columns are not there at all
    if rows:
        missing = [c for c in ("chain_i", "chain_j", "iLIS") if c not in rows[0]]
        if missing:
            eprint(f"ERROR: AFM-LIS output for {job_dir.name} is missing expected columns "
                   f"{missing}. Present columns: {sorted(rows[0].keys())}. "
                   f"Upstream lis.py has probably changed its output schema -- pin the tool ref.")
        if "cLIR_indices_i" not in rows[0]:
            eprint(f"WARNING: AFM-LIS output for {job_dir.name} has no 'cLIR_indices_i' column. "
                   f"interface_consistency and vLRR_fraction cannot be computed.")
    return rows


def choose_lis_row(rows: List[Dict[str, str]], cif: Path, job_dir: Path,
                   nlr_chain: str, eff_chain: str) -> Optional[Dict[str, str]]:
    if not rows:
        return None
    rel = cif.relative_to(job_dir).as_posix()
    pair_rows = [r for r in rows if {r.get("chain_i"), r.get("chain_j")} == {nlr_chain, eff_chain}]
    exact = [r for r in pair_rows if str(r.get("structure_file", "")).replace("\\", "/").endswith(rel)]
    if len(exact) == 1:
        return exact[0]
    base = [r for r in pair_rows if Path(str(r.get("structure_file", ""))).name == cif.name]
    if len(base) == 1:
        return base[0]
    mm = re.search(r"_model_(\d+)\.cif$", cif.name)
    if mm:
        num = mm.group(1)
        cand = [r for r in pair_rows if r.get("model", "") == num]
        if len(cand) == 1:
            return cand[0]
    # [v2] used to return None silently
    eprint(f"WARNING: could not match an AFM-LIS row for {cif.name} "
           f"(chains {nlr_chain}/{eff_chain}; {len(pair_rows)} candidate rows). "
           f"iLIS/cLIR/pLDDT will be blank for this model.")
    return None


def parse_lis_row(row: Optional[Dict[str, str]], nlr_chain: str, eff_chain: str) -> Dict[str, Any]:
    if not row:
        return {}

    # [v8] Integrity checks against the official AFM-LIS definitions.
    # This catches a parser/schema mismatch instead of letting a wrong column
    # silently enter the ranking.
    lis_v = as_float(row.get("LIS"))
    clis_v = as_float(row.get("cLIS"))
    ilis_v = as_float(row.get("iLIS"))
    lia_v = as_float(row.get("LIA"))
    clia_v = as_float(row.get("cLIA"))
    ilia_v = as_float(row.get("iLIA"))
    if lis_v is None or clis_v is None or ilis_v is None:
        eprint("ERROR: official AFM-LIS row is missing LIS/cLIS/iLIS; refusing this row.")
        return {}
    expected_ilis = math.sqrt(max(0.0, lis_v * clis_v))
    if abs(ilis_v - expected_ilis) > 5e-4:  # official CSV is rounded to 4 decimals
        eprint(f"ERROR: AFM-LIS iLIS integrity check failed: output={ilis_v}, "
               f"sqrt(LIS*cLIS)={expected_ilis}. Refusing this row.")
        return {}
    if lia_v is not None and clia_v is not None and ilia_v is not None:
        expected_ilia = math.sqrt(max(0.0, lia_v * clia_v))
        if abs(ilia_v - expected_ilia) > 0.2:  # official iLIA CSV rounded to 0.1
            eprint(f"ERROR: AFM-LIS iLIA integrity check failed: output={ilia_v}, "
                   f"sqrt(LIA*cLIA)={expected_ilia}. Refusing this row.")
            return {}
    nlr_is_i = row.get("chain_i") == nlr_chain

    def pick(a: str, b: str) -> str:
        return row.get(a if nlr_is_i else b, "")

    return {
        "iLIS": as_float(row.get("iLIS")),
        "iLIA": as_float(row.get("iLIA")),
        "iLISA": as_float(row.get("iLISA")),
        "afm_lis_ipSAE_crosscheck": as_float(row.get("ipSAE")),
        "LIS": as_float(row.get("LIS")),
        "cLIS": as_float(row.get("cLIS")),
        "LIA": as_float(row.get("LIA")),
        "cLIA": as_float(row.get("cLIA")),
        "nlr_pLDDT": as_float(pick("pLDDT_i", "pLDDT_j")),
        "effector_pLDDT": as_float(pick("pLDDT_j", "pLDDT_i")),
        "LIR_nlr": as_float(pick("LIR_i", "LIR_j")),
        "LIR_effector": as_float(pick("LIR_j", "LIR_i")),
        "cLIR_nlr": as_float(pick("cLIR_i", "cLIR_j")),
        "cLIR_effector": as_float(pick("cLIR_j", "cLIR_i")),
        "LIpLDDT_nlr": as_float(pick("LIpLDDT_i", "LIpLDDT_j")),
        "LIpLDDT_effector": as_float(pick("LIpLDDT_j", "LIpLDDT_i")),
        "cLIpLDDT_nlr": as_float(pick("cLIpLDDT_i", "cLIpLDDT_j")),
        "cLIpLDDT_effector": as_float(pick("cLIpLDDT_j", "cLIpLDDT_i")),
        "LIR_indices_nlr": pick("LIR_indices_i", "LIR_indices_j"),
        "LIR_indices_effector": pick("LIR_indices_j", "LIR_indices_i"),
        "cLIR_indices_nlr": pick("cLIR_indices_i", "cLIR_indices_j"),
        "cLIR_indices_effector": pick("cLIR_indices_j", "cLIR_indices_i"),
        "afm_lis_pairwise_ipTM": as_float(row.get("ipTM")),
        "afm_lis_pTM": as_float(row.get("pTM")),
    }


def hardlink_or_copy(src: Path, dst: Path) -> None:
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def parse_ipsae_txt(path: Path, nlr_chain: str, eff_chain: str) -> Dict[str, Any]:
    """
    [v2] Parse the ipsae.py .txt output BY HEADER NAME.

    The original parsed by hard-coded column index (parts[5], parts[7],
    parts[19] ...). Those indices are correct for the current ipsae.py, but
    combined with auto-downloading the latest version, any upstream column
    insertion would silently make every ipSAE value a different quantity --
    and it would still look like a plausible 0-1 number. Header parsing makes
    that failure mode impossible; the positional layout is only a fallback.
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    col: Optional[Dict[str, int]] = None
    for line in lines:
        parts = line.split()
        if len(parts) >= 6 and parts[0].strip().lower() in {"chn1", "chain1"} \
                and any(p.strip().lower() == "ipsae" for p in parts):
            col = {name: i for i, name in enumerate(parts)}
            break
    header_found = col is not None
    if col is None:
        # [v8] Fail closed.  A positional fallback can silently parse the wrong
        # quantity after an upstream schema change while still producing a
        # plausible-looking 0-1 number.
        eprint(f"ERROR: no recognizable header line found in official ipSAE output {path.name}. "
               f"Refusing positional parsing; pin/check the official tool version.")
        return {"ipsae_header_parsed": False}

    def gf(parts: List[str], name: str) -> Optional[float]:
        i = col.get(name)
        if i is None or i >= len(parts):
            return None
        try:
            return float(parts[i])
        except ValueError:
            return None

    def gi(parts: List[str], name: str) -> Optional[int]:
        v = gf(parts, name)
        return int(v) if v is not None else None

    type_idx = col.get("Type", 4)
    c1_idx = col.get("Chn1", 0)
    c2_idx = col.get("Chn2", 1)

    rows: List[Dict[str, Any]] = []
    for line in lines:
        parts = line.split()
        if len(parts) <= max(type_idx, c1_idx, c2_idx):
            continue
        if parts[type_idx] not in {"asym", "max"}:
            continue
        c1, c2 = parts[c1_idx], parts[c2_idx]
        if {c1, c2} != {nlr_chain, eff_chain}:
            continue
        ipsae = gf(parts, "ipSAE")
        if ipsae is None:
            continue
        rows.append({
            "c1": c1, "c2": c2, "type": parts[type_idx],
            "ipSAE": ipsae,
            "ipSAE_d0chn": gf(parts, "ipSAE_d0chn"),
            "ipSAE_d0dom": gf(parts, "ipSAE_d0dom"),
            "pDockQ": gf(parts, "pDockQ"),
            "pDockQ2": gf(parts, "pDockQ2"),
            "LIS_ipsae_script": gf(parts, "LIS"),
            "nres1": gi(parts, "nres1"), "nres2": gi(parts, "nres2"),
            "dist1": gi(parts, "dist1"), "dist2": gi(parts, "dist2"),
        })

    out: Dict[str, Any] = {"ipsae_header_parsed": header_found}
    for r in rows:
        if r["type"] == "max":
            out.update({
                "ipSAE_official_max": r["ipSAE"],
                "ipSAE_official_d0chn_max": r["ipSAE_d0chn"],
                "ipSAE_official_d0dom_max": r["ipSAE_d0dom"],
                "pDockQ_official": r["pDockQ"],
                "pDockQ2_official": r["pDockQ2"],
                "LIS_ipsae_script": r["LIS_ipsae_script"],
                "ipSAE_nres1_max": r["nres1"], "ipSAE_nres2_max": r["nres2"],
                "ipSAE_dist1_max": r["dist1"], "ipSAE_dist2_max": r["dist2"],
            })
        elif r["c1"] == nlr_chain and r["c2"] == eff_chain:
            out["ipSAE_official_nlr_to_effector"] = r["ipSAE"]
            out["ipSAE_official_d0dom_nlr_to_effector"] = r["ipSAE_d0dom"]
        elif r["c1"] == eff_chain and r["c2"] == nlr_chain:
            out["ipSAE_official_effector_to_nlr"] = r["ipSAE"]
            out["ipSAE_official_d0dom_effector_to_nlr"] = r["ipSAE_d0dom"]
    return out


def run_official_ipsae(ipsae_script: Path, full_path: Path, cif_path: Path, summary_path: Path,
                       nlr_chain: str, eff_chain: str) -> Dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="ipsae_run_") as td_s:
        td = Path(td_s)
        full_tmp = td / full_path.name
        cif_tmp = td / cif_path.name
        summary_tmp = td / summary_path.name
        hardlink_or_copy(full_path, full_tmp)
        hardlink_or_copy(cif_path, cif_tmp)
        hardlink_or_copy(summary_path, summary_tmp)
        cmd = [sys.executable, str(ipsae_script), str(full_tmp), str(cif_tmp),
               str(IPSAE_PAE_CUTOFF), str(IPSAE_DIST_CUTOFF)]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              encoding="utf-8", errors="replace", env=child_env())  # [v5]
        if proc.returncode != 0:
            eprint(f"WARNING: ipSAE failed for {cif_path.name}: "
                   f"{(proc.stderr or '')[-1000:]} {(proc.stdout or '')[-1000:]}")
            return {}
        ps = str(int(IPSAE_PAE_CUTOFF)).zfill(2) if IPSAE_PAE_CUTOFF < 10 else str(int(IPSAE_PAE_CUTOFF))
        ds = str(int(IPSAE_DIST_CUTOFF)).zfill(2) if IPSAE_DIST_CUTOFF < 10 else str(int(IPSAE_DIST_CUTOFF))
        out_txt = Path(str(cif_tmp.with_suffix("")) + f"_{ps}_{ds}.txt")
        if not out_txt.exists():
            # [v2] be forgiving about the exact naming convention
            cands = sorted(td.glob("*.txt"))
            cands = [c for c in cands if "byres" not in c.name.lower()]
            if len(cands) == 1:
                out_txt = cands[0]
            else:
                eprint(f"WARNING: ipSAE output missing for {cif_path.name} "
                       f"(looked for {out_txt.name}; found {[c.name for c in cands]})")
                return {}
        return parse_ipsae_txt(out_txt, nlr_chain, eff_chain)


def load_biopdb_structure(cif_path: Path):
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(cif_path.stem, str(cif_path))


def protein_chain_ids(structure) -> List[str]:
    ids: List[str] = []
    model = next(structure.get_models())
    for chain in model:
        if any(is_aa(res, standard=False) for res in chain.get_residues()):
            ids.append(chain.id)
    return ids


def make_subset_structure(structure, chain_ids: Sequence[str]):
    new_s = Structure("subset")
    new_m = Model(0)
    new_s.add(new_m)
    model = next(structure.get_models())
    wanted = set(chain_ids)
    for chain in model:
        if chain.id in wanted:
            new_m.add(copy.deepcopy(chain))
    return new_s


# ------------------------- [v3] benchmark-identical SASA ---------------------

def _fibonacci_sphere(n: int = SASA_SPHERE_POINTS) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


_SASA_SPHERE = _fibonacci_sphere()


def _sasa_per_atom(coords: np.ndarray, radii: np.ndarray,
                   which: Optional[Sequence[int]] = None,
                   probe: float = SASA_PROBE) -> np.ndarray:
    """
    Shrake-Rupley per-atom SASA, numerically identical to the loop in
    Rx_benchmarking_table.py but vectorised over neighbours.

    `which` restricts the computation to a subset of atom indices (their
    neighbour environment is still the full `coords` array). Returns an array
    the same length as `which` (or as `coords` when which is None).
    """
    r = radii + probe
    tree = cKDTree(coords)
    rmax = float(r.max())
    idx = range(len(coords)) if which is None else list(which)
    out = np.zeros(len(idx))
    npts = len(_SASA_SPHERE)
    for k, i in enumerate(idx):
        neigh = [j for j in tree.query_ball_point(coords[i], r[i] + rmax) if j != i]
        pts = coords[i] + _SASA_SPHERE * r[i]
        if neigh:
            nc = coords[neigh]
            nr2 = r[neigh] ** 2
            d2 = ((pts[:, None, :] - nc[None, :, :]) ** 2).sum(-1)
            acc = np.all(d2 > nr2[None, :], axis=1)
        else:
            acc = np.ones(npts, dtype=bool)
        out[k] = 4 * np.pi * r[i] ** 2 * acc.sum() / npts
    return out


def structure_heavy_atoms(structure) -> Tuple[np.ndarray, np.ndarray]:
    """[v3] Heavy-atom coordinates + Bondi radii from a Bio.PDB structure."""
    coords: List[Any] = []
    radii: List[float] = []
    for atom in structure.get_atoms():
        el = (getattr(atom, "element", "") or "").upper()
        if el in ("H", "D"):
            continue
        coords.append(atom.coord)
        radii.append(BONDI_VDW.get(el, SASA_DEFAULT_RADIUS))
    if not coords:
        return np.empty((0, 3)), np.empty((0,))
    return np.asarray(coords, dtype=float), np.asarray(radii, dtype=float)


def buried_surface_area_benchmark(struct_a, struct_b) -> Optional[float]:
    """
    [v3] TOTAL buried surface area (SASA_A + SASA_B - SASA_AB), using exactly
    the Rx_benchmarking_table.py definition. Divide by 2 for the per-side
    number that the benchmark table reports.

    Only atoms whose accessibility can actually change are recomputed: an atom
    i of chain A can only be occluded by an atom j of chain B if
    |c_i - c_j| <= (r_i + probe) + (r_j + probe). Atoms failing that test have
    identical SASA in the isolated and complexed states and contribute exactly
    zero, so restricting to them is exact, not an approximation.
    """
    xa, ra = structure_heavy_atoms(struct_a)
    xb, rb = structure_heavy_atoms(struct_b)
    if xa.size == 0 or xb.size == 0:
        return None

    Ra, Rb = ra + SASA_PROBE, rb + SASA_PROBE
    tree_b, tree_a = cKDTree(xb), cKDTree(xa)
    rb_max, ra_max = float(Rb.max()), float(Ra.max())

    aff_a = [i for i in range(len(xa))
             if any(np.linalg.norm(xa[i] - xb[j]) <= Ra[i] + Rb[j]
                    for j in tree_b.query_ball_point(xa[i], Ra[i] + rb_max))]
    aff_b = [j for j in range(len(xb))
             if any(np.linalg.norm(xb[j] - xa[i]) <= Rb[j] + Ra[i]
                    for i in tree_a.query_ball_point(xb[j], Rb[j] + ra_max))]
    if not aff_a and not aff_b:
        return 0.0

    both_x = np.vstack([xa, xb])
    both_r = np.concatenate([ra, rb])
    na = len(xa)

    iso_a = _sasa_per_atom(xa, ra, aff_a).sum() if aff_a else 0.0
    iso_b = _sasa_per_atom(xb, rb, aff_b).sum() if aff_b else 0.0
    cpx_a = _sasa_per_atom(both_x, both_r, aff_a).sum() if aff_a else 0.0
    cpx_b = _sasa_per_atom(both_x, both_r, [na + j for j in aff_b]).sum() if aff_b else 0.0

    return float(max(0.0, (iso_a - cpx_a) + (iso_b - cpx_b)))


def sasa_area(structure, backend: str = "benchmark") -> Tuple[Optional[float], str]:
    """Total SASA of a structure. Only used by the freesasa/biopdb backends."""
    if backend == "freesasa":
        try:
            import freesasa  # type: ignore
            result, _classes = freesasa.calcBioPDB(structure)
            return float(result.totalArea()), "FreeSASA"
        except Exception:
            return None, "unavailable"
    if backend == "biopdb":
        try:
            from Bio.PDB.SASA import ShrakeRupley
            st = copy.deepcopy(structure)
            sr = ShrakeRupley()
            sr.compute(st, level="S")
            return float(st.sasa), "Bio.PDB.ShrakeRupley"
        except Exception:
            return None, "unavailable"
    x, r = structure_heavy_atoms(structure)
    if x.size == 0:
        return None, "unavailable"
    return float(_sasa_per_atom(x, r).sum()), "benchmark_ShrakeRupley_256pt_Bondi"


def geometry_metrics(cif_path: Path, requested_nlr: str, requested_eff: str,
                     sasa_backend: str = "benchmark") -> Tuple[Dict[str, Any], str, str]:
    """Return clash + BSA; also resolved structure chain IDs."""
    try:
        structure = load_biopdb_structure(cif_path)
    except Exception as exc:
        eprint(f"WARNING: could not parse mmCIF {cif_path.name}: {exc}")
        return {}, requested_nlr, requested_eff

    model = next(structure.get_models())
    struct_ids = [c.id for c in model]
    nlr_chain, eff_chain = requested_nlr, requested_eff
    if nlr_chain not in struct_ids or eff_chain not in struct_ids:
        pids = protein_chain_ids(structure)
        if len(pids) >= 2:
            eprint(f"WARNING: chains {requested_nlr}/{requested_eff} not in {cif_path.name}; "
                   f"falling back to first two protein chains {pids[0]}/{pids[1]}")
            nlr_chain, eff_chain = pids[0], pids[1]

    def coords_for(chain_id: str) -> np.ndarray:
        try:
            ch = model[chain_id]
        except KeyError:
            return np.empty((0, 3))
        coords = []
        for atom in ch.get_atoms():
            element = (getattr(atom, "element", "") or "").upper()
            if element == "H" or element == "D":
                continue
            coords.append(atom.coord)
        return np.asarray(coords, dtype=float) if coords else np.empty((0, 3))

    ca = coords_for(nlr_chain)
    cb = coords_for(eff_chain)
    clash_count: Optional[int] = None
    min_dist: Optional[float] = None
    if ca.size and cb.size:
        ta, tb = cKDTree(ca), cKDTree(cb)
        neigh = ta.query_ball_tree(tb, r=HEAVY_ATOM_CLASH_CUTOFF)
        clash_count = int(sum(len(x) for x in neigh))
        d, _ = ta.query(cb, k=1)
        min_dist = float(np.min(d)) if len(d) else None

    a = make_subset_structure(structure, [nlr_chain])
    b = make_subset_structure(structure, [eff_chain])
    sasa_a = sasa_b = sasa_ab = None
    bsa_total = bsa_half = bsa_frac = None

    if sasa_backend == "benchmark":
        # [v3] Same definition as Rx_benchmarking_table.py.
        bsa_total = buried_surface_area_benchmark(a, b)
        method = "benchmark_ShrakeRupley_256pt_Bondi"
        if bsa_total is None:
            method = "unavailable"
            eprint(f"WARNING: BSA unavailable for {cif_path.name} (empty chain)")
    else:
        ab = make_subset_structure(structure, [nlr_chain, eff_chain])
        sasa_a, method_a = sasa_area(a, sasa_backend)
        sasa_b, method_b = sasa_area(b, sasa_backend)
        sasa_ab, method_ab = sasa_area(ab, sasa_backend)
        # [v2] The three SASA calls can each fail independently. Subtracting a
        # ShrakeRupley area from a FreeSASA area is meaningless (different
        # radii and probe defaults), so refuse rather than emit a
        # plausible-looking BSA.
        methods = {method_a, method_b, method_ab}
        if len(methods) == 1 and "unavailable" not in methods:
            method = method_ab
            if sasa_a is not None and sasa_b is not None and sasa_ab is not None:
                bsa_total = max(0.0, sasa_a + sasa_b - sasa_ab)
        else:
            method = "MIXED_OR_UNAVAILABLE:" + "/".join(sorted(methods))
            eprint(f"WARNING: inconsistent SASA backends for {cif_path.name} ({method}); "
                   f"BSA set to blank for this model")

    if bsa_total is not None:
        bsa_half = bsa_total / 2.0
        if sasa_a is not None and sasa_b is not None:
            denom = sasa_a + sasa_b
            bsa_frac = bsa_total / denom if denom > 0 else None

    # [v4] geometric interface residue sets, author numbering
    iface_n, iface_e, n_pairs = interface_residue_sets(model, nlr_chain, eff_chain)

    return {
        "iface_indices_nlr": indices_to_str(iface_n),
        "iface_indices_effector": indices_to_str(iface_e),
        "iface_res_nlr_n": len(iface_n) or None,
        "iface_res_effector_n": len(iface_e) or None,
        "geometric_contact_pairs_lt4p5A": n_pairs,
        "AB_heavy_atom_clash_pairs_lt1p5A": clash_count,
        "AB_min_heavy_atom_distance_A": min_dist,
        "SASA_nlr_A2": sasa_a,
        "SASA_effector_A2": sasa_b,
        "SASA_AB_A2": sasa_ab,
        "BSA_total_A2": bsa_total,
        "BSA_interface_area_A2": bsa_half,
        "BSA_fraction_combined": bsa_frac,
        "BSA_method": method,
    }, nlr_chain, eff_chain


# --------------- [v4] geometric interface residues (no lis.py) ---------------

GEOM_CONTACT_CUTOFF = 4.5   # same as Rx_benchmarking_table.py


def interface_residue_sets(model, chain_a: str, chain_b: str,
                           cutoff: float = GEOM_CONTACT_CUTOFF):
    """
    [v4] Interface residues straight from the mmCIF: any residue with a
    heavy atom within `cutoff` of the partner chain.

    Why this exists: vLRR_fraction and interface consistency used to depend on
    AFM-LIS's cLIR_indices columns, which (a) vanish silently if upstream
    renames a column and (b) have an unspecified index base, so intersecting
    them with author residue ranges was never safe. Bio.PDB gives author
    numbering directly -- the same numbering your vLRR ranges are written in --
    so the ambiguity disappears.

    Also returns the residue-pair count, which is exactly the
    "Geometric contacts < 4.5 A" column of the Rx benchmark table.
    """
    def heavy(chain_id: str):
        out = []
        try:
            ch = model[chain_id]
        except KeyError:
            return out
        for res in ch:
            for atom in res:
                el = (getattr(atom, "element", "") or "").upper()
                if el in ("H", "D"):
                    continue
                out.append((res.id[1], atom.coord))
        return out

    ha, hb = heavy(chain_a), heavy(chain_b)
    if not ha or not hb:
        return set(), set(), None
    xb = np.asarray([c for _, c in hb], dtype=float)
    tree = cKDTree(xb)
    sa: Set[int] = set()
    sb: Set[int] = set()
    pairs: Set[Tuple[int, int]] = set()
    for rn, xyz in ha:
        for j in tree.query_ball_point(xyz, cutoff):
            rb = hb[j][0]
            sa.add(int(rn))
            sb.add(int(rb))
            pairs.add((int(rn), int(rb)))
    return sa, sb, len(pairs)


def indices_to_str(s: Set[int]) -> str:
    return ",".join(str(x) for x in sorted(s))


def parse_index_set(value: Any) -> Set[int]:
    if value is None:
        return set()
    s = str(value).strip().strip('"').strip()
    if not s:
        return set()
    s = s.strip("[]")
    if not s:
        return set()
    out: Set[int] = set()
    # [v2] also tolerate space- and semicolon-separated lists
    for part in re.split(r"[,;\s]+", s):
        part = part.strip()
        if not part:
            continue
        # strip a possible "A:123" chain prefix
        if ":" in part:
            part = part.split(":", 1)[1]
        if "-" in part.lstrip("-"):
            a, _, b = part.partition("-")
            try:
                aa, bb = int(a), int(b)
                out.update(range(min(aa, bb), max(aa, bb) + 1))
            except ValueError:
                pass
        else:
            try:
                out.add(int(part))
            except ValueError:
                pass
    return out


def mean_pairwise_jaccard(sets: List[Set[int]]) -> Optional[float]:
    """
    [v2] Returns None (not 0.0) when there is not enough data.

    The original returned 0.0 whenever fewer than two non-empty sets were
    available. If the upstream lis.py column were renamed, every candidate
    would silently get interface_consistency = 0.0, which reads as
    'all seeds disagree' rather than 'this was never computed'.
    """
    nonempty = [s for s in sets if s]
    if len(nonempty) < 2:
        return None
    vals: List[float] = []
    for i in range(len(nonempty)):
        for j in range(i + 1, len(nonempty)):
            u = nonempty[i] | nonempty[j]
            vals.append(len(nonempty[i] & nonempty[j]) / len(u) if u else 0.0)
    return float(np.mean(vals)) if vals else None


def vlrr_bin_label(frac: Optional[float], lo: float, hi: float) -> str:
    """[v2] Continuous vLRR fraction -> the same Yes/Partially/No vocabulary."""
    if frac is None:
        return ""
    if frac >= hi:
        return "Yes"
    if frac >= lo:
        return "Partially"
    return "No"


VLRR_ORDER = {"Yes": 0, "Partially": 1, "No": 2}


def summarize_rows(raw_rows: List[Dict[str, Any]],
                   vlrr_ranges: Dict[str, Set[int]],
                   vlrr_lo: float, vlrr_hi: float,
                   interface_source: str = "geometric",
                lis_extra_args: Optional[Sequence[str]] = None) -> List[Dict[str, Any]]:
    groups: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for r in raw_rows:
        groups[(r["backbone"], r["candidate_id"], r["effector"])].append(r)

    summary: List[Dict[str, Any]] = []
    for (backbone, candidate, effector), rows in groups.items():
        first = rows[0]
        s: Dict[str, Any] = {
            "backbone": backbone,
            "candidate_id": candidate,
            "effector": effector,
            "vLRR_Interface": first.get("vLRR_Interface", ""),
            "n_models": len(rows),
        }
        for col in SUMMARY_NUMERIC_COLUMNS:
            s[f"{col}_mean"] = safe_mean(r.get(col) for r in rows)

        # [v9] worst single model, for the magnitude escape hatch below
        _clash = [as_float(r.get("AB_heavy_atom_clash_pairs_lt1p5A")) for r in rows]
        _clash = [x for x in _clash if x is not None]
        s["AB_heavy_atom_clash_pairs_lt1p5A_max"] = max(_clash) if _clash else None

        # [v2] medians + SDs for everything the ranking actually uses
        for col in MEDIAN_SD_COLUMNS:
            s[f"{col}_median"] = safe_median(r.get(col) for r in rows)
            s[f"{col}_sd"] = safe_sd(r.get(col) for r in rows)

        # [v4] consistency uses the same interface definition as vLRR_fraction
        nkey = "iface_indices_nlr" if interface_source == "geometric" else "cLIR_indices_nlr"
        ekey = "iface_indices_effector" if interface_source == "geometric" else "cLIR_indices_effector"
        n_sets = [parse_index_set(r.get(nkey)) for r in rows]
        e_sets = [parse_index_set(r.get(ekey)) for r in rows]
        if not any(n_sets):   # selected source empty -> try the other one
            n_sets = [parse_index_set(r.get("cLIR_indices_nlr" if interface_source == "geometric"
                                            else "iface_indices_nlr")) for r in rows]
            e_sets = [parse_index_set(r.get("cLIR_indices_effector" if interface_source == "geometric"
                                            else "iface_indices_effector")) for r in rows]
        jac_n = mean_pairwise_jaccard(n_sets)
        jac_e = mean_pairwise_jaccard(e_sets)
        s["cLIR_consistency_nlr_jaccard"] = jac_n
        s["cLIR_consistency_effector_jaccard"] = jac_e
        # [v2] None-safe average
        parts = [x for x in (jac_n, jac_e) if x is not None]
        s["interface_consistency_mean"] = float(np.mean(parts)) if parts else None

        s["models_with_cLIR"] = sum(1 for a, b in zip(n_sets, e_sets) if a and b)
        # [v10] Jaccard below is calculated only among non-empty interface sets.
        # Record how many of the expected models actually contain an interface,
        # so sparse 2/5 support cannot masquerade as robust 5/5 consistency.
        s["interface_model_support_fraction"] = (
            s["models_with_cLIR"] / len(rows) if rows else None
        )
        s["models_with_iLIS_gt0"] = sum(1 for r in rows if (as_float(r.get("iLIS")) or 0) > 0)
        s["af3_has_clash_n"] = sum(1 for r in rows if bool(r.get("af3_has_clash")))
        s["n_models_with_clash"] = sum(                                   # [v9]
            1 for r in rows
            if (as_float(r.get("AB_heavy_atom_clash_pairs_lt1p5A")) or 0) > 0)

        # [v2] PAE and contact columns are now part of the critical check
        def missing(r: Dict[str, Any]) -> bool:
            return (r.get("iLIS") is None
                    or r.get("ipSAE_official_max") is None
                    or r.get("nlr_pLDDT") is None
                    or r.get("effector_pLDDT") is None
                    or r.get("PAE_nlr_to_effector_mean") is None
                    or r.get("contact_prob_max") is None)
        s["critical_metric_missing_n"] = sum(1 for r in rows if missing(r))

        # [v2] computed vLRR class + agreement with the manual annotation
        s["vLRR_class_computed"] = vlrr_bin_label(s.get("vLRR_fraction_median"), vlrr_lo, vlrr_hi)
        manual = str(s.get("vLRR_Interface", "")).strip().title()
        comp = s["vLRR_class_computed"]
        if manual and comp:
            s["vLRR_class_agreement"] = "AGREE" if manual == comp else f"DIFFER({manual}->{comp})"
        else:
            s["vLRR_class_agreement"] = ""

        s["ranking_mode"] = "computed_vLRR_fraction" if backbone in vlrr_ranges else "manual_vLRR_class"
        # [v2] flag the case where the computed key was requested but unusable
        s["vLRR_fallback_to_manual"] = bool(
            backbone in vlrr_ranges and not s["vLRR_class_computed"]
        )
        s["annotation_warning"] = first.get("annotation_warning", "")
        s["ipsae_header_parsed"] = bool(first.get("ipsae_header_parsed", False))
        s["BSA_method"] = first.get("BSA_method", "")
        summary.append(s)

    def desc(v: Any) -> float:
        x = as_float(v)
        return -x if x is not None else float("inf")

    def sort_key(s: Dict[str, Any]):
        """
        [v2] Ranking key.
          1. vLRR localisation  -- computed fraction (binned) when available,
             otherwise the manual annotation. Binned rather than raw so that a
             continuous value does not lexicographically dominate everything.
          2. iLIS median        -- PAE channel
          3. contact pairs>0.3  -- distogram channel, INDEPENDENT of PAE
          4. ipSAE median       -- PAE channel, confirmatory
          5. fraction of models that reproduce a non-empty interface
          6. interface-residue Jaccard among interface-bearing models
        Medians throughout: one catastrophic seed should not move a candidate.
        """
        manual = VLRR_ORDER.get(str(s.get("vLRR_Interface", "")).strip().title(), 3)
        if s.get("ranking_mode") == "computed_vLRR_fraction":
            comp = s.get("vLRR_class_computed", "")
            # [v2] If --vlrr-ranges was given but the fraction could not be
            # computed for this candidate (missing cLIR indices, or residue
            # numbering that does not overlap the supplied range), fall back
            # to the manual class for THIS row rather than letting every
            # candidate collapse to the same primary key -- which would
            # silently turn the ranking into "iLIS only".
            primary = VLRR_ORDER.get(comp, manual) if comp else manual
        else:
            primary = manual
        return (
            primary,
            desc(s.get("iLIS_median")),
            desc(s.get("n_contact_pairs_gt0.3_median")),
            desc(s.get("ipSAE_official_max_median")),
            desc(s.get("interface_model_support_fraction")),
            desc(s.get("interface_consistency_mean")),
        )

    # rank within backbone
    by_backbone: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for s in summary:
        by_backbone[s["backbone"]].append(s)
    for _backbone, items in by_backbone.items():
        items.sort(key=sort_key)
        for rank, s in enumerate(items, 1):
            s["suggested_rank_within_backbone"] = rank

    # [v2] rank within (backbone, effector). Comparing absolute iLIS across
    # different effectors is confounded by effector size/sequence; the
    # within-effector rank is the defensible one for a 5-effector screen.
    by_be: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for s in summary:
        by_be[(s["backbone"], norm_effector(s["effector"]))].append(s)
    for (_bb, _eff), items in by_be.items():
        items.sort(key=sort_key)
        n = len(items)
        for rank, s in enumerate(items, 1):
            s["rank_within_effector"] = rank
            s["n_in_effector_group"] = n

    summary.sort(key=lambda s: (s["backbone"], s.get("suggested_rank_within_backbone", 9999)))
    return summary


def make_qc_table(summary_rows: List[Dict[str, Any]], expected_models: int,
                  nlr_plddt_warn: Optional[float], effector_plddt_warn: Optional[float]) -> List[Dict[str, Any]]:
    qc: List[Dict[str, Any]] = []
    for s in summary_rows:
        reasons: List[str] = []
        n_models = int(s.get("n_models", 0))
        if n_models != expected_models:
            reasons.append(f"MODEL_COUNT_{n_models}_OF_{expected_models}")
        if int(s.get("critical_metric_missing_n", 0) or 0) > 0:
            reasons.append("MISSING_CRITICAL_METRIC")
        if int(s.get("af3_has_clash_n", 0) or 0) > 0:
            reasons.append("AF3_HAS_CLASH")
        # [v9] Clash QC: reproducibility, not single-occurrence.
        #
        # The old rule flagged a candidate if ANY of its models contained even
        # one heavy-atom pair below 1.5 A. On this dataset that flagged 33/61,
        # 25 of which had a clash in only 1-2 of 5 models -- i.e. sampling
        # noise, not a property of the model. A clash present in the MEDIAN
        # model (equivalently >=3 of 5) is a stable feature worth inspecting.
        #
        # Note the 1.5 A distance itself is not strict: MolProbity's serious
        # clash is ~0.4 A of vdW overlap (~3.0 A for two carbons) and even
        # short hydrogen bonds are ~2.5 A. AF3's native has_clash is far more
        # permissive still (>100 clashing atoms in a chain, or >50% of a
        # chain) because it is a whole-prediction failure detector, not an
        # interface geometry check -- it is blind to a four-atom overlap at
        # an interface, so it cannot serve as the reference standard here.
        clash_mean = as_float(s.get("AB_heavy_atom_clash_pairs_lt1p5A_mean"))
        clash_median = as_float(s.get("AB_heavy_atom_clash_pairs_lt1p5A_median"))
        clash_worst = as_float(s.get("AB_heavy_atom_clash_pairs_lt1p5A_max"))
        if clash_median is not None and clash_median > 0:
            reasons.append("AB_CLASH_SEVERE" if (clash_mean or 0) >= 3.0
                           else "AB_CLASH_REPEATED")
        elif clash_worst is not None and clash_worst >= 10:
            # Inert on the current data (worst single model among the
            # 1-2/5 candidates is 4 pairs); kept so a future run cannot hide
            # a catastrophic single model behind a zero median.
            reasons.append("AB_CLASH_SINGLE_MODEL_SEVERE")
        nlr_plddt = as_float(s.get("nlr_pLDDT_mean"))
        eff_plddt = as_float(s.get("effector_pLDDT_mean"))
        if nlr_plddt_warn is not None and nlr_plddt is not None and nlr_plddt < nlr_plddt_warn:
            reasons.append("LOW_NLR_pLDDT")
        if effector_plddt_warn is not None and eff_plddt is not None and eff_plddt < effector_plddt_warn:
            reasons.append("LOW_EFFECTOR_pLDDT")
        if not s.get("vLRR_Interface") and not s.get("vLRR_class_computed"):
            reasons.append("MISSING_vLRR_CLASS")
        # [v2] surface things that used to be invisible
        if s.get("annotation_warning"):
            reasons.append(str(s["annotation_warning"]))
        if s.get("interface_consistency_mean") is None:
            reasons.append("NO_INTERFACE_CONSISTENCY")
        if not s.get("ipsae_header_parsed", False):
            reasons.append("IPSAE_POSITIONAL_FALLBACK")
        bsa_method = str(s.get("BSA_method", ""))
        if bsa_method.startswith("MIXED_OR_UNAVAILABLE"):
            reasons.append("BSA_BACKEND_INCONSISTENT")
        if s.get("vLRR_class_agreement", "").startswith("DIFFER"):
            reasons.append(str(s["vLRR_class_agreement"]))
        if s.get("vLRR_fallback_to_manual"):
            reasons.append("vLRR_FRACTION_UNCOMPUTABLE_USED_MANUAL")

        qc.append({
            "backbone": s["backbone"],
            "candidate_id": s["candidate_id"],
            "effector": s["effector"],
            "vLRR_Interface_manual": s.get("vLRR_Interface", ""),
            "vLRR_class_computed": s.get("vLRR_class_computed", ""),
            "vLRR_fraction_median": s.get("vLRR_fraction_median"),
            "vLRR_class_agreement": s.get("vLRR_class_agreement", ""),
            "n_models_found": n_models,
            "n_models_expected": expected_models,
            "QC_status": "CHECK" if reasons else "PASS",
            "QC_notes": "; ".join(reasons),
            "nlr_pLDDT_mean": s.get("nlr_pLDDT_mean"),
            "effector_pLDDT_mean": s.get("effector_pLDDT_mean"),
            "nlr_intra_PAE_mean": s.get("nlr_intra_PAE_mean_mean"),
            "effector_intra_PAE_mean": s.get("effector_intra_PAE_mean_mean"),
            "af3_has_clash_n": s.get("af3_has_clash_n"),
            "AB_heavy_atom_clash_pairs_mean": s.get("AB_heavy_atom_clash_pairs_lt1p5A_mean"),
            "AB_heavy_atom_clash_pairs_median": s.get("AB_heavy_atom_clash_pairs_lt1p5A_median"),
            "AB_heavy_atom_clash_pairs_worst_model": s.get("AB_heavy_atom_clash_pairs_lt1p5A_max"),
            "n_models_with_clash": s.get("n_models_with_clash"),
            "AB_min_heavy_atom_distance_mean_A": s.get("AB_min_heavy_atom_distance_A_mean"),
            "models_with_cLIR": s.get("models_with_cLIR"),
            "interface_model_support_fraction": s.get("interface_model_support_fraction"),
            "interface_consistency_mean": s.get("interface_consistency_mean"),
            "iLIS_median": s.get("iLIS_median"),
            "iLIS_sd": s.get("iLIS_sd"),
            "contact_prob_max_median": s.get("contact_prob_max_median"),
            "n_contact_pairs_gt0.3_median": s.get("n_contact_pairs_gt0.3_median"),
            "ipSAE_official_max_median": s.get("ipSAE_official_max_median"),
            "n_tokens_nlr_mean": s.get("n_tokens_nlr_mean"),
            "n_tokens_effector_mean": s.get("n_tokens_effector_mean"),
            "n_tokens_atp_mean": s.get("n_tokens_atp_mean"),
            "BSA_method": s.get("BSA_method", ""),
        })
    return qc


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: Optional[List[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = []
        seen: Set[str] = set()
        for r in rows:
            for k in r.keys():
                if k not in seen:
                    seen.add(k)
                    fieldnames.append(k)
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_xlsx(path: Path, sheets: Sequence[Tuple[str, List[Dict[str, Any]]]]) -> bool:
    try:
        import xlsxwriter
    except ImportError:
        eprint("xlsxwriter not installed; CSV outputs were still created.")
        return False
    wb = xlsxwriter.Workbook(str(path))
    header = wb.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1, "text_wrap": True, "valign": "top"})
    numfmt = wb.add_format({"num_format": "0.0000"})
    textfmt = wb.add_format({"valign": "top"})
    for sheet_name, rows in sheets:
        ws = wb.add_worksheet(sheet_name[:31])
        if not rows:
            ws.write(0, 0, "No data")
            continue
        fields: List[str] = []
        seen: Set[str] = set()
        for r in rows:
            for k in r.keys():
                if k not in seen:
                    seen.add(k)
                    fields.append(k)
        for c, f in enumerate(fields):
            ws.write(0, c, f, header)
        for rr, row in enumerate(rows, start=1):
            for cc, f in enumerate(fields):
                v = row.get(f)
                if isinstance(v, (int, float)) and not isinstance(v, bool) and v is not None:
                    ws.write_number(rr, cc, float(v), numfmt)
                elif isinstance(v, bool):
                    ws.write(rr, cc, str(v))
                elif v is not None:
                    ws.write(rr, cc, str(v), textfmt)
        ws.freeze_panes(1, 0)
        ws.autofilter(0, 0, len(rows), len(fields) - 1)
        for c, f in enumerate(fields):
            width = min(36, max(10, len(f) + 2))
            if "indices" in f or "file" in f or "notes" in f:
                width = 30
            ws.set_column(c, c, width)
    wb.close()
    return True


def process_job(job: Dict[str, Any], lookup: Dict[Tuple[str, str], Dict[str, str]],
                lis_script: Path, ipsae_script: Path, tool_out: Path, lis_workers: int,
                vlrr_ranges: Dict[str, Set[int]], sasa_backend: str = "benchmark",
                interface_source: str = "geometric",
                lis_extra_args: Optional[Sequence[str]] = None) -> List[Dict[str, Any]]:
    job_dir: Path = job["job_dir"]
    backbone = job["backbone"]
    candidate = job["candidate_id"]
    effector = job["effector"]

    ann = lookup.get((backbone, candidate), {})
    vlrr_class = ann.get("vLRR_Interface", "")
    expected_eff = ann.get("effector", "")
    annotation_warning = ""
    if expected_eff and norm_effector(expected_eff) != norm_effector(effector):
        annotation_warning = f"LOOKUP_EFFECTOR_MISMATCH:{expected_eff}"

    models = find_model_files(job_dir)
    if not models:
        eprint(f"WARNING: no complete AF3 models found in {job_dir}")
        return []

    lis_csv = tool_out / "lis" / f"{backbone}_{candidate}_{norm_effector(effector)}.csv"
    lis_rows = run_lis_for_job(lis_script, job_dir, lis_csv, lis_workers, lis_extra_args)
    if not lis_rows:
        raise RuntimeError(
            f"Official AFM-LIS produced no usable rows for {job_dir.name}. "
            f"Stopping rather than ranking with missing iLIS."
        )

    vlrr_set = vlrr_ranges.get(backbone)

    raw: List[Dict[str, Any]] = []
    for idx, model in enumerate(models, start=1):
        cif, full_path, summary_path = model["cif"], model["full"], model["summary"]
        label = f"{job['folder_name']}/{model['model_id']}"
        try:
            full = read_json(full_path)
            summary = read_json(summary_path)
        except Exception as exc:
            eprint(f"WARNING: failed to read JSON for {cif}: {exc}")
            continue

        chain_order = chain_order_from_json(summary, full)
        if len(chain_order) < 2:
            eprint(f"ERROR: {label}: fewer than two distinct chains resolved "
                   f"({chain_order}); skipping this model.")
            continue
        if len(set(chain_order)) != len(chain_order):
            eprint(f"ERROR: {label}: duplicate chain IDs in {chain_order}; skipping.")
            continue
        nlr_chain, eff_chain = chain_order[0], chain_order[1]
        atp_chain = chain_order[2] if len(chain_order) >= 3 else ""

        # [v2] chain-assignment sanity check. Assignment is positional, so a
        # job submitted with a different chain order would silently invert
        # NLR and effector. The NLR should be by far the largest chain.
        counts = chain_token_counts(full)
        n_tok_nlr = counts.get(nlr_chain)
        n_tok_eff = counts.get(eff_chain)
        n_tok_atp = counts.get(atp_chain) if atp_chain else None
        if idx == 1:   # [v5] print the chain inventory once per job
            inv = ", ".join(f"{c}={counts.get(c, 0)} tok" for c in chain_order)
            eprint(f"  chains: {inv}  ->  NLR={nlr_chain}, effector={eff_chain}"
                   + (f", ATP={atp_chain}" if atp_chain else ""))
        chain_warning = ""
        if n_tok_nlr is not None and n_tok_eff is not None and n_tok_nlr <= n_tok_eff:
            chain_warning = f"CHAIN_ORDER_SUSPECT:nlr={n_tok_nlr}<=eff={n_tok_eff}"
            eprint(f"WARNING: {label}: chain1 ({nlr_chain}, {n_tok_nlr} tokens) is not larger "
                   f"than chain2 ({eff_chain}, {n_tok_eff} tokens). "
                   f"Check that chain1=NLR really holds for this job.")

        geom, struct_nlr, struct_eff = geometry_metrics(cif, nlr_chain, eff_chain, sasa_backend)
        lis_row = choose_lis_row(lis_rows, cif, job_dir, nlr_chain, eff_chain)
        lis_metrics = parse_lis_row(lis_row, nlr_chain, eff_chain)
        if lis_metrics.get("iLIS") is None:
            raise RuntimeError(
                f"Official AFM-LIS iLIS could not be matched/validated for {label}. "
                f"Stopping rather than using an approximate replacement."
            )

        ipsae_metrics = run_official_ipsae(
            ipsae_script, full_path, cif, summary_path, nlr_chain, eff_chain
        )
        if (not ipsae_metrics.get("ipsae_header_parsed", False)
                or ipsae_metrics.get("ipSAE_official_max") is None):
            raise RuntimeError(
                f"Official ipSAE output could not be parsed/validated for {label}. "
                f"Stopping rather than using a positional or approximate replacement."
            )

        pae_metrics = compute_pae_qc(full, nlr_chain, eff_chain, label=label)
        contact_metrics = compute_contact_metrics(full, nlr_chain, eff_chain, label=label)  # [v2]

        matrix_order = matrix_chain_order_from_summary(summary, full)
        try:
            i_n = matrix_order.index(nlr_chain)
            i_e = matrix_order.index(eff_chain)
        except ValueError:
            raise RuntimeError(
                f"Could not map NLR/effector chains {nlr_chain}/{eff_chain} onto "
                f"AF3 chain-pair matrix order {matrix_order} for {label}."
            )

        row: Dict[str, Any] = {
            "backbone": backbone,
            "candidate_id": candidate,
            "effector": effector,
            "folder_name": job["folder_name"],
            "model_id": str(model["model_id"]),
            "model_ordinal": idx,
            "vLRR_Interface": vlrr_class,
            "annotation_warning": "; ".join(x for x in (annotation_warning, chain_warning) if x),
            "nlr_chain": nlr_chain,
            "effector_chain": eff_chain,
            "atp_chain": atp_chain,
            "n_tokens_nlr": n_tok_nlr,
            "n_tokens_effector": n_tok_eff,
            "n_tokens_atp": n_tok_atp,
            "structure_nlr_chain_resolved": struct_nlr,
            "structure_effector_chain_resolved": struct_eff,
            "structure_file": str(cif),
            "full_data_file": str(full_path),
            "summary_file": str(summary_path),
            "af3_global_iptm": as_float(summary.get("iptm")),
            "af3_global_ptm": as_float(summary.get("ptm")),
            "af3_ranking_score": as_float(summary.get("ranking_score")),
            "af3_has_clash": bool(summary.get("has_clash", False)),
            "pairwise_iptm_nlr_effector": matrix_value(summary, "chain_pair_iptm", i_n, i_e),
            "pairwise_pae_min_nlr_to_effector": matrix_value(summary, "chain_pair_pae_min", i_n, i_e),
            "pairwise_pae_min_effector_to_nlr": matrix_value(summary, "chain_pair_pae_min", i_e, i_n),
        }
        row.update(lis_metrics)
        row.update(ipsae_metrics)
        row.update(pae_metrics)
        row.update(contact_metrics)  # [v2]
        row.update(geom)

        # [v4] vLRR localisation from whichever interface definition is selected.
        clir_nlr = parse_index_set(row.get("cLIR_indices_nlr"))
        geom_nlr = parse_index_set(row.get("iface_indices_nlr"))
        row["cLIR_index_min"] = min(clir_nlr) if clir_nlr else None
        row["cLIR_index_max"] = max(clir_nlr) if clir_nlr else None
        row["iface_index_min"] = min(geom_nlr) if geom_nlr else None
        row["iface_index_max"] = max(geom_nlr) if geom_nlr else None

        if interface_source == "geometric":
            primary_set, used = geom_nlr, "geometric_4.5A"
            if not primary_set and clir_nlr:
                primary_set, used = clir_nlr, "cLIR_fallback"
        else:
            primary_set, used = clir_nlr, "cLIR"
            if not primary_set and geom_nlr:
                primary_set, used = geom_nlr, "geometric_fallback"
        row["interface_source_used"] = used if primary_set else ""

        if vlrr_set and primary_set:
            inside = len(primary_set & vlrr_set)
            row["vLRR_cLIR_n_in_range"] = inside
            row["vLRR_fraction"] = inside / len(primary_set)
        else:
            row["vLRR_cLIR_n_in_range"] = None
            row["vLRR_fraction"] = None

        raw.append(row)
        eprint(f"  {job['folder_name']}: {idx}/{len(models)} done")
    return raw


def report_clir_numbering(raw_rows: List[Dict[str, Any]]) -> List[str]:
    """
    [v2] Diagnostic: print the observed cLIR index range per backbone.

    CRITICAL SANITY CHECK -- read this before trusting vLRR_fraction.
    If lis.py emits 0-based token indices rather than 1-based author residue
    numbers, intersecting them with your vLRR residue range is meaningless.
    Compare the printed min/max against the real length of each NLR chain.
    """
    lines: List[str] = []
    by_bb: Dict[str, List[int]] = defaultdict(list)
    for r in raw_rows:
        lo, hi = r.get("cLIR_index_min"), r.get("cLIR_index_max")
        if lo is not None:
            by_bb[r["backbone"]].append(int(lo))
        if hi is not None:
            by_bb[r["backbone"]].append(int(hi))
    for bb, vals in sorted(by_bb.items()):
        if not vals:
            continue
        toks = [as_float(r.get("n_tokens_nlr")) for r in raw_rows if r["backbone"] == bb]
        toks = [t for t in toks if t is not None]
        ntok = int(max(toks)) if toks else None
        lines.append(f"  {bb}: cLIR indices span {min(vals)}..{max(vals)}"
                     + (f"  (NLR chain has {ntok} tokens)" if ntok else ""))
    return lines


def main() -> int:
    global IPSAE_PAE_CUTOFF  # [v3] --ipsae-pae-cutoff overrides the module default
    _ipsae_default = IPSAE_PAE_CUTOFF
    ap = argparse.ArgumentParser(
        description="Verified pipeline: author-official AFM-LIS + ipSAE, native AF3 fields, and labelled project-specific derived metrics.")
    ap.add_argument("--roots", nargs="+", required=True,
                    help="One or more directories containing rx*/sr35* candidate folders.")
    ap.add_argument("--interface-lookup", required=True,
                    help="CSV containing backbone,candidate_id,effector,vLRR_Interface.")
    ap.add_argument("--out", required=True, help="Output analysis directory.")
    ap.add_argument("--expected-models", type=int, default=EXPECTED_MODELS_DEFAULT)
    ap.add_argument("--lis-workers", type=int, default=1,
                    help="Workers passed to official lis.py for each candidate folder.")
    ap.add_argument("--refresh-tools", action="store_true")
    ap.add_argument("--nlr-plddt-warn", type=float, default=None,
                    help="Optional QC warning threshold; no pLDDT threshold is imposed by default.")
    ap.add_argument("--effector-plddt-warn", type=float, default=None,
                    help="Optional QC warning threshold; no pLDDT threshold is imposed by default.")
    # [v2] new options
    ap.add_argument("--vlrr-ranges", nargs="*", default=None, metavar="BACKBONE=RANGES",
                    help='Compute vLRR localisation objectively instead of using the manual '
                         'annotation, e.g. --vlrr-ranges "Rx=489-593" "Sr35=500-620,700-740". '
                         'Residue numbering must match whatever lis.py reports in cLIR_indices '
                         '-- check the cLIR numbering diagnostic printed at the end of the run.')
    ap.add_argument("--vlrr-bin-low", type=float, default=VLRR_BIN_LOW_DEFAULT,
                    help="vLRR fraction below this -> 'No' (default 0.20).")
    ap.add_argument("--vlrr-bin-high", type=float, default=VLRR_BIN_HIGH_DEFAULT,
                    help="vLRR fraction at or above this -> 'Yes' (default 0.60).")
    ap.add_argument("--lis-extra-args", nargs="*", default=None,
                    help="Extra arguments passed through to lis.py, e.g. "
                         "--lis-extra-args --platform af3")
    ap.add_argument("--interface-source", choices=["geometric", "clir"], default="geometric",
                    help="Which interface residue set feeds vLRR_fraction and interface "
                         "consistency. 'geometric' (default) reads heavy-atom contacts <4.5 A "
                         "straight from the mmCIF in author numbering -- no dependence on "
                         "AFM-LIS output columns and no index-base ambiguity. 'clir' uses "
                         "AFM-LIS cLIR_indices instead.")
    ap.add_argument("--sasa-backend", choices=["benchmark", "freesasa", "biopdb"],
                    default="benchmark",
                    help="BSA backend. 'benchmark' (default) reproduces "
                         "Rx_benchmarking_table.py exactly: Shrake-Rupley, 1.4 A probe, "
                         "256-point sphere, Bondi radii, heavy atoms only. The other two "
                         "are faster but NOT numerically comparable with the benchmark table.")
    ap.add_argument("--ipsae-pae-cutoff", type=float, default=_ipsae_default,
                    help="ipSAE PAE cutoff in A (default 10, matching the Rx benchmark table). "
                         "Your earlier ipSAE_15 analysis used 15 -- the two are not comparable.")
    ap.add_argument("--afm-lis-ref", default=AFM_LIS_REF_DEFAULT,
                    help="Git ref for lis.py. PIN THIS to a commit SHA for reproducibility.")
    ap.add_argument("--ipsae-ref", default=IPSAE_REF_DEFAULT,
                    help="Git ref for ipsae.py. PIN THIS to a commit SHA for reproducibility.")
    ap.add_argument("--expect-lis-sha256", default=None,
                    help="Abort if the downloaded lis.py does not match this SHA256.")
    ap.add_argument("--expect-ipsae-sha256", default=None,
                    help="Abort if the downloaded ipsae.py does not match this SHA256.")
    args = ap.parse_args()

    IPSAE_PAE_CUTOFF = float(args.ipsae_pae_cutoff)  # [v3]

    roots = [Path(x) for x in args.roots]
    lookup_path = Path(args.interface_lookup)
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    tool_dir = out / "_tools"
    tool_out = out / "_tool_outputs"

    if not lookup_path.exists():
        raise SystemExit(f"Lookup CSV not found: {lookup_path}")
    lookup = load_interface_lookup(lookup_path)

    vlrr_ranges = parse_vlrr_ranges(args.vlrr_ranges)
    if vlrr_ranges:
        for bb, rs in sorted(vlrr_ranges.items()):
            eprint(f"vLRR range for {bb}: {len(rs)} residues ({min(rs)}..{max(rs)})")
    else:
        eprint("NOTE: no --vlrr-ranges given; ranking will fall back to the MANUAL "
               "vLRR_Interface annotation as its primary key.")

    lis_url = AFM_LIS_URL_TMPL.format(ref=args.afm_lis_ref)
    ipsae_url = IPSAE_URL_TMPL.format(ref=args.ipsae_ref)
    if args.afm_lis_ref == "main" or args.ipsae_ref == "main":
        eprint("WARNING: analysis tools are pinned to 'main'. Re-running later may silently "
               "produce different numbers. Record the SHA256s below and pass them back via "
               "--expect-lis-sha256 / --expect-ipsae-sha256, or pin --afm-lis-ref/--ipsae-ref "
               "to a commit SHA.")

    lis_script = ensure_tool(lis_url, tool_dir / "lis.py", args.refresh_tools, args.expect_lis_sha256)
    ipsae_script = ensure_tool(ipsae_url, tool_dir / "ipsae.py", args.refresh_tools, args.expect_ipsae_sha256)

    provenance = [
        "=== AUTHOR/NATIVE METRICS ===",
        "AFM-LIS metrics are obtained by executing the authors' official lis.py; "
        "this wrapper does not reimplement iLIS.",
        f"AFM-LIS URL: {lis_url}",
        f"AFM-LIS ref: {args.afm_lis_ref}",
        f"AFM-LIS SHA256: {sha256_file(lis_script)}",
        f"IPSAE URL: {ipsae_url}",
        f"IPSAE ref: {args.ipsae_ref}",
        f"IPSAE SHA256: {sha256_file(ipsae_script)}",
        f"AFM-LIS parameters passed to official lis.py: PAE cutoff={LIS_PAE_CUTOFF} A, Cbeta cutoff={LIS_CB_CUTOFF} A. "
        "Current lis.py implementation applies strict < cutoffs internally.",
        f"ipSAE parameters: PAE cutoff={IPSAE_PAE_CUTOFF} A, distance argument={IPSAE_DIST_CUTOFF} A",
        f"contact probability thresholds: {list(CONTACT_PROB_THRESHOLDS)}",
        f"SASA backend: {args.sasa_backend}",
        f"Interface residue source: {args.interface_source} "
        f"(geometric = heavy-atom contacts < {GEOM_CONTACT_CUTOFF} A, author numbering)",
        f"AB heavy-atom clash count cutoff: <{HEAVY_ATOM_CLASH_CUTOFF} A",
        "Clash QC rule: flagged only when the MEDIAN model has >=1 clashing pair "
        "(equivalently >=3 of 5 models); AB_CLASH_SEVERE when the mean is >=3 "
        "pairs/model, otherwise AB_CLASH_REPEATED. A single model with >=10 pairs "
        "is flagged separately. Clash counts never enter the ranking.",
        "Pairwise analysis: chain1=NLR, chain2=effector; chain3=ATP excluded from local interface metrics.",
        "",
        "=== DERIVED / PROJECT-SPECIFIC ANALYSES (NOT OFFICIAL AF3/AFM-LIS/IPSAE METRICS) ===",
        "contact_prob_max/mean and counts above 0.3/0.5: summaries derived from AF3's native contact_probs matrix.",
        "intra-/inter-chain PAE mean/median/p90: summaries derived from AF3's native PAE matrix.",
        "heavy-atom clash-pair count (<1.5 A): project-specific geometry QC; AF3 has_clash is the native AF3 clash flag.",
        "geometric interface residues/contact pairs (<4.5 A): project-specific structural definition.",
        "BSA: project-specific Shrake-Rupley implementation chosen for consistency with the Rx benchmark.",
        "vLRR_fraction, Yes/Partially/No binning, Jaccard interface consistency and final ranking: project-specific analyses.",
        "",
        "BSA convention:",
        "  BSA_total_A2          = SASA(NLR) + SASA(effector) - SASA(NLR+effector)  [total buried]",
        "  BSA_interface_area_A2 = BSA_total_A2 / 2                                  [per-side]",
        "  Rx_benchmarking_table.py reports the PER-SIDE number, i.e. BSA_interface_area_A2.",
        "  Quote that column, not BSA_total_A2, or the two tables are not comparable.",
        "",
        "Ranking key (all medians across models):",
        "  1. vLRR localisation (computed fraction binned, else manual class)",
        "  2. iLIS median                 [PAE channel]",
        "  3. n_contact_pairs_gt0.3       [distogram channel, independent of PAE]",
        "  4. ipSAE median                [PAE channel, confirmatory]",
        "  5. interface model-support fraction (non-empty interface models / all models)",
        "  6. interface consistency (Jaccard among interface-bearing models)",
        f"vLRR bins: >= {args.vlrr_bin_high} -> Yes; >= {args.vlrr_bin_low} -> Partially; else No",
        "Interpretation note: rank_within_effector is the primary cross-candidate comparison. "
        "suggested_rank_within_backbone mixes different effectors and is exploratory because "
        "raw contact-pair counts and confidence calibration can be effector-dependent.",
    ]
    if vlrr_ranges:
        for bb, rs in sorted(vlrr_ranges.items()):
            provenance.append(f"vLRR range {bb}: {min(rs)}..{max(rs)} ({len(rs)} residues)")

    jobs = discover_jobs(roots)
    if not jobs:
        raise SystemExit("No candidate folders matching rx*_and_*_atp or sr35*_and_*_atp were found.")
    eprint(f"Found {len(jobs)} candidate folders.")

    raw_rows: List[Dict[str, Any]] = []
    for n, job in enumerate(jobs, start=1):
        eprint(f"[{n}/{len(jobs)}] {job['folder_name']}")
        raw_rows.extend(process_job(job, lookup, lis_script, ipsae_script,
                                    tool_out, args.lis_workers, vlrr_ranges,
                                    args.sasa_backend, args.interface_source,
                                    args.lis_extra_args))

    if not raw_rows:
        raise SystemExit("No model rows were successfully analysed.")

    summary_rows = summarize_rows(raw_rows, vlrr_ranges, args.vlrr_bin_low,
                                  args.vlrr_bin_high, args.interface_source)
    qc_rows = make_qc_table(summary_rows, args.expected_models,
                            args.nlr_plddt_warn, args.effector_plddt_warn)

    qc_key = {(q["backbone"], q["candidate_id"], q["effector"]): q for q in qc_rows}
    for s in summary_rows:
        q = qc_key.get((s["backbone"], s["candidate_id"], s["effector"]), {})
        s["QC_status"] = q.get("QC_status", "")
        s["QC_notes"] = q.get("QC_notes", "")

    # [v2] cLIR numbering diagnostic
    numbering_lines = report_clir_numbering(raw_rows)
    if numbering_lines:
        provenance.append("")
        provenance.append("cLIR index ranges observed (check against real chain lengths "
                          "before trusting vLRR_fraction):")
        provenance.extend(numbering_lines)

    (out / "tool_provenance.txt").write_text("\n".join(provenance) + "\n", encoding="utf-8")

    raw_csv = out / "raw_data.csv"
    summary_csv = out / "summary_mean.csv"
    qc_csv = out / "qc.csv"
    write_csv(raw_csv, raw_rows)
    write_csv(summary_csv, summary_rows)
    write_csv(qc_csv, qc_rows)
    write_xlsx(out / "AF3_candidate_metrics.xlsx", [
        ("Raw_Data", raw_rows),
        ("Summary_Mean", summary_rows),
        ("QC", qc_rows),
    ])

    # ------------------------------ run report ------------------------------
    n_contact = sum(1 for r in raw_rows if r.get("contact_prob_max") is not None)
    n_vlrr = sum(1 for r in raw_rows if r.get("vLRR_fraction") is not None)
    n_consistency = sum(1 for s in summary_rows if s.get("interface_consistency_mean") is not None)
    n_check = sum(1 for q in qc_rows if q["QC_status"] == "CHECK")

    print(f"Done. Analysed {len(raw_rows)} models across {len(summary_rows)} candidate/backbone pairs.")
    print(f"  contact_probs available:      {n_contact}/{len(raw_rows)} models")
    print(f"  vLRR_fraction computed:       {n_vlrr}/{len(raw_rows)} models")
    print(f"  interface consistency usable: {n_consistency}/{len(summary_rows)} candidates")
    print(f"  QC status CHECK:              {n_check}/{len(qc_rows)} candidates")
    if numbering_lines:
        print("cLIR index ranges observed (verify numbering before trusting vLRR_fraction):")
        for line in numbering_lines:
            print(line)
    print(f"Raw:        {raw_csv}")
    print(f"Summary:    {summary_csv}")
    print(f"QC:         {qc_csv}")
    print(f"Excel:      {out / 'AF3_candidate_metrics.xlsx'}")
    print(f"Provenance: {out / 'tool_provenance.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
