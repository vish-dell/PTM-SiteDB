from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
try:
    from scipy.stats import fisher_exact, mannwhitneyu
except Exception:
    fisher_exact = None
    mannwhitneyu = None

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
PANEL_PATH = PROC / "benchmark_panel.csv"
DEFAULT_PTM_MAPPED = PROC / "ptm_sites_mapped.tsv"
DEFAULT_PTM = PROC / "ptm_sites.tsv"
DEFAULT_OPPORTUNITIES = PROC / "context_opportunities.tsv"
DEFAULT_STRUCTURE = PROC / "alphafold_ca_coordinates.tsv"
RESULTS = PROC / "cohort_results"
RESULTS.mkdir(parents=True, exist_ok=True)

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def normalize_chrom(s: object) -> str:
    value = str(s).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if value.upper() in {"M", "MT"}:
        return "MT"
    return value


def normalize_gene(s: object) -> str:
    return str(s).strip().upper()


def context96(left: str, ref: str, alt: str, right: str) -> str | None:
    left, ref, alt, right = left.upper(), ref.upper(), alt.upper(), right.upper()
    if any(len(x) != 1 for x in (left, ref, alt, right)) or any(x not in "ACGT" for x in (left, ref, alt, right)) or ref == alt:
        return None
    if ref in "AG":
        left, ref, alt, right = right.translate(COMP), ref.translate(COMP), alt.translate(COMP), left.translate(COMP)
    return f"{left}[{ref}>{alt}]{right}"


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    nums = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = nums.dropna()
    n = len(valid)
    if not n:
        return out
    arr = valid.to_numpy(dtype=float)
    order = np.argsort(arr)
    q = np.empty(n, dtype=float)
    prev = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        value = min(prev, arr[idx] * n / rank)
        q[idx] = min(1.0, value)
        prev = q[idx]
    out.loc[valid.index] = q
    return out


def pick_col(df: pd.DataFrame, choices: list[str], required: bool = True) -> str | None:
    def norm(x: object) -> str:
        return re.sub(r"[^a-z0-9]", "", str(x).lower())
    look = {norm(c): c for c in df.columns}
    for choice in choices:
        if norm(choice) in look:
            return look[norm(choice)]
    for c in df.columns:
        cc = norm(c)
        if any(norm(x) in cc for x in choices):
            return c
    if required:
        raise RuntimeError(f"Cannot find one of columns {choices}. Actual columns: {list(df.columns)}")
    return None


def resolve_fasta_contig(fasta, chrom: str) -> str:
    refs = set(fasta.references)
    n = str(chrom)
    options = [n, n.removeprefix("chr"), "chr" + n.removeprefix("chr")]
    if n.upper() in {"M", "MT", "CHRM", "CHRMT"}:
        options.extend(["MT", "M", "chrM", "chrMT"])
    for x in options:
        if x in refs:
            return x
    raise KeyError(n)


def load_coordinates(path: Path, min_plddt: float | None) -> tuple[dict[tuple[str, int], tuple[float, float, float]], dict]:
    d = pd.read_csv(path, sep=None, engine="python")
    required = {"gene", "position", "x", "y", "z"}
    if not required.issubset(d.columns):
        raise RuntimeError(f"Coordinate file must contain {sorted(required)}")
    d["gene"] = d["gene"].map(normalize_gene)
    for c in ("position", "x", "y", "z"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    n_before = len(d)
    d = d.dropna(subset=["position", "x", "y", "z"]).copy()
    d["position"] = d["position"].astype(int)
    n_low = 0
    if min_plddt is not None and "plddt" in d.columns:
        d["plddt"] = pd.to_numeric(d["plddt"], errors="coerce")
        n_low = int((d["plddt"] < min_plddt).fillna(False).sum())
        d = d[d["plddt"].isna() | (d["plddt"] >= min_plddt)]
    # Keep one coordinate per gene-position; duplicates make structural numbering ambiguous.
    duplicated = d.duplicated(["gene", "position"], keep=False)
    if duplicated.any():
        raise RuntimeError("Coordinate file has duplicate gene-position rows. Resolve isoform/structure numbering before analysis.")
    coords = CoordMap({(r.gene, int(r.position)): (float(r.x), float(r.y), float(r.z)) for r in d.itertuples()})
    return coords, {"coordinate_rows_read": n_before, "coordinate_rows_valid": len(coords), "low_confidence_rows_excluded": n_low}


def read_maf(path: Path, fasta_path: Path, opportunities: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    try:
        from fasta_utils import IndexedFasta
    except ImportError:
        try:
            from scripts.fasta_utils import IndexedFasta
        except ImportError:
            try:
                from research.scripts.fasta_utils import IndexedFasta
            except ImportError as exc:
                raise RuntimeError("Included portable FASTA indexer could not be imported; verify the project folder is intact.") from exc
    if not Path(str(fasta_path) + ".fai").exists():
        # Constructing IndexedFasta will build an index, but a missing index at this point
        # is reported explicitly so inputs can be checked before analysis.
        try:
            with IndexedFasta(fasta_path):
                pass
        except Exception as exc:
            raise RuntimeError(f"Could not create/read FASTA index for {fasta_path}: {exc}") from exc
    maf = pd.read_csv(path, sep="\t", comment="#", low_memory=False, compression="infer")
    gene_c = pick_col(maf, ["Hugo_Symbol", "gene", "gene_symbol"])
    chr_c = pick_col(maf, ["Chromosome", "chromosome", "chr"])
    pos_c = pick_col(maf, ["Start_Position", "start_position", "genomic_position"])
    ref_c = pick_col(maf, ["Reference_Allele", "ref"])
    alt_c = pick_col(maf, ["Tumor_Seq_Allele2", "Tumor_Seq_Allele1", "alt"])
    class_c = pick_col(maf, ["Variant_Classification", "consequence", "classification"], required=False)
    sample_c = pick_col(maf, ["Tumor_Sample_Barcode", "sample", "sample_id"], required=False)

    maf_rows = []
    skipped = Counter()
    eligible_genes = set(opportunities["gene"].astype(str).map(normalize_gene)) if "gene" in opportunities.columns else set()
    with IndexedFasta(str(fasta_path)) as fasta:
        for r in maf.itertuples(index=False, name=None):
            vals = dict(zip(maf.columns, r))
            gene = normalize_gene(vals.get(gene_c, ""))
            if eligible_genes and gene not in eligible_genes:
                skipped["outside_panel_or_opportunity_gene"] += 1
                continue
            chrom = normalize_chrom(vals.get(chr_c, ""))
            try:
                pos = int(vals.get(pos_c))
            except Exception:
                skipped["invalid_start_position"] += 1
                continue
            ref = str(vals.get(ref_c, "")).upper()
            alt = str(vals.get(alt_c, "")).upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                skipped["non_snv_or_invalid_allele"] += 1
                continue
            classification = str(vals.get(class_c, "")).strip().lower() if class_c else ""
            if classification in {"nan", "none", "na", "n/a", "."}:
                classification = ""
            if classification and not ("missense" in classification or "nonsynonymous" in classification):
                skipped["non_missense_maf_classification"] += 1
                continue
            try:
                contig = resolve_fasta_contig(fasta, chrom)
                triplet = fasta.fetch(contig, pos - 2, pos + 1).upper()
            except Exception:
                skipped["fasta_contig_or_flank_unavailable"] += 1
                continue
            if len(triplet) != 3 or triplet[1] != ref:
                skipped["reference_mismatch_or_edge"] += 1
                continue
            ctx = context96(triplet[0], ref, alt, triplet[2])
            maf_rows.append({
                "gene": gene, "chromosome": chrom, "genomic_position": pos,
                "ref": ref, "alt": alt, "context96": ctx,
                "sample": str(vals.get(sample_c, "")) if sample_c else "",
            })

    variants = pd.DataFrame(maf_rows)
    if variants.empty:
        raise RuntimeError("No usable single-nucleotide missense variants remained after MAF/FASTA validation.")

    opp = opportunities.copy()
    required = {"gene", "protein_pos", "context96", "chromosome", "genomic_position", "ref", "alt"}
    if not required.issubset(opp.columns):
        raise RuntimeError(f"Opportunity map must contain {sorted(required)}; found {list(opp.columns)}")
    opp["gene"] = opp["gene"].map(normalize_gene)
    opp["chromosome"] = opp["chromosome"].map(normalize_chrom)
    opp["genomic_position"] = pd.to_numeric(opp["genomic_position"], errors="coerce")
    opp["protein_pos"] = pd.to_numeric(opp["protein_pos"], errors="coerce")
    opp["ref"] = opp["ref"].astype(str).str.upper()
    opp["alt"] = opp["alt"].astype(str).str.upper()
    opp["context96"] = opp["context96"].astype(str)
    opp = opp.dropna(subset=["genomic_position", "protein_pos"]).copy()
    opp["genomic_position"] = opp["genomic_position"].astype(int)
    opp["protein_pos"] = opp["protein_pos"].astype(int)
    key = ["gene", "chromosome", "genomic_position", "ref", "alt"]
    if opp.duplicated(key).any():
        raise RuntimeError("Opportunity map contains duplicate genomic variant keys. Ensure one MANE/selected transcript per gene.")
    joined = variants.merge(opp[key + ["protein_pos", "context96", "transcript_id"] if "transcript_id" in opp.columns else key + ["protein_pos", "context96"]], on=key, how="left", suffixes=("_maf", "_opp"), indicator=True)
    # The independent genomic context calculated from MAF + FASTA must agree with the enumerated opportunity map.
    joined["context_ok"] = joined["context96_maf"].eq(joined["context96_opp"])
    mapped = joined[(joined["_merge"] == "both") & joined["context_ok"]].copy()
    mapped = mapped.rename(columns={"context96_maf": "context96"})
    mapped = mapped.drop(columns=[c for c in ["_merge", "context_ok", "context96_opp"] if c in mapped.columns])
    unmapped = joined[joined["_merge"] != "both"]
    ctx_mismatch = joined[(joined["_merge"] == "both") & ~joined["context_ok"]]
    skipped["maf_missense_not_mapped_to_selected_transcript_opportunity"] += int(len(unmapped))
    skipped["opportunity_context_mismatch"] += int(len(ctx_mismatch))
    candidate_by_gene = variants.groupby("gene").size().to_dict()
    unmapped_by_gene = unmapped.groupby("gene").size().to_dict()
    context_mismatch_by_gene = ctx_mismatch.groupby("gene").size().to_dict()
    meta = {
        "maf_rows": int(len(maf)),
        "candidate_maf_missense_snvs": int(len(variants)),
        "variants_mapped_to_selected_transcript": int(len(mapped)),
        "skipped": dict(skipped),
        "candidate_missense_by_gene": {str(k): int(v) for k, v in candidate_by_gene.items()},
        "unmapped_by_gene": {str(k): int(v) for k, v in unmapped_by_gene.items()},
        "context_mismatch_by_gene": {str(k): int(v) for k, v in context_mismatch_by_gene.items()},
        "mapping_fraction": float(len(mapped) / len(variants)) if len(variants) else 0.0,
        "reference_fasta": str(fasta_path.resolve()),
    }
    return mapped, meta


class CoordMap(dict):
    """dict[(gene, position)] -> (x, y, z) that also caches a per-gene array index."""
    _index = None


def gene_arrays(coords: dict, gene: str) -> tuple[np.ndarray, np.ndarray]:
    index = getattr(coords, "_index", None)
    if index is None:
        grouped = defaultdict(list)
        for (g, p), xyz in coords.items():
            grouped[g].append((int(p), xyz))
        index = {}
        for g, rows in grouped.items():
            rows.sort()
            index[g] = (np.array([r[0] for r in rows], dtype=int), np.array([r[1] for r in rows], dtype=float))
        if isinstance(coords, CoordMap):
            coords._index = index
    return index.get(gene, (np.array([], dtype=int), np.zeros((0, 3), dtype=float)))


def near_matrix(P: np.ndarray, X: np.ndarray, sep: int, dist: float) -> np.ndarray:
    """Boolean n x n matrix: residues i, j are >= sep apart in sequence and within dist (CA-CA) in 3D."""
    n = len(P)
    NM = np.zeros((n, n), dtype=bool)
    d2max = float(dist) * float(dist)
    step = 256
    for i in range(0, n, step):
        blk = X[i:i + step]
        d2 = ((blk[:, None, :] - X[None, :, :]) ** 2).sum(axis=2)
        NM[i:i + step] = (d2 <= d2max) & (np.abs(P[i:i + step, None] - P[None, :]) >= sep)
    return NM


def near_mask(positions: np.ndarray, ptm_positions: np.ndarray, coords: dict, gene: str, sep: int, dist: float) -> np.ndarray:
    positions = np.asarray(positions, dtype=int)
    out = np.zeros(len(positions), dtype=bool)
    anchors = [(int(p), coords[(gene, int(p))]) for p in ptm_positions if (gene, int(p)) in coords]
    if not anchors or len(positions) == 0:
        return out
    uniq, inv = np.unique(positions, return_inverse=True)
    have = np.array([(gene, int(u)) in coords for u in uniq], dtype=bool)
    flag = np.zeros(len(uniq), dtype=bool)
    if have.any():
        U = np.array([coords[(gene, int(u))] for u in uniq[have]], dtype=float)
        ap = np.array([a[0] for a in anchors], dtype=int)
        A = np.array([a[1] for a in anchors], dtype=float)
        d2 = ((U[:, None, :] - A[None, :, :]) ** 2).sum(axis=2)
        far = np.abs(uniq[have][:, None] - ap[None, :]) >= sep
        flag[have] = ((d2 <= float(dist) ** 2) & far).any(axis=1)
    return flag[inv]


def collapse_recurrence(obs: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Recurrent hotspots must not count as independent events.

    distinct_residue (primary): each mutated residue counts once, however many tumours carry it
    or how many different substitutions it has. all_events: legacy per-variant counting (sensitivity only).
    """
    if mode == "all_events" or obs.empty:
        return obs
    if mode != "distinct_residue":
        raise ValueError(f"Unknown recurrence mode: {mode}")
    sort_cols = [c for c in ("protein_pos", "genomic_position", "alt") if c in obs.columns]
    return obs.sort_values(sort_cols).drop_duplicates("protein_pos", keep="first")


def null_probabilities(obs: pd.DataFrame, opp: pd.DataFrame, near_opp: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[str]]:
    counts = Counter(obs["context96"].astype(str))
    obs_contexts = sorted(counts)
    opp_contexts = opp["context96"].astype(str).to_numpy()
    if not obs_contexts:
        return np.array([], dtype=int), np.array([], dtype=float), []
    context_counts = np.array([counts[c] for c in obs_contexts], dtype=int)
    probs = []
    for ctx in obs_contexts:
        mask = opp_contexts == ctx
        if int(mask.sum()) == 0:
            raise RuntimeError(f"No gene-specific coding opportunity exists for observed SBS96 context {ctx}.")
        probs.append(float(np.mean(near_opp[mask])))
    return context_counts, np.asarray(probs, dtype=float), obs_contexts


def sample_null(n_perm: int, context_counts: np.ndarray, probs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(context_counts) == 0:
        return np.zeros(n_perm, dtype=int)
    draw = rng.binomial(n=context_counts.reshape(1, -1), p=probs.reshape(1, -1), size=(n_perm, len(context_counts)))
    return draw.sum(axis=1)


def decoy_calibration(gene: str, anchor_positions: np.ndarray, eligible_opp: pd.DataFrame, obs: pd.DataFrame,
                      context_counts: np.ndarray, contexts: list[str], coords: dict, sep: int, dist: float,
                      n_decoy: int, seed: int) -> dict:
    """PTM-specificity control.

    The real PTM anchors are compared with the same number of NON-PTM residues of the same amino-acid
    identity (S->S, T->T, Y->Y), drawn from structured residues of the same gene. For every anchor set
    the statistic is  excess = observed proximal mutated residues - exact context-null expectation.
    A PTM-specific signal must exceed what matched non-PTM residues give in the same gene, which removes
    the effect of 'mutations in this gene simply cluster near S/T/Y-rich or functional regions'.
    """
    base = {"n_decoy": 0, "decoy_status": "not_run"}
    if "ref_aa" not in eligible_opp.columns:
        return {**base, "decoy_status": "unavailable_no_ref_aa_in_opportunity_map"}
    P, X = gene_arrays(coords, gene)
    if len(P) == 0:
        return {**base, "decoy_status": "unavailable_no_gene_coordinates"}
    aa_by_pos = (eligible_opp.drop_duplicates("protein_pos").set_index("protein_pos")["ref_aa"].astype(str).str.upper().to_dict())
    anchors = [int(p) for p in anchor_positions]
    anchor_aa = [aa_by_pos.get(p) for p in anchors]
    if any(a is None for a in anchor_aa):
        return {**base, "decoy_status": "unavailable_anchor_residue_identity_unknown"}
    need = Counter(anchor_aa)
    anchor_set = set(anchors)
    pools = {aa: np.array([i for i, p in enumerate(P) if int(p) not in anchor_set and aa_by_pos.get(int(p)) == aa], dtype=int)
             for aa in need}
    ratio = min(len(pools[aa]) / k for aa, k in need.items())
    # Require >= 2x as many eligible decoy residues as anchors so decoy sets are not near-copies of the real set.
    if any(len(pools[aa]) < 2 * k for aa, k in need.items()):
        return {**base, "decoy_status": "unavailable_decoy_pool_too_small", "decoy_pool_ratio": float(ratio)}

    NM = near_matrix(P, X, sep, dist)
    idx_opp = np.searchsorted(P, eligible_opp["protein_pos"].astype(int).to_numpy())
    idx_obs = np.searchsorted(P, obs["protein_pos"].astype(int).to_numpy())
    K = len(contexts)
    ctx_index = {c: i for i, c in enumerate(contexts)}
    opp_ctx = np.array([ctx_index.get(c, -1) for c in eligible_opp["context96"].astype(str)], dtype=int)
    valid = opp_ctx >= 0
    tot = np.bincount(opp_ctx[valid], minlength=K).astype(float)

    def excess(anchor_idx: np.ndarray) -> tuple[float, int, float]:
        near_pos = NM[:, anchor_idx].any(axis=1)
        observed = int(near_pos[idx_obs].sum()) if len(idx_obs) else 0
        if K == 0:
            return float(observed), observed, 0.0
        sums = np.bincount(opp_ctx[valid], weights=near_pos[idx_opp][valid].astype(float), minlength=K)
        expected = float((context_counts * (sums / tot)).sum())
        return observed - expected, observed, expected

    real_idx = np.searchsorted(P, np.array(anchors, dtype=int))
    real_excess, _, _ = excess(real_idx)
    rng = np.random.default_rng(seed + 7919 + sum(ord(c) for c in gene))
    decoy_excess = np.empty(int(n_decoy), dtype=float)
    for d in range(int(n_decoy)):
        chosen = np.concatenate([rng.choice(pools[aa], size=k, replace=False) for aa, k in need.items()])
        decoy_excess[d], _, _ = excess(chosen)
    mean_d = float(decoy_excess.mean())
    sd_d = float(decoy_excess.std(ddof=1)) if len(decoy_excess) > 1 else 0.0
    if sd_d > 1e-12:
        z = float((real_excess - mean_d) / sd_d)
    else:
        z = 0.0 if abs(real_excess - mean_d) < 1e-12 else float("nan")
    return {
        "n_decoy": int(n_decoy), "decoy_status": "ok", "decoy_pool_ratio": float(ratio),
        "excess_over_context_expectation": float(real_excess),
        "decoy_mean_excess": mean_d, "decoy_sd_excess": sd_d, "decoy_z": z,
        "p_decoy": float((1 + np.sum(decoy_excess >= real_excess - 1e-12)) / (len(decoy_excess) + 1)),
        "decoy_control": "same-gene non-PTM residues matched on amino-acid identity and anchor count",
    }


def setting_gene_result(gene: str, ptm_positions: list[int], mut: pd.DataFrame, opp: pd.DataFrame,
                        coords: dict, sep: int, dist: float, n_perm: int, seed: int,
                        structure_mapping_ok: bool, recurrence: str = "distinct_residue",
                        n_decoy: int = 0) -> dict:
    ptm_positions = sorted(set(int(x) for x in ptm_positions))
    mutations = mut[mut["gene"] == gene].copy()
    if not ptm_positions:
        return {"gene": gene, "status": "untestable_no_mapped_ptm_anchors", "primary_test_eligible": False}

    # Genes with zero observed missense events remain testable when their transcript,
    # PTM anchors, coordinates and opportunity map pass QC (observed 0, p = 1).
    obs_has_coords = mutations["protein_pos"].map(lambda x: (gene, int(x)) in coords) if not mutations.empty else pd.Series([], dtype=bool, index=mutations.index)
    obs_all = mutations[obs_has_coords].copy()
    missing_structure = int((~obs_has_coords).sum())
    ptm_coords = [p for p in ptm_positions if (gene, p) in coords]
    if not ptm_coords:
        return {"gene": gene, "status": "untestable_no_ptm_anchor_coordinates", "primary_test_eligible": False,
                "n_missense_mapped": int(len(mutations)), "n_mutations_missing_coordinates": missing_structure,
                "n_ptm_anchors_source": len(ptm_positions), "n_ptm_anchors_with_coordinates": 0}
    if len(mutations) > 0 and obs_all.empty:
        return {"gene": gene, "status": "untestable_no_mutation_coordinates", "primary_test_eligible": False,
                "n_missense_mapped": int(len(mutations)), "n_mutations_missing_coordinates": missing_structure,
                "n_ptm_anchors_source": len(ptm_positions), "n_ptm_anchors_with_coordinates": len(ptm_coords)}

    eligible_opp = opp[opp["protein_pos"].map(lambda x: (gene, int(x)) in coords)].copy()
    if eligible_opp.empty:
        return {"gene": gene, "status": "untestable_no_structured_opportunities", "primary_test_eligible": False}
    anchor_positions = np.asarray(ptm_coords, dtype=int)
    opp_pos = eligible_opp["protein_pos"].astype(int).to_numpy()
    near_opp = near_mask(opp_pos, anchor_positions, coords, gene, sep, dist)

    obs = collapse_recurrence(obs_all, recurrence)
    near_obs = near_mask(obs["protein_pos"].astype(int).to_numpy(), anchor_positions, coords, gene, sep, dist)
    observed = int(near_obs.sum())
    try:
        context_counts, probs, contexts = null_probabilities(obs, eligible_opp, near_opp)
    except RuntimeError as exc:
        return {"gene": gene, "status": "untestable_unmatched_context_opportunity", "primary_test_eligible": False,
                "error": str(exc), "n_missense_mapped": int(len(mutations)),
                "n_mutations_missing_coordinates": missing_structure,
                "n_ptm_anchors_source": len(ptm_positions), "n_ptm_anchors_with_coordinates": len(ptm_coords)}
    seed_gene = int(seed + sep * 1009 + round(dist * 101) + sum(ord(c) for c in gene))
    rng = np.random.default_rng(seed_gene)
    null_counts = sample_null(n_perm, context_counts, probs, rng)
    p = float((1 + np.sum(null_counts >= observed)) / (n_perm + 1))
    expected = float(np.mean(null_counts))

    # Sensitivity: legacy per-variant counting (recurrent hotspots counted repeatedly).
    p_all = p
    observed_all = observed
    if recurrence == "distinct_residue" and len(obs_all) != len(obs):
        near_all = near_mask(obs_all["protein_pos"].astype(int).to_numpy(), anchor_positions, coords, gene, sep, dist)
        observed_all = int(near_all.sum())
        cc_a, pr_a, _ = null_probabilities(obs_all, eligible_opp, near_opp)
        null_a = sample_null(n_perm, cc_a, pr_a, np.random.default_rng(seed_gene + 1))
        p_all = float((1 + np.sum(null_a >= observed_all)) / (n_perm + 1))

    decoy = {}
    if n_decoy and n_decoy > 0:
        decoy = decoy_calibration(gene, anchor_positions, eligible_opp, obs, context_counts, contexts,
                                  coords, sep, dist, n_decoy, seed_gene)

    total_anchor_coverage = len(ptm_coords) / len(ptm_positions) if ptm_positions else 0.0
    status = "ok_no_observed_missense_events" if len(mutations) == 0 else "ok"
    testable = True
    if missing_structure > 0:
        status, testable = "partial_structure_mapping", False
    elif not structure_mapping_ok:
        status, testable = "partial_genomic_transcript_mapping", False
    elif total_anchor_coverage < 1.0:
        status, testable = "partial_ptm_anchor_mapping", False
    out = {
        "gene": gene,
        "status": status,
        "primary_test_eligible": bool(testable),
        "sequence_separation_min": sep,
        "distance_max_A": dist,
        "ptm_sites_mapped": len(ptm_positions),
        "ptm_anchors_with_coordinates": len(ptm_coords),
        "ptm_anchor_coordinate_fraction": total_anchor_coverage,
        "missense_mutations_mapped": int(len(mutations)),
        "missense_mutations_with_coordinates": int(len(obs_all)),
        "n_distinct_mutated_residues": int(obs_all["protein_pos"].nunique()) if len(obs_all) else 0,
        "mutations_without_coordinates": missing_structure,
        "recurrence_mode": recurrence,
        "observed_proximal_mutation_events": observed,
        "observed_proximal_events_all_variants": observed_all,
        "expected_proximal_events_context_null": expected,
        "fold_enrichment": float(observed / expected) if expected > 0 else (None if observed == 0 else math.inf),
        "p_context_null": p,
        "p_context_null_all_variants": p_all,
        "q_context_null": np.nan,
        "permutations": n_perm,
        "n_distinct_observed_contexts": len(contexts),
        "null_model": "gene-specific SBS96-matched coding missense opportunity permutation",
        "coordinate_source": "AlphaFold-derived or user-supplied CA coordinates",
    }
    out.update(decoy)
    return out


def site_results(primary: pd.DataFrame, ptm: pd.DataFrame, mut: pd.DataFrame, opportunities: pd.DataFrame,
                 coords: dict, n_perm: int, seed: int) -> pd.DataFrame:
    rows = []
    for gene, anchors in ptm.groupby("gene"):
        gene_mut_all = mut[mut["gene"] == gene]
        gene_opp = opportunities[opportunities["gene"] == gene]
        for anchor in sorted(set(pd.to_numeric(anchors["position"], errors="coerce").dropna().astype(int))):
            if (gene, int(anchor)) not in coords or gene_mut_all.empty:
                rows.append({"gene": gene, "position": int(anchor), "status": "untestable_missing_anchor_or_mutations"})
                continue
            obs_has = gene_mut_all["protein_pos"].map(lambda x: (gene, int(x)) in coords)
            obs = collapse_recurrence(gene_mut_all[obs_has], "distinct_residue")
            opp = gene_opp[gene_opp["protein_pos"].map(lambda x: (gene, int(x)) in coords)]
            if obs.empty or opp.empty:
                rows.append({"gene": gene, "position": int(anchor), "status": "untestable_no_structured_events"})
                continue
            # Site test evaluates one PTM anchor at a time, with the same pre-specified separation/distance threshold.
            op_pos = opp["protein_pos"].astype(int).to_numpy()
            ob_pos = obs["protein_pos"].astype(int).to_numpy()
            coords_gene = coords
            op_near = np.zeros(len(op_pos), dtype=bool)
            op_xyz = coords_gene[(gene, int(anchor))]
            for i, p in enumerate(op_pos):
                if abs(int(p) - int(anchor)) < 8:
                    continue
                xyz = coords_gene.get((gene, int(p)))
                if xyz is None:
                    continue
                d2 = sum((xyz[j]-op_xyz[j])**2 for j in range(3))
                op_near[i] = d2 <= 100.0
            ob_near = np.zeros(len(ob_pos), dtype=bool)
            for i, p in enumerate(ob_pos):
                if abs(int(p) - int(anchor)) < 8:
                    continue
                xyz = coords_gene.get((gene, int(p)))
                if xyz is None:
                    continue
                ob_near[i] = sum((xyz[j]-op_xyz[j])**2 for j in range(3)) <= 100.0
            observed = int(ob_near.sum())
            try:
                counts, probs, _ = null_probabilities(obs, opp, op_near)
            except RuntimeError:
                rows.append({"gene": gene, "position": int(anchor), "status": "untestable_unmatched_context_opportunity"})
                continue
            rng = np.random.default_rng(seed + 991 * int(anchor) + sum(ord(c) for c in gene))
            null = sample_null(n_perm, counts, probs, rng)
            pval = float((1 + np.sum(null >= observed)) / (n_perm + 1))
            rows.append({"gene": gene, "position": int(anchor), "ptm_type": ";".join(sorted(set(anchors["ptm_type"].astype(str)))),
                         "n_missense_events": int(len(obs)), "observed_near_events": observed,
                         "expected_near_events": float(np.mean(null)), "p_context_null": pval,
                         "q_context_null": np.nan, "status": "ok", "note": "site-level secondary analysis; not the primary protein-level endpoint"})
    out = pd.DataFrame(rows)
    if not out.empty and "p_context_null" in out.columns:
        out["q_context_null"] = benjamini_hochberg(out["p_context_null"])
    return out


def make_strata(df: pd.DataFrame, burden_col: str = "n_distinct_mutated_residues",
                anchor_col: str = "ptm_anchors_with_coordinates", bins: int | None = None) -> np.ndarray:
    """Strata from mutation burden x PTM-anchor count, so groups are compared like-for-like.

    bins=None picks tertiles for >= 60 genes, halves for 20-59, and burden-only halves below 20.
    """
    n = len(df)
    if bins is None:
        bins = 3 if n >= 60 else 2 if n >= 20 else 1
    def tert(col):
        v = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        r = v.rank(method="first")
        return pd.qcut(r, q=min(bins, max(1, len(v))), labels=False, duplicates="drop").astype(int).to_numpy()
    return tert(burden_col) * 100 + tert(anchor_col)


def stratified_label_permutation_test(scores: np.ndarray, is_driver: np.ndarray, strata: np.ndarray,
                                      n_perm: int = 10000, seed: int = 1) -> dict:
    """One-sided test that drivers score higher, permuting driver labels only WITHIN strata.

    Statistic: sum over driver genes of (score - stratum mean). Under within-stratum exchangeability this
    removes any difference explained by mutation burden or PTM-anchor count.
    """
    scores = np.asarray(scores, dtype=float)
    is_driver = np.asarray(is_driver, dtype=bool)
    strata = np.asarray(strata)
    if len(scores) == 0 or is_driver.sum() == 0 or (~is_driver).sum() == 0:
        return {"status": "not_estimable_missing_group"}
    groups, centered = [], np.empty_like(scores)
    informative = 0
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        centered[idx] = scores[idx] - scores[idx].mean()
        groups.append(idx)
        if 0 < is_driver[idx].sum() < len(idx):
            informative += 1
    if informative == 0:
        return {"status": "not_estimable_no_stratum_contains_both_groups"}
    obs_stat = float(centered[is_driver].sum())
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(int(n_perm)):
        t = 0.0
        for idx in groups:
            k = int(is_driver[idx].sum())
            if 0 < k < len(idx):
                t += centered[rng.choice(idx, size=k, replace=False)].sum()
        count += t >= obs_stat - 1e-12
    return {"status": "ok", "statistic": obs_stat, "one_sided_p_value": float((1 + count) / (n_perm + 1)),
            "n_strata": int(len(groups)), "n_informative_strata": int(informative), "n_permutations": int(n_perm)}


def primary_group_comparison(result: pd.DataFrame, panel: pd.DataFrame, n_strat_perm: int = 10000) -> dict:
    """Cohort-level driver-vs-not-listed comparisons.

    v1.0 test (kept, now SECONDARY): one-sided Mann-Whitney U on -log10(p_context_null).
    v1.1 PRIMARY: PTM-specific score (decoy_z) compared between groups with a stratified permutation test
    (strata = mutation burden x PTM-anchor count). A burden-stratified test on the raw score is also reported.
    """
    required_result = {"gene", "primary_test_eligible", "p_context_null", "q_context_null", "status"}
    required_panel = {"gene", "intogen_driver"}
    if result.empty or panel.empty or not required_result.issubset(result.columns) or not required_panel.issubset(panel.columns):
        return {"status": "not_estimable_missing_inputs", "test": "one-sided Mann–Whitney U on -log10(p_context_null)"}
    labels = panel[["gene", "intogen_driver"]].copy()
    labels["gene"] = labels["gene"].map(normalize_gene)
    data = labels.merge(result, on="gene", how="left", suffixes=("_panel", ""))
    if "intogen_driver_panel" in data.columns and "intogen_driver" not in data.columns:
        data["intogen_driver"] = data["intogen_driver_panel"]
    eligible = data["primary_test_eligible"].fillna(False).astype(bool)
    eligible &= data["status"].astype(str).isin({"ok", "ok_no_observed_missense_events"})
    pvals = pd.to_numeric(data["p_context_null"], errors="coerce")
    qvals = pd.to_numeric(data["q_context_null"], errors="coerce")
    eligible &= pvals.notna() & qvals.notna()
    tested = data.loc[eligible].copy()
    pvals = pd.to_numeric(tested["p_context_null"], errors="coerce").clip(lower=1e-300, upper=1.0)
    tested["evidence_score"] = -np.log10(pvals)
    is_driver = tested["intogen_driver"].astype(str).str.lower().isin({"true", "1", "yes"})
    driver_scores = tested.loc[is_driver, "evidence_score"].to_numpy(dtype=float)
    not_listed_scores = tested.loc[~is_driver, "evidence_score"].to_numpy(dtype=float)
    output = {
        "status": "ok" if len(driver_scores) and len(not_listed_scores) else "not_estimable_missing_group",
        "test": "SECONDARY (v1.0): one-sided Mann–Whitney U on per-gene -log10(p_context_null); drivers > not-listed",
        "eligible_gene_count": int(len(tested)),
        "driver_genes_tested": int(len(driver_scores)),
        "not_listed_genes_tested": int(len(not_listed_scores)),
        "driver_median_evidence_score": float(np.median(driver_scores)) if len(driver_scores) else None,
        "not_listed_median_evidence_score": float(np.median(not_listed_scores)) if len(not_listed_scores) else None,
        "median_evidence_score_difference_driver_minus_not_listed": (float(np.median(driver_scores) - np.median(not_listed_scores)) if len(driver_scores) and len(not_listed_scores) else None),
        "mann_whitney_u": None, "one_sided_p_value": None, "cliffs_delta": None,
        "secondary_hit_fraction_fisher_p_greater": None,
        "secondary_driver_q_lt_005_hits": None, "secondary_not_listed_q_lt_005_hits": None,
        "burden_stratified_raw_score_test": {"status": "not_run"},
        "primary_v1_1_decoy_calibrated_stratified_test": {"status": "not_run"},
    }
    if len(driver_scores) and len(not_listed_scores):
        if mannwhitneyu is not None:
            test = mannwhitneyu(driver_scores, not_listed_scores, alternative="greater", method="auto")
            output["mann_whitney_u"] = float(test.statistic)
            output["one_sided_p_value"] = float(test.pvalue)
        pairwise = driver_scores[:, None] - not_listed_scores[None, :]
        output["cliffs_delta"] = float((np.sum(pairwise > 0) - np.sum(pairwise < 0)) / pairwise.size)
        hits = pd.to_numeric(tested["q_context_null"], errors="coerce") < 0.05
        driver_hits = int((hits & is_driver).sum())
        other_hits = int((hits & ~is_driver).sum())
        output["secondary_driver_q_lt_005_hits"] = driver_hits
        output["secondary_not_listed_q_lt_005_hits"] = other_hits
        if fisher_exact is not None:
            fisher = fisher_exact([[driver_hits, len(driver_scores)-driver_hits], [other_hits, len(not_listed_scores)-other_hits]], alternative="greater")
            f_p = float(fisher.pvalue)
            output["secondary_hit_fraction_fisher_p_greater"] = f_p if math.isfinite(f_p) else None
        if {"n_distinct_mutated_residues", "ptm_anchors_with_coordinates"}.issubset(tested.columns):
            strata = make_strata(tested)
            output["burden_stratified_raw_score_test"] = stratified_label_permutation_test(
                tested["evidence_score"].to_numpy(float), is_driver.to_numpy(), strata, n_strat_perm, seed=1)
            if "decoy_z" in tested.columns:
                dz = pd.to_numeric(tested["decoy_z"], errors="coerce")
                keep = dz.notna() & np.isfinite(dz)
                sub = tested.loc[keep]
                if len(sub):
                    res = stratified_label_permutation_test(
                        dz[keep].to_numpy(float), is_driver[keep].to_numpy(), make_strata(sub), n_strat_perm, seed=1)
                    res.update({"score": "decoy_z (PTM anchors vs matched non-PTM residues)",
                                "driver_genes_with_decoy_score": int(is_driver[keep].sum()),
                                "not_listed_genes_with_decoy_score": int((~is_driver[keep]).sum()),
                                "driver_median_decoy_z": float(dz[keep][is_driver[keep]].median()) if is_driver[keep].any() else None,
                                "not_listed_median_decoy_z": float(dz[keep][~is_driver[keep]].median()) if (~is_driver[keep]).any() else None})
                    output["primary_v1_1_decoy_calibrated_stratified_test"] = res
    return output


def sha256_file(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _source_manifest_dict() -> dict:
    path = PROC / "source_manifest.json"
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    except Exception:
        return {}


def _intogen_release_from_manifest() -> str | None:
    source = _source_manifest_dict().get("sources", {}).get("IntOGen", {})
    return source.get("release_declared_inside_archive")


def _ptm_source_metadata_from_manifest() -> dict:
    source = _source_manifest_dict().get("sources", {}).get("dbPTM_phosphorylation", {})
    return {key: source.get(key) for key in ("url", "path", "sha256", "extracted_table", "retrieval_date_utc", "status") if source.get(key) is not None}


def main() -> int:
    ap = argparse.ArgumentParser(description="PTM-anchor 3D analysis with gene-specific SBS96 opportunity matching.")
    ap.add_argument("maf", type=Path)
    ap.add_argument("fasta", type=Path)
    ap.add_argument("--cohort", default="cohort")
    ap.add_argument("--opportunities", type=Path, default=DEFAULT_OPPORTUNITIES)
    ap.add_argument("--structure", type=Path, default=DEFAULT_STRUCTURE)
    ap.add_argument("--ptm", type=Path, default=None)
    ap.add_argument("--panel", type=Path, default=PANEL_PATH)
    ap.add_argument("--sep", type=int, default=8)
    ap.add_argument("--dist", type=float, default=10.0)
    ap.add_argument("--permutations", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-plddt", type=float, default=None, help="Optional pLDDT threshold applied to coordinate rows")
    ap.add_argument("--decoys", type=int, default=1000, help="Matched non-PTM decoy anchor sets per gene (PTM-specificity control); 0 disables")
    ap.add_argument("--recurrence", choices=["distinct_residue", "all_events"], default="distinct_residue",
                    help="distinct_residue (primary): a recurrently mutated residue counts once. all_events: legacy sensitivity mode")
    ap.add_argument("--robustness", action="store_true", help="Also run the full predeclared 3 x 4 secondary threshold matrix")
    args = ap.parse_args()

    for p in [args.opportunities, args.structure, args.panel]:
        if not p.exists():
            raise SystemExit(f"Required real analysis input missing: {p}")
    ptm_path = args.ptm or (DEFAULT_PTM_MAPPED if DEFAULT_PTM_MAPPED.exists() else DEFAULT_PTM)
    if not ptm_path.exists():
        raise SystemExit(f"PTM source missing: {ptm_path}; run bootstrap first.")

    ptm = pd.read_csv(ptm_path, sep="\t", low_memory=False)
    panel = pd.read_csv(args.panel)
    opportunities = pd.read_csv(args.opportunities, sep="\t", low_memory=False)
    if not {"gene", "position", "ptm_type"}.issubset(ptm.columns):
        raise SystemExit("Mapped PTM file must have gene, position, ptm_type columns.")
    ptm["gene"] = ptm["gene"].map(normalize_gene)
    ptm["position"] = pd.to_numeric(ptm["position"], errors="coerce")
    ptm = ptm.dropna(subset=["position"]).copy()
    ptm["position"] = ptm["position"].astype(int)
    ptm["ptm_type"] = ptm["ptm_type"].astype(str)
    ptm = ptm[ptm["gene"].isin(set(panel["gene"].map(normalize_gene)))].copy()
    opportunities["gene"] = opportunities["gene"].map(normalize_gene)
    opportunities["context96"] = opportunities["context96"].astype(str)
    opportunities["protein_pos"] = pd.to_numeric(opportunities["protein_pos"], errors="coerce")
    opportunities = opportunities.dropna(subset=["protein_pos", "context96"]).copy()
    opportunities["protein_pos"] = opportunities["protein_pos"].astype(int)
    mut, maf_qc = read_maf(args.maf, args.fasta, opportunities)
    mut["gene"] = mut["gene"].map(normalize_gene)
    coords, coord_qc = load_coordinates(args.structure, args.min_plddt)

    panel["gene"] = panel["gene"].map(normalize_gene)
    if "intogen_driver" not in panel.columns:
        panel["intogen_driver"] = False
    panel["intogen_driver"] = panel["intogen_driver"].astype(str).str.lower().isin({"true", "1", "yes"})
    ptm_positions = {g: sorted(set(group["position"].astype(int))) for g, group in ptm.groupby("gene")}

    result_rows = []
    for gene in sorted(panel["gene"].unique()):
        gene_mut = mut[mut["gene"] == gene]
        gene_candidate_n = int(maf_qc["candidate_missense_by_gene"].get(gene, 0))
        gene_unmapped_n = int(maf_qc["unmapped_by_gene"].get(gene, 0))
        gene_context_mismatch_n = int(maf_qc["context_mismatch_by_gene"].get(gene, 0))
        transcript_mapping_ok = (gene_unmapped_n == 0 and gene_context_mismatch_n == 0)
        result = setting_gene_result(gene, ptm_positions.get(gene, []), mut, opportunities, coords,
                                     args.sep, args.dist, args.permutations, args.seed,
                                     structure_mapping_ok=transcript_mapping_ok,
                                     recurrence=args.recurrence, n_decoy=args.decoys)
        if not transcript_mapping_ok:
            result["status"] = "partial_genomic_transcript_mapping"
            result["primary_test_eligible"] = False
        result["intogen_driver"] = bool(panel.loc[panel["gene"] == gene, "intogen_driver"].iloc[0])
        result["n_maf_candidate_missense_events"] = gene_candidate_n
        result["n_maf_mapped_events"] = int(len(gene_mut))
        result["n_missense_not_mapped_to_selected_transcript"] = gene_unmapped_n
        result["n_context_mismatch_events"] = gene_context_mismatch_n
        result["maf_mapping_fraction_gene"] = float(len(gene_mut) / gene_candidate_n) if gene_candidate_n else 0.0
        result["maf_mapping_fraction_global"] = maf_qc["mapping_fraction"]
        result["source_ptm_file"] = str(ptm_path.relative_to(ROOT)) if ptm_path.is_relative_to(ROOT) else str(ptm_path)
        result_rows.append(result)
    result = pd.DataFrame(result_rows)
    if "p_context_null" not in result.columns:
        result["p_context_null"] = np.nan
    result["q_context_null"] = np.nan
    eligible = result["primary_test_eligible"].fillna(False).astype(bool) if "primary_test_eligible" in result.columns else pd.Series(False, index=result.index)
    result.loc[eligible, "q_context_null"] = benjamini_hochberg(result.loc[eligible, "p_context_null"]).values
    result["q_decoy"] = np.nan
    if "p_decoy" in result.columns:
        dmask = eligible & pd.to_numeric(result["p_decoy"], errors="coerce").notna()
        if dmask.any():
            result.loc[dmask, "q_decoy"] = benjamini_hochberg(result.loc[dmask, "p_decoy"]).values
    result["is_primary_hit"] = eligible & (pd.to_numeric(result["q_context_null"], errors="coerce") < 0.05)
    result = result.sort_values(["primary_test_eligible", "q_context_null", "gene"], ascending=[False, True, True])

    safe = re.sub(r"[^A-Za-z0-9_-]", "_", args.cohort).strip("_") or "cohort"
    primary_out = RESULTS / f"{safe}_ptm_proximity.tsv"
    result.to_csv(primary_out, sep="\t", index=False)

    # Site-level evidence is secondary and is calculated only from actual loaded PTM anchors.
    sites = site_results(result, ptm, mut, opportunities, coords, args.permutations, args.seed)
    site_out = RESULTS / f"{safe}_ptm_sites.tsv"
    sites.to_csv(site_out, sep="\t", index=False)

    robustness_summary = []
    if args.robustness:
        grid = []
        for sep in (8, 12, 20):
            for dist in (6.0, 8.0, 10.0, 12.0):
                setting_rows = []
                for gene in sorted(panel["gene"].unique()):
                    gene_unmapped_n = int(maf_qc["unmapped_by_gene"].get(gene, 0))
                    gene_context_mismatch_n = int(maf_qc["context_mismatch_by_gene"].get(gene, 0))
                    transcript_mapping_ok = (gene_unmapped_n == 0 and gene_context_mismatch_n == 0)
                    r = setting_gene_result(gene, ptm_positions.get(gene, []), mut, opportunities, coords,
                                            sep, dist, args.permutations, args.seed, structure_mapping_ok=transcript_mapping_ok,
                                            recurrence=args.recurrence, n_decoy=0)
                    if not transcript_mapping_ok:
                        r["status"] = "partial_genomic_transcript_mapping"
                        r["primary_test_eligible"] = False
                    r["intogen_driver"] = bool(panel.loc[panel["gene"] == gene, "intogen_driver"].iloc[0])
                    setting_rows.append(r)
                sdf = pd.DataFrame(setting_rows)
                sdf["q_context_null"] = np.nan
                sdf["q_family_all_12_settings"] = np.nan
                e = sdf["primary_test_eligible"].fillna(False).astype(bool) if "primary_test_eligible" in sdf.columns else pd.Series(False, index=sdf.index)
                sdf.loc[e, "q_context_null"] = benjamini_hochberg(sdf.loc[e, "p_context_null"]).values
                grid.append((sep, dist, sdf, e))
        # One BH family over every eligible gene x setting test in this cohort (12 settings together).
        fam_p = pd.concat([g[2].loc[g[3], "p_context_null"] for g in grid], ignore_index=True)
        fam_q = benjamini_hochberg(fam_p).to_numpy() if len(fam_p) else np.array([])
        offset = 0
        for sep, dist, sdf, e in grid:
            n = int(e.sum())
            if n:
                sdf.loc[e, "q_family_all_12_settings"] = fam_q[offset:offset + n]
            offset += n
            hits = e & (pd.to_numeric(sdf["q_context_null"], errors="coerce") < 0.05)
            hits_family = e & (pd.to_numeric(sdf["q_family_all_12_settings"], errors="coerce") < 0.05)
            robust_file = RESULTS / f"{safe}_robustness_sep{sep}_dist{int(dist)}.tsv"
            sdf.to_csv(robust_file, sep="\t", index=False)
            robustness_summary.append({
                "cohort": safe, "sequence_separation_min": sep, "distance_max_A": dist,
                "is_primary_setting": bool(sep == 8 and dist == 10.0),
                "genes_in_panel": int(len(sdf)), "genes_tested": int(e.sum()),
                "genes_q_lt_005_per_setting": int(hits.sum()),
                "genes_q_lt_005_family_corrected": int(hits_family.sum()),
                "median_fold_enrichment": float(pd.to_numeric(sdf.loc[e, "fold_enrichment"], errors="coerce").median()) if e.any() else np.nan,
                "result_file": robust_file.name,
            })
        robustness_out = RESULTS / f"{safe}_robustness_summary.tsv"
        pd.DataFrame(robustness_summary).to_csv(robustness_out, sep="\t", index=False)
    else:
        robustness_out = None

    summary = {
        "cohort": safe,
        "primary_hypothesis": "Among testable PTM-bearing genes, the distribution of -log10(gene-specific SBS96-null p-values) is shifted higher for IntOGen-listed driver genes than for PTM-bearing genes not listed in the selected IntOGen release.",
        "primary_group_test": "One-sided Mann–Whitney U comparison of per-gene -log10(p_context_null) scores; the benchmark page reports the cohort-level result. Gene-level BH q-values and hit fractions are complementary gene-level summaries.",
        "primary_group_comparison": primary_group_comparison(result, panel),
        "primary_settings": {"recurrence_mode": args.recurrence, "decoy_sets_per_gene": args.decoys, "sequence_separation_min": args.sep, "CA_distance_max_angstrom": args.dist,
                             "permutations": args.permutations, "seed": args.seed, "FDR": 0.05},
        "maf_qc": maf_qc,
        "structure_qc": coord_qc,
        "ptm_records_in_analysis": int(len(ptm)),
        "ptm_bearing_genes_in_panel": int(ptm["gene"].nunique()),
        "genes_in_panel": int(panel["gene"].nunique()),
        "genes_tested_primary": int(eligible.sum()),
        "genes_q_lt_005_primary": int((eligible & (pd.to_numeric(result["q_context_null"], errors="coerce") < 0.05)).sum()),
        "primary_outcome_gene_level": ("No primary inference possible: zero genes passed all structural/transcript/context QC." if int(eligible.sum()) == 0 else ("No individual gene had q<0.05 under the primary setting." if int((eligible & (pd.to_numeric(result["q_context_null"], errors="coerce") < 0.05)).sum()) == 0 else "One or more genes passed the gene-level BH q<0.05 threshold; check the cohort-level driver/background test and independent validation.")),
        "null_model": "For each gene, randomly sample candidate missense SNV opportunities with the same observed SBS96 context counts, preserving the trinucleotide spectrum and mapping to the selected MANE/Ensembl canonical transcript.",
        "limitations": [
            "This empirical null conditions on gene-specific coding missense opportunities and observed SBS96 counts; it is not a fitted COSMIC signature mixture and does not model regional replication timing or gene-specific callable coverage.",
            "Distance is CA-CA only (AlphaFold CA coordinates); side-chain/functional-group distances are not modelled.",
            "Decoy anchors control for PTM specificity within a gene (same amino acid, same anchor count) but not for unmeasured regional mutation-rate variation.",
            "Only single-nucleotide missense variants are included; indels and multi-nucleotide variants are excluded.",
            "A gene with missing structure positions, transcript mapping gaps, or unavailable context opportunities is untestable/partial, not negative.",
            "IntOGen absence means not listed in the selected release, not a validated biological non-driver label.",
            "AlphaFold structures are predictions; use experimental structures if available and suitable, and report mapping QC."
        ],
        "primary_result_file": str(primary_out.resolve()),
        "site_level_secondary_file": str(site_out.resolve()),
        "robustness_summary_file": str(robustness_out.resolve()) if robustness_out else None,
        "input_provenance": {
            "maf_path": str(args.maf.resolve()),
            "maf_sha256": sha256_file(args.maf),
            "reference_fasta_path": str(args.fasta.resolve()),
            "reference_fasta_sha256": sha256_file(args.fasta),
            "opportunity_map_path": str(args.opportunities.resolve()),
            "opportunity_map_sha256": sha256_file(args.opportunities),
            "coordinate_map_path": str(args.structure.resolve()),
            "coordinate_map_sha256": sha256_file(args.structure),
            "ptm_input_path": str(ptm_path.resolve()),
            "ptm_input_sha256": sha256_file(ptm_path),
            "panel_path": str(args.panel.resolve()),
            "panel_sha256": sha256_file(args.panel),
        },
        "source_versions": {
            "IntOGen": "Release 2024.09.20 when that official source archive is used",
            "PTM": "dbPTM source release/data retrieval date recorded in source_manifest.json",
        },
    }
    summary_out = RESULTS / f"{safe}_analysis_summary.json"
    summary_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("PRIMARY RESULTS:", primary_out)
    print("SITE-LEVEL SECONDARY:", site_out)
    if robustness_out:
        print("ROBUSTNESS SUMMARY:", robustness_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
