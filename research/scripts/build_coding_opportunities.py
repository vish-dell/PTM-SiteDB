from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Standard genetic code; stop codons are '*'.
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def reverse_complement(s: str) -> str:
    return s.translate(COMP)[::-1].upper()


def context96(left: str, ref: str, alt: str, right: str) -> str | None:
    left, ref, alt, right = left.upper(), ref.upper(), alt.upper(), right.upper()
    if len(left) != 1 or len(ref) != 1 or len(alt) != 1 or len(right) != 1:
        return None
    if any(x not in "ACGT" for x in (left, ref, alt, right)) or ref == alt:
        return None
    if ref in "AG":
        left, ref, alt, right = right.translate(COMP), ref.translate(COMP), alt.translate(COMP), left.translate(COMP)
    return f"{left}[{ref}>{alt}]{right}"


def parse_attrs(s: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for m in re.finditer(r'([A-Za-z0-9_]+) "([^"]*)"', s):
        out[m.group(1)] = m.group(2)
    return out


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def choose_transcripts(gtf: Path, wanted: set[str], seen_out: set[str] | None = None):
    tags: dict[str, set[str]] = defaultdict(set)
    meta: dict[str, dict[str, str]] = {}
    cds: dict[str, list[dict]] = defaultdict(list)
    with open_text(gtf) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, phase, attr_s = parts
            attr = parse_attrs(attr_s)
            tx = attr.get("transcript_id")
            gene = (attr.get("gene_name") or attr.get("gene_symbol") or "").upper()
            if not tx or not gene or (wanted and gene not in wanted):
                continue
            meta.setdefault(tx, {
                "gene": gene,
                "transcript_id": tx,
                "gene_id": attr.get("gene_id", ""),
                "transcript_type": attr.get("transcript_type", attr.get("transcript_biotype", "")),
            })
            if "tag" in attr:
                tags[tx].add(attr["tag"])
            # Attributes can contain repeated tag fields; capture all occurrences explicitly.
            tags[tx].update(re.findall(r'tag "([^"]+)"', attr_s))
            if feature == "CDS":
                try:
                    cds[tx].append({"chrom": chrom, "start": int(start), "end": int(end), "strand": strand, "phase": phase})
                except ValueError:
                    pass

    if seen_out is not None:
        seen_out.update(m["gene"] for m in meta.values())
    grouped: dict[str, list[str]] = defaultdict(list)
    for tx, segments in cds.items():
        m = meta.get(tx, {})
        if not segments or not m.get("gene"):
            continue
        tx_type = (m.get("transcript_type") or "").lower()
        if tx_type and "protein_coding" not in tx_type:
            continue
        grouped[m["gene"]].append(tx)

    selected = {}
    for gene, txs in grouped.items():
        def rank(tx):
            t = tags[tx]
            length = sum(x["end"] - x["start"] + 1 for x in cds[tx])
            priority = 3 if "MANE_Select" in t else 2 if "Ensembl_canonical" in t else 1
            return (priority, length, tx)
        tx = sorted(txs, key=rank, reverse=True)[0]
        selected[gene] = {**meta[tx], "tags": sorted(tags[tx]), "segments": cds[tx]}
    return selected


def resolve_contig(fasta, name: str) -> str:
    refs = set(fasta.references)
    n = str(name)
    candidates = [n, n.removeprefix("chr"), "chr" + n.removeprefix("chr")]
    if n in {"M", "MT", "chrM", "chrMT"}:
        candidates += ["MT", "M", "chrM", "chrMT"]
    for c in candidates:
        if c in refs:
            return c
    raise KeyError(f"Contig {name!r} is absent from FASTA index (examples: {list(fasta.references[:8])})")


def write_gene_loss_report(panel: pd.DataFrame, wanted: set[str], loss: dict[str, str], path: Path,
                           retained: set[str] | None = None) -> dict:
    """Write one row per lost panel gene and print a plain accounting. Lost genes are excluded, never counted as negatives."""
    panel = panel.copy()
    panel["gene"] = panel["gene"].astype(str).str.upper()
    if "intogen_driver" in panel.columns:
        panel["_driver"] = panel["intogen_driver"].astype(str).str.lower().isin({"true", "1", "yes"})
    else:
        panel["_driver"] = False
    group = dict(zip(panel["gene"], panel["_driver"]))
    rows = [{"gene": g, "intogen_driver": group.get(g, ""), "reason": r} for g, r in sorted(loss.items()) if g in wanted]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["gene", "intogen_driver", "reason"]).to_csv(path, sep="\t", index=False)
    total = len(wanted)
    lost = len(rows)
    by_reason: dict[str, int] = defaultdict(int)
    for r in rows:
        by_reason[r["reason"]] += 1
    d_total = sum(1 for g in wanted if group.get(g) is True)
    n_total = total - d_total
    d_lost = sum(1 for r in rows if r["intogen_driver"] is True)
    n_lost = lost - d_lost
    acc = {
        "panel_genes_requested": total, "panel_genes_lost": lost, "panel_genes_retained": total - lost,
        "lost_by_reason": dict(by_reason),
        "driver_genes_retained": f"{d_total - d_lost}/{d_total}", "not_listed_genes_retained": f"{n_total - n_lost}/{n_total}",
        "note": "Lost genes are excluded from the analysis; they are not counted as negative findings.",
    }
    print(f"\nGene accounting: {total - lost}/{total} panel genes retained; {lost} lost -> {path}")
    for reason, n in sorted(by_reason.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {reason}")
    print(f"  drivers retained {acc['driver_genes_retained']}, not-listed retained {acc['not_listed_genes_retained']}")
    if total and lost / total > 0.2:
        print("WARNING: more than 20% of the panel was lost. If most are 'gene_symbol_not_found_in_gtf', the GTF probably uses older "
              "gene symbols than the panel; use a current GENCODE release for your assembly. Unequal loss between drivers and "
              "not-listed genes can bias the comparison - check the two retained counts above.")
    return acc


def main() -> int:
    ap = argparse.ArgumentParser(description="Enumerate every possible missense single-nucleotide opportunity on one transcript per target gene.")
    ap.add_argument("--gtf", required=True, type=Path, help="GTF/GTF.GZ matching the FASTA assembly")
    ap.add_argument("--fasta", required=True, type=Path, help="Reference FASTA with .fai index")
    ap.add_argument("--panel", type=Path, default=Path("data/processed/benchmark_panel.csv"))
    ap.add_argument("--out", type=Path, default=Path("data/processed/context_opportunities.tsv"))
    ap.add_argument("--proteins-out", type=Path, default=Path("data/processed/canonical_proteins.fasta"))
    ap.add_argument("--manifest-out", type=Path, default=Path("data/processed/opportunity_manifest.json"))
    ap.add_argument("--loss-out", type=Path, default=None, help="Gene-loss report (default: gene_loss_report.tsv next to --out)")
    ap.add_argument("--max-genes", type=int, default=0, help="Optional development cap; 0 means all panel genes")
    args = ap.parse_args()

    if not Path(str(args.fasta) + ".fai").exists():
        raise SystemExit(f"FASTA index required. IndexedFasta can create it; check write access for {args.fasta} and rerun PREPARE_GENOMIC_INPUTS.ps1.")
    if not args.panel.exists():
        raise SystemExit(f"Panel not found: {args.panel}; run scripts/build_200gene_panel.py first.")
    try:
        from fasta_utils import IndexedFasta
    except ImportError:
        try:
            from scripts.fasta_utils import IndexedFasta
        except ImportError:
            try:
                from research.scripts.fasta_utils import IndexedFasta
            except ImportError as exc:
                raise SystemExit("Included portable FASTA indexer could not be imported; verify the project folder is intact.") from exc

    panel = pd.read_csv(args.panel)
    if "gene" not in panel.columns:
        raise SystemExit("Panel must have a gene column")
    wanted = set(panel["gene"].astype(str).str.upper())
    if args.max_genes > 0:
        wanted = set(sorted(wanted)[:args.max_genes])
    seen_genes: set[str] = set()
    txs = choose_transcripts(args.gtf, wanted, seen_out=seen_genes)
    # Every panel gene is accounted for: retained, or listed with the reason it was lost.
    loss: dict[str, str] = {}
    for g in sorted(wanted):
        if g not in seen_genes:
            loss[g] = "gene_symbol_not_found_in_gtf"
        elif g not in txs:
            loss[g] = "no_protein_coding_transcript_with_cds"
    loss_path = args.loss_out or (args.out.parent / "gene_loss_report.tsv")
    if not txs:
        write_gene_loss_report(panel, wanted, loss, loss_path)
        raise SystemExit("No protein-coding panel transcripts could be selected. Confirm GTF assembly and gene symbols. "
                         f"See {loss_path} for the per-gene reasons.")

    fasta = IndexedFasta(str(args.fasta))
    out_rows = []
    protein_path = args.proteins_out
    protein_path.parent.mkdir(parents=True, exist_ok=True)
    skipped = defaultdict(int)
    with protein_path.open("w", encoding="utf-8") as pf:
        for ix, (gene, txinfo) in enumerate(sorted(txs.items()), 1):
            segments = txinfo["segments"]
            strand = segments[0]["strand"]
            if any(s["strand"] != strand for s in segments):
                skipped["inconsistent_strand"] += 1
                loss[gene] = "inconsistent_strand"
                continue
            ordered = sorted(segments, key=lambda s: s["start"], reverse=(strand == "-"))
            seq_parts = []
            genomic_positions = []
            chroms = []
            try:
                for seg in ordered:
                    chrom = resolve_contig(fasta, seg["chrom"])
                    raw = fasta.fetch(chrom, seg["start"] - 1, seg["end"]).upper()
                    coords = list(range(seg["start"], seg["end"] + 1))
                    if strand == "-":
                        raw = reverse_complement(raw)
                        coords.reverse()
                    seq_parts.append(raw)
                    genomic_positions.extend(coords)
                    chroms.extend([chrom] * len(coords))
            except Exception as exc:
                skipped["unresolvable_or_fetch_error"] += 1
                loss[gene] = "unresolvable_or_fetch_error"
                print(f"SKIP {gene}: {exc}")
                continue
            cds_seq = "".join(seq_parts)
            if len(cds_seq) < 3 or len(cds_seq) % 3 != 0:
                skipped["non_triplet_complete_cds"] += 1
                loss[gene] = "non_triplet_complete_cds"
                continue
            aa = []
            valid = True
            for i in range(0, len(cds_seq), 3):
                codon = cds_seq[i:i+3]
                if any(x not in "ACGT" for x in codon) or codon not in CODON_TABLE:
                    valid = False
                    break
                aa.append(CODON_TABLE[codon])
            if not valid:
                skipped["ambiguous_cds"] += 1
                loss[gene] = "ambiguous_cds"
                continue
            protein = "".join(aa)
            # GTF CDS may include a terminal stop codon; remove it from the protein FASTA.
            protein_for_fasta = protein[:-1] if protein.endswith("*") else protein
            protein_id = f"{gene}|{txinfo['transcript_id']}"
            pf.write(f">{protein_id}\n")
            for j in range(0, len(protein_for_fasta), 60):
                pf.write(protein_for_fasta[j:j+60] + "\n")

            for nt_idx, genome_pos in enumerate(genomic_positions):
                aa_pos0 = nt_idx // 3
                if aa_pos0 >= len(protein) or protein[aa_pos0] == "*":
                    continue
                codon_start = aa_pos0 * 3
                local_idx = nt_idx - codon_start
                codon = cds_seq[codon_start:codon_start+3]
                if len(codon) != 3 or any(x not in "ACGT" for x in codon):
                    continue
                chrom = chroms[nt_idx]
                try:
                    triplet = fasta.fetch(chrom, genome_pos - 2, genome_pos + 1).upper()
                except Exception:
                    continue
                if len(triplet) != 3 or any(x not in "ACGT" for x in triplet):
                    continue
                genomic_ref = triplet[1]
                cds_ref = cds_seq[nt_idx]
                # Strand consistency check; reference in the CDS should match genome orientation.
                expected_cds_ref = genomic_ref if strand == "+" else genomic_ref.translate(COMP)
                if cds_ref != expected_cds_ref:
                    skipped["strand_reference_mismatch"] += 1
                    continue
                for genomic_alt in "ACGT":
                    if genomic_alt == genomic_ref:
                        continue
                    cds_alt = genomic_alt if strand == "+" else genomic_alt.translate(COMP)
                    alt_codon = codon[:local_idx] + cds_alt + codon[local_idx+1:]
                    alt_aa = CODON_TABLE.get(alt_codon)
                    ref_aa = protein[aa_pos0]
                    if not alt_aa or alt_aa == "*" or alt_aa == ref_aa:
                        continue
                    ctx = context96(triplet[0], genomic_ref, genomic_alt, triplet[2])
                    if ctx is None:
                        continue
                    out_rows.append({
                        "gene": gene,
                        "protein_pos": aa_pos0 + 1,
                        "context96": ctx,
                        "chromosome": chrom,
                        "genomic_position": genome_pos,
                        "ref": genomic_ref,
                        "alt": genomic_alt,
                        "transcript_id": txinfo["transcript_id"],
                        "transcript_selection": "MANE_Select" if "MANE_Select" in txinfo["tags"] else "Ensembl_canonical" if "Ensembl_canonical" in txinfo["tags"] else "longest_CDS_fallback",
                        "ref_aa": ref_aa,
                        "alt_aa": alt_aa,
                    })
            if ix % 50 == 0:
                print(f"Processed transcript {ix}/{len(txs)}")

    fasta.close()
    out = pd.DataFrame(out_rows)
    if out.empty:
        raise SystemExit("No missense SNV opportunities generated; verify GTF/FASTA assembly and panel genes.")
    out.to_csv(args.out, sep="\t", index=False)
    retained = set(out["gene"].astype(str))
    for g in txs:
        if g not in retained and g not in loss:
            loss[g] = "no_missense_opportunities_generated"
    accounting = write_gene_loss_report(panel, wanted, loss, loss_path, retained=retained)
    manifest = {
        "gtf": str(args.gtf.resolve()),
        "reference_fasta": str(args.fasta.resolve()),
        "panel_file": str(args.panel.resolve()),
        "transcripts_selected": int(len(txs)),
        "genes_with_opportunities": int(out["gene"].nunique()),
        "missense_snv_opportunities": int(len(out)),
        "canonical_protein_fasta": str(protein_path.resolve()),
        "opportunity_file": str(args.out.resolve()),
        "transcript_priority": ["MANE_Select", "Ensembl_canonical", "longest CDS fallback"],
        "skipped_reasons": dict(skipped),
        "gene_accounting": accounting,
        "gene_loss_report": str(loss_path.resolve()),
        "context_definition": "SBS96, normalized to C/T reference substitutions",
        "notes": [
            "Each row is one possible coding single-nucleotide substitution that changes an amino acid without creating a stop codon.",
            "This map is based on one selected transcript per gene and is only valid for the exact FASTA/GTF assembly pair.",
            "Indels and multi-nucleotide substitutions are excluded from the SBS96 missense analysis."
        ],
    }
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in manifest.items() if k != "skipped_reasons"}, indent=2))
    print("Skip counters:", dict(skipped))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
