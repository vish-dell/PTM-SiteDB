from __future__ import annotations

import argparse
import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import requests
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def read_fasta(path: Path) -> dict[str, str]:
    out = {}
    key = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if key is not None:
                out[key] = "".join(chunks).upper()
            key = line[1:].split()[0]
            chunks = []
        elif key is not None:
            chunks.append(line.strip())
    if key is not None:
        out[key] = "".join(chunks).upper()
    return out


def parse_fasta_sequence(text: str) -> str:
    return "".join(x.strip() for x in text.splitlines() if x and not x.startswith(">")).upper()


def uniprot_sequence(session: requests.Session, accession: str) -> str:
    r = session.get(UNIPROT_FASTA.format(accession=accession), timeout=45)
    r.raise_for_status()
    seq = parse_fasta_sequence(r.text)
    if not seq:
        raise RuntimeError(f"No UniProt sequence returned for {accession}")
    return seq


def fetch_pdb_url(session: requests.Session, accession: str) -> tuple[str, str]:
    r = session.get(AF_API + accession, timeout=30)
    r.raise_for_status()
    body = r.json()
    models = body if isinstance(body, list) else [body]
    for model in models:
        if model.get("pdbUrl"):
            return model["pdbUrl"], "AlphaFold DB API"
    for version in ("v6", "v5", "v4"):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_{version}.pdb"
        probe = session.get(url, timeout=45)
        if probe.ok and probe.text.startswith(("HEADER", "TITLE", "REMARK", "ATOM")):
            return url, f"AlphaFold DB fallback {version}"
    raise RuntimeError(f"No AlphaFold PDB URL returned for {accession}")


def parse_pdb_ca(text: str) -> dict[int, dict]:
    # AlphaFold PDB B-factor column is pLDDT. Keep one CA per residue.
    out = {}
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or line[12:16].strip() != "CA":
            continue
        altloc = line[16:17]
        if altloc not in (" ", "A"):
            continue
        try:
            chain = line[21:22].strip() or "A"
            if chain != "A" and any(x.get("chain") == "A" for x in out.values()):
                continue
            resno = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            b = float(line[60:66]) if line[60:66].strip() else float("nan")
            aa3 = line[17:20].strip().upper()
            aa = seq1(aa3, custom_map={"UNK": "X"})
        except Exception:
            continue
        if resno not in out or altloc == "A":
            out[resno] = {"x": x, "y": y, "z": z, "plddt": b, "aa": aa, "chain": chain}
    return out


def alignment_map(canonical: str, structure_seq: str) -> tuple[dict[int, int], float, float]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5
    aln = aligner.align(canonical, structure_seq)[0]
    coords = aln.coordinates
    mapping = {}
    exact = aligned = 0
    for i in range(coords.shape[1] - 1):
        a0, a1 = int(coords[0, i]), int(coords[0, i+1])
        b0, b1 = int(coords[1, i]), int(coords[1, i+1])
        da, db = a1-a0, b1-b0
        if da <= 0 or db <= 0:
            continue
        n = min(da, db)
        for k in range(n):
            ca_idx, st_idx = a0+k, b0+k
            aligned += 1
            if canonical[ca_idx] == structure_seq[st_idx]:
                exact += 1
                mapping[ca_idx+1] = st_idx+1
    identity = exact / aligned if aligned else 0.0
    coverage = len(mapping) / len(canonical) if canonical else 0.0
    return mapping, identity, coverage


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch AlphaFold CA coordinates and map them to selected MANE/Ensembl transcript protein numbering.")
    ap.add_argument("--proteins", type=Path, default=PROC / "canonical_proteins.fasta")
    ap.add_argument("--ptm", type=Path, default=PROC / "ptm_sites.tsv")
    ap.add_argument("--coordinates-out", type=Path, default=PROC / "alphafold_ca_coordinates.tsv")
    ap.add_argument("--ptm-out", type=Path, default=PROC / "ptm_sites_mapped.tsv")
    ap.add_argument("--manifest-out", type=Path, default=PROC / "structure_manifest.json")
    ap.add_argument("--max-genes", type=int, default=0, help="Optional development cap; 0 means all genes")
    ap.add_argument("--sleep", type=float, default=0.05)
    args = ap.parse_args()

    if not args.proteins.exists() or not args.ptm.exists():
        raise SystemExit("Run scripts/build_coding_opportunities.py first and ensure real PTM data is loaded.")
    proteins = read_fasta(args.proteins)
    protein_by_gene = {key.split("|", 1)[0].upper(): seq for key, seq in proteins.items() if "|" in key}
    ptm = pd.read_csv(args.ptm, sep="\t", dtype={"gene": str, "uniprot": str})
    needed = {"gene", "uniprot", "position", "ptm_type"}
    if not needed.issubset(ptm.columns):
        raise SystemExit(f"PTM file needs columns: {sorted(needed)}")
    ptm["gene"] = ptm["gene"].astype(str).str.upper()
    ptm["uniprot"] = ptm["uniprot"].astype(str).str.strip()
    ptm["position"] = pd.to_numeric(ptm["position"], errors="coerce")
    ptm = ptm.dropna(subset=["position"])
    ptm["position"] = ptm["position"].astype(int)

    sessions = requests.Session()
    sessions.headers.update({"User-Agent": "PTM-SiteDB-Integrated-Explorer/3.0 (research-use)"})
    coords_rows, mapped_ptm_rows, manifests = [], [], {}
    accession_counts = defaultdict(Counter)
    for r in ptm.itertuples(index=False):
        # Prefer primary UniProt accessions so the selected AF model is for the same accession.
        if "-" not in str(r.uniprot) and str(r.uniprot) not in {"", "nan", "None"}:
            accession_counts[r.gene][r.uniprot] += 1
    genes = sorted(set(protein_by_gene) & set(accession_counts))
    if args.max_genes > 0:
        genes = genes[:args.max_genes]

    for ix, gene in enumerate(genes, 1):
        canonical = protein_by_gene[gene]
        gene_ptm = ptm[ptm["gene"] == gene]
        candidate_accessions = [a for a, _ in accession_counts[gene].most_common()]
        success = False
        errors = []
        for acc in candidate_accessions:
            try:
                uni_seq = uniprot_sequence(sessions, acc)
                pdb_url, api_source = fetch_pdb_url(sessions, acc)
                pdb = sessions.get(pdb_url, timeout=90)
                pdb.raise_for_status()
                ca = parse_pdb_ca(pdb.text)
                structure_resnos = sorted(ca)
                structure_seq = "".join(ca[k]["aa"] for k in structure_resnos)
                if not ca or not structure_seq:
                    raise RuntimeError("PDB did not contain CA coordinates")
                # Coordinate index in the PDB should follow its residue order. Align the selected UniProt sequence to the model sequence;
                # only identical amino acids are mapped into canonical transcript numbering below.
                uni_to_struct, uni_identity, uni_coverage = alignment_map(uni_seq, structure_seq)
                can_to_uni, identity, coverage = alignment_map(canonical, uni_seq)
                if identity < 0.80 or coverage < 0.75 or uni_identity < 0.90 or uni_coverage < 0.75:
                    raise RuntimeError(f"sequence mapping QC failed (canonical identity={identity:.3f}, coverage={coverage:.3f}; model identity={uni_identity:.3f}, coverage={uni_coverage:.3f})")

                canonical_to_struct = {}
                for can_pos, uni_pos in can_to_uni.items():
                    st_seq_pos = uni_to_struct.get(uni_pos)
                    st_pos = structure_resnos[st_seq_pos-1] if st_seq_pos and 1 <= st_seq_pos <= len(structure_resnos) else None
                    if st_pos and st_pos in ca and canonical[can_pos-1] == uni_seq[uni_pos-1] == ca[st_pos]["aa"]:
                        canonical_to_struct[can_pos] = st_pos
                        cr = ca[st_pos]
                        coords_rows.append({
                            "gene": gene, "position": can_pos,
                            "x": cr["x"], "y": cr["y"], "z": cr["z"], "plddt": cr["plddt"],
                            "uniprot": acc, "structure_residue_number": st_pos,
                            "coordinate_source": api_source, "sequence_identity": identity,
                        })

                # Keep only PTM sites reported on the exact model UniProt accession, then map them by aligned identity.
                # Site records on other isoforms/accessions are counted in the QC manifest and are not silently projected.
                same_acc = gene_ptm[gene_ptm["uniprot"] == acc]
                gene_site_rows = []
                uni_to_can = {uni_pos: can_pos for can_pos, uni_pos in can_to_uni.items() if canonical[can_pos-1] == uni_seq[uni_pos-1]}
                for r in same_acc.itertuples(index=False):
                    can_pos = uni_to_can.get(int(r.position))
                    if can_pos:
                        gene_site_rows.append({
                            "gene": gene, "uniprot": acc, "position": can_pos,
                            "source_position": int(r.position), "ptm_type": r.ptm_type,
                            "flanking_sequence": getattr(r, "flanking_sequence", ""),
                            "source": getattr(r, "source", "dbPTM"),
                            "mapping_note": "sequence-mapped by exact amino-acid identity; CA availability is checked separately by primary QC",
                        })
                mapped_ptm_rows.extend(gene_site_rows)
                manifests[gene] = {
                    "status": "mapped", "uniprot": acc, "coordinate_rows": len(canonical_to_struct),
                    "source_ptm_records_for_gene": int(len(gene_ptm)),
                    "ptm_records_on_selected_accession": int(len(same_acc)),
                    "records_on_other_accessions_not_projected": int(len(gene_ptm) - len(same_acc)),
                    "sequence_mapped_ptm_records": int(len(gene_site_rows)),
                    "sequence_mapping_fraction_selected_accession": float(len(gene_site_rows) / len(same_acc)) if len(same_acc) else None,
                    "ptm_sites_without_exact_sequence_mapping": int(len(same_acc) - len(gene_site_rows)),
                    "canonical_to_uniprot_identity": identity, "canonical_coverage": coverage,
                    "uniprot_to_model_identity": uni_identity, "model_coverage": uni_coverage,
                    "pdb_url": pdb_url,
                }
                success = True
                break
            except Exception as exc:
                errors.append(f"{acc}: {exc}")
        if not success:
            manifests[gene] = {
                "status": "unavailable_or_mapping_failed",
                "source_ptm_records_for_gene": int(len(gene_ptm)),
                "ptm_accessions": sorted(set(gene_ptm["uniprot"].astype(str))),
                "errors": errors[-4:],
            }
        if ix % 25 == 0 or ix == len(genes):
            print(f"Genes attempted {ix}/{len(genes)}; mapped={sum(x.get('status') == 'mapped' for x in manifests.values())}")
        if args.sleep:
            time.sleep(args.sleep)

    coords_df = pd.DataFrame(coords_rows)
    ptm_df = pd.DataFrame(mapped_ptm_rows)
    args.coordinates_out.parent.mkdir(parents=True, exist_ok=True)
    coords_df.to_csv(args.coordinates_out, sep="\t", index=False)
    ptm_df.to_csv(args.ptm_out, sep="\t", index=False)
    mapping_qc = pd.DataFrame([{"gene": gene, **vals} for gene, vals in sorted(manifests.items())])
    mapping_qc_path = PROC / "ptm_site_mapping_qc.tsv"
    mapping_qc_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_qc.to_csv(mapping_qc_path, sep="\t", index=False)
    manifest = {
        "genes_attempted": len(genes), "genes_mapped": sum(x.get("status") == "mapped" for x in manifests.values()),
        "mapping_qc_file": str(mapping_qc_path.resolve()),
        "coordinate_rows": int(len(coords_df)), "mapped_experimental_ptm_records": int(len(ptm_df)),
        "coordinate_file": str(args.coordinates_out.resolve()), "mapped_ptm_file": str(args.ptm_out.resolve()),
        "gene_results": manifests,
        "important_limitations": [
            "AlphaFold coordinates are predictions, not experimental structures.",
            "Residue coordinates are mapped only where amino acids align identically between the selected transcript, UniProt accession and AlphaFold model.",
            "PTM sites from non-selected UniProt accessions/isoforms are not silently merged into the selected canonical numbering.",
            "A missing coordinate or unmapped PTM makes that position untestable, not a negative biological result."
        ]
    }
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in manifest.items() if k != "gene_results"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
