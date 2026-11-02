"""Build a discovery/replication summary from REAL batch CSV outputs.

Usage:
  python replication_report.py results/proximity_brca_tcga_sep8_dist10.csv results/proximity_luad_tcga_sep8_dist10.csv --out results/replication_sep8_dist10.csv

This script never invents a hit. A candidate is labelled replicated only when the same
(gene, UniProt) row is present and passes the user-specified criteria in BOTH files.
"""
import argparse, csv
from pathlib import Path


def read(path):
    with open(path, newline='', encoding='utf-8') as f:
        return {(r.get('gene',''), r.get('uniprot','')): r for r in csv.DictReader(f)}


def f(r, k):
    try: return float(r[k]) if r.get(k) not in (None, '', 'None') else None
    except (TypeError, ValueError): return None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('discovery')
    ap.add_argument('replication')
    ap.add_argument('--out', default='results/replication_summary.csv')
    ap.add_argument('--q-cutoff', type=float, default=0.05)
    ap.add_argument('--require-permutation', action='store_true')
    a=ap.parse_args()
    d=read(a.discovery); r=read(a.replication)
    rows=[]
    for key in sorted(set(d)|set(r)):
        dr, rr=d.get(key), r.get(key)
        dq=f(dr,'q_perm') if dr else None
        rq=f(rr,'q_perm') if rr else None
        de=f(dr,'q_bh') if dr else None
        re=f(rr,'q_bh') if rr else None
        dpass=(dq is not None and dq <= a.q_cutoff) if a.require_permutation else (de is not None and de <= a.q_cutoff)
        rpass=(rq is not None and rq <= a.q_cutoff) if a.require_permutation else (re is not None and re <= a.q_cutoff)
        if dr is None: status='replication_only'
        elif rr is None: status='discovery_only'
        elif dpass and rpass: status='replicated'
        else: status='shared_not_replicated'
        rows.append({'gene':key[0],'uniprot':key[1],'status':status,
                     'discovery_q_bh':de,'replication_q_bh':re,
                     'discovery_q_perm':dq,'replication_q_perm':rq,
                     'discovery_observed':dr.get('observed') if dr else None,
                     'replication_observed':rr.get('observed') if rr else None})
    out=Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ['gene','uniprot','status'])
        w.writeheader(); w.writerows(rows)
    print(f'Wrote {out} from real CSV inputs; no biological values were generated.')

if __name__=='__main__': main()
