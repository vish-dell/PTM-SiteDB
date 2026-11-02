"""Summarize distance/sequence-separation sensitivity from real batch CSV outputs.

Example:
  python sensitivity_report.py --study brca_tcga --seps 8 12 20 --dists 6 8 10 12

The script only reads existing CSV results. It never creates biological values.
"""
import argparse, csv
from pathlib import Path


def read(path):
    if not path.exists():
        return []
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def num(r, k):
    try:
        return float(r.get(k)) if r.get(k) not in ('', None, 'None') else None
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', required=True)
    ap.add_argument('--seps', type=int, nargs='+', default=[8, 12, 20])
    ap.add_argument('--dists', type=float, nargs='+', default=[6, 8, 10, 12])
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    rows = []
    for sep in args.seps:
        for dist in args.dists:
            # batch_proximity.py writes one file per (study, sep, dist); the distance is in the name.
            candidates = [Path('results') / f'proximity_{args.study}_sep{sep}_dist{dist:g}.csv']
            path = next((p for p in candidates if p.exists()), None)
            if path is None:
                rows.append({'study': args.study, 'sep': sep, 'dist': dist, 'status': 'missing_results'})
                continue
            data = read(path)
            tested = [r for r in data if num(r, 'p_value') is not None]
            perm = [r for r in data if num(r, 'p_perm') is not None]
            sig_exact = [r for r in data if (num(r, 'q_bh') is not None and num(r, 'q_bh') <= .05)]
            sig_perm = [r for r in data if (num(r, 'q_perm') is not None and num(r, 'q_perm') <= .05)]
            rows.append({
                'study': args.study, 'sep': sep, 'dist': dist,
                'status': 'ok', 'proteins_tested': len(tested),
                'permutation_tested': len(perm),
                'exact_q05': len(sig_exact), 'permutation_q05': len(sig_perm),
                'exact_genes': ';'.join(sorted({r.get('gene','') for r in sig_exact})),
                'permutation_genes': ';'.join(sorted({r.get('gene','') for r in sig_perm})),
                'source_file': str(path),
            })
    out = Path(args.out or f'results/sensitivity_{args.study}.csv')
    out.parent.mkdir(exist_ok=True)
    cols = ['study','sep','dist','status','proteins_tested','permutation_tested','exact_q05','permutation_q05','exact_genes','permutation_genes','source_file']
    with out.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f'Wrote {out}; rows are summaries of existing real CSV results only.')

if __name__ == '__main__':
    main()
