"""Compare project results with a user-supplied, cited benchmark CSV and compute real agreement metrics.

No benchmark dataset is bundled: reference calls must come from a real, cited source.
Benchmark columns: gene, uniprot (optional), site (optional), reference_result, reference_source.
  reference_result must be 'positive' (expected enriched) or 'negative' (expected not enriched).
Project input: a proximity (protein-level) or sites (site-level) CSV from batch_proximity.py.
A benchmark row is matched on (gene, uniprot, site); blank benchmark fields are ignored when matching.

Outputs:
  <out>              per-row comparison (our call vs reference call; 'not_tested' is never counted as negative)
  <out>.summary.csv  TP/FP/FN/TN, precision, recall, specificity over rows that were actually tested
"""
import argparse, csv
from pathlib import Path


def read(path):
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def num(v):
    try:
        return float(v) if v not in (None, '', 'None') else None
    except (TypeError, ValueError):
        return None


def project_q(row, metric):
    cols = {'perm': ('q_perm_all_sites', 'q_perm'), 'exact': ('q_bh_all_sites', 'q_bh')}[metric]
    for c in cols:
        if num(row.get(c)) is not None:
            return num(row[c])
    return None


def find_match(project_rows, b):
    g, u, s = (b.get(k, '').strip() for k in ('gene', 'uniprot', 'site'))
    hits = [r for r in project_rows
            if r.get('gene', '').strip() == g
            and (not u or r.get('uniprot', '').strip() == u)
            and (not s or r.get('site', '').strip() == s)]
    return hits[0] if hits else None


def compare(project_rows, benchmark_rows, q_cutoff=0.05, metric='perm'):
    rows = []
    for b in benchmark_rows:
        ref = (b.get('reference_result') or '').strip().lower()
        p = find_match(project_rows, b)
        q = project_q(p, metric) if p else None
        if ref not in ('positive', 'negative'):
            call, outcome = ('not_tested' if q is None else ('positive' if q <= q_cutoff else 'negative')), 'invalid_reference'
        elif p is None or q is None:
            call, outcome = 'not_tested', 'not_comparable'
        else:
            call = 'positive' if q <= q_cutoff else 'negative'
            outcome = {('positive', 'positive'): 'TP', ('positive', 'negative'): 'FN',
                       ('negative', 'positive'): 'FP', ('negative', 'negative'): 'TN'}[(ref, call)]
        rows.append({'gene': b.get('gene', ''), 'uniprot': b.get('uniprot', ''), 'site': b.get('site', ''),
                     'our_q': '' if q is None else q, 'our_call': call, 'reference_result': ref,
                     'reference_source': b.get('reference_source', ''), 'outcome': outcome})
    return rows


def summarize(rows):
    c = {k: sum(r['outcome'] == k for r in rows) for k in ('TP', 'FP', 'FN', 'TN')}
    div = lambda a, b: round(a / b, 4) if b else None
    return {**c, 'n_not_comparable': sum(r['outcome'] == 'not_comparable' for r in rows),
            'n_tested': sum(r['outcome'] in ('TP','FP','FN','TN') for r in rows),
            'n_invalid_reference': sum(r['outcome'] == 'invalid_reference' for r in rows),
            'precision': div(c['TP'], c['TP'] + c['FP']), 'recall': div(c['TP'], c['TP'] + c['FN']),
            'specificity': div(c['TN'], c['TN'] + c['FP'])}


def write_csv(path, rows, cols):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('project_csv'); ap.add_argument('benchmark_csv')
    ap.add_argument('--out', default='results/benchmark_comparison.csv')
    ap.add_argument('--q-cutoff', type=float, default=0.05)
    ap.add_argument('--metric', choices=['perm', 'exact'], default='perm',
                    help='Which q-value defines our positive call (default: composition-matched permutation).')
    a = ap.parse_args()
    rows = compare(read(a.project_csv), read(a.benchmark_csv), a.q_cutoff, a.metric)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    write_csv(out, rows, ['gene', 'uniprot', 'site', 'our_q', 'our_call', 'reference_result', 'reference_source', 'outcome'])
    summ = summarize(rows)
    write_csv(out.with_suffix('.summary.csv'), [summ], list(summ))
    print(f'Wrote {out} and its summary; {summ}. Untested rows are never counted as negatives.')


if __name__ == '__main__':
    main()
