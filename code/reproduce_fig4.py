from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from result_paths import reference_dir
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import t as student_t
import reproduce_fig1_fig2 as base
SSW_YEARS = np.array([1949, 1951, 1952, 1954, 1956, 1957, 1959, 1962, 1965, 1967, 1968, 1969, 1970, 1972, 1976, 1978, 1979, 1980, 1983, 1984, 1986, 1987, 1988, 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009, 2012, 2017, 2018, 2020, 2022], dtype=int)
MODE_NAMES = [f'S{i}' for i in range(2, 7)] + [f'T{i}' for i in range(1, 8)]

def load_inputs(folder: Path):
    sn = np.load(reference_dir(folder) / 'S_basis_and_correlations.npz')
    tn = np.load(folder / 'T_basis_and_correlations.npz')
    np.testing.assert_array_equal(sn['dates'], base.DATES)
    np.testing.assert_array_equal(tn['dates'], base.DATES)
    s = np.asarray(sn['scores'][:, :6], dtype=float)
    t = np.asarray(tn['scores'][:, :7], dtype=float)
    sden = np.sqrt(np.sum(s[base.BASIS_MASK] ** 2, axis=0))
    tden = np.sqrt(np.sum(t[base.BASIS_MASK] ** 2, axis=0))
    (s, t) = (s / sden, t / tden)
    y = np.full(len(base.DATES), np.nan)
    with open(reference_dir(folder) / 'spv_daily.csv', encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            i = int((np.datetime64(row['date']) - base.START) / np.timedelta64(1, 'D'))
            y[i] = float(row['spv_standardized'])
    return (s, t, y)

def design_at_h(s, t, y, h):
    dst = np.flatnonzero(base.TARGET_MASK & np.isfinite(y))
    start = dst - h
    hist = start[:, None] - np.arange(60)[None, :]
    inside = (hist >= 0).all(axis=1)
    inside &= base.WINTER[hist[:, -1]] == base.WINTER[dst]
    inside &= np.isfinite(s[hist, 0]).all(axis=1)
    inside &= np.isfinite(s[start, 1:]).all(axis=1)
    inside &= np.isfinite(t[start, :]).all(axis=1)
    (dst, start, hist) = (dst[inside], start[inside], hist[inside])
    X = np.column_stack((s[hist, 0], s[start, 1:6], t[start, 0:7]))
    return (X, y[dst], base.WINTER[dst])

def ols_cr1(X, y, winter):
    (n, k) = X.shape
    p = k + 1
    groups = np.unique(winter)
    G = len(groups)
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=0)
    if np.any(sd == 0):
        raise ValueError('constant design column')
    Z = (X - mu) / sd
    yy = y - y.mean()
    gram = Z.T @ Z
    cf = cho_factor(gram, check_finite=True)
    gamma = cho_solve(cf, Z.T @ yy)
    beta = gamma / sd
    intercept = float(y.mean() - mu @ beta)
    e = yy - Z @ gamma
    inv = cho_solve(cf, np.eye(k))
    A = np.zeros((p, p), dtype=float)
    A[0, 0] = 1 / n
    A[1:, 1:] = inv
    meat = np.zeros((p, p), dtype=float)
    for g in groups:
        take = winter == g
        score = np.r_[np.sum(e[take]), Z[take].T @ e[take]]
        meat += np.outer(score, score)
    V = G / (G - 1) * ((n - 1) / (n - p)) * (A @ meat @ A)
    se = np.sqrt(np.diag(V)[1:]) / sd
    pv = 2 * student_t.sf(np.abs(beta / se), G - 1)
    return (intercept, beta, se, pv, np.sqrt(np.linalg.cond(gram)), k)

def plot_fig4(output: Path, beta, pv):
    from publication_style import plot_fig4 as draw_publication
    return draw_publication(output, beta, pv)

def main():
    pa = argparse.ArgumentParser(description='Reproduce report Fig. 4 from the All75 source-field EMT/SPV outputs.\n\nFirst run reproduce_fig1_fig2.py; this script reads its correct 241-day basis,\nthe signed 7-Apr projections, and standardized SPV. It fits All/SSW/non-SSW\nfull-sample descriptive regressions, not leave-one-winter predictions.\n')
    pa.add_argument('--input', type=Path, default=Path(__file__).with_name('fig1_fig2_results'))
    pa.add_argument('--output', type=Path, default=Path(__file__).with_name('fig3_fig4_results'))
    args = pa.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (s, t, y) = load_inputs(args.input)
    beta = np.full((3, 60, 12), np.nan)
    se = beta.copy()
    pv = beta.copy()
    ns = np.zeros((3, 60), int)
    conditions = np.zeros((3, 60))
    ranks = np.zeros((3, 60), int)
    for h in range(1, 61):
        (X, response, winters) = design_at_h(s, t, y, h)
        for (j, mask) in enumerate((np.ones(len(winters), dtype=bool), np.isin(winters, SSW_YEARS), ~np.isin(winters, SSW_YEARS))):
            (xj, yj, wj) = (X[mask], response[mask], winters[mask])
            (intercept, bj, sej, pj, cond, rank) = ols_cr1(xj, yj, wj)
            beta[j, h - 1] = bj[60:]
            se[j, h - 1] = sej[60:]
            pv[j, h - 1] = pj[60:]
            ns[j, h - 1] = len(yj)
            conditions[j, h - 1] = cond
            ranks[j, h - 1] = rank
        if h % 10 == 0:
            print(f'Fig4: h={h}/60', flush=True)
    np.savez_compressed(args.output / 'fig4_coefficients_CR1.npz', beta=beta, se=se, p=pv, n=ns, standardized_condition=conditions, slope_rank=ranks, h=np.arange(1, 61), names=np.array(MODE_NAMES), groups=np.array(['All', 'SSW', 'non-SSW']))
    with open(args.output / 'fig4_coefficients_CR1.csv', 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(['group', 'h', 'mode', 'b', 'SE_CR1', 'p_two_sided', 'n'])
        for (j, group) in enumerate(('All', 'SSW', 'non-SSW')):
            for h in range(60):
                for (k, mode) in enumerate(MODE_NAMES):
                    w.writerow([group, h + 1, mode, beta[j, h, k], se[j, h, k], pv[j, h, k], ns[j, h]])
    plot_fig4(args.output / 'Fig4', beta, pv)
    qa = {'all_h7_S6_beta': float(beta[0, 6, 4]), 'all_h7_S6_p': float(pv[0, 6, 4]), 'max_abs_beta': np.max(abs(beta), axis=(1, 2)).tolist(), 'p_le_0p05_cells': np.sum(pv <= 0.05, axis=(1, 2)).tolist(), 'min_p': np.min(pv, axis=(1, 2)).tolist(), 'max_standardized_condition': float(np.max(conditions)), 'all_slope_ranks': np.unique(ranks).tolist(), 'n_h1': ns[:, 0].tolist(), 'n_h7': ns[:, 6].tolist(), 'n_h15': ns[:, 14].tolist(), 'n_h30': ns[:, 29].tolist(), 'n_h60': ns[:, 59].tolist()}
    (args.output / 'fig4_QA.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
    print(json.dumps(qa, indent=2), flush=True)
if __name__ == '__main__':
    main()
