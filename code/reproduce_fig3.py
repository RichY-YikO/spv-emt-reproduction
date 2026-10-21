from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_name, '4')
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve
import reproduce_fig1_fig2 as base
from reproduce_fig4 import SSW_YEARS

def all75_anomalies(directory: Path, level: str):
    raw = base.read_annual_fields(directory, level)
    (anomaly, clim) = base.anomaly_and_climatology(raw)
    del raw
    valid_count = np.array([np.isfinite(anomaly[base.CAL_KEY == k, 0]).sum() for k in range(366)])
    return (anomaly, clim, valid_count)

def leaveout_delta(anomaly: np.ndarray, valid_count: np.ndarray, w: int):
    delta = np.zeros((366, 1584), dtype=np.float64)
    idx = np.flatnonzero(base.WINTER == w)
    for i in idx:
        k = base.CAL_KEY[i]
        if np.isfinite(anomaly[i]).all():
            if valid_count[k] < 2:
                raise ValueError(f'Too few training observations on date key {k}')
            delta[k] = anomaly[i] / (valid_count[k] - 1)
    return delta

def fold_basis(anomaly: np.ndarray, all_gram: np.ndarray, w: int, n_keep: int, level: str):
    block = anomaly[base.BASIS_MASK & (base.WINTER == w)]
    if block.shape != (241, 1584) or not np.isfinite(block).all():
        raise ValueError(f'{level} winter {w}: invalid 241-day basis block')
    gram = all_gram - 75.0 / 74.0 * (block.T @ block)
    (lam, U) = eigh(gram, subset_by_index=(1584 - n_keep, 1583), driver='evr', check_finite=False, overwrite_a=True)
    U = U[:, ::-1]
    if level == 'S' and np.dot(U[:, 0], np.broadcast_to(np.cos(np.deg2rad(base.LAT))[:, None], (11, 144)).ravel()) < 0:
        U[:, 0] *= -1
    for j in range(1 if level == 'S' else 0, n_keep):
        if U[np.argmax(abs(U[:, j])), j] < 0:
            U[:, j] *= -1
    return U

def project_fold(anomaly: np.ndarray, delta: np.ndarray, U: np.ndarray):
    scores = np.full((len(base.DATES), U.shape[1]), np.nan, dtype=np.float64)
    idx = np.flatnonzero(base.PROJECT_MASK & np.isfinite(anomaly).all(axis=1))
    correction = delta @ U
    scores[idx] = anomaly[idx] @ U + correction[base.CAL_KEY[idx]]
    return scores

def fold_spv(polar_all: np.ndarray, delta_s: np.ndarray, w: int):
    weights = np.broadcast_to(np.cos(np.deg2rad(base.LAT))[:, None], (11, 144)).ravel()
    weights = weights / weights.sum()
    polar = polar_all + (delta_s @ weights)[base.CAL_KEY]
    residual = np.full(len(base.DATES), np.nan)
    x = ((base.DATES - np.datetime64('1985-01-01')) / np.timedelta64(1, 'D')).astype(float) / 365.2425
    train = base.TARGET_MASK & (base.WINTER != w) & np.isfinite(polar)
    for month in (11, 12, 1, 2, 3, 4):
        take = train & (base.CAL_MONTH == month)
        D = np.column_stack((np.ones(take.sum()), x[take]))
        beta = np.linalg.lstsq(D, polar[take], rcond=None)[0]
        all_month = base.TARGET_MASK & (base.CAL_MONTH == month) & np.isfinite(polar)
        residual[all_month] = polar[all_month] - beta[0] - beta[1] * x[all_month]
    mu = residual[train].mean()
    sd = residual[train].std(ddof=0)
    return ((residual - mu) / sd, mu, sd)

def make_design(s: np.ndarray, t: np.ndarray, y: np.ndarray, h: int):
    dst = np.flatnonzero(base.TARGET_MASK & np.isfinite(y))
    start = dst - h
    hist = start[:, None] - np.arange(60)[None, :]
    valid = (hist >= 0).all(axis=1)
    valid &= base.WINTER[hist[:, -1]] == base.WINTER[dst]
    valid &= np.isfinite(s[hist, 0]).all(axis=1)
    valid &= np.isfinite(s[start, 1:6]).all(axis=1)
    valid &= np.isfinite(t[start, :7]).all(axis=1)
    (dst, start, hist) = (dst[valid], start[valid], hist[valid])
    X = np.column_stack((s[hist, 0], s[start, 1:6], t[start, :7]))
    return (dst, X, y[dst])

def stats(X, y):
    n = len(y)
    return {'n': n, 'sx': X.sum(axis=0), 'sxx': X.T @ X, 'sy': y.sum(), 'sxy': X.T @ y}

def subtract(a, b):
    return {k: a[k] - b[k] for k in a}

def combine(parts):
    first = parts[0]
    return {k: sum((p[k] for p in parts)) for k in first}

def fit_from_stats(st, ncol: int):
    n = st['n']
    mean = st['sx'][:ncol] / n
    ym = st['sy'] / n
    centered = st['sxx'][:ncol, :ncol] - n * np.outer(mean, mean)
    variance = np.diag(centered) / n
    if np.any(variance <= 0):
        raise ValueError('Non-positive training feature variance')
    sd = np.sqrt(variance)
    gram = centered / sd[:, None] / sd[None, :]
    rhs = (st['sxy'][:ncol] - n * mean * ym) / sd
    try:
        beta_scaled = cho_solve(cho_factor(gram, check_finite=False), rhs, check_finite=False)
    except np.linalg.LinAlgError:
        beta_scaled = np.linalg.lstsq(gram, rhs, rcond=None)[0]
    beta = beta_scaled / sd
    intercept = ym - np.dot(mean, beta)
    return (intercept, beta)

def score_sums(y, pred):
    return np.array([len(y), y.sum(), pred.sum(), np.dot(y, y), np.dot(pred, pred), np.dot(y, pred), np.dot(pred - y, pred - y)], float)

def summarize(arr):
    (n, sy, sp, sy2, sp2, syp, se2) = np.moveaxis(arr, -1, 0)
    cov = syp - sy * sp / n
    vy = sy2 - sy * sy / n
    vp = sp2 - sp * sp / n
    r = np.divide(cov, np.sqrt(vy * vp), out=np.full_like(cov, np.nan), where=(n >= 2) & (vy > 0) & (vp > 0))
    q = np.divide(se2, sy2, out=np.full_like(se2, np.nan), where=sy2 > 0)
    return (r, q)

def draw_weights():
    all_years = np.arange(1948, 2023)
    non = all_years[~np.isin(all_years, SSW_YEARS)]
    ss = SSW_YEARS.copy()
    rng = np.random.Generator(np.random.PCG64(2026091301))
    all_w = rng.multinomial(75, np.full(75, 1 / 75), size=2000)
    non_w = rng.multinomial(36, np.full(36, 1 / 36), size=2000)
    rng = np.random.Generator(np.random.PCG64(2026091202))
    ss_w = rng.multinomial(39, np.full(39, 1 / 39), size=2000)
    return [(all_years, all_w), (ss, ss_w), (non, non_w)]

def plot_fig3(out: Path, metrics, intervals):
    from publication_style import plot_fig3 as draw_publication
    return draw_publication(out, metrics, intervals)

def main():
    pa = argparse.ArgumentParser(description='Reproduce report Fig. 3: 75-fold leave-one-winter-out direct SPV forecasts.\n\nThis script reads the original annual 10/500 hPa files. It never uses the\n241-day merged files or full-sample U as prediction inputs. Run after\nreproduce_fig1_fig2.py only for shared calendar constants and QA comparisons.\n')
    pa.add_argument('--s-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '10hpa')
    pa.add_argument('--t-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '500hpa')
    pa.add_argument('--output', type=Path, default=Path(__file__).with_name('fig3_fig4_results'))
    pa.add_argument('--first-winter', type=int, default=1948, help='Diagnostic subset start; final figure requires 1948')
    pa.add_argument('--last-winter', type=int, default=2022, help='Diagnostic subset end; final figure requires 2022')
    args = pa.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    print('Reading annual source fields and computing All75 daily anomalies', flush=True)
    (anom_s, _, count_s) = all75_anomalies(args.s_dir, 'S')
    (anom_t, _, count_t) = all75_anomalies(args.t_dir, 'T')
    basis_s = anom_s[base.BASIS_MASK]
    basis_t = anom_t[base.BASIS_MASK]
    (all_gram_s, all_gram_t) = (basis_s.T @ basis_s, basis_t.T @ basis_t)
    del basis_s, basis_t
    weights = np.broadcast_to(np.cos(np.deg2rad(base.LAT))[:, None], (11, 144)).ravel()
    weights = weights / weights.sum()
    polar_all = anom_s @ weights
    winter_years = np.arange(1948, 2023)
    labels = np.isin(winter_years, SSW_YEARS)
    all_sums = np.zeros((75, 60, 2, 7))
    cond_sums = np.zeros_like(all_sums)
    pred_all = np.full((len(base.DATES), 60, 2), np.nan)
    pred_cond = np.full_like(pred_all, np.nan)
    observed = np.full(len(base.DATES), np.nan)
    fold_scales = np.full((75, 2), np.nan)
    for w in range(args.first_winter, args.last_winter + 1):
        wi = w - 1948
        ds = leaveout_delta(anom_s, count_s, w)
        dt = leaveout_delta(anom_t, count_t, w)
        Us = fold_basis(anom_s, all_gram_s, w, 6, 'S')
        Ut = fold_basis(anom_t, all_gram_t, w, 7, 'T')
        ss = project_fold(anom_s, ds, Us)
        tt = project_fold(anom_t, dt, Ut)
        (yy, mu_y, sd_y) = fold_spv(polar_all, ds, w)
        fold_scales[wi] = [mu_y, sd_y]
        take_obs = (base.WINTER == w) & base.TARGET_MASK & np.isfinite(yy)
        observed[take_obs] = yy[take_obs]
        for h in range(1, 61):
            (dst, X, y) = make_design(ss, tt, yy, h)
            ww = base.WINTER[dst]
            block = []
            for year in winter_years:
                take = ww == year
                block.append(stats(X[take], y[take]))
            allstat = combine(block)
            ssstat = combine([block[i] for i in range(75) if labels[i]])
            nonstat = combine([block[i] for i in range(75) if not labels[i]])
            target = ww == w
            if target.sum() == 0:
                continue
            test_stat = block[wi]
            all_train = subtract(allstat, test_stat)
            group_train = subtract(ssstat if labels[wi] else nonstat, test_stat)
            for (branch, st, pred, sums) in ((0, all_train, pred_all, all_sums), (1, group_train, pred_cond, cond_sums)):
                for (model, ncol) in enumerate((60, 72)):
                    (alpha, beta) = fit_from_stats(st, ncol)
                    yh = alpha + X[target, :ncol] @ beta
                    pred[dst[target], h - 1, model] = yh
                    sums[wi, h - 1, model] = score_sums(y[target], yh)
        if (w - args.first_winter) % 5 == 0 or w == args.last_winter:
            print(f'Completed winter {w}/{args.last_winter}', flush=True)
    if args.first_winter != 1948 or args.last_winter != 2022:
        np.savez_compressed(args.output / 'fig3_partial_diagnostic.npz', pred_all=pred_all, pred_cond=pred_cond, observed=observed, all_sums=all_sums, cond_sums=cond_sums)
        print('Partial diagnostic saved; no final Fig3 scored', flush=True)
        return
    grouped = [all_sums, cond_sums[labels], cond_sums[~labels]]
    metric = np.empty((3, 2, 2, 60))
    bands = np.empty((3, 2, 2, 60, 2))
    count = np.empty((3, 60), int)
    for (j, (annual, (_, weights_boot))) in enumerate(zip(grouped, draw_weights())):
        summed = annual.sum(axis=0)
        (rr, qq) = summarize(summed)
        metric[j, 0] = rr.T
        metric[j, 1] = qq.T
        count[j] = summed[:, 0, 0].astype(int)
        boot_sums = np.einsum('bw,whmk->bhmk', weights_boot, annual, optimize=True)
        (br, bq) = summarize(boot_sums)
        for (row, bmetric) in enumerate((br, bq)):
            bounds = np.quantile(bmetric, [0.025, 0.975], axis=0).transpose(2, 1, 0)
            bands[j, row] = bounds
    np.savez_compressed(args.output / 'fig3_frozen_predictions.npz', dates=base.DATES, winter=base.WINTER, observed=observed, pred_all=pred_all, pred_conditional=pred_cond, all_winter_sums=all_sums, conditional_winter_sums=cond_sums, fold_spv_mean_scale=fold_scales, metrics=metric, intervals=bands, counts=count, h=np.arange(1, 61))
    with open(args.output / 'fig3_scores.csv', 'w', encoding='utf-8-sig', newline='') as f:
        out = csv.writer(f)
        out.writerow(['group', 'h', 'n', 'r_B', 'r_M', 'Q_B', 'Q_M', 'r_B_lo', 'r_B_hi', 'r_M_lo', 'r_M_hi', 'Q_B_lo', 'Q_B_hi', 'Q_M_lo', 'Q_M_hi'])
        for (j, name) in enumerate(('All', 'SSW', 'non-SSW')):
            for h in range(60):
                out.writerow([name, h + 1, count[j, h], metric[j, 0, 0, h], metric[j, 0, 1, h], metric[j, 1, 0, h], metric[j, 1, 1, h], *bands[j, 0, 0, h], *bands[j, 0, 1, h], *bands[j, 1, 0, h], *bands[j, 1, 1, h]])
    plot_fig3(args.output / 'Fig3', metric, bands)
    checks = {'h': [7, 15, 30, 60], 'n': count[:, [6, 14, 29, 59]].tolist(), 'r_B': metric[:, 0, 0, [6, 14, 29, 59]].tolist(), 'r_M': metric[:, 0, 1, [6, 14, 29, 59]].tolist(), 'Q_B': metric[:, 1, 0, [6, 14, 29, 59]].tolist(), 'Q_M': metric[:, 1, 1, [6, 14, 29, 59]].tolist()}
    (args.output / 'fig3_QA.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    print(json.dumps(checks, indent=2), flush=True)
if __name__ == '__main__':
    main()
