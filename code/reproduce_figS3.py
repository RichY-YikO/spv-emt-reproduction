from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_name, '4')
import numpy as np
from scipy.linalg import cho_factor, cho_solve
import reproduce_fig1_fig2 as base
import reproduce_fig3 as fig3
from reproduce_fig4 import SSW_YEARS
NAMES = ['B', 'BS', 'BT']
COLS = [np.arange(60), np.arange(65), np.r_[np.arange(60), np.arange(65, 72)]]

def fit_columns(st, cols):
    n = st['n']
    mean = st['sx'][cols] / n
    ym = st['sy'] / n
    centered = st['sxx'][np.ix_(cols, cols)] - n * np.outer(mean, mean)
    sd = np.sqrt(np.diag(centered) / n)
    if not np.isfinite(sd).all() or (sd <= 0).any():
        raise ValueError('Invalid training predictor variance')
    gram = centered / sd[:, None] / sd[None, :]
    rhs = (st['sxy'][cols] - n * mean * ym) / sd
    try:
        scaled = cho_solve(cho_factor(gram, check_finite=False), rhs, check_finite=False)
    except np.linalg.LinAlgError:
        scaled = np.linalg.lstsq(gram, rhs, rcond=None)[0]
    beta = scaled / sd
    return (ym - mean @ beta, beta)

def score_groups(all_sums, cond_sums):
    winter_years = np.arange(1948, 2023)
    ssw = np.isin(winter_years, SSW_YEARS)
    grouped = [all_sums, cond_sums[ssw], cond_sums[~ssw]]
    metrics = np.empty((3, 2, 3, 60))
    bands = np.empty((3, 2, 3, 60, 2))
    counts = np.empty((3, 60), dtype=int)
    for (group, (annual, (_, weights))) in enumerate(zip(grouped, fig3.draw_weights())):
        (r, q) = fig3.summarize(annual.sum(axis=0))
        (metrics[group, 0], metrics[group, 1]) = (r.T, q.T)
        counts[group] = annual.sum(axis=0)[:, 0, 0].astype(int)
        boot = np.einsum('bw,whmk->bhmk', weights, annual, optimize=True)
        (br, bq) = fig3.summarize(boot)
        for (metric, value) in enumerate((br, bq)):
            bands[group, metric] = np.quantile(value, [0.025, 0.975], axis=0).transpose(2, 1, 0)
    return (metrics, bands, counts)

def plot(path, metrics, bands):
    from publication_style import plot_s3 as draw_publication
    return draw_publication(path, metrics, bands)

def write_csv(path, metrics, bands, counts):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['group', 'h', 'model', 'n_pairs', 'r', 'Q', 'r_lo', 'r_hi', 'Q_lo', 'Q_hi'])
        for (g, group) in enumerate(('All', 'SSW', 'non-SSW')):
            for hi in range(60):
                for (m, name) in enumerate(NAMES):
                    writer.writerow([group, hi + 1, name, int(counts[g, hi]), metrics[g, 0, m, hi], metrics[g, 1, m, hi], *bands[g, 0, m, hi], *bands[g, 1, m, hi]])

def main():
    p = argparse.ArgumentParser(description='Reproduce report Fig. S3: B versus separate S and T extensions.\n\nThe 75-fold 10/500 hPa EMT, SPV target and pair mask follow reproduce_fig3.py.\nVerified Fig. 3 B predictions are reused; BS/BT are independently refit for\neach held-out winter, group branch and lead. No 100 hPa heat flux is used.\n')
    p.add_argument('--s-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '10hpa')
    p.add_argument('--t-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '500hpa')
    p.add_argument('--fig3-archive', type=Path, default=Path(__file__).with_name('fig3_fig4_results') / 'fig3_frozen_predictions.npz')
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('figS3_S9_results'))
    p.add_argument('--first-winter', type=int, default=1948)
    p.add_argument('--last-winter', type=int, default=2022)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.fig3_archive) as z:
        np.testing.assert_array_equal(z['dates'], base.DATES)
        observed = z['observed'].copy()
        old_all = z['all_winter_sums'].copy()
        old_cond = z['conditional_winter_sums'].copy()
        old_pred_all = z['pred_all'].copy()
        old_pred_cond = z['pred_conditional'].copy()
        old_counts = z['counts'].copy()
    all_sums = np.full((75, 60, 3, 7), np.nan)
    cond_sums = np.full_like(all_sums, np.nan)
    all_sums[:, :, 0] = old_all[:, :, 0]
    cond_sums[:, :, 0] = old_cond[:, :, 0]
    pred_all = np.full((len(base.DATES), 60, 2), np.nan)
    pred_cond = np.full_like(pred_all, np.nan)
    print('Reading annual 10/500 hPa source fields', flush=True)
    (anom_s, _, count_s) = fig3.all75_anomalies(args.s_dir, 'S')
    (anom_t, _, count_t) = fig3.all75_anomalies(args.t_dir, 'T')
    gram_s = anom_s[base.BASIS_MASK].T @ anom_s[base.BASIS_MASK]
    gram_t = anom_t[base.BASIS_MASK].T @ anom_t[base.BASIS_MASK]
    cap_w = np.broadcast_to(np.cos(np.deg2rad(base.LAT))[:, None], (11, 144)).ravel()
    cap_w /= cap_w.sum()
    polar_all = anom_s @ cap_w
    years = np.arange(1948, 2023)
    labels = np.isin(years, SSW_YEARS)
    for w in range(args.first_winter, args.last_winter + 1):
        wi = w - 1948
        ds = fig3.leaveout_delta(anom_s, count_s, w)
        dt = fig3.leaveout_delta(anom_t, count_t, w)
        Us = fig3.fold_basis(anom_s, gram_s, w, 6, 'S')
        Ut = fig3.fold_basis(anom_t, gram_t, w, 7, 'T')
        ss = fig3.project_fold(anom_s, ds, Us)
        tt = fig3.project_fold(anom_t, dt, Ut)
        (yy, _, _) = fig3.fold_spv(polar_all, ds, w)
        own_target = (base.WINTER == w) & base.TARGET_MASK & np.isfinite(yy)
        if not np.allclose(yy[own_target], observed[own_target], rtol=0, atol=1e-10):
            raise ValueError(f'Winter {w}: SPV differs from Fig3 archive')
        for h in range(1, 61):
            (dst, X, y) = fig3.make_design(ss, tt, yy, h)
            winter = base.WINTER[dst]
            target = winter == w
            if not target.any():
                continue
            test_dates = dst[target]
            if target.sum() != old_all[wi, h - 1, 0, 0]:
                raise ValueError(f'Winter {w}, h={h}: Fig3 sample count differs')
            if not (np.isfinite(old_pred_all[test_dates, h - 1, 0]).all() and np.isfinite(old_pred_cond[test_dates, h - 1, 0]).all()):
                raise ValueError(f'Winter {w}, h={h}: Fig3 pair dates differ')
            total_st = fig3.stats(X, y)
            same_ssw = labels[winter - 1948]
            ssw_st = fig3.stats(X[same_ssw], y[same_ssw])
            non_st = fig3.stats(X[~same_ssw], y[~same_ssw])
            held_st = fig3.stats(X[target], y[target])
            all_train = fig3.subtract(total_st, held_st)
            cond_train = fig3.subtract(ssw_st if labels[wi] else non_st, held_st)
            for (train, sums, predictions) in ((all_train, all_sums, pred_all), (cond_train, cond_sums, pred_cond)):
                for m in (1, 2):
                    (alpha, beta) = fit_columns(train, COLS[m])
                    yh = alpha + X[target][:, COLS[m]] @ beta
                    predictions[test_dates, h - 1, m - 1] = yh
                    sums[wi, h - 1, m] = fig3.score_sums(y[target], yh)
        if (w - args.first_winter) % 5 == 0 or w == args.last_winter:
            print(f'FigS3 completed held-out winter {w}', flush=True)
    if args.first_winter != 1948 or args.last_winter != 2022:
        np.savez_compressed(args.output / 'figS3_partial_diagnostic.npz', all_sums=all_sums, conditional_sums=cond_sums)
        return
    if not (np.isfinite(all_sums).all() and np.isfinite(cond_sums).all()):
        raise ValueError('Missing a winter/lead/model score')
    (metrics, bands, counts) = score_groups(all_sums, cond_sums)
    np.testing.assert_array_equal(counts, old_counts)
    with np.load(args.fig3_archive) as old:
        np.testing.assert_allclose(metrics[:, :, 0], old['metrics'][:, :, 0], atol=1e-12, rtol=0)
        np.testing.assert_allclose(bands[:, :, 0], old['intervals'][:, :, 0], atol=1e-12, rtol=0)
    np.testing.assert_allclose(metrics[0, 0, 1:, 29], [0.28571986, 0.29039979], atol=5e-08, rtol=0)
    np.testing.assert_allclose(metrics[0, 1, 1:, 29], [0.91889913, 0.91624888], atol=5e-08, rtol=0)
    plot(args.output / 'FigS3', metrics, bands)
    write_csv(args.output / 'FigS3_scores.csv', metrics, bands, counts)
    np.savez_compressed(args.output / 'FigS3_scores.npz', h=np.arange(1, 61), model_names=np.array(NAMES), metrics=metrics, intervals=bands, counts=counts, all_winter_sums=all_sums, conditional_winter_sums=cond_sums, dates=base.DATES, pred_all_BS_BT=pred_all, pred_conditional_BS_BT=pred_cond)
    qa = {'all_h30_BS_r_Q': metrics[0, :, 1, 29].tolist(), 'all_h30_BT_r_Q': metrics[0, :, 2, 29].tolist(), 'all_h30_B_r_Q': metrics[0, :, 0, 29].tolist(), 'h30_pair_counts_All_SSW_nonSSW': counts[:, 29].tolist()}
    (args.output / 'FigS3_QA.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
    print(json.dumps(qa, indent=2), flush=True)
if __name__ == '__main__':
    main()
