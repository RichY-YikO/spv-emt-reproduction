from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_name, '4')
import numpy as np
from scipy.linalg import eigh
import reproduce_fig1_fig2 as base
import reproduce_fig3 as fig3
from reproduce_fig4 import SSW_YEARS
DJF_MASK = base.VALID_WINTER & ((base.CAL_MONTH == 12) | (base.CAL_MONTH == 1) | (base.CAL_MONTH == 2)) & (base.CAL_KEY != 59)
if int(DJF_MASK.sum()) != 6750:
    raise ValueError('DJF basis must contain 75 winters x 90 days')

def djf_basis(anomaly, gram_all, w, n_keep, reference):
    block = anomaly[DJF_MASK & (base.WINTER == w)]
    if block.shape != (90, 1584) or not np.isfinite(block).all():
        raise ValueError(f'Winter {w}: invalid 90-day DJF basis block')
    gram = gram_all - 75.0 / 74.0 * (block.T @ block)
    (_, U) = eigh(gram, subset_by_index=(1584 - n_keep, 1583), driver='evr', check_finite=False, overwrite_a=True)
    U = U[:, ::-1]
    overlap_before = np.sum(U * reference, axis=0)
    U *= np.where(overlap_before < 0, -1.0, 1.0)[None, :]
    overlap_after = np.sum(U * reference, axis=0)
    if (overlap_after < -1e-13).any():
        raise ValueError('Rankwise sign alignment failed')
    return (U, overlap_after)

def score_all(djf_all, djf_cond, default_all, default_cond):
    years = np.arange(1948, 2023)
    ssw = np.isin(years, SSW_YEARS)
    grouped = [(default_all, djf_all), (default_cond[ssw], djf_cond[ssw]), (default_cond[~ssw], djf_cond[~ssw])]
    metrics = np.empty((3, 2, 2, 2, 60))
    bands = np.empty((3, 2, 2, 2, 60, 2))
    counts = np.empty((3, 60), dtype=int)
    for (g, ((sep, djf), (_, weights))) in enumerate(zip(grouped, fig3.draw_weights())):
        np.testing.assert_array_equal(sep[:, :, 0, 0], djf[:, :, 0, 0])
        for (scheme, annual) in enumerate((sep, djf)):
            (r, q) = fig3.summarize(annual.sum(axis=0))
            metrics[g, scheme, 0] = r.T
            metrics[g, scheme, 1] = q.T
            counts[g] = annual.sum(axis=0)[:, 0, 0].astype(int)
            boot = np.einsum('bw,whmk->bhmk', weights, annual, optimize=True)
            (br, bq) = fig3.summarize(boot)
            for (metric, value) in enumerate((br, bq)):
                bands[g, scheme, metric] = np.quantile(value, [0.025, 0.975], axis=0).transpose(2, 1, 0)
    return (metrics, bands, counts)

def plot(path, metrics, bands):
    from publication_style import plot_s9 as draw_publication
    return draw_publication(path, metrics, bands)

def main():
    p = argparse.ArgumentParser(description='Reproduce report Fig. S9: default Sep-Apr versus DJF EMT basis window.\n\nThe only experimental change is the season used to estimate spatial U. Every\nheld-out winter has All74 climate/SPV preprocessing and the same dates in both\nschemes; all July-April fields, including the S1 history, are reprojected to\ntheir corresponding basis. B/M are independently refit for each scheme.\n')
    p.add_argument('--s-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '10hpa')
    p.add_argument('--t-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '500hpa')
    p.add_argument('--fig3-archive', type=Path, default=Path(__file__).with_name('fig3_fig4_results') / 'fig3_frozen_predictions.npz')
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('figS3_S9_results'))
    p.add_argument('--first-winter', type=int, default=1948)
    p.add_argument('--last-winter', type=int, default=2022)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.fig3_archive) as old:
        np.testing.assert_array_equal(old['dates'], base.DATES)
        obs = old['observed'].copy()
        sep_all = old['all_winter_sums'].copy()
        sep_cond = old['conditional_winter_sums'].copy()
        sep_metrics = old['metrics'].copy()
        sep_bands = old['intervals'].copy()
        old_pred_all = old['pred_all'].copy()
        old_pred_cond = old['pred_conditional'].copy()
    djf_all = np.full_like(sep_all, np.nan)
    djf_cond = np.full_like(sep_cond, np.nan)
    pred_all = np.full((len(base.DATES), 60, 2), np.nan)
    pred_cond = np.full_like(pred_all, np.nan)
    overlaps = np.full((75, 13), np.nan)
    print('Reading annual 10/500 hPa source fields', flush=True)
    (anom_s, _, count_s) = fig3.all75_anomalies(args.s_dir, 'S')
    (anom_t, _, count_t) = fig3.all75_anomalies(args.t_dir, 'T')
    sep_gram_s = anom_s[base.BASIS_MASK].T @ anom_s[base.BASIS_MASK]
    sep_gram_t = anom_t[base.BASIS_MASK].T @ anom_t[base.BASIS_MASK]
    djf_gram_s = anom_s[DJF_MASK].T @ anom_s[DJF_MASK]
    djf_gram_t = anom_t[DJF_MASK].T @ anom_t[DJF_MASK]
    cap_w = np.broadcast_to(np.cos(np.deg2rad(base.LAT))[:, None], (11, 144)).ravel()
    cap_w /= cap_w.sum()
    polar_all = anom_s @ cap_w
    years = np.arange(1948, 2023)
    labels = np.isin(years, SSW_YEARS)
    for w in range(args.first_winter, args.last_winter + 1):
        wi = w - 1948
        ds = fig3.leaveout_delta(anom_s, count_s, w)
        dt = fig3.leaveout_delta(anom_t, count_t, w)
        sep_s = fig3.fold_basis(anom_s, sep_gram_s, w, 6, 'S')
        sep_t = fig3.fold_basis(anom_t, sep_gram_t, w, 7, 'T')
        (U_s, overlaps[wi, :6]) = djf_basis(anom_s, djf_gram_s, w, 6, sep_s)
        (U_t, overlaps[wi, 6:]) = djf_basis(anom_t, djf_gram_t, w, 7, sep_t)
        s_djf = fig3.project_fold(anom_s, ds, U_s)
        t_djf = fig3.project_fold(anom_t, dt, U_t)
        (yy, _, _) = fig3.fold_spv(polar_all, ds, w)
        own = (base.WINTER == w) & base.TARGET_MASK & np.isfinite(yy)
        if not np.allclose(yy[own], obs[own], atol=1e-10, rtol=0):
            raise ValueError(f'Winter {w}: SPV target changed')
        for h in range(1, 61):
            (dst, X, y) = fig3.make_design(s_djf, t_djf, yy, h)
            winter = base.WINTER[dst]
            target = winter == w
            if not target.any():
                continue
            held_dates = dst[target]
            if target.sum() != sep_all[wi, h - 1, 0, 0]:
                raise ValueError(f'Winter {w}, h={h}: DJF changed date mask')
            if not (np.isfinite(old_pred_all[held_dates, h - 1]).all() and np.isfinite(old_pred_cond[held_dates, h - 1]).all()):
                raise ValueError(f'Winter {w}, h={h}: Sep-Apr date mismatch')
            total_st = fig3.stats(X, y)
            ssw_take = labels[winter - 1948]
            ssw_st = fig3.stats(X[ssw_take], y[ssw_take])
            non_st = fig3.stats(X[~ssw_take], y[~ssw_take])
            held_st = fig3.stats(X[target], y[target])
            all_train = fig3.subtract(total_st, held_st)
            cond_train = fig3.subtract(ssw_st if labels[wi] else non_st, held_st)
            for (train, sums, predictions) in ((all_train, djf_all, pred_all), (cond_train, djf_cond, pred_cond)):
                for (model, ncol) in enumerate((60, 72)):
                    (alpha, beta) = fig3.fit_from_stats(train, ncol)
                    yh = alpha + X[target, :ncol] @ beta
                    predictions[held_dates, h - 1, model] = yh
                    sums[wi, h - 1, model] = fig3.score_sums(y[target], yh)
        if (w - args.first_winter) % 5 == 0 or w == args.last_winter:
            print(f'FigS9 completed held-out winter {w}', flush=True)
    if args.first_winter != 1948 or args.last_winter != 2022:
        np.savez_compressed(args.output / 'FigS9_partial_diagnostic.npz', all_sums=djf_all, conditional_sums=djf_cond, rankwise_signed_overlaps=overlaps)
        return
    if not (np.isfinite(djf_all).all() and np.isfinite(djf_cond).all() and np.isfinite(overlaps).all()):
        raise ValueError('Some DJF scores or rankwise overlaps missing')
    (metrics, bands, counts) = score_all(djf_all, djf_cond, sep_all, sep_cond)
    np.testing.assert_allclose(metrics[:, 0], sep_metrics, atol=1e-12, rtol=0)
    np.testing.assert_allclose(bands[:, 0], sep_bands, atol=1e-12, rtol=0)
    expected = np.array([[[0.33304605, 0.33389748], [0.88992776, 0.88935802]], [[0.40755255, 0.40832539], [0.81986196, 0.81923101]], [[0.31065527, 0.31082662], [0.87419217, 0.87410539]]])
    actual = metrics[:, :, :, 1, 29].transpose(0, 2, 1)
    np.testing.assert_allclose(actual, expected, atol=5e-08, rtol=0)
    plot(args.output / 'FigS9', metrics, bands)
    np.savez_compressed(args.output / 'FigS9_scores.npz', h=np.arange(1, 61), schemes=np.array(['Sep-Apr', 'DJF']), metrics=metrics, intervals=bands, counts=counts, DJF_all_winter_sums=djf_all, DJF_conditional_winter_sums=djf_cond, rankwise_signed_overlaps=overlaps, dates=base.DATES, DJF_pred_all=pred_all, DJF_pred_conditional=pred_cond)
    with (args.output / 'FigS9_scores.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['group', 'basis_window', 'h', 'model', 'n_pairs', 'r', 'Q', 'r_lo', 'r_hi', 'Q_lo', 'Q_hi'])
        for (g, group) in enumerate(('All', 'SSW', 'non-SSW')):
            for (scheme, name) in enumerate(('Sep-Apr', 'DJF')):
                for hi in range(60):
                    for (model, model_name) in enumerate(('B', 'M')):
                        writer.writerow([group, name, hi + 1, model_name, int(counts[g, hi]), metrics[g, scheme, 0, model, hi], metrics[g, scheme, 1, model, hi], *bands[g, scheme, 0, model, hi], *bands[g, scheme, 1, model, hi]])
    qa = {'DJF_basis_days_per_winter': 90, 'min_rankwise_aligned_dot_product': float(overlaps.min()), 'h30_M_r_Q_SepApr_DJF': {group: {'r': actual[g, 0].tolist(), 'Q': actual[g, 1].tolist()} for (g, group) in enumerate(('All', 'SSW', 'non-SSW'))}, 'h30_pair_counts_All_SSW_nonSSW': counts[:, 29].tolist()}
    (args.output / 'FigS9_QA.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
    print(json.dumps(qa, indent=2), flush=True)
if __name__ == '__main__':
    main()
