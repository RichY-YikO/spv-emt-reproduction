from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import reproduce_fig1_fig2 as base
import reproduce_fig3 as fig3
from reproduce_fig4 import SSW_YEARS
MONTHS = np.array([11, 12, 1, 2, 3, 4])
MONTH_LABELS = ['November', 'December', 'January', 'February', 'March', 'April']
MODEL_LABELS = ['B', 'M']

def monthly_winter_sums(observed, prediction):
    out = np.zeros((75, 6, 60, 2, 7), dtype=np.float64)
    winter_idx = base.WINTER - 1948
    for (mi, month) in enumerate(MONTHS):
        date_ok = (base.CAL_MONTH == month) & np.isfinite(observed)
        for hi in range(60):
            take = date_ok & np.isfinite(prediction[:, hi, 0]) & np.isfinite(prediction[:, hi, 1])
            wi = winter_idx[take]
            y = observed[take]
            for model in range(2):
                pred = prediction[take, hi, model]
                out[:, mi, hi, model] = np.stack((np.bincount(wi, minlength=75), np.bincount(wi, weights=y, minlength=75), np.bincount(wi, weights=pred, minlength=75), np.bincount(wi, weights=y * y, minlength=75), np.bincount(wi, weights=pred * pred, minlength=75), np.bincount(wi, weights=y * pred, minlength=75), np.bincount(wi, weights=(pred - y) ** 2, minlength=75)), axis=1)
    return out

def score(monthly, boot_weights):
    summed = monthly.sum(axis=0)
    (r, q) = fig3.summarize(summed)
    metrics = np.stack((r.transpose(0, 2, 1), q.transpose(0, 2, 1)), axis=1)
    boot_sum = np.einsum('bw,wmhsk->bmhsk', boot_weights, monthly, optimize=True)
    (br, bq) = fig3.summarize(boot_sum)
    bands = np.empty((6, 2, 2, 60, 2), dtype=np.float64)
    for (metric, value) in enumerate((br, bq)):
        bands[:, metric] = np.quantile(value, [0.025, 0.975], axis=0).transpose(1, 3, 2, 0)
    counts = summed[:, :, 0, 0].astype(int)
    n_winters = (monthly[:, :, :, 0, 0] > 0).sum(axis=0)
    return (metrics, bands, counts, n_winters)

def plot(path: Path, metrics, bands, q_upper: float):
    from publication_style import plot_monthly as draw_publication
    return draw_publication(path, metrics, bands, q_upper)

def write_csv(path, metrics, bands, counts, n_winters):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['target_month', 'h', 'model', 'n_pairs', 'n_winters', 'r', 'Q', 'r_lo', 'r_hi', 'Q_lo', 'Q_hi'])
        for (mi, month) in enumerate(MONTHS):
            for hi in range(60):
                for (model, name) in enumerate(MODEL_LABELS):
                    writer.writerow([int(month), hi + 1, name, int(counts[mi, hi]), int(n_winters[mi, hi]), metrics[mi, 0, model, hi], metrics[mi, 1, model, hi], *bands[mi, 0, model, hi], *bands[mi, 1, model, hi]])

def main():
    p = argparse.ArgumentParser(description="Reproduce report Fig. S4/S5/S6 by target-month scoring of frozen Fig. 3 predictions.\n\nNo EMT, SPV preprocessing or regression is refit by month. The 2000 paired\nwinter bootstrap uses the report's same PCG64 winter weights as Fig. 3.\n")
    p.add_argument('--fig3-archive', type=Path, default=Path(__file__).with_name('fig3_fig4_results') / 'fig3_frozen_predictions.npz')
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('figS3_S9_results'))
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.fig3_archive) as z:
        np.testing.assert_array_equal(z['dates'], base.DATES)
        observed = z['observed'].copy()
        pred_all = z['pred_all'].copy()
        pred_cond = z['pred_conditional'].copy()
        expected_all = z['all_winter_sums'].copy()
        expected_cond = z['conditional_winter_sums'].copy()
    mask_ssw = np.isin(np.arange(1948, 2023), SSW_YEARS)
    boot = fig3.draw_weights()
    monthly_all = monthly_winter_sums(observed, pred_all)
    np.testing.assert_allclose(monthly_all.sum(axis=1), expected_all, rtol=0, atol=1e-08)
    monthly_cond = monthly_winter_sums(observed, pred_cond)
    np.testing.assert_allclose(monthly_cond.sum(axis=1), expected_cond, rtol=0, atol=1e-08)
    qa = {}
    for (fig_name, data, weights, expected_first, march_range, april_range, winter_range, q_upper) in (('FigS4', monthly_all, boot[0][1], [2250, 2325, 2325, 2100], [1767, 2325], [1710, 2214], [57, 75], 1.6), ('FigS5', monthly_cond[mask_ssw], boot[1][1], [1170, 1209, 1209, 1092], [930, 1209], [900, 1152], [30, 39], 1.85), ('FigS6', monthly_cond[~mask_ssw], boot[2][1], [1080, 1116, 1116, 1008], [837, 1116], [810, 1062], [27, 36], 1.6)):
        (metrics, bands, counts, n_winters) = score(data, weights)
        np.testing.assert_array_equal(counts[:4].min(axis=1), expected_first)
        np.testing.assert_array_equal(counts[:4].max(axis=1), expected_first)
        np.testing.assert_array_equal([counts[4].min(), counts[4].max()], march_range)
        np.testing.assert_array_equal([counts[5].min(), counts[5].max()], april_range)
        np.testing.assert_array_equal(n_winters[:4].min(axis=1), [data.shape[0]] * 4)
        np.testing.assert_array_equal([n_winters[4:].min(), n_winters[4:].max()], winter_range)
        if np.nanmax(bands[:, 1, :, :, 1]) > q_upper + 1e-12:
            raise ValueError(f'{fig_name}: Q interval would be clipped at {q_upper}')
        plot(args.output / fig_name, metrics, bands, q_upper)
        write_csv(args.output / f'{fig_name}_monthly_scores.csv', metrics, bands, counts, n_winters)
        np.savez_compressed(args.output / f'{fig_name}_monthly_scores.npz', months=MONTHS, h=np.arange(1, 61), metrics=metrics, intervals=bands, counts=counts, n_winters=n_winters, winter_sums=data)
        qa[fig_name] = {'first_four_month_counts': expected_first, 'march_min_max': [int(counts[4].min()), int(counts[4].max())], 'april_min_max': [int(counts[5].min()), int(counts[5].max())], 'winter_count_min_max': [int(n_winters.min()), int(n_winters.max())], 'Q_interval_max': float(np.nanmax(bands[:, 1, :, :, 1]))}
        print(f'{fig_name} complete: {qa[fig_name]}', flush=True)
    (args.output / 'figS4_S5_S6_QA.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
if __name__ == '__main__':
    main()
