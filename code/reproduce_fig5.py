from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_name, '4')
import numpy as np
from netCDF4 import Dataset
from scipy.linalg import cho_factor, cho_solve
import reproduce_fig1_fig2 as base
import reproduce_fig3 as fig3
from reproduce_fig4 import SSW_YEARS, MODE_NAMES
MODEL_NAMES = ['B', 'BH', 'M', 'MH'] + [f'M_minus_{name}' for name in MODE_NAMES] + [f'M_minus_{name}_H' for name in MODE_NAMES]
COLS = [np.arange(60), np.r_[np.arange(60), 72], np.arange(72), np.arange(73)]
for j in range(12):
    COLS.append(np.r_[np.arange(60), np.arange(60, 72)[np.arange(12) != j]])
for j in range(12):
    COLS.append(np.r_[np.arange(60), np.arange(60, 72)[np.arange(12) != j], 72])

def daily_heat_flux(v_dir: Path, t_dir: Path, cache: Path, refresh=False):
    if cache.is_file() and (not refresh):
        with np.load(cache) as z:
            np.testing.assert_array_equal(z['dates'], base.DATES)
            return z['H'].copy()
    H = np.full(len(base.DATES), np.nan, dtype=np.float64)
    weights = np.cos(np.deg2rad(np.arange(75.0, 44.0, -1.0)))
    weights /= weights.sum()
    epoch = np.datetime64('1970-01-01T00:00:00', 's')
    for year in range(1948, 2024):
        vp = v_dir / f'era5_v_100hPa_{year}.nc'
        tp = t_dir / f'era5_temperature_100hPa_{year}.nc'
        with Dataset(vp) as vd, Dataset(tp) as td:
            vt = np.asarray(vd.variables['valid_time'][:], np.int64)
            tt = np.asarray(td.variables['valid_time'][:], np.int64)
            np.testing.assert_array_equal(vt, tt)
            np.testing.assert_array_equal(vd.variables['latitude'][:], td.variables['latitude'][:])
            np.testing.assert_array_equal(vd.variables['longitude'][:], td.variables['longitude'][:])
            np.testing.assert_array_equal(vd.variables['pressure_level'][:], [100])
            np.testing.assert_array_equal(td.variables['pressure_level'][:], [100])
            if len(vt) % 4 or not np.array_equal(vt.reshape(-1, 4) % 86400, np.tile([0, 21600, 43200, 64800], (len(vt) // 4, 1))):
                raise ValueError(f'{year}: wrong UTC sequence')
            file_days = (epoch + vt[::4].astype('timedelta64[s]')).astype('datetime64[D]')
            pos = ((file_days - base.START) / np.timedelta64(1, 'D')).astype(int)
            if (pos < 0).any() or (pos >= len(H)).any():
                raise ValueError(f'{year}: date outside study axis')
            daily = np.empty(len(file_days), dtype=np.float64)
            for lo_day in range(0, len(daily), 32):
                hi_day = min(lo_day + 32, len(daily))
                (lo, hi) = (4 * lo_day, 4 * hi_day)
                v = np.asarray(np.ma.filled(vd.variables['v'][lo:hi, 0], np.nan), np.float64)
                t = np.asarray(np.ma.filled(td.variables['t'][lo:hi, 0], np.nan), np.float64)
                if v.shape != t.shape or v.shape[1:] != (31, 360):
                    raise ValueError(f'{year}: mismatched v/T grid')
                if not (np.isfinite(v).all() and np.isfinite(t).all()):
                    raise ValueError(f'{year}: missing source grid values')
                v -= v.mean(axis=-1, keepdims=True)
                t -= t.mean(axis=-1, keepdims=True)
                by_lat = np.mean(v * t, axis=-1)
                by_utc = by_lat @ weights
                daily[lo_day:hi_day] = by_utc.reshape(-1, 4).mean(axis=1)
            H[pos] = daily
        if year % 10 == 0 or year in (1948, 2023):
            print(f'100 hPa heat flux complete through {year}', flush=True)
    if int(np.isfinite(H).sum()) != 25143:
        raise ValueError(f'Expected 25143 daily fluxes; got {np.isfinite(H).sum()}')
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, dates=base.DATES, H=H, units='K m s-1')
    with (cache.parent / 'h100_daily.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['date', 'H_K_m_s-1'])
        for (d, value) in zip(base.DATES, H):
            if np.isfinite(value):
                writer.writerow([str(d), f'{value:.15g}'])
    return H

def fold_h40(H: np.ndarray, w: int, sums: np.ndarray, counts: np.ndarray):
    held = (base.WINTER == w) & np.isfinite(H)
    held_sum = np.bincount(base.CAL_KEY[held], weights=H[held], minlength=366)
    held_n = np.bincount(base.CAL_KEY[held], minlength=366)
    n = counts - held_n
    clim = np.divide(sums - held_sum, n, out=np.full(366, np.nan), where=n > 0)
    anomaly = H - clim[base.CAL_KEY]
    good = np.isfinite(anomaly)
    x = np.where(good, anomaly, 0.0)
    cs = np.r_[0.0, np.cumsum(x)]
    cn = np.r_[0, np.cumsum(good.astype(np.int32))]
    h40 = np.full(len(H), np.nan)
    idx = np.arange(39, len(H))
    complete = cn[idx + 1] - cn[idx - 39] == 40
    complete &= base.WINTER[idx] == base.WINTER[idx - 39]
    h40[idx[complete]] = (cs[idx[complete] + 1] - cs[idx[complete] - 39]) / 40.0
    return h40

def fit_columns(st: dict, cols: np.ndarray):
    n = st['n']
    mean = st['sx'][cols] / n
    ym = st['sy'] / n
    centered = st['sxx'][np.ix_(cols, cols)] - n * np.outer(mean, mean)
    variance = np.diag(centered) / n
    if np.any(variance <= 0):
        raise ValueError('Nonpositive predictor variance')
    sd = np.sqrt(variance)
    gram = centered / sd[:, None] / sd[None, :]
    rhs = (st['sxy'][cols] - n * mean * ym) / sd
    try:
        beta_scaled = cho_solve(cho_factor(gram, check_finite=False), rhs, check_finite=False)
    except np.linalg.LinAlgError:
        beta_scaled = np.linalg.lstsq(gram, rhs, rcond=None)[0]
    beta = beta_scaled / sd
    return (ym - mean @ beta, beta)

def make_fig5_design(s, t, y, H40, h):
    (dst, X0, response) = fig3.make_design(s, t, y, h)
    start = dst - h
    good = np.isfinite(H40[start])
    if not good.all():
        raise ValueError(f'h={h}: {np.count_nonzero(~good)} Fig3 pairs lack H40')
    X = np.column_stack((X0, H40[start]))
    return (dst, X, response)

def q_from_sums(arr):
    return arr[..., 6] / arr[..., 3]

def score_fig5(all_sums, cond_sums):
    years = np.arange(1948, 2023)
    ssw_mask = np.isin(years, SSW_YEARS)
    grouped = [all_sums, cond_sums[ssw_mask], cond_sums[~ssw_mask]]
    metrics = np.empty((3, 3, 60))
    intervals = np.empty((3, 3, 60, 2))
    C = np.empty((3, 12, 60))
    counts = np.empty((3, 60), dtype=int)
    for (g, (annual, (_, boot_w))) in enumerate(zip(grouped, fig3.draw_weights())):
        totals = annual.sum(axis=0)
        q = q_from_sums(totals)
        counts[g] = totals[:, 0, 0].astype(int)
        metrics[g] = np.stack((q[:, 0] - q[:, 1], q[:, 0] - q[:, 2], q[:, 2] - q[:, 3]), axis=0)
        C[g] = (q[:, 4:16] - q[:, 16:28] - (q[:, 2] - q[:, 3])[:, None]).T
        errors = annual[:, :, [0, 1, 2, 3], 6]
        denominator = annual[:, :, 0, 3]
        boot_e = np.einsum('bw,whm->bhm', boot_w, errors, optimize=True)
        boot_d = boot_w @ denominator
        bq = boot_e / boot_d[:, :, None]
        bgain = np.stack((bq[:, :, 0] - bq[:, :, 1], bq[:, :, 0] - bq[:, :, 2], bq[:, :, 2] - bq[:, :, 3]), axis=1)
        for m in range(3):
            intervals[g, m] = np.quantile(bgain[:, m], [0.025, 0.975], axis=0).T
    return (metrics, intervals, C, counts)

def plot_fig5(path: Path, gains, bands, C):
    from publication_style import plot_fig5 as draw_publication
    return draw_publication(path, gains, bands, C)

def main():
    p = argparse.ArgumentParser(description='Reproduce report Fig. 5 from 100 hPa v/T and annual 10/500 hPa fields.\n\nRequires the verified Fig. 3 archive beside this script. Usage:\n    python reproduce_fig5.py\n    python reproduce_fig5.py --heat-only\n')
    p.add_argument('--v-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / 'v100')
    p.add_argument('--heat-t-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / 't100')
    p.add_argument('--s-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '10hpa')
    p.add_argument('--t-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '500hpa')
    p.add_argument('--fig3-archive', type=Path, default=Path(__file__).with_name('fig3_fig4_results') / 'fig3_frozen_predictions.npz')
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('fig5_results'))
    p.add_argument('--heat-only', action='store_true')
    p.add_argument('--plot-only', action='store_true', help='Redraw Fig5 from saved fig5_scores.npz without recomputing data')
    p.add_argument('--refresh-heat', action='store_true')
    p.add_argument('--first-winter', type=int, default=1948)
    p.add_argument('--last-winter', type=int, default=2022)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.plot_only:
        with np.load(args.output / 'fig5_scores.npz') as saved:
            plot_fig5(args.output / 'Fig5', saved['gains'], saved['intervals'], saved['C'])
        return
    H = daily_heat_flux(args.v_dir, args.heat_t_dir, args.output / 'h100_daily.npz', args.refresh_heat)
    if args.heat_only:
        return
    with np.load(args.fig3_archive) as old:
        np.testing.assert_array_equal(old['dates'], base.DATES)
        old_observed = old['observed'].copy()
        old_all = old['all_winter_sums'].copy()
        old_cond = old['conditional_winter_sums'].copy()
        old_pred_all = old['pred_all'].copy()
        old_pred_cond = old['pred_conditional'].copy()
    all_sums = np.full((75, 60, len(MODEL_NAMES), 7), np.nan)
    cond_sums = np.full_like(all_sums, np.nan)
    for (out_sums, saved) in ((all_sums, old_all), (cond_sums, old_cond)):
        out_sums[:, :, 0] = saved[:, :, 0]
        out_sums[:, :, 2] = saved[:, :, 1]
    del old_all, old_cond
    valid_h = np.isfinite(H)
    totals_h = np.bincount(base.CAL_KEY[valid_h], weights=H[valid_h], minlength=366)
    counts_h = np.bincount(base.CAL_KEY[valid_h], minlength=366)
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
    models_new = [1, 3] + list(range(4, 28))
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
        if not np.allclose(yy[own_target], old_observed[own_target], rtol=0, atol=1e-10):
            raise ValueError(f'Winter {w}: SPV differs from Fig3 archive')
        h40 = fold_h40(H, w, totals_h, counts_h)
        for h in range(1, 61):
            (dst, X, y) = make_fig5_design(ss, tt, yy, h40, h)
            winter = base.WINTER[dst]
            target = winter == w
            if target.sum() == 0:
                continue
            hold_dst = dst[target]
            if not (np.isfinite(old_pred_all[hold_dst, h - 1, :]).all() and np.isfinite(old_pred_cond[hold_dst, h - 1, :]).all()):
                raise ValueError(f'Winter {w}, h={h}: Fig3 mask mismatch')
            old_all_n = int(all_sums[wi, h - 1, 0, 0])
            if int(target.sum()) != old_all_n:
                raise ValueError(f'Winter {w}, h={h}: pair count {target.sum()} != {old_all_n}')
            total_st = fig3.stats(X, y)
            ssw_st = fig3.stats(X[labels[winter - 1948]], y[labels[winter - 1948]])
            non_st = fig3.stats(X[~labels[winter - 1948]], y[~labels[winter - 1948]])
            held_st = fig3.stats(X[target], y[target])
            train_all = fig3.subtract(total_st, held_st)
            train_cond = fig3.subtract(ssw_st if labels[wi] else non_st, held_st)
            for (st, dest) in ((train_all, all_sums), (train_cond, cond_sums)):
                for m in models_new:
                    (alpha, beta) = fit_columns(st, COLS[m])
                    pred = alpha + X[target][:, COLS[m]] @ beta
                    dest[wi, h - 1, m] = fig3.score_sums(y[target], pred)
        if (w - args.first_winter) % 5 == 0 or w == args.last_winter:
            print(f'Fig5 models complete through held-out winter {w}', flush=True)
    if args.first_winter != 1948 or args.last_winter != 2022:
        np.savez_compressed(args.output / 'fig5_partial_diagnostic.npz', all_sums=all_sums, conditional_sums=cond_sums, model_names=np.array(MODEL_NAMES))
        return
    if not (np.isfinite(all_sums).all() and np.isfinite(cond_sums).all()):
        raise ValueError('Some winter/lead/model scores remain missing')
    (metrics, bands, C, counts) = score_fig5(all_sums, cond_sums)
    expected_h15 = np.array([[0.01308967, 0.06727554, 0.00336721], [0.02101925, 0.07607823, 0.00546155], [-0.00126533, 0.05081174, -0.00249221]])
    expected_c_s6_h10 = np.array([0.0108935834, 0.01868539, 0.0023668219])
    np.testing.assert_allclose(metrics[:, :, 14], expected_h15, rtol=0, atol=5e-07)
    np.testing.assert_allclose(C[:, 4, 9], expected_c_s6_h10, rtol=0, atol=5e-07)
    np.testing.assert_array_equal(counts[:, 14], [12729, 6636, 6093])
    plot_fig5(args.output / 'Fig5', metrics, bands, C)
    np.savez_compressed(args.output / 'fig5_scores.npz', h=np.arange(1, 61), model_names=np.array(MODEL_NAMES), all_winter_sums=all_sums, conditional_winter_sums=cond_sums, gains=metrics, intervals=bands, C=C, counts=counts)
    with (args.output / 'fig5_gains.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['group', 'h', 'n', 'delta_H', 'delta_M', 'delta_H_given_M', 'delta_H_lo', 'delta_H_hi', 'delta_M_lo', 'delta_M_hi', 'delta_H_given_M_lo', 'delta_H_given_M_hi'])
        for (g, name) in enumerate(('All', 'SSW', 'non-SSW')):
            for hi in range(60):
                writer.writerow([name, hi + 1, counts[g, hi], *metrics[g, :, hi], *bands[g, 0, hi], *bands[g, 1, hi], *bands[g, 2, hi]])
    with (args.output / 'fig5_C.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['group', 'removed_mode', 'h', 'C'])
        for (g, name) in enumerate(('All', 'SSW', 'non-SSW')):
            for (j, mode) in enumerate(MODE_NAMES):
                for hi in range(60):
                    writer.writerow([name, mode, hi + 1, C[g, j, hi]])
    checkpoints = {'h15_gains': {name: metrics[g, :, 14].tolist() for (g, name) in enumerate(('All', 'SSW', 'non-SSW'))}, 'h10_C_S6': {name: float(C[g, 4, 9]) for (g, name) in enumerate(('All', 'SSW', 'non-SSW'))}, 'h15_counts': counts[:, 14].tolist(), 'daily_heat_count': int(np.isfinite(H).sum())}
    (args.output / 'fig5_QA.json').write_text(json.dumps(checkpoints, indent=2), encoding='utf-8')
    print(json.dumps(checkpoints, indent=2), flush=True)
if __name__ == '__main__':
    main()
