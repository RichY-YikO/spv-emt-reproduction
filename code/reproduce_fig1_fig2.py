from __future__ import annotations
from result_paths import reference_dir, validation_dir, publish_figure
import argparse
import csv
import json
import os
from pathlib import Path
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(_name, '4')
import numpy as np
import xarray as xr
from scipy.linalg import eigh
START = np.datetime64('1948-07-01', 'D')
STOP = np.datetime64('2023-07-01', 'D')
DATES = np.arange(START, STOP, np.timedelta64(1, 'D'))
YEARS = np.arange(1948, 2024)
LAT = np.arange(90.0, 64.9, -2.5)
LON = np.arange(-180.0, 180.0, 2.5)
ANCHORS = {'S': [(87.5, -125, -1), (70, -117.5, -1), (67.5, 130, +1), (65, 5, -1), (65, 50, +1), (65, -55, -1)], 'T': [(87.5, 140, +1), (75, 57.5, +1), (72.5, -7.5, -1), (70, 95, +1), (72.5, -170, -1), (67.5, 60, -1), (65, -172.5, +1)]}
TARGET_FRAC = {'S': [71.2290375, 13.53716076, 9.21706744, 1.81732623, 1.59358208, 1.48726614], 'T': [34.618625, 14.11098934, 10.33605102, 7.0138291, 6.6897466, 5.46338466, 2.81400795]}
TARGET_LAG = {'S': [0, 21, 28, 25, 14, 7], 'T': [0, 40, 10, 5, 55, 30, 9]}
TARGET_R = {'S': [-0.9793716982, -0.2126080489, 0.1467663583, -0.0593467101, 0.0443956955, -0.2068961306], 'T': [0.1719116161, 0.1212370352, -0.1490875841, -0.0486052585, -0.1041460274, -0.0934302017, -0.1339798666]}

def date_parts(dates: np.ndarray):
    days = dates.astype('datetime64[D]')
    years = days.astype('datetime64[Y]').astype(int) + 1970
    months = days.astype('datetime64[M]').astype(int) % 12 + 1
    keys = np.array([int((np.datetime64(f'2000-{str(d)[5:]}') - np.datetime64('2000-01-01')) / np.timedelta64(1, 'D')) for d in days], dtype=np.int16)
    winter = years - (months < 7)
    return (years, months, keys, winter)
(CAL_YEAR, CAL_MONTH, CAL_KEY, WINTER) = date_parts(DATES)
VALID_WINTER = (WINTER >= 1948) & (WINTER <= 2022)
BASIS_MASK = VALID_WINTER & ((CAL_MONTH >= 9) | (CAL_MONTH <= 4)) & ~((CAL_MONTH == 9) & (np.array([int(str(d)[8:]) for d in DATES]) == 1)) & ~(CAL_KEY == 59)
TARGET_MASK = VALID_WINTER & ((CAL_MONTH >= 11) | (CAL_MONTH <= 4)) & ~(CAL_KEY == 59)
PROJECT_MASK = VALID_WINTER & ((CAL_MONTH >= 7) | (CAL_MONTH <= 4))

def read_annual_fields(folder: Path, level: str) -> np.ndarray:
    out = np.full((len(DATES), len(LAT) * len(LON)), np.nan, dtype=np.float64)
    for year in YEARS:
        path = folder / f'{year}.nc'
        if not path.is_file():
            raise FileNotFoundError(path)
        with xr.open_dataset(path, mask_and_scale=True, decode_times=level == 'T') as ds:
            if 'z' not in ds:
                raise ValueError(f'{path}: no variable z')
            lat_name = 'latitude' if level == 'S' else 'lat'
            lon_name = 'longitude' if level == 'S' else 'lon'
            np.testing.assert_allclose(ds[lat_name].values, LAT, rtol=0, atol=1e-09)
            np.testing.assert_allclose(ds[lon_name].values, LON, rtol=0, atol=1e-09)
            field = np.asarray(ds['z'].values, dtype=np.float64)
            if field.shape[1:] != (11, 144):
                raise ValueError(f'{path}: unexpected grid {field.shape}')
            if level == 'S':
                if field.shape[0] != 365 or not np.array_equal(ds.time.values, np.arange(1, 366)):
                    raise ValueError(f'{path}: expected anonymous days 1..365')
                file_dates = np.arange(np.datetime64(f'{year}-01-01'), np.datetime64(f'{year + 1}-01-01'), np.timedelta64(1, 'D'))
                file_dates = file_dates[np.array([str(d)[5:] != '02-29' for d in file_dates])]
            else:
                file_dates = ds.time.values.astype('datetime64[D]')
                if len(np.unique(file_dates)) != len(file_dates):
                    raise ValueError(f'{path}: duplicate calendar dates')
                if field.shape[0] != len(file_dates):
                    raise ValueError(f'{path}: time mismatch')
            use = (file_dates >= START) & (file_dates < STOP)
            pos = (file_dates[use] - START).astype('timedelta64[D]').astype(int)
            out[pos] = field[use].reshape(-1, 1584)
    if np.isinf(out).any():
        raise ValueError(f'{level}: infinite source values')
    return out

def anomaly_and_climatology(raw: np.ndarray):
    clim = np.full((366, raw.shape[1]), np.nan, dtype=np.float64)
    for key in np.unique(CAL_KEY):
        sub = raw[CAL_KEY == key]
        count = np.isfinite(sub).sum(axis=0)
        total = np.nansum(sub, axis=0)
        clim[key] = np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)
    return (raw - clim[CAL_KEY], clim)

def make_spv(anom_s: np.ndarray):
    weights = np.broadcast_to(np.cos(np.deg2rad(LAT))[:, None], (11, 144)).ravel()
    if np.isnan(anom_s[TARGET_MASK]).any():
        raise ValueError('Missing 10 hPa grid values on a nominal SPV target date')
    polar = anom_s @ weights / weights.sum()
    residual = np.full(len(DATES), np.nan, dtype=np.float64)
    slopes = {}
    origin = np.datetime64('1985-01-01', 'D')
    x = ((DATES - origin) / np.timedelta64(1, 'D')).astype(float) / 365.2425
    for month in (11, 12, 1, 2, 3, 4):
        use = TARGET_MASK & (CAL_MONTH == month) & np.isfinite(polar)
        if not use.any():
            raise ValueError(f'No SPV values in month {month}')
        beta = np.linalg.lstsq(np.column_stack((np.ones(use.sum()), x[use])), polar[use], rcond=None)[0]
        residual[use] = polar[use] - (beta[0] + beta[1] * x[use])
        slopes[str(month)] = beta.tolist()
    use = TARGET_MASK & np.isfinite(residual)
    (mean, scale) = (residual[use].mean(), residual[use].std(ddof=0))
    Y = (residual - mean) / scale
    return (Y, polar, residual, mean, scale, slopes)

def fit_basis(anom: np.ndarray, level: str):
    X = anom[BASIS_MASK]
    if X.shape != (18075, 1584) or not np.isfinite(X).all():
        raise ValueError(f'{level}: bad basis matrix {X.shape}; finite={np.isfinite(X).all()}')
    mu = X.mean(axis=0)
    centered = X - mu
    gram = centered.T @ centered
    trace = float(np.trace(gram))
    (lam, U) = eigh(gram, subset_by_index=(1484, 1583), driver='evr', check_finite=True, overwrite_a=True)
    (lam, U) = (lam[::-1], U[:, ::-1])
    fractions = 100.0 * lam / trace
    for (j, (lat, lon, wanted)) in enumerate(ANCHORS[level]):
        ilat = int(round((90.0 - lat) / 2.5))
        ilon = int(round((lon + 180.0) / 2.5))
        if U[ilat * 144 + ilon, j] * wanted < 0:
            U[:, j] *= -1
    for j in range(len(ANCHORS[level]), U.shape[1]):
        idx = np.argmax(np.abs(U[:, j]))
        if U[idx, j] < 0:
            U[:, j] *= -1
    scores = np.full((len(DATES), 100), np.nan, dtype=np.float64)
    good = PROJECT_MASK & np.isfinite(anom).all(axis=1)
    scores[good] = (anom[good] - mu) @ U
    if not np.allclose(U.T @ U, np.eye(100), atol=1e-10):
        raise AssertionError(f'{level}: spatial vectors are not orthogonal')
    return (mu, lam / (18075 - 1), U, fractions, scores)

def lag_pairs(scores: np.ndarray, target_idx: np.ndarray, lag: int):
    src_idx = target_idx - lag
    inside = (src_idx >= 0) & (src_idx < len(DATES))
    valid = inside.copy()
    valid[inside] &= WINTER[src_idx[inside]] == WINTER[target_idx[inside]]
    valid[inside] &= np.isfinite(scores[src_idx[inside], 0])
    return (src_idx[valid], valid)

def corr_columns(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    Xc = X - X.mean(axis=0)
    yc = y - y.mean()
    return Xc.T @ yc / np.sqrt(np.sum(Xc * Xc, axis=0) * np.sum(yc * yc))

def make_surrogates(y: np.ndarray, target_dates: np.ndarray, count: int):
    if count == 0:
        return None
    rng = np.random.Generator(np.random.PCG64(2026091206))
    winter = WINTER[(target_dates - START).astype('timedelta64[D]').astype(int)]
    jump = np.flatnonzero((np.diff(target_dates) != np.timedelta64(1, 'D')) | (np.diff(winter) != 0)) + 1
    edges = np.r_[0, jump, len(y)]
    if len(edges) - 1 != 93:
        raise ValueError(f'Expected 93 continuous SPV segments, got {len(edges) - 1}')
    result = np.empty((len(y), count), dtype=np.float64)
    for (a, b) in zip(edges[:-1], edges[1:]):
        n = b - a
        phase = rng.uniform(0.0, 2 * np.pi, size=(n // 2 + 1, count))
        phase[0] = 0.0
        if n % 2 == 0:
            phase[-1] = 0.0
        fourier = np.fft.rfft(y[a:b])[:, None] * np.exp(1j * phase)
        result[a:b] = np.fft.irfft(fourier, n=n, axis=0)
    return result

def observed_and_null(scores, Y, y_sur, level):
    target_idx = np.flatnonzero(TARGET_MASK & np.isfinite(Y))
    observed = np.full((121, 100), np.nan)
    fixed_lo = np.full((121, len(ANCHORS[level])), np.nan)
    fixed_hi = fixed_lo.copy()
    extreme_sur = None if y_sur is None else np.zeros((100, y_sur.shape[1]))
    counts = np.zeros(121, dtype=int)
    for lag in range(-60, 61):
        (source, valid) = lag_pairs(scores, target_idx, lag)
        X = scores[source]
        y = Y[target_idx[valid]]
        observed[lag + 60] = corr_columns(X, y)
        counts[lag + 60] = len(y)
        if y_sur is not None:
            Z = y_sur[valid]
            ncol = 100 if lag >= 0 else len(ANCHORS[level])
            Xc = X[:, :ncol] - X[:, :ncol].mean(axis=0)
            Zc = Z - Z.mean(axis=0)
            R = Xc.T @ Zc / np.sqrt(np.sum(Xc * Xc, axis=0)[:, None] * np.sum(Zc * Zc, axis=0)[None, :])
            fixed_lo[lag + 60] = np.quantile(R[:len(ANCHORS[level])], 0.025, axis=1)
            fixed_hi[lag + 60] = np.quantile(R[:len(ANCHORS[level])], 0.975, axis=1)
            if lag >= 0:
                replace = np.abs(R) > np.abs(extreme_sur)
                extreme_sur[replace] = R[replace]
    positive = observed[60:121]
    chosen = np.argmax(np.abs(positive), axis=0)
    extreme_obs = positive[chosen, np.arange(100)]
    extreme_lo = extreme_hi = None
    if extreme_sur is not None:
        extreme_lo = np.quantile(extreme_sur, 0.025, axis=1)
        extreme_hi = np.quantile(extreme_sur, 0.975, axis=1)
    return (observed, counts, chosen, extreme_obs, fixed_lo, fixed_hi, extreme_lo, extreme_hi)

def save_csv(path, names, rows):
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(names)
        w.writerows(rows)

def plot_figure(path_base: Path, level: str, U, fractions, observed, chosen, extreme_obs, fixed_lo, fixed_hi, extreme_lo, extreme_hi):
    from publication_style import plot_maps as draw_publication
    return draw_publication(path_base, level, U, fractions, observed, chosen, extreme_obs, fixed_lo, fixed_hi, extreme_lo, extreme_hi)

def main():
    p = argparse.ArgumentParser(description="Reproduce SPV, Fig. 1 and Fig. 2 from the annual polar-cap NetCDF fields.\n\nRun with the Python environment that has numpy, scipy, xarray, netCDF4,\nmatplotlib and cartopy.  The optional --surrogates 0 mode is a quick numeric\ncheck; use the default 999 for the report's final null envelopes.\n")
    p.add_argument('--s-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '10hpa')
    p.add_argument('--t-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '500hpa')
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('fig1_fig2_results'))
    p.add_argument('--surrogates', type=int, default=999, help='999 for final report; 0 for rapid numeric and layout check')
    args = p.parse_args()
    if args.surrogates not in (0, 999):
        raise ValueError('Use 0 for a quick check or 999 for the report method')
    args.output.mkdir(parents=True, exist_ok=True)
    reference = reference_dir(args.output)
    reference.mkdir(parents=True, exist_ok=True)
    diagnostics = validation_dir(args.output) / 'fig1_fig2'
    diagnostics.mkdir(parents=True, exist_ok=True)
    render = validation_dir(args.output) / '_render'
    report = {'surrogates': args.surrogates, 'source_s': str(args.s_dir), 'source_t': str(args.t_dir), 'software': {'numpy': np.__version__, 'xarray': xr.__version__}}
    print('Reading 10 hPa annual source fields', flush=True)
    raw_s = read_annual_fields(args.s_dir, 'S')
    missing_leap = np.flatnonzero((CAL_KEY == 59) & np.isnan(raw_s).all(axis=1))
    if len(missing_leap) != 18:
        raise ValueError(f'Expected 18 missing 10 hPa Feb 29 dates, got {len(missing_leap)}')
    (anom_s, clim_s) = anomaly_and_climatology(raw_s)
    del raw_s
    (Y, polar, residual, mean, scale, slopes) = make_spv(anom_s)
    target_idx = np.flatnonzero(TARGET_MASK & np.isfinite(Y))
    if len(target_idx) != 13575:
        raise ValueError(f'Expected 13575 SPV days, got {len(target_idx)}')
    report['spv'] = {'n': len(target_idx), 'residual_mean': float(mean), 'residual_std_ddof0': float(scale), 'report_std': 4247.457942402572, 'std_error': float(scale - 4247.457942402572), 'monthly_intercept_slope': slopes}
    save_csv(reference / 'spv_daily.csv', ['date', 'winter_start_year', 'polar_anomaly', 'detrended_residual', 'spv_standardized'], ((str(DATES[i]), int(WINTER[i]), polar[i], residual[i], Y[i]) for i in target_idx))
    y_sur = make_surrogates(Y[target_idx], DATES[target_idx], args.surrogates)
    for (level, anom) in (('S', anom_s), ('T', None)):
        if level == 'T':
            print('Reading 500 hPa annual source fields', flush=True)
            raw_t = read_annual_fields(args.t_dir, 'T')
            (anom, clim_t) = anomaly_and_climatology(raw_t)
            del raw_t
        print(f'Computing {level} top-100 EMT basis', flush=True)
        (mu, lam, U, frac, scores) = fit_basis(anom, level)
        if level == 'S':
            del anom_s
        else:
            del anom
        print(f'Computing {level} correlations and null envelopes', flush=True)
        (obs, counts, lag, r, lo, hi, elo, ehi) = observed_and_null(scores, Y, y_sur, level)
        np.savez_compressed((reference if level == 'S' else args.output) / f'{level}_basis_and_correlations.npz', U=U, eigenvalue=lam, fraction_percent=frac, mu=mu, scores=scores, correlation_lags=np.arange(-60, 61), correlation=obs, pair_counts=counts, maxabs_lag_0_60=lag, signed_maxabs_r=r, fixed_null_lo=lo, fixed_null_hi=hi, extreme_null_lo=np.array([]) if elo is None else elo, extreme_null_hi=np.array([]) if ehi is None else ehi, dates=DATES)
        if level == 'T':
            save_csv(args.output / f'{level}_top100.csv', ['mode', 'variance_percent', 'maxabs_lag_0_60', 'signed_r', 'null_2p5', 'null_97p5'], ((f'{level}{j + 1}', frac[j], int(lag[j]), r[j], '' if elo is None else elo[j], '' if ehi is None else ehi[j]) for j in range(100)))
        report[level] = {'basis_days': int(BASIS_MASK.sum()), 'target_fractions': TARGET_FRAC[level], 'computed_fractions': frac[:len(ANCHORS[level])].tolist(), 'fraction_errors': (frac[:len(ANCHORS[level])] - TARGET_FRAC[level]).tolist(), 'target_lags': TARGET_LAG[level], 'computed_lags': lag[:len(ANCHORS[level])].tolist(), 'target_r': TARGET_R[level], 'computed_r': r[:len(ANCHORS[level])].tolist(), 'r_errors': (r[:len(ANCHORS[level])] - TARGET_R[level]).tolist(), 'counts_m60_m30_0_30_60': counts[[0, 30, 60, 90, 120]].tolist(), 'anchors': []}
        for (j, (lat, lon, wanted)) in enumerate(ANCHORS[level]):
            (ilat, ilon) = (int(round((90 - lat) / 2.5)), int(round((lon + 180) / 2.5)))
            value = float(U[ilat * 144 + ilon, j])
            report[level]['anchors'].append({'mode': f'{level}{j + 1}', 'lat': lat, 'lon': lon, 'required_sign': wanted, 'computed_value': value})
        if level == 'S':
            spectrum = np.abs(np.fft.rfft(U[:, :6].T.reshape(6, 11, 144), axis=2) / 144) ** 2
            spectrum[:, :, 1:72] *= 2
            power = spectrum.sum(axis=1)
            fractions_wave = 100 * np.column_stack((power[:, 0], power[:, 1], power[:, 2], power[:, 3:].sum(axis=1))) / power.sum(axis=1)[:, None]
            save_csv(args.output / 'table_S1.csv', ['mode', 'k0_percent', 'k1_percent', 'k2_percent', 'k_ge_3_percent'], ((f'S{j + 1}', *fractions_wave[j]) for j in range(6)))
            report['table_S1'] = fractions_wave.tolist()
        if level == 'S':
            from reproduce_fig1a_no_orientation import run, draw_full
            final_folder = args.output
            run(args.s_dir, final_folder)
            draw_full(render, archive=final_folder / 'S_final_pearson.npz')
            publish_figure(render / 'Fig1', args.output)
        else:
            plot_figure(render / 'Fig2', level, U, frac, obs, lag, r, lo, hi, elo, ehi)
            publish_figure(render / 'Fig2', args.output)
        print(f'{level}: fraction% {np.round(frac[:len(ANCHORS[level])], 8)}', flush=True)
        print(f'{level}: selected lags {lag[:len(ANCHORS[level])].tolist()}', flush=True)
        print(f'{level}: signed r {np.round(r[:len(ANCHORS[level])], 10)}', flush=True)
        del scores, U
    with open(diagnostics / 'qa_report.json', 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print(f'Finished: {args.output}', flush=True)
if __name__ == '__main__':
    main()
