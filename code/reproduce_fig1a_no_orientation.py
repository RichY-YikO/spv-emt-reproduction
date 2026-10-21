from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import reproduce_fig1_fig2 as base
import numpy as np
from result_paths import reference_dir, validation_dir, pearson_reports, publish_figure
ROOT = Path(__file__).resolve().parent
STATISTIC_KIND = 'pearson_negated_raw_scores_v1'

def require_pearson_cache(path):
    with np.load(path, allow_pickle=False) as z:
        if 'statistic_kind' not in z or str(z['statistic_kind'].item()) != STATISTIC_KIND:
            raise ValueError('Final Fig1/FigS1 require a Pearson cache; run --from-cache first.')
        if 'none' not in str(z['orientation'].item()) or int(z['surrogates']) != 999:
            raise ValueError('Expected unoriented raw basis and 999 surrogates.')

def load_spv(results=None):
    Y = np.full(len(base.DATES), np.nan)
    with (reference_dir(results) / 'spv_daily.csv').open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            i = int((np.datetime64(row['date'], 'D') - base.START) / np.timedelta64(1, 'D'))
            Y[i] = float(row['spv_standardized'])
    return Y

def recompute_cached(source, output):
    (source, output) = (Path(source), Path(output))
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    with np.load(source, allow_pickle=False) as z:
        data = {key: z[key].copy() for key in ['U_raw', 'scores_raw', 'mu', 'eigenvalue', 'fraction_percent', 'dates', 'orientation']}
        old_r = z['correlation'].copy()
        old_counts = z['pair_counts'].copy()
    if 'none' not in str(data['orientation'].item()):
        raise ValueError('Cache must contain the final unoriented raw eigensolver output.')
    np.testing.assert_array_equal(data['dates'], base.DATES)
    scores = -data['scores_raw']
    Y = load_spv(source.parent)
    target = np.flatnonzero(base.TARGET_MASK & np.isfinite(Y))
    assert len(target) == 13575
    surrogate = base.make_surrogates(Y[target], base.DATES[target], 999)
    print('Recomputing 100 Pearson lag curves and 999 phase surrogates.', flush=True)
    (obs, counts, lags, signed, lo, hi, elo, ehi) = base.observed_and_null(scores, Y, surrogate, 'S')
    np.testing.assert_array_equal(counts, old_counts)
    data.update(scores_negated=scores, correlation_lags=np.arange(-60, 61), correlation=obs, pair_counts=counts, maxabs_lag_0_60=lags, signed_maxabs_r=signed, fixed_null_lo=lo, fixed_null_hi=hi, extreme_null_lo=elo, extreme_null_hi=ehi, statistic=np.asarray('Pearson correlation; negated scores centered by their own paired mean'), statistic_kind=np.asarray(STATISTIC_KIND), surrogates=np.asarray(999), seed=np.asarray(2026091206))
    np.savez_compressed(output / 'S_final_pearson.npz', **data)
    base.save_csv(output / 'S_top100.csv', ['mode', 'variance_percent', 'maxabs_lag_0_60', 'signed_r', 'null_2p5', 'null_97p5'], ((f'S{j + 1}', data['fraction_percent'][j], int(lags[j]), signed[j], elo[j], ehi[j]) for j in range(100)))
    report = {'statistic_kind': STATISTIC_KIND, 'orientation_applied': False, 'score_multiplier': -1, 'mean': 'own mean of NEGATED scores on each lag paired sample', 'surrogates': 999, 'seed': 2026091206, 'source_cache_sha256': source_hash, 'decomposition_recomputed': False, 'first_six_r': signed[:6].tolist(), 'first_six_lags': lags[:6].tolist(), 'max_change_from_previous': float(np.max(np.abs(obs - old_r))), 'Fig1_envelope': 'signed maximum absolute correlation reselected over lags 0..60 for each surrogate', 'FigS1_envelope': 'fixed-lag 2.5/97.5 percentiles, no lag selection', 'elapsed_seconds': time.monotonic() - started}
    (pearson_reports(output) / 'calculation_and_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)

def run(source, output):
    output.mkdir(parents=True, exist_ok=True)
    protected = list(ROOT.glob('*_results/*.npz')) + list(ROOT.glob('*_results/*.pdf')) + list(ROOT.glob('*_results/*.png'))
    protected += [p for p in (ROOT / 'JGR_Atmospheres_vector_figures').rglob('*') if p.suffix in ('.pdf', '.png')]
    protected += list(reference_dir(output).glob('*.npz'))
    protected = [p for p in protected if p.resolve() != (output / 'S_final_pearson.npz').resolve()]
    hashes = {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    started = time.monotonic()
    print('Reading original 10 hPa annual fields; reference files stay unchanged.', flush=True)
    raw = base.read_annual_fields(source, 'S')
    (anom, clim) = base.anomaly_and_climatology(raw)
    del raw, clim
    X = anom[base.BASIS_MASK]
    assert X.shape == (18075, 1584) and np.isfinite(X).all()
    mu = X.mean(axis=0)
    centered = X - mu
    gram = centered.T @ centered
    trace = float(np.trace(gram))
    del X, centered
    print('Computing the same eigh/evr top-100 eigenvectors; BOTH orientation loops omitted.', flush=True)
    (lam, U) = base.eigh(gram, subset_by_index=(1484, 1583), driver='evr', check_finite=True, overwrite_a=True)
    (lam, U) = (lam[::-1], U[:, ::-1])
    raw_scores = np.full((len(base.DATES), 100), np.nan)
    good = base.PROJECT_MASK & np.isfinite(anom).all(axis=1)
    raw_scores[good] = (anom[good] - mu) @ U
    del anom, gram
    scores = -raw_scores
    Y = load_spv(output)
    target = np.flatnonzero(base.TARGET_MASK & np.isfinite(Y))
    assert len(target) == 13575
    surrogate = base.make_surrogates(Y[target], base.DATES[target], 999)
    (obs, counts, lags, signed, lo, hi, elo, ehi) = base.observed_and_null(scores, Y, surrogate, 'S')
    fractions = 100 * lam / trace
    np.savez_compressed(output / 'S_final_pearson.npz', U_raw=U, scores_raw=raw_scores, scores_negated=scores, mu=mu, eigenvalue=lam / (18075 - 1), fraction_percent=fractions, dates=base.DATES, correlation_lags=np.arange(-60, 61), correlation=obs, pair_counts=counts, maxabs_lag_0_60=lags, signed_maxabs_r=signed, fixed_null_lo=lo, fixed_null_hi=hi, extreme_null_lo=elo, extreme_null_hi=ehi, orientation=np.asarray('none; unmodified scipy eigh evr output'), statistic=np.asarray('Pearson correlation; negated scores centered by their own paired mean'), statistic_kind=np.asarray(STATISTIC_KIND), surrogates=np.asarray(999), seed=np.asarray(2026091206))
    base.save_csv(output / 'S_top100.csv', ['mode', 'variance_percent', 'maxabs_lag_0_60', 'signed_r', 'null_2p5', 'null_97p5'], ((f'S{j + 1}', fractions[j], int(lags[j]), signed[j], elo[j], ehi[j]) for j in range(100)))
    with np.load(reference_dir(output) / 'S_basis_and_correlations.npz') as old:
        overlap = np.sum(U * old['U'], axis=0)
        np.testing.assert_allclose(np.abs(overlap), 1, rtol=0, atol=1e-10)
        np.testing.assert_allclose(lam / (18075 - 1), old['eigenvalue'], rtol=1e-12)
        np.testing.assert_allclose(fractions, old['fraction_percent'], rtol=1e-12)
        expected_scores = old['scores'] * np.sign(overlap)[None, :]
        np.testing.assert_allclose(raw_scores, expected_scores, rtol=1e-10, atol=1e-06, equal_nan=True)
        expected = -old['correlation'] * np.sign(overlap)[None, :]
        np.testing.assert_allclose(obs, expected, rtol=0, atol=1e-11)
    np.testing.assert_allclose(U.T @ U, np.eye(100), rtol=0, atol=1e-10)
    report = {'orientation_applied': False, 'eigensolver': 'scipy.linalg.eigh, driver=evr', 'score_multiplier': -1, 'mean': 'own mean of NEGATED scores on each lag paired sample', 'statistic_kind': STATISTIC_KIND, 'surrogates': 999, 'seed': 2026091206, 'first_six_r': signed[:6].tolist(), 'first_six_lags': lags[:6].tolist(), 'first_six_raw_vs_oriented_sign': np.sign(overlap[:6]).astype(int).tolist(), 'original_files_unchanged': len(hashes), 'elapsed_seconds': time.monotonic() - started}
    for (rel, h) in hashes.items():
        assert hashlib.sha256(Path(rel).read_bytes()).hexdigest() == h, rel
    (pearson_reports(output) / 'calculation_and_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)

def draw(output, signed, low, high):
    import publication_style as style
    from publication_revision import glossy_points, GlossMarker, HandlerGlossMarker
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FormatStrFormatter
    style.setup()
    fig = style.plt.figure(figsize=(7.2, 2.7))
    ax = fig.add_axes([0.66 / 7.2, 0.5 / 2.7, 6.32 / 7.2, 2.0 / 2.7])
    ranks = np.arange(1, 101)
    ax.fill_between(ranks, low, high, color='.83', lw=0)
    ax.axvspan(0.5, 6.5, facecolor='none', edgecolor='.78', hatch='////', lw=0, zorder=0)
    ax.plot(ranks, signed, color='.25', lw=0.85)
    glossy_points(ax, ranks[6:], signed[6:])
    glossy_points(ax, ranks[:6], signed[:6], red=True)
    ax.axhline(0, color='.5', lw=0.7)
    ax.set(xlim=(0.5, 100.5), ylim=(-0.25, 1), yticks=np.round(np.arange(-0.2, 1.01, 0.2), 8), xlabel='Rank of microstates', ylabel='Identified correlation')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    style.curve_axes(ax)
    style.panel(ax, 0)
    style.boxed_legend(ax, handles=[Line2D([], [], color='.83', lw=5), GlossMarker()], labels=['95% CI', 'Selected states'], loc='upper right', bbox_to_anchor=(0.985, 0.985), fontsize=10.5, handlelength=2.3, handler_map={GlossMarker: HandlerGlossMarker()})
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), signed)
    style.save(fig, output / 'Fig1a_no_orientation')

def draw_full(output, archive=None):
    from publication_revision import plot_maps
    output = Path(output)
    archive = Path(archive) if archive is not None else ROOT / 'fig1_fig2_results' / 'S_final_pearson.npz'
    require_pearson_cache(archive)
    with np.load(archive, allow_pickle=False) as z:
        plot_maps(output / 'Fig1', 'S', -z['U_raw'], z['fraction_percent'], z['correlation'], z['maxabs_lag_0_60'], z['signed_maxabs_r'], z['fixed_null_lo'], z['fixed_null_hi'], z['extreme_null_lo'], z['extreme_null_hi'], correlation_ylim=(-0.25, 1))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Standalone Fig1(a): eigensolver output without either sign-orientation rule.\n\nFinal Fig1/FigS1 use ordinary Pearson correlation of negated raw scores.\nEach lag is centered using the mean of its own valid, NEGATED paired scores.\nUse --from-cache to recompute inference without reading fields or decomposing.\n')
    parser.add_argument('--s-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / 'raw' / '10hpa')
    parser.add_argument('--output', type=Path, default=ROOT / 'fig1_fig2_results')
    parser.add_argument('--from-cache', type=Path, help='Reuse raw eigenvectors/scores; recompute ONLY Pearson inference.')
    args = parser.parse_args()
    if args.from_cache:
        recompute_cached(args.from_cache, args.output)
    else:
        run(args.s_dir, args.output)
        render = validation_dir(args.output) / '_render'
        draw_full(render, archive=args.output / 'S_final_pearson.npz')
        publish_figure(render / 'Fig1', args.output)
