from pathlib import Path
import argparse
import json
import reproduce_fig1_fig2 as base
import numpy as np
from result_paths import reference_dir, pearson_reports
from reproduce_fig1a_no_orientation import require_pearson_cache, load_spv
ROOT = Path(__file__).resolve().parent

def verify(previous=None):
    path = ROOT / 'fig1_fig2_results' / 'S_final_pearson.npz'
    require_pearson_cache(path)
    with np.load(path) as a:
        z = {k: a[k] for k in a.files}
    raw = z['scores_raw']
    x = -raw
    np.testing.assert_array_equal(x, z['scores_negated'])
    np.testing.assert_array_equal(z['dates'], base.DATES)
    np.testing.assert_allclose(z['U_raw'].T @ z['U_raw'], np.eye(100), atol=1e-10, rtol=0)
    if previous:
        with np.load(previous) as old:
            for key in ['U_raw', 'scores_raw', 'scores_negated', 'mu', 'eigenvalue', 'fraction_percent', 'dates', 'pair_counts']:
                np.testing.assert_array_equal(z[key], old[key], err_msg=key)
    Y = load_spv()
    target = np.flatnonzero(base.TARGET_MASK & np.isfinite(Y))
    check = np.empty_like(z['correlation'])
    raw_check = check.copy()
    independent_null_checks = 0
    sur = base.make_surrogates(Y[target], base.DATES[target], 999)
    for (k, h) in enumerate(range(-60, 61)):
        (source, valid) = base.lag_pairs(x, target, h)
        yy = Y[target[valid]]
        yc = yy - yy.mean()
        assert len(yy) == z['pair_counts'][k]
        for j in range(100):
            xx = x[source, j]
            xc = xx - xx.mean()
            rr = raw[source, j]
            rc = rr - rr.mean()
            denom = np.sqrt(np.dot(xc, xc) * np.dot(yc, yc))
            check[k, j] = np.dot(xc, yc) / denom
            raw_check[k, j] = np.dot(rc, yc) / np.sqrt(np.dot(rc, rc) * np.dot(yc, yc))
        if h in [-60, -30, 0, 17, 60]:
            yy_sur = sur[valid]
            yc_sur = yy_sur - yy_sur.mean(axis=0)
            for j in [0, 3, 5]:
                xx = x[source, j]
                xc = xx - xx.mean()
                rr = (xc[:, None] * yc_sur).sum(axis=0) / np.sqrt((xc * xc).sum() * (yc_sur * yc_sur).sum(axis=0))
                np.testing.assert_allclose(np.quantile(rr, [0.025, 0.975]), [z['fixed_null_lo'][k, j], z['fixed_null_hi'][k, j]], atol=2e-12, rtol=0)
                independent_null_checks += 1
    np.testing.assert_allclose(z['correlation'], check, atol=2e-12, rtol=0)
    np.testing.assert_allclose(z['correlation'], -raw_check, atol=2e-12, rtol=0)
    chosen = np.argmax(np.abs(check[60:]), axis=0)
    np.testing.assert_array_equal(chosen, z['maxabs_lag_0_60'])
    np.testing.assert_allclose(check[60:][chosen, np.arange(100)], z['signed_maxabs_r'], atol=2e-12, rtol=0)
    with np.load(reference_dir() / 'S_basis_and_correlations.npz') as old:
        overlaps = np.sum(z['U_raw'] * old['U'], axis=0)
        np.testing.assert_allclose(np.abs(overlaps), 1, atol=1e-10, rtol=0)
        sign = -np.sign(overlaps)
        np.testing.assert_allclose(z['correlation'], old['correlation'] * sign, atol=2e-11, rtol=0)
        for (prefix, n) in [('fixed', 6), ('extreme', 100)]:
            a = old[prefix + '_null_lo']
            b = old[prefix + '_null_hi']
            expected_lo = np.where(sign[:n] > 0, a, -b)
            expected_hi = np.where(sign[:n] > 0, b, -a)
            np.testing.assert_allclose(z[prefix + '_null_lo'], expected_lo, atol=2e-11, rtol=0)
            np.testing.assert_allclose(z[prefix + '_null_hi'], expected_hi, atol=2e-11, rtol=0)
    result = {'statistic': 'standard Pearson of negated raw scores', 'all_12100_correlations_verified': True, 'negative_of_raw_Pearson': True, 'paired_dates_and_counts_unchanged': True, 'raw_basis_scores_and_spatial_directions_unchanged': bool(previous), 'all_fixed_and_maxabs_envelopes_match_independent_Pearson_reference': True, 'independent_999_surrogate_fixed_lag_spot_checks': independent_null_checks, 'max_direct_formula_error': float(np.max(np.abs(check - z['correlation'])))}
    (pearson_reports() / 'Pearson_validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Focused numerical audit of final Fig1/FigS1; no plotting or decomposition.\n\nOptional --previous CACHE checks unchanged raw arrays against the backed-up\ncustom-statistic cache. Uses independently expressed Pearson formulas.\n')
    p.add_argument('--previous', type=Path)
    verify(p.parse_args().previous)
