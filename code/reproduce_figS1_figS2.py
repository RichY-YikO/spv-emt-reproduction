from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from result_paths import validation_dir

def load_archive(path: Path, prefix: str, n_modes: int):
    with np.load(path) as z:
        lags = np.asarray(z['correlation_lags'], dtype=int)
        dates = z['dates']
        r = np.asarray(z['correlation'][:, :n_modes], dtype=float)
        lo = np.asarray(z['fixed_null_lo'], dtype=float)
        hi = np.asarray(z['fixed_null_hi'], dtype=float)
        n = np.asarray(z['pair_counts'], dtype=int)
    np.testing.assert_array_equal(lags, np.arange(-60, 61))
    if r.shape != (121, n_modes) or lo.shape != r.shape or hi.shape != r.shape:
        raise ValueError(f'{prefix}: incomplete observed/null lag curves')
    if not (np.isfinite(r).all() and np.isfinite(lo).all() and np.isfinite(hi).all()):
        raise ValueError(f'{prefix}: non-finite correlation or null envelope')
    return (lags, r, lo, hi, n, dates)

def figure_archive_path(input_folder: Path, prefix: str):
    input_folder = Path(input_folder)
    if prefix == 'S':
        folder = input_folder
        path = folder / 'S_final_pearson.npz'
        if not path.is_file():
            raise FileNotFoundError('FigS1 requires final no-orientation Fig1 cache: ' + str(path))
        with np.load(path) as z:
            assert 'none' in str(z['orientation'].item())
        from reproduce_fig1a_no_orientation import require_pearson_cache
        require_pearson_cache(path)
        return path
    return input_folder / 'T_basis_and_correlations.npz'

def load_figure_archive(input_folder: Path, prefix: str, n_modes: int):
    return load_archive(figure_archive_path(input_folder, prefix), prefix, n_modes)

def draw(path: Path, prefix: str, lags, r, lo, hi):
    from publication_style import plot_lags as draw_publication
    return draw_publication(path, prefix, lags, r, lo, hi)

def export_csv(path: Path, prefix: str, lags, r, lo, hi, n):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['mode', 'lag_days', 'n_pairs', 'pearson_r', 'null_2p5', 'null_97p5'])
        for j in range(r.shape[1]):
            for (k, lag) in enumerate(lags):
                writer.writerow([f'{prefix}{j + 1}', int(lag), int(n[k]), r[k, j], lo[k, j], hi[k, j]])

def main():
    p = argparse.ArgumentParser(description='Plot S1 from the final no-orientation Fig1 Pearson cache.\n\nS1 uses negated raw scores centered by their own paired mean (Pearson).\nS2 retains the original T Pearson correlations. Both gray bands are saved\nfixed-lag phase-surrogate null envelopes, not lag-selection envelopes.\nNo preprocessing, decomposition or surrogate inference is recomputed here.\n')
    p.add_argument('--input', type=Path, default=Path(__file__).with_name('fig1_fig2_results'))
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('figS1_figS2_results'))
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    qa_path = validation_dir(args.input) / 'fig1_fig2' / 'qa_report.json'
    if qa_path.is_file() and json.loads(qa_path.read_text(encoding='utf-8'))['surrogates'] != 999:
        raise ValueError('Fig. S1/S2 require the 999-surrogate All75 Fig. 1/2 archive')
    report = {}
    dates_ref = None
    for (prefix, n_modes, expected) in (('S', 6, (9057, 13575)), ('T', 7, (9075, 13575))):
        (lags, r, lo, hi, n, dates) = load_figure_archive(args.input, prefix, n_modes)
        if dates_ref is not None:
            np.testing.assert_array_equal(dates_ref, dates)
        dates_ref = dates
        np.testing.assert_array_equal([n.min(), n.max()], expected)
        draw(args.output / f"Fig{('S1' if prefix == 'S' else 'S2')}", prefix, lags, r, lo, hi)
        export_csv(args.output / f"fig{('S1' if prefix == 'S' else 'S2')}_correlations.csv", prefix, lags, r, lo, hi, n)
        report[prefix] = {'n_pairs_min_max': [int(n.min()), int(n.max())], 'max_signed_r_lag': [int(lags[np.argmax(r[:, j])]) for j in range(n_modes)], 'min_signed_r_lag': [int(lags[np.argmin(r[:, j])]) for j in range(n_modes)], 'r_at_zero': r[60].tolist(), 'statistic': 'Pearson correlation', 'gray_band': 'fixed-lag phase-surrogate null envelope', 'source': str(figure_archive_path(args.input, prefix))}
    with np.load(figure_archive_path(args.input, 'S')) as s:
        report['S1_r_lag0'] = float(s['correlation'][60, 0])
    with np.load(args.input / 'T_basis_and_correlations.npz') as t:
        report['T3_r_lag10'] = float(t['correlation'][70, 2])
        report['T2_r_lag40'] = float(t['correlation'][100, 1])
    np.testing.assert_allclose(report['S1_r_lag0'], 0.9793716982451013, atol=1e-11, rtol=0)
    np.testing.assert_allclose(report['T3_r_lag10'], -0.1490875841, atol=1e-09, rtol=0)
    np.testing.assert_allclose(report['T2_r_lag40'], 0.1212370352, atol=1e-09, rtol=0)
    (args.output / 'figS1_figS2_QA.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
if __name__ == '__main__':
    main()
