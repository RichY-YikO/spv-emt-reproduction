from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import reproduce_fig1_fig2 as f12
import reproduce_fig3 as f3
import reproduce_fig4 as f4
import reproduce_fig5 as f5
import reproduce_figS1_figS2 as fs12
import reproduce_figS3 as fs3
import reproduce_figS4_S5_S6 as fm
import reproduce_figS7_S8 as fs78
import reproduce_figS9 as fs9
ROOT = Path(__file__).resolve().parent
FIGURES = [f'Fig{i}' for i in range(1, 6)] + [f'FigS{i}' for i in range(1, 10)]

def replot(destination, figures=None):
    destination = Path(destination).resolve()
    selected = set(FIGURES if figures is None else figures)
    if not selected or not selected.issubset(FIGURES):
        raise ValueError('Invalid figure selection')
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in list(ROOT.glob('*_results/*.npz')) + list(ROOT.glob('reference_inputs/*.npz'))}
    outputs = []

    def target(folder, name):
        p = destination / folder / name
        outputs.append(str(p.relative_to(destination)))
        print('Drawing', folder, name, flush=True)
        return p
    for (level, name) in [('S', 'Fig1'), ('T', 'Fig2')]:
        if name in selected:
            if level == 'S':
                from reproduce_fig1a_no_orientation import draw_full
                draw_full(target('fig1_fig2_results', name).parent)
            else:
                with np.load(ROOT / 'fig1_fig2_results' / 'T_basis_and_correlations.npz') as z:
                    f12.plot_figure(target('fig1_fig2_results', name), level, z['U'], z['fraction_percent'], z['correlation'], z['maxabs_lag_0_60'], z['signed_maxabs_r'], z['fixed_null_lo'], z['fixed_null_hi'], z['extreme_null_lo'], z['extreme_null_hi'])
        if ('FigS1' if level == 'S' else 'FigS2') in selected:
            (lags, r, lo, hi, counts, _) = fs12.load_figure_archive(ROOT / 'fig1_fig2_results', level, 6 if level == 'S' else 7)
            base = target('figS1_figS2_results', 'FigS1' if level == 'S' else 'FigS2')
            fs12.draw(base, level, lags, r, lo, hi)
            fs12.export_csv(base.parent / ('figS1_correlations.csv' if level == 'S' else 'figS2_correlations.csv'), level, lags, r, lo, hi, counts)
    if 'Fig3' in selected:
        with np.load(ROOT / 'fig3_fig4_results' / 'fig3_frozen_predictions.npz') as z:
            f3.plot_fig3(target('fig3_fig4_results', 'Fig3'), z['metrics'], z['intervals'])
    if 'Fig4' in selected:
        with np.load(ROOT / 'fig3_fig4_results' / 'fig4_coefficients_CR1.npz') as z:
            f4.plot_fig4(target('fig3_fig4_results', 'Fig4'), z['beta'], z['p'])
    if 'Fig5' in selected:
        with np.load(ROOT / 'fig5_results' / 'fig5_scores.npz') as z:
            f5.plot_fig5(target('fig5_results', 'Fig5'), z['gains'], z['intervals'], z['C'])
    for folder in ['figS3_S9_results']:
        if 'FigS3' in selected:
            with np.load(ROOT / folder / 'FigS3_scores.npz') as z:
                fs3.plot(target(folder, 'FigS3'), z['metrics'], z['intervals'])
        for name in ['FigS4', 'FigS5', 'FigS6']:
            archive = ROOT / folder / f'{name}_monthly_scores.npz'
            if name in selected and archive.exists():
                with np.load(archive) as z:
                    fm.plot(target(folder, name), z['metrics'], z['intervals'], float(np.nanmax(z['intervals'][:, 1, ..., 1])) + 0.05)
    folder = 'figS3_S9_results'
    if 'FigS7' in selected:
        with np.load(ROOT / folder / 'FigS7_scores.npz') as z:
            np.testing.assert_array_equal(z['leads'], fs78.LEADS)
            np.testing.assert_array_equal(z['tau'], fs78.TAU)
            fs78.plot_s7(target(folder, 'FigS7'), z['metrics'], z['intervals'])
    if 'FigS8' in selected:
        with np.load(ROOT / folder / 'FigS8_scores.npz') as z:
            np.testing.assert_array_equal(z['leads'], fs78.LEADS)
            fs78.plot_s8(target(folder, 'FigS8'), z['midpoint_year'], z['metrics'])
    if 'FigS9' in selected:
        with np.load(ROOT / folder / 'FigS9_scores.npz') as z:
            fs9.plot(target(folder, 'FigS9'), z['metrics'], z['intervals'])
    after = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in list(ROOT.glob('*_results/*.npz')) + list(ROOT.glob('reference_inputs/*.npz'))}
    assert hashes == after, 'Scientific data archives changed'
    manifest_path = destination / 'replot_manifest.json'
    previous = json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.exists() else {}
    all_outputs = list(dict.fromkeys(previous.get('figures', []) + outputs))
    manifest_path.write_text(json.dumps({'figures': all_outputs, 'last_redraw': outputs, 'selection': sorted(selected), 'data_sha256': hashes, 'archives_unchanged': True}, indent=2), encoding='utf-8')
    collection = destination / 'JGR_Atmospheres_vector_figures'
    (collection / 'PNG_preview').mkdir(parents=True, exist_ok=True)
    for rel in outputs:
        if 'figS3_S4_S5_results' in rel:
            continue
        base = destination / rel
        shutil.copy2(base.with_suffix('.pdf'), collection / base.with_suffix('.pdf').name)
        shutil.copy2(base.with_suffix('.png'), collection / 'PNG_preview' / base.with_suffix('.png').name)
    print('Completed:', len(outputs), 'figure files; all cached data unchanged.', flush=True)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Replot cached results without scientific recomputation.\n\npython replot_publication.py --plot-only [--output-root DIRECTORY]\nDefault: write intermediate redraws under outputs/replot_workspace.\n')
    parser.add_argument('--plot-only', action='store_true', help='Cached drawing only (always enabled).')
    parser.add_argument('--output-root', type=Path, default=ROOT.parent / 'outputs' / 'replot_workspace')
    parser.add_argument('--figures', nargs='+', choices=FIGURES, help='Only redraw these figure types.')
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    replot(args.output_root, args.figures)
