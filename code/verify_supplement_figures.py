from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
import supplement_style as style
import reproduce_figS1_figS2 as lag_source
ROOT = Path(__file__).resolve().parent
figures = []

def capture(fig, path):
    fig.canvas.draw()
    figures.append(fig)

def band(ax, index, x, low, high):
    vertices = ax.collections[index].get_paths()[0].vertices
    for (t, lo, hi) in zip(x, low, high):
        y = vertices[vertices[:, 0] == t, 1]
        np.testing.assert_array_equal([y.min(), y.max()], [lo, hi])
    assert ax.get_ylim()[0] <= np.min(low)
    assert ax.get_ylim()[1] >= np.max(high)

def geometry(fig, count):
    assert len(fig.axes) == count
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        panel = next((t for t in ax.texts if t.get_gid() == 'panel-label'))
        assert panel.get_fontsize() == 14
        assert panel.get_window_extent(renderer).y0 >= ax.get_window_extent(renderer).y1
        titles = [t for t in ax.texts if t.get_gid() == 'group-title']
        for title in titles:
            assert title.get_position() == (0.5, 1.06)
            assert title.get_verticalalignment() == panel.get_verticalalignment()
            assert not title.get_window_extent(renderer).overlaps(panel.get_window_extent(renderer))
        legend = ax.get_legend()
        if legend:
            box = legend.get_window_extent(renderer)
            axesbox = ax.get_window_extent(renderer)
            assert axesbox.contains(box.x0, box.y0) and axesbox.contains(box.x1, box.y1)
            np.testing.assert_allclose(legend.get_frame().get_facecolor(), [1, 1, 1, 1])
    if fig._supylabel is not None:
        ybox = fig._supylabel.get_window_extent(renderer)
        for ax in fig.axes:
            for tick in ax.get_yticklabels():
                if tick.get_visible():
                    assert not ybox.overlaps(tick.get_window_extent(renderer)), 'Y label overlaps ticks'
    assert fig._supxlabel is not None

def verify():
    style.save = capture
    result = {}
    for (prefix, name, n) in [('S', 'FigS1', 6), ('T', 'FigS2', 7)]:
        (lags, r, lo, hi, counts, dates) = lag_source.load_figure_archive(ROOT / 'fig1_fig2_results', prefix, n)
        if prefix == 'S':
            with np.load(ROOT / 'fig1_fig2_results' / 'S_final_pearson.npz') as z:
                np.testing.assert_array_equal(r, z['correlation'][:, :6])
                np.testing.assert_array_equal(lo, z['fixed_null_lo'])
                np.testing.assert_array_equal(hi, z['fixed_null_hi'])
        style.plot_lags(name, prefix, lags, r, lo, hi)
        fig = figures.pop()
        for (j, ax) in enumerate(fig.axes):
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), r[:, j])
            np.testing.assert_array_equal(ax.get_xticks(), [-60, -30, 0, 30, 60])
            band(ax, 0, lags, lo[:, j], hi[:, j])
        geometry(fig, n)
        plt.close(fig)
        result[name] = {'cache_values_and_fixed_lag_envelope_verified': True, 'panels': n}
    folder = ROOT / 'figS3_S9_results'
    x = np.arange(1, 61)
    with np.load(folder / 'FigS3_scores.npz') as z:
        style.plot_s3('FigS3', z['metrics'], z['intervals'])
        fig = figures.pop()
        for row in range(2):
            for g in range(3):
                ax = fig.axes[row * 3 + g]
                for m in range(3):
                    np.testing.assert_array_equal(ax.lines[m].get_ydata(), z['metrics'][g, row, m])
                    band(ax, m, x, z['intervals'][g, row, m, :, 0], z['intervals'][g, row, m, :, 1])
        geometry(fig, 6)
        plt.close(fig)
        result['FigS3'] = {'all_curves_and_bands_unchanged': True}
    for name in ['FigS4', 'FigS5', 'FigS6']:
        with np.load(folder / (name + '_monthly_scores.npz')) as z:
            style.plot_monthly(name, z['metrics'], z['intervals'], 0)
            fig = figures.pop()
            for mi in range(6):
                for metric in range(2):
                    row = (0 if mi < 3 else 2) + metric
                    col = mi % 3
                    ax = fig.axes[row * 3 + col]
                    for m in range(2):
                        np.testing.assert_array_equal(ax.lines[m].get_ydata(), z['metrics'][mi, metric, m])
                        band(ax, m, x, z['intervals'][mi, metric, m, :, 0], z['intervals'][mi, metric, m, :, 1])
            geometry(fig, 12)
            plt.close(fig)
            result[name] = {'all_curves_and_bands_unchanged': True}
    with np.load(folder / 'FigS7_scores.npz') as z:
        style.plot_s7('FigS7', z['metrics'], z['intervals'])
        fig = figures.pop()
        for row in range(2):
            for li in range(3):
                ax = fig.axes[row * 3 + li]
                for m in range(2):
                    np.testing.assert_array_equal(ax.lines[m].get_ydata(), z['metrics'][li, row, m])
                    band(ax, m, z['tau'], z['intervals'][li, row, m, :, 0], z['intervals'][li, row, m, :, 1])
        geometry(fig, 6)
        plt.close(fig)
        result['FigS7'] = {'all_curves_and_bands_unchanged': True}
    with np.load(folder / 'FigS8_scores.npz') as z:
        style.plot_s8('FigS8', z['midpoint_year'], z['metrics'])
        fig = figures.pop()
        for row in range(2):
            for li in range(3):
                for m in range(2):
                    ax = fig.axes[row * 3 + li]
                    np.testing.assert_array_equal(ax.lines[m].get_ydata(), z['metrics'][li, row, m])
        geometry(fig, 6)
        plt.close(fig)
        result['FigS8'] = {'all_curves_unchanged': True}
    with np.load(folder / 'FigS9_scores.npz') as z:
        style.plot_s9('FigS9', z['metrics'], z['intervals'])
        fig = figures.pop()
        for row in range(2):
            for g in range(3):
                ax = fig.axes[row * 3 + g]
                for scheme in range(2):
                    for m in range(2):
                        i = scheme * 2 + m
                        np.testing.assert_array_equal(ax.lines[i].get_ydata(), z['metrics'][g, scheme, row, m])
                        band(ax, i, x, z['intervals'][g, scheme, row, m, :, 0], z['intervals'][g, scheme, row, m, :, 1])
                        if scheme:
                            assert ax.lines[i].get_linestyle() == '--' and ax.lines[i].get_path_effects()
        geometry(fig, 6)
        plt.close(fig)
        result['FigS9'] = {'all_curves_and_bands_unchanged': True, 'dashed_halo_verified': True}
    report = ROOT.parent / 'validation' / 'si_artists.json'
    report.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    verify()
