from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import ScaledTranslation
GROUPS = ('All winters', 'SSW winters', 'Normal winters')
COLORS = ('#2868a1', '#c83b3b')
MODES = [f'S{k}' for k in range(2, 7)] + [f'T{k}' for k in range(1, 8)]

def upright_modes():
    return [f'$\\mathrm{{{m[0]}}}_{{{m[1:]}}}$' for m in MODES]

def aligned_panel(ax, index, title='', size=11):
    baseline = 1.055
    ax.text(0, baseline, f'({chr(97 + index)})', transform=ax.transAxes, ha='left', va='baseline', fontsize=size, fontweight='bold', color='black', clip_on=False, gid='panel-label')
    if title:
        ax.text(0.5, baseline, title, transform=ax.transAxes, ha='center', va='baseline', fontsize=11, clip_on=False, gid='group-title')

def boxed_legend(ax, fontsize=9, **kwargs):
    config = dict(loc='best', frameon=True, facecolor='white', framealpha=1, edgecolor='.6', fontsize=fontsize, borderpad=0.25, borderaxespad=0.3, handlelength=1.9, labelspacing=0.16, handletextpad=0.4)
    config.update(kwargs)
    return ax.legend(**config)

def setup():
    plt.rcParams.update({'font.family': 'Times New Roman', 'mathtext.fontset': 'stix', 'font.size': 10.5, 'axes.labelsize': 11, 'axes.titlesize': 11, 'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5, 'legend.fontsize': 9, 'axes.linewidth': 0.8, 'xtick.major.size': 4, 'ytick.major.size': 4, 'xtick.major.width': 0.9, 'ytick.major.width': 0.9, 'xtick.direction': 'out', 'ytick.direction': 'out', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white', 'hatch.linewidth': 0.35, 'axes.axisbelow': True, 'path.simplify': False})

def panel(ax, index, title=''):
    ax.text(0, 1.045, f'({chr(97 + index)})', transform=ax.transAxes, ha='left', va='bottom', fontsize=11, fontweight='bold', color='black', clip_on=False, gid='panel-label')
    if title:
        ax.set_title(title, pad=12, x=0.57)

def legend(ax, handles=None, labels=None, loc='best', fontsize=9, handlelength=1.9):
    kw = dict(loc=loc, frameon=True, facecolor='white', framealpha=0.83, edgecolor='none', borderpad=0.28, borderaxespad=0.3, handlelength=handlelength, labelspacing=0.22, handletextpad=0.4, fontsize=fontsize)
    return ax.legend(handles, labels, **kw) if handles is not None else ax.legend(**kw)

def curve_axes(ax, xlabel=None, ylabel=None):
    ax.grid(color='.88', linewidth=0.4, alpha=0.65)
    ax.tick_params(colors='.12', pad=3)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

def limits(values, bands=None, reference=None):
    arr = np.asarray(values).ravel()
    if bands is not None:
        arr = np.r_[arr, np.asarray(bands).ravel()]
    arr = arr[np.isfinite(arr)]
    if reference is not None:
        arr = np.r_[arr, reference]
    (lo, hi) = (float(arr.min()), float(arr.max()))
    pad = max(0.025, (hi - lo) * 0.07)
    return (lo - pad, hi + pad)

def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for artist in fig.findobj():
        if hasattr(artist, 'set_rasterized'):
            artist.set_rasterized(False)
    fig.canvas.draw()
    for ax in fig.axes:
        for collection in ax.collections:
            collection.set_rasterized(False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    checks = []
    for ax in fig.axes:
        box = ax.get_window_extent(renderer)
        leg = ax.get_legend()
        if leg:
            lb = leg.get_window_extent(renderer)
            assert box.contains(lb.x0, lb.y0) and box.contains(lb.x1, lb.y1), f'Legend outside axes: axes={box.bounds}, legend={lb.bounds}'
        for text in ax.texts:
            if text.get_gid() == 'panel-label':
                tb = text.get_window_extent(renderer)
                assert tb.y0 >= box.y1, 'Panel label must be outside axes'
                checks.append(text.get_text())
    assert not fig.legends, 'Figure-level legends are prohibited'
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight', pad_inches=0.12)
    fig.savefig(path.with_suffix('.png'), dpi=600, bbox_inches='tight', pad_inches=0.12)
    path.with_suffix('.layout.json').write_text(json.dumps({'panels': checks, 'legends_inside': True, 'rasterization_disabled': True, 'font': 'Times New Roman', 'width_inches': fig.get_figwidth()}, indent=2), encoding='utf-8')
    plt.close(fig)

def score_grid(path, metrics, bands, colors=COLORS, labels=('B', 'M')):
    setup()
    (fig, axes) = plt.subplots(2, 3, figsize=(7.2, 4.8), layout='constrained', sharex=True, sharey='row')
    x = np.arange(1, 61)
    for row in range(2):
        yrange = limits(metrics[:, row], bands[:, row], 0 if row == 0 else 1)
        for g in range(3):
            ax = axes[row, g]
            for (m, label) in enumerate(labels):
                ax.fill_between(x, bands[g, row, m, :, 0], bands[g, row, m, :, 1], color=colors[m], alpha=0.16, lw=0)
                ax.plot(x, metrics[g, row, m], color=colors[m], lw=1.5, label=label)
            ax.axhline(0 if row == 0 else 1, color='.5', lw=0.7, ls='--')
            ax.set(xlim=(1, 60), ylim=yrange, xticks=[1, 15, 30, 45, 60])
            panel(ax, row * 3 + g, GROUPS[g] if row == 0 else '')
            curve_axes(ax, 'Lead (days)' if row == 1 else None, ('Correlation ($r$)' if row == 0 else 'MSE ratio ($Q$)') if g == 0 else None)
    legend(axes[0, 0])
    save(fig, path)

def plot_fig3(path, metrics, intervals):
    setup()
    (fig, axes) = plt.subplots(2, 3, figsize=(7.2, 4.2), layout='constrained', sharex=True, sharey='row')
    x = np.arange(1, 61)
    for row in range(2):
        yrange = (-0.1, 1.02) if row == 0 else (-0.02, 1.1)
        assert np.nanmin(intervals[:, row, ..., 0]) >= yrange[0]
        assert np.nanmax(intervals[:, row, ..., 1]) <= yrange[1]
        for g in range(3):
            ax = axes[row, g]
            for (m, label) in enumerate(('Baseline', 'Full model')):
                ax.fill_between(x, intervals[g, row, m, :, 0], intervals[g, row, m, :, 1], color=COLORS[m], alpha=0.16, lw=0)
                ax.plot(x, metrics[g, row, m], color=COLORS[m], lw=1.5, label=label)
            ax.axhline(0 if row == 0 else 1, color='.5', lw=0.7, ls='--')
            ax.set(xlim=(1, 60), ylim=yrange, xticks=[1, 10, 20, 30, 40, 50, 60])
            ax.tick_params(labelbottom=True)
            aligned_panel(ax, row * 3 + g, GROUPS[g] if row == 0 else '')
            curve_axes(ax, ylabel=('Correlation ($r$)' if row == 0 else 'MSE ratio ($Q$)') if g == 0 else None)
    boxed_legend(axes[0, 0])
    fig.supxlabel('Forecast lead (days)', fontsize=11)
    save(fig, path)

def plot_s3(path, metrics, bands):
    from supplement_style import plot_s3 as draw_si
    return draw_si(path, metrics, bands)

def plot_fig4(path, beta, pv):
    from matplotlib.patches import Rectangle
    setup()
    (fig, axes) = plt.subplots(2, 3, figsize=(7.2, 5.2), layout='constrained', sharex=True, sharey=True)
    (xe, ye) = (np.arange(0.5, 13.5), np.arange(0.5, 61.5))
    for g in range(3):
        (top, bottom) = axes[:, g]
        mesh1 = top.pcolormesh(xe, ye, np.abs(beta[g]), shading='flat', cmap='RdPu', vmin=0, vmax=30, rasterized=False)
        mesh2 = bottom.pcolormesh(xe, ye, np.log10(pv[g]), shading='flat', cmap='plasma', vmin=-14, vmax=0, rasterized=False)
        for (h, k) in np.argwhere(pv[g] <= 0.05):
            bottom.add_patch(Rectangle((k + 0.5, h + 0.5), 1, 1, facecolor='none', edgecolor=(1, 1, 1, 0.45), lw=0, hatch='////', zorder=3))
        for (row, ax) in enumerate((top, bottom)):
            ax.set(xlim=(0.5, 12.5), ylim=(0.5, 60.5), yticks=[1, 15, 30, 45, 60])
            aligned_panel(ax, row * 3 + g, GROUPS[g] if row == 0 else '', size=14)
            if g == 0:
                ax.set_ylabel('Forecast lead (days)')
        bottom.set_xticks(np.arange(1, 13), upright_modes(), rotation=0, ha='center', fontsize=9.5)
    cb1 = fig.colorbar(mesh1, ax=axes[0, :], pad=0.02, fraction=0.035)
    cb1.set_label('$|\\beta|$', labelpad=4)
    cb1.solids.set_edgecolor('face')
    cb2 = fig.colorbar(mesh2, ax=axes[1, :], pad=0.02, fraction=0.035)
    cb2.set_label('$\\log_{10}(p)$', labelpad=4)
    cb2.solids.set_edgecolor('face')
    fig.supxlabel('Predictor eigen microstate', fontsize=11)
    save(fig, path)

def plot_fig5(path, gains, bands, C):
    setup()
    colors = ('#2868a1', '#c83b3b', '#743d9a')
    labels = ('$Q(B)-Q(BH)$', '$Q(B)-Q(M)$', '$Q(M)-Q(MH)$')
    fig = plt.figure(figsize=(7.2, 4.4))
    axes = np.empty((2, 3), dtype=object)
    width = 1.8
    height = width / 1.3
    for (row, bottom) in enumerate((2.645, 0.56)):
        for (col, left) in enumerate((0.65, 2.7, 4.75)):
            axes[row, col] = fig.add_axes([left / 7.2, bottom / 4.4, width / 7.2, height / 4.4], sharey=axes[row, 0] if col else None)
            if col:
                axes[row, col].tick_params(labelleft=False)
    x = np.arange(1, 61)
    for g in range(3):
        (top, bottom) = axes[:, g]
        top.set_box_aspect(1 / 1.3)
        bottom.set_box_aspect(1 / 1.3)
        for m in range(3):
            top.fill_between(x, bands[g, m, :, 0], bands[g, m, :, 1], color=colors[m], alpha=0.16, lw=0)
            top.plot(x, gains[g, m], color=colors[m], lw=1.5, label=labels[m])
        top.axhline(0, color='.5', lw=0.7, ls='--')
        top.set(xlim=(1, 60), ylim=limits(gains, bands, 0), xticks=[1, 15, 30, 45, 60])
        top.set_yticks([-0.05, 0, 0.05, 0.1, 0.15])
        curve_axes(top, ylabel='MSE ratio reduction' if g == 0 else None)
        aligned_panel(top, g, GROUPS[g])
        mesh = bottom.pcolormesh(np.arange(13), np.arange(0.5, 61.5), C[g].T, cmap='RdBu_r', vmin=-0.02, vmax=0.02, shading='flat', rasterized=False)
        bottom.axvline(5, color='.5', lw=0.7)
        bottom.set(xlim=(0, 12), ylim=(0.5, 60.5), yticks=[1, 10, 20, 30, 40, 50, 60])
        bottom.set_xticks(np.arange(12) + 0.5, upright_modes(), rotation=0, ha='center', fontsize=9.5)
        if g == 0:
            bottom.set_ylabel('Forecast lead (days)')
        aligned_panel(bottom, g + 3)
    boxed_legend(axes[0, 0], fontsize=8)
    cax = fig.add_axes([6.7 / 7.2, 0.56 / 4.4, 0.08 / 7.2, height / 4.4])
    colorbar = fig.colorbar(mesh, cax=cax, label='$C_j$', ticks=[-0.02, -0.01, 0, 0.01, 0.02])
    colorbar.solids.set_edgecolor('face')
    fig.text((0.65 + 6.55) / 2 / 7.2, 2.18 / 4.4, 'Forecast lead (days)', ha='center', va='center', fontsize=11, gid='top-shared-xlabel')
    fig.supxlabel('Removed microstate', fontsize=11, y=0.015)
    save(fig, path)

def plot_monthly(path, metrics, bands, q_upper):
    from supplement_style import plot_monthly as draw_si
    return draw_si(path, metrics, bands, q_upper)

def plot_s7(path, metrics, bands):
    from supplement_style import plot_s7 as draw_si
    return draw_si(path, metrics, bands)

def plot_s8(path, years, metrics):
    from supplement_style import plot_s8 as draw_si
    return draw_si(path, years, metrics)

def plot_s9(path, metrics, bands):
    from supplement_style import plot_s9 as draw_si
    return draw_si(path, metrics, bands)

def plot_lags(path, prefix, lags, r, lo, hi):
    from supplement_style import plot_lags as draw_si
    return draw_si(path, prefix, lags, r, lo, hi)

def plot_maps(path, level, U, fractions, observed, chosen, extreme_obs, fixed_lo, fixed_hi, extreme_lo, extreme_hi):
    from publication_revision import plot_maps as draw_revision
    return draw_revision(path, level, U, fractions, observed, chosen, extreme_obs, fixed_lo, fixed_hi, extreme_lo, extreme_hi)
