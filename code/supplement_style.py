import numpy as np
from publication_style import setup, save, boxed_legend, curve_axes, limits, COLORS, GROUPS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
LEADS = [1, 10, 20, 30, 40, 50, 60]

def finish(fig, path):
    fig.canvas.draw()
    fig.set_layout_engine(None)
    save(fig, path)

def si_panel(ax, index, title=''):
    ax.text(0, 1.06, f'({chr(97 + index)})', transform=ax.transAxes, fontsize=14, fontweight='bold', ha='left', va='baseline', clip_on=False, gid='panel-label')
    if title:
        ax.text(0.5, 1.06, title, transform=ax.transAxes, fontsize=11, ha='center', va='baseline', clip_on=False, gid='group-title')

def si_legend(ax, **kwargs):
    options = dict(fontsize=8, handlelength=2.15)
    options.update(kwargs)
    boxed_legend(ax, **options)

def six_grid():
    setup()
    (fig, axes) = plt.subplots(2, 3, figsize=(7.2, 4.4), layout='constrained', sharex=True, sharey='row')
    for ax in axes.flat:
        ax.set_box_aspect(1 / 1.3)
        ax.tick_params(labelbottom=True)
    return (fig, axes)

def score_grid(path, metrics, bands, colors, labels):
    (fig, axes) = six_grid()
    x = np.arange(1, 61)
    for row in range(2):
        yrange = limits(metrics[:, row], bands[:, row], 0 if row == 0 else 1)
        for g in range(3):
            ax = axes[row, g]
            for (m, label) in enumerate(labels):
                ax.fill_between(x, bands[g, row, m, :, 0], bands[g, row, m, :, 1], color=colors[m], alpha=0.16, lw=0)
                ax.plot(x, metrics[g, row, m], color=colors[m], lw=1.5, label=label)
            ax.axhline(0 if row == 0 else 1, color='.5', lw=0.7, ls='--')
            ax.set(xlim=(1, 60), ylim=yrange, xticks=LEADS)
            si_panel(ax, row * 3 + g, GROUPS[g] if row == 0 else '')
            curve_axes(ax, ylabel=('Correlation ($r$)' if row == 0 else 'MSE ratio ($Q$)') if g == 0 else None)
    si_legend(axes[0, 0])
    fig.supxlabel('Forecast lead (days)', fontsize=11)
    finish(fig, path)

def plot_s3(path, metrics, bands):
    score_grid(path, metrics, bands, ('#2868a1', '#238447', '#e27b18'), ('Baseline', 'Baseline + S', 'Baseline + T'))

def plot_monthly(path, metrics, bands, q_upper):
    setup()
    months = ('November', 'December', 'January', 'February', 'March', 'April')
    (fig, axes) = plt.subplots(4, 3, figsize=(7.2, 8.0), layout='constrained', sharex=True)
    x = np.arange(1, 61)
    for (mi, month) in enumerate(months):
        col = mi % 3
        for metric in range(2):
            row = (0 if mi < 3 else 2) + metric
            ax = axes[row, col]
            ax.set_box_aspect(1 / 1.3)
            for m in range(2):
                ax.fill_between(x, bands[mi, metric, m, :, 0], bands[mi, metric, m, :, 1], color=COLORS[m], alpha=0.16, lw=0)
                ax.plot(x, metrics[mi, metric, m], color=COLORS[m], lw=1.5, label=('Baseline', 'Full model')[m])
            ax.axhline(0 if metric == 0 else 1, color='.5', lw=0.7, ls='--')
            ax.set(xlim=(1, 60), ylim=limits(metrics[:, metric], bands[:, metric], 0 if metric == 0 else 1), xticks=LEADS)
            ax.tick_params(labelbottom=True)
            si_panel(ax, row * 3 + col, month if metric == 0 else '')
            curve_axes(ax, ylabel=('Correlation ($r$)' if metric == 0 else 'MSE ratio ($Q$)') if col == 0 else None)
    si_legend(axes[0, 0])
    fig.supxlabel('Forecast lead (days)', fontsize=11)
    finish(fig, path)

def plot_s7(path, metrics, bands):
    (fig, axes) = six_grid()
    x = np.arange(-30, 31)
    for row in range(2):
        for (li, h) in enumerate((7, 15, 30)):
            ax = axes[row, li]
            for m in range(2):
                ax.fill_between(x, bands[li, row, m, :, 0], bands[li, row, m, :, 1], color=COLORS[m], alpha=0.16, lw=0)
                ax.plot(x, metrics[li, row, m], color=COLORS[m], lw=1.5, label=('Baseline', 'Full model')[m])
            ax.axhline(0 if row == 0 else 1, color='.5', ls='--', lw=0.7)
            ax.axvline(0, color='.4', ls=':', lw=0.8)
            ax.axvline(h, color='.4', ls='--', lw=0.8)
            ax.set(xlim=(-30, 30), xticks=[-30, -15, 0, 15, 30], ylim=limits(metrics[:, row], bands[:, row], 0 if row == 0 else 1))
            si_panel(ax, row * 3 + li, f'{h}-day lead' if row == 0 else '')
            curve_axes(ax, ylabel=('Correlation ($r$)' if row == 0 else 'MSE ratio ($Q$)') if li == 0 else None)
    si_legend(axes[0, 0])
    fig.supxlabel('Target offset $\\tau$ (days)', fontsize=11)
    finish(fig, path)

def plot_s8(path, years, metrics):
    (fig, axes) = six_grid()
    for row in range(2):
        for (li, h) in enumerate((7, 15, 30)):
            ax = axes[row, li]
            for m in range(2):
                ax.plot(years, metrics[li, row, m], color=COLORS[m], lw=1.5, label=('Baseline', 'Full model')[m])
            ax.axhline(0 if row == 0 else 1, color='.5', ls='--', lw=0.7)
            ax.axvline(1979, color='.4', ls=':', lw=0.8)
            ax.set(xlim=(float(np.min(years)) - 0.5, float(np.max(years)) + 0.5), xticks=[1960, 1979, 2000, 2018], ylim=limits(metrics[:, row], reference=0 if row == 0 else 1))
            si_panel(ax, row * 3 + li, f'{h}-day lead' if row == 0 else '')
            curve_axes(ax, ylabel=('Correlation ($r$)' if row == 0 else 'MSE ratio ($Q$)') if li == 0 else None)
    si_legend(axes[0, 0])
    fig.supxlabel('Window midpoint year', fontsize=11)
    finish(fig, path)

def plot_s9(path, metrics, bands):
    import matplotlib.patheffects as pe
    (fig, axes) = six_grid()
    x = np.arange(1, 61)
    for row in range(2):
        for g in range(3):
            ax = axes[row, g]
            for scheme in range(2):
                for m in range(2):
                    ax.fill_between(x, bands[g, scheme, row, m, :, 0], bands[g, scheme, row, m, :, 1], color=COLORS[m], alpha=0.1 if scheme == 0 else 0.06, lw=0, zorder=1)
            for scheme in range(2):
                for m in range(2):
                    (line,) = ax.plot(x, metrics[g, scheme, row, m], color=COLORS[m] if scheme == 0 else ('#123f75', '#941b28')[m], lw=1.6 if scheme == 0 else 1.4, ls='-' if scheme == 0 else (0, (6, 3.5)), zorder=3 + scheme)
                    if scheme == 1:
                        line.set_path_effects([pe.Stroke(linewidth=2.1, foreground='white'), pe.Normal()])
            ax.axhline(0 if row == 0 else 1, color='.5', ls=':', lw=0.7, zorder=0)
            ax.set(xlim=(1, 60), xticks=LEADS, ylim=limits(metrics[:, :, row], bands[:, :, row], 0 if row == 0 else 1))
            si_panel(ax, row * 3 + g, GROUPS[g] if row == 0 else '')
            curve_axes(ax, ylabel=('Correlation ($r$)' if row == 0 else 'MSE ratio ($Q$)') if g == 0 else None)
    handles = [Line2D([], [], color=COLORS[0], lw=1.6), Line2D([], [], color=COLORS[1], lw=1.6), Line2D([], [], color='.25', lw=1.6), Line2D([], [], color='.25', lw=1.4, ls=(0, (6, 3.5)))]
    si_legend(axes[0, 0], handles=handles, labels=['Baseline', 'Full model', 'Sep-Apr EMT', 'DJF EMT'], handlelength=2.15)
    fig.supxlabel('Forecast lead (days)', fontsize=11)
    finish(fig, path)

def plot_lags(path, prefix, lags, r, lo, hi):
    setup()
    n = r.shape[1]
    cols = 3 if n == 6 else 4
    (fw, fh) = (7.2, 4.6 if n == 6 else 3.5)
    fig = plt.figure(figsize=(fw, fh))
    (width, gap) = (1.88, 0.36) if n == 6 else (1.25, 0.46)
    axes = []
    for j in range(n):
        (row, col) = divmod(j, cols)
        left = 0.72 + col * (width + gap) + ((width + gap) / 2 if n == 7 and row == 1 else 0)
        bottom = (2.55 if row == 0 else 0.44) if n == 6 else 1.96 if row == 0 else 0.44
        ax = fig.add_axes([left / fw, bottom / fh, width / fw, width / 1.25 / fh])
        axes.append(ax)
        ax.fill_between(lags, lo[:, j], hi[:, j], color='.85', lw=0)
        ax.plot(lags, r[:, j], color='.15', lw=1.3)
        ax.axvline(lags[np.argmax(r[:, j])], color='#c83b3b', ls=':', lw=1)
        ax.axvline(lags[np.argmin(r[:, j])], color='#2868a1', ls=':', lw=1)
        ax.axhline(0, color='.5', ls=':', lw=0.7)
        ax.set(xlim=(-60, 60), xticks=[-60, -30, 0, 30, 60], ylim=limits(r[:, j], np.stack((lo[:, j], hi[:, j]))))
        si_panel(ax, j, f'$\\mathrm{{{prefix}}}_{{{j + 1}}}$')
        curve_axes(ax)
    fig.supxlabel('Lag (days)', fontsize=11)
    fig.supylabel('Correlation ($r$)', fontsize=11)
    si_legend(axes[0] if n == 6 else axes[4], handles=[Line2D([], [], color='.15', lw=1.3), Patch(facecolor='.85')], labels=['$r$', '95% null envelope'], loc='best' if n == 6 else 'upper right')
    finish(fig, path)
