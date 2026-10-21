import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib.transforms import ScaledTranslation
from matplotlib.path import Path as MplPath
from publication_style import setup, curve_axes, panel, save, boxed_legend

class GlossMarker:
    pass

class HandlerGlossMarker(HandlerBase):

    def create_artists(self, legend, handle, xdescent, ydescent, width, height, fontsize, trans):
        x = width / 2 - xdescent
        y = height / 2 - ydescent
        artists = []
        for (diameter, color, dx, dy) in [(4, '#8c1425', 0, 0), (3.2, '#e5263d', 0, 0), (1.63, '#ffe4df', -0.6, 0.6)]:
            artists.append(Line2D([x + dx], [y + dy], marker='o', markersize=diameter, markerfacecolor=color, markeredgecolor='none', linestyle='none', transform=trans))
        return artists

def glossy_points(ax, x, y, red=False):
    scale = 4 / 9
    (edge, body, shine) = ('#8c1425', '#e5263d', '#ffe4df') if red else ('#080808', '#333333', '#dedede')
    ax.scatter(x, y, s=36 * scale, facecolor=edge, edgecolor=edge, lw=0.3, zorder=4)
    ax.scatter(x, y, s=23 * scale, facecolor=body, edgecolor='none', zorder=5)
    transform = ax.transData + ScaledTranslation(-0.6 / 72, 0.6 / 72, ax.figure.dpi_scale_trans)
    ax.scatter(x, y, s=6 * scale, facecolor=shine, edgecolor='none', transform=transform, zorder=6)

def plot_maps(path, level, U, fractions, observed, chosen, extreme_obs, fixed_lo, fixed_hi, extreme_lo, extreme_hi, *, correlation_ylim=None):
    import cartopy.crs as ccrs
    from scipy.interpolate import RectBivariateSpline
    setup()
    n = 6 if level == 'S' else 7
    cols = 3 if level == 'S' else 4
    (fw, fh) = (7.2, 7.0)
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([0.66 / fw, (fh - 2.2) / fh, 6.32 / fw, 2.0 / fh])
    ranks = np.arange(1, 101)
    shown = np.asarray(extreme_obs)
    if extreme_lo is not None:
        low = np.asarray(extreme_lo)
        high = np.asarray(extreme_hi)
        assert np.all(low <= high)
        ax.fill_between(ranks, low, high, color='.83', lw=0)
    ax.axvspan(0.5, n + 0.5, facecolor='none', edgecolor='.78', hatch='////', lw=0, zorder=0)
    ax.plot(ranks, shown, color='.25', lw=0.85)
    glossy_points(ax, ranks[n:], shown[n:])
    glossy_points(ax, ranks[:n], shown[:n], red=True)
    ax.axhline(0, color='.5', lw=0.7)
    ticks = np.round(np.arange(-0.2, 1.01, 0.2) if level == 'S' else np.arange(-0.18, 0.181, 0.06), 8)
    ax.set(xlim=(0.5, 100.5), ylim=correlation_ylim if correlation_ylim is not None else (-0.2, 1) if level == 'S' else (-0.18, 0.18), yticks=ticks, xlabel='Rank of microstates', ylabel='Identified correlation')
    from matplotlib.ticker import FormatStrFormatter
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f' if level == 'S' else '%.2f'))
    curve_axes(ax)
    panel(ax, 0)
    boxed_legend(ax, handles=[Line2D([], [], color='.83', lw=5), GlossMarker()], labels=['95% CI', 'Selected states'], loc='upper right', bbox_to_anchor=(0.985, 0.985), fontsize=10.5, handlelength=2.3, handler_map={GlossMarker: HandlerGlossMarker()})
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    bounds = np.linspace(-0.02, 0.02, 11)
    cmap = plt.get_cmap('Spectral_r')
    norm = mcolors.BoundaryNorm(bounds, cmap.N, extend='both')
    theta = np.linspace(0, 2 * np.pi, 200)
    circle = MplPath(np.vstack((np.sin(theta), np.cos(theta))).T * 0.5 + 0.5)
    lat = np.arange(90, 64.9, -2.5)
    lon = np.arange(-180, 180, 2.5)
    latfine = np.linspace(65, 90, 101)
    lonfine = np.linspace(-180, 180, 361)
    lonext = np.r_[lon[-4:] - 360, lon, lon[:4] + 360]
    diameter = 1.56
    gap = 0.18
    start = (fw - (cols * diameter + (cols - 1) * gap)) / 2
    import matplotlib.patheffects as pe
    for j in range(n):
        (row, col) = divmod(j, cols)
        left = start + col * (diameter + gap) + ((diameter + gap) / 2 if n == 7 and row == 1 else 0)
        bottom = 2.4 if row == 0 else 0.45
        a = fig.add_axes([left / fw, bottom / fh, diameter / fw, diameter / fh], projection=proj)
        native = U[:, j].reshape(11, 144)[::-1]
        periodic = np.concatenate((native[:, -4:], native, native[:, :4]), axis=1)
        smooth = RectBivariateSpline(lat[::-1], lonext, periodic, kx=3, ky=3)(latfine, lonfine)
        mesh = a.contourf(lonfine, latfine, smooth, levels=bounds, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, extend='both')
        a.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
        a.set_boundary(circle, transform=a.transAxes)
        a.coastlines(resolution='50m', lw=0.45, color='.2')
        a.gridlines(crs=ccrs.PlateCarree(), xlocs=np.arange(-180, 180, 30), ylocs=np.arange(65, 90, 5), linewidth=0.35, color='.5', alpha=0.6, linestyle=(0, (3, 3)))
        for (longitude, text) in [(-120, '120°W'), (-60, '60°W'), (120, '120°E'), (60, '60°E'), (180, '180°'), (0, '0°')]:
            angle = np.deg2rad(longitude)
            radius = 0.5 + 5.5 / (diameter * 72)
            (px, py) = (0.5 + radius * np.sin(angle), 0.5 - radius * np.cos(angle))
            rotation = (longitude + 90) % 180 - 90
            annotation = a.annotate(text, (px, py), xycoords='axes fraction', ha='center', va='center', rotation=rotation, rotation_mode='anchor', fontsize=7.5, clip_on=False, gid='longitude-label')
            annotation.set_path_effects([pe.Stroke(linewidth=1.2, foreground='white'), pe.Normal()])
        for latitude in (70, 80, 90):
            (px, py) = proj.transform_point(0, latitude, ccrs.PlateCarree())
            text = a.annotate(f'{latitude}°N', (px, py), xycoords='data', xytext=(0, 2), textcoords='offset points', fontsize=6.5, ha='center', va='bottom', color='.18', clip_on=True, gid='latitude-label')
            text.set_path_effects([pe.Stroke(linewidth=1, foreground='white', alpha=0.65), pe.Normal()])
        panel(a, j + 1)
    cax = fig.add_axes([0.66 / fw, 0.12 / fh, 6.32 / fw, 0.07 / fh])
    fig.colorbar(mesh, cax=cax, orientation='horizontal', ticks=bounds[::2], extend='both')
    save(fig, path)
