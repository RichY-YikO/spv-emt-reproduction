from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import openpyxl
import reproduce_fig1_fig2 as base
import reproduce_fig3 as fig3
from reproduce_fig4 import SSW_YEARS
EVENT_STRINGS = ('1950-03-05', '1952-02-22', '1952-11-19', '1954-12-18', '1955-01-26', '1957-02-04', '1958-02-01', '1960-01-17', '1963-01-27', '1965-12-16', '1966-02-22', '1968-01-07', '1968-11-28', '1969-03-13', '1970-01-02', '1971-01-18', '1971-03-20', '1973-01-31', '1977-01-09', '1979-02-22', '1980-02-29', '1981-03-04', '1984-02-24', '1985-01-01', '1987-01-23', '1987-12-08', '1988-03-14', '1989-02-21', '1998-12-15', '1999-02-26', '2000-03-20', '2001-02-11', '2001-12-30', '2002-02-17', '2003-01-18', '2004-01-05', '2006-01-21', '2007-02-24', '2008-02-22', '2009-01-24', '2010-02-09', '2010-03-24', '2013-01-06', '2018-02-12', '2019-01-01', '2021-01-05', '2023-02-16')
EXPECTED_EVENT_SHA256 = 'f7ed7bea677d7d13b63b5a7a8501e5604f9a3bef36cf7732a46a3d6236c0d238'
DEFAULT_EVENT_FILE = Path(__file__).resolve().parents[1] / 'data' / 'raw' / 'SSWevent.xlsx'
LEADS = np.array([7, 15, 30])
TAU = np.arange(-30, 31)

def events_and_winters(event_file):
    digest = hashlib.sha256(event_file.read_bytes()).hexdigest()
    if digest != EXPECTED_EVENT_SHA256:
        raise ValueError(f'SSWevent.xlsx SHA256 mismatch: {digest}')
    workbook = openpyxl.load_workbook(event_file, read_only=True, data_only=True)
    try:
        sheet = workbook['Sheet1']
        raw = [sheet.cell(row, 2).value for row in range(2, 53)]
    finally:
        workbook.close()
    if len(raw) != 51 or any((not hasattr(value, 'date') for value in raw)):
        raise ValueError('Sheet1!B2:B52 must contain 51 Excel dates')
    raw_dates = np.array([value.date().isoformat() for value in raw], dtype='datetime64[D]')
    in_period = (raw_dates >= np.datetime64('1948-07-01')) & (raw_dates <= np.datetime64('2023-06-30'))
    events = raw_dates[in_period]
    expected = np.array(EVENT_STRINGS, dtype='datetime64[D]')
    np.testing.assert_array_equal(events, expected)
    source_rows = np.arange(2, 53)[in_period]
    winter = base.WINTER[((events - base.START) / np.timedelta64(1, 'D')).astype(int)]
    np.testing.assert_array_equal(np.unique(winter), SSW_YEARS)
    if len(events) != 47 or len(np.unique(winter)) != 39:
        raise ValueError('Appendix A must have 47 events in 39 SSW winters')
    audit = {'source_file': str(event_file), 'sha256': digest, 'sha256_matches_report': True, 'sheet_range': 'Sheet1!B2:B52', 'source_rows': len(raw_dates), 'excluded_dates': [str(d) for d in raw_dates[~in_period]], 'study_period_rows': len(events), 'matches_report_appendix_A_dates': True}
    return (events, winter, source_rows, audit)

def event_sums(observed, pred_cond, events, event_winters):
    sums = np.zeros((39, 61, 3, 2, 7), dtype=np.float64)
    for (event, winter) in zip(events, event_winters):
        wi = int(np.searchsorted(SSW_YEARS, winter))
        for (ti, tau) in enumerate(TAU):
            date = event + np.timedelta64(int(tau), 'D')
            pos = int((date - base.START) / np.timedelta64(1, 'D'))
            if pos < 0 or pos >= len(base.DATES) or base.WINTER[pos] != winter:
                continue
            y = observed[pos]
            for (li, h) in enumerate(LEADS):
                pred = pred_cond[pos, h - 1, :]
                if not (np.isfinite(y) and np.isfinite(pred).all()):
                    continue
                for model in range(2):
                    a = float(pred[model])
                    sums[wi, ti, li, model] += [1.0, y, a, y * y, a * a, y * a, (a - y) ** 2]
    return (events, event_winters, sums)

def event_metrics(sums, boot_weights):
    total = sums.sum(axis=0)
    (r, q) = fig3.summarize(total)
    metrics = np.stack((r.transpose(1, 2, 0), q.transpose(1, 2, 0)), axis=1)
    counts = total[:, :, 0, 0].T.astype(int)
    n_winter = (sums[:, :, :, 0, 0] > 0).sum(axis=0).T
    boot = np.einsum('bw,wtlmk->btlmk', boot_weights, sums, optimize=True)
    (br, bq) = fig3.summarize(boot)
    bands = np.empty((3, 2, 2, 61, 2))
    for (metric, values) in enumerate((br, bq)):
        bands[:, metric] = np.quantile(values, [0.025, 0.975], axis=0).transpose(2, 3, 1, 0)
    return (metrics, bands, counts, n_winter)

def rolling_metrics(all_winter_sums):
    all_years = np.arange(1948, 2023)
    centers = np.array([(all_years[i] + 1 + all_years[i + 9] + 1) / 2 for i in range(66)], dtype=float)
    out = np.empty((66, 3, 2, 7))
    for wi in range(66):
        out[wi] = all_winter_sums[wi:wi + 10, LEADS - 1].sum(axis=0)
    (r, q) = fig3.summarize(out)
    metrics = np.stack((r.transpose(1, 2, 0), q.transpose(1, 2, 0)), axis=1)
    counts = out[:, :, 0, 0].T.astype(int)
    return (centers, metrics, counts, out)

def plot_s7(path, metrics, bands):
    from publication_style import plot_s7 as draw_publication
    return draw_publication(path, metrics, bands)

def plot_s8(path, years, metrics):
    from publication_style import plot_s8 as draw_publication
    return draw_publication(path, years, metrics)

def main():
    p = argparse.ArgumentParser(description='Reproduce report Fig. S7 event-aligned and Fig. S8 rolling-winter scores.\n\nBoth figures score frozen leave-one-winter-out Fig. 3 predictions. Fig. S7\nreads the original SSWevent.xlsx and verifies its hash and Appendix A dates.\n')
    p.add_argument('--fig3-archive', type=Path, default=Path(__file__).with_name('fig3_fig4_results') / 'fig3_frozen_predictions.npz')
    p.add_argument('--event-file', type=Path, default=DEFAULT_EVENT_FILE)
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('figS3_S9_results'))
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with np.load(args.fig3_archive) as z:
        np.testing.assert_array_equal(z['dates'], base.DATES)
        observed = z['observed'].copy()
        pred_cond = z['pred_conditional'].copy()
        all_winter_sums = z['all_winter_sums'].copy()
    (events, event_winters, source_rows, source_audit) = events_and_winters(args.event_file)
    (events, event_winters, e_sums) = event_sums(observed, pred_cond, events, event_winters)
    (e_metrics, e_bands, e_count, e_winters) = event_metrics(e_sums, fig3.draw_weights()[1][1])
    for li in (0, 1):
        np.testing.assert_array_equal([e_count[li].min(), e_count[li].max()], [41, 47])
        np.testing.assert_array_equal([e_winters[li].min(), e_winters[li].max()], [34, 39])
    np.testing.assert_array_equal([e_count[2].min(), e_count[2].max()], [44, 47])
    np.testing.assert_array_equal([e_winters[2].min(), e_winters[2].max()], [37, 39])
    np.testing.assert_array_equal([e_count[2, 30], e_winters[2, 30]], [46, 38])
    np.testing.assert_allclose(e_metrics[2, :, 1, 30], [0.23006, 0.74859], atol=5e-06, rtol=0)
    plot_s7(args.output / 'FigS7', e_metrics, e_bands)
    np.savez_compressed(args.output / 'FigS7_scores.npz', events=events, event_winters=event_winters, tau=TAU, leads=LEADS, metrics=e_metrics, intervals=e_bands, n_events=e_count, n_winters=e_winters, winter_sums=e_sums)
    with (args.output / 'FigS7_scores.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['h', 'tau', 'model', 'n_events', 'n_winters', 'r', 'Q', 'r_lo', 'r_hi', 'Q_lo', 'Q_hi'])
        for (li, h) in enumerate(LEADS):
            for (ti, tau) in enumerate(TAU):
                for (m, name) in enumerate(('B', 'M')):
                    writer.writerow([int(h), int(tau), name, int(e_count[li, ti]), int(e_winters[li, ti]), e_metrics[li, 0, m, ti], e_metrics[li, 1, m, ti], *e_bands[li, 0, m, ti], *e_bands[li, 1, m, ti]])
    with (args.output / 'FigS7_events_from_report.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['event_center_date', 'winter_start_year', 'source'])
        writer.writerows(((str(d), int(w), 'report Appendix A; verified against SSWevent.xlsx') for (d, w) in zip(events, event_winters)))
    with (args.output / 'FigS7_events_from_source.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['source_row', 'event_center_date', 'winter_start_year'])
        writer.writerows(((int(row), str(d), int(w)) for (row, d, w) in zip(source_rows, events, event_winters)))
    (centers, roll_metrics, roll_counts, roll_sums) = rolling_metrics(all_winter_sums)
    if len(centers) != 66:
        raise ValueError('Expected 66 ten-winter windows')
    np.testing.assert_allclose([centers[0], centers[-1]], [1953.5, 2018.5], atol=0)
    np.testing.assert_equal(roll_counts[0, 0], 1700)
    np.testing.assert_allclose(roll_metrics[0, :, 0, 0], [0.819250792, 0.3375227152], atol=1e-09, rtol=0)
    plot_s8(args.output / 'FigS8', centers, roll_metrics)
    np.savez_compressed(args.output / 'FigS8_scores.npz', leads=LEADS, midpoint_year=centers, metrics=roll_metrics, counts=roll_counts, window_sums=roll_sums)
    with (args.output / 'FigS8_scores.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['h', 'first_winter', 'last_winter', 'midpoint_year', 'model', 'n_pairs', 'r', 'Q'])
        for (li, h) in enumerate(LEADS):
            for (wi, center) in enumerate(centers):
                for (m, name) in enumerate(('B', 'M')):
                    writer.writerow([int(h), 1948 + wi, 1957 + wi, center, name, int(roll_counts[li, wi]), roll_metrics[li, 0, m, wi], roll_metrics[li, 1, m, wi]])
    qa = {'S7': {'source': 'SSWevent.xlsx Sheet1!B2:B52; verified against report Appendix A', 'source_audit': source_audit, 'n_events': len(events), 'n_ssw_winters': len(np.unique(event_winters)), 'ranges_by_h': {str(int(h)): {'events': [int(e_count[i].min()), int(e_count[i].max())], 'winters': [int(e_winters[i].min()), int(e_winters[i].max())]} for (i, h) in enumerate(LEADS)}, 'h30_tau0_n_events_winters': [int(e_count[2, 30]), int(e_winters[2, 30])], 'h30_tau0_M_r_Q': e_metrics[2, :, 1, 30].tolist()}, 'S8': {'n_windows': len(centers), 'midpoint_first_last': [centers[0], centers[-1]], 'first_window_h7_B_n': int(roll_counts[0, 0]), 'first_window_h7_B_r_Q': roll_metrics[0, :, 0, 0].tolist()}}
    (args.output / 'FigS7_S8_QA.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
    print(json.dumps(qa, indent=2), flush=True)
if __name__ == '__main__':
    main()
