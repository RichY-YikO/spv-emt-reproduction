from pathlib import Path
import shutil

CODE = Path(__file__).resolve().parent
RESULTS = CODE / 'fig1_fig2_results'

def reference_dir(results=None):
    results = RESULTS if results is None else Path(results)
    return results.resolve().parent / 'reference_inputs'

def validation_dir(results=None):
    results = RESULTS if results is None else Path(results)
    return CODE.parent / 'validation' if results.resolve() == RESULTS.resolve() else results.resolve().parent / 'validation'

def pearson_reports(results=None):
    folder = validation_dir(results) / 'fig1_pearson'
    folder.mkdir(parents=True, exist_ok=True)
    return folder

def publish_figure(base, results=None):
    base = Path(base)
    validation = validation_dir(results)
    root = validation.parent
    for suffix, folder in [('.pdf', root / 'figures' / 'pdf'), ('.png', root / 'figures' / 'png'), ('.layout.json', validation / 'layout')]:
        folder.mkdir(parents=True, exist_ok=True)
        target = folder / (base.name + suffix)
        source = base.with_suffix(suffix)
        if source.resolve() != target.resolve():
            shutil.move(str(source), str(target))
