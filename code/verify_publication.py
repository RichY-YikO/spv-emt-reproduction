from pathlib import Path
import argparse
import json
from pypdf import PdfReader

def audit_pdf(path, layout_dir=None):
    path = Path(path)
    if layout_dir is not None:
        layout_path = Path(layout_dir) / (path.stem + '.layout.json')
    elif path.with_suffix('.layout.json').is_file():
        layout_path = path.with_suffix('.layout.json')
    elif path.parent.resolve() == (Path(__file__).resolve().parents[1] / 'figures' / 'pdf').resolve():
        layout_path = Path(__file__).resolve().parents[1] / 'validation' / 'layout' / (path.stem + '.layout.json')
    else:
        layout_path = path.parent.parent / 'validation' / 'layout' / (path.stem + '.layout.json')
    reader = PdfReader(path)
    assert len(reader.pages) == 1, path
    images = []
    fonts = {}
    visited = set()

    def resources(res):
        if res is None:
            return
        res = res.get_object()
        for (name, ref) in res.get('/XObject', {}).items():
            obj = ref.get_object()
            key = (getattr(ref, 'idnum', None), getattr(ref, 'generation', None))
            if key in visited:
                continue
            visited.add(key)
            if obj.get('/Subtype') == '/Image':
                images.append(str(name))
            if obj.get('/Subtype') == '/Form':
                resources(obj.get('/Resources'))
        for (name, ref) in res.get('/Font', {}).items():
            font = ref.get_object()
            assert font.get('/Subtype') != '/Type3', f'Type3 font: {path}'
            descendants = font.get('/DescendantFonts', [font])
            for d in descendants:
                d = d.get_object()
                descriptor = d.get('/FontDescriptor')
                assert descriptor is not None, f'No font descriptor: {path}'
                descriptor = descriptor.get_object()
                assert any((k in descriptor for k in ['/FontFile', '/FontFile2', '/FontFile3'])), f'Unembedded font: {path}'
            fonts[str(name)] = str(font.get('/BaseFont'))
    for page in reader.pages:
        resources(page.get('/Resources'))
    assert not images, f'Raster images found: {path}: {images}'
    text = reader.pages[0].extract_text()
    assert 'non-SSW' not in text and 'non-ssw' not in text, path
    layout = json.loads(layout_path.read_text(encoding='utf-8'))
    for label in layout['panels']:
        assert label in text, f'Missing panel {label}: {path}'
    content = reader.pages[0].get_contents().get_data()
    assert b' m' in content and b' l' in content and (b'BT' in content), f'Missing vector paths/text: {path}'
    return {'file': str(path), 'pages': 1, 'raster_image_objects': 0, 'embedded_fonts': fonts, 'panel_count': len(layout['panels']), 'legends_inside': layout['legends_inside']}
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Audit vector PDF resources, fonts and panel labels without historical backups.')
    parser.add_argument('--output-root', type=Path, default=Path(__file__).resolve().parents[1] / 'figures' / 'pdf')
    parser.add_argument('--layout-dir', type=Path, help='Directory containing the figure layout validation JSON files.')
    args = parser.parse_args()
    paths = sorted(args.output_root.glob('Fig*.pdf'))
    if not paths:
        raise FileNotFoundError('No figure PDFs in ' + str(args.output_root))
    print(json.dumps([audit_pdf(p, args.layout_dir) for p in paths], indent=2))
